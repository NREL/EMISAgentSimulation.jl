"""
This function creates parallel workers for making price predictions.
"""
function create_parallel_workers(case::CaseDefinition, hpc::Bool)
    data_dir = get_data_dir(case)
    dir_name = joinpath(data_dir, "investors")
    investor_names = readdir(dir_name)

    investor_scenarios = AxisArrays.AxisArray(zeros(length(investor_names)), investor_names)

    for investor in investor_names
        investor_dir = joinpath(dir_name, investor)
        forecast_type = get_forecast_type(case)

        if forecast_type == "perfect"
            num_scenarios = 1
        elseif forecast_type == "imperfect"
            if get_uncertainty(case)
                file_name = joinpath(investor_dir, "markets_data", "scenario_data.csv")
                @assert isfile(file_name)
                scenario_df = read_data(file_name)
                num_scenarios = DataFrames.nrow(scenario_df)
            else
                num_scenarios = 1
            end
        end

        investor_scenarios[investor] = num_scenarios
    end

    parallelize_investors = get_parallel_investors(case)
    parallelize_scenarios = get_parallel_scenarios(case)

    if parallelize_investors
        if parallelize_scenarios
            num_workers_required = sum(investor_scenarios)
        else
            num_workers_required = length(investor_scenarios)
        end
    else
        if parallelize_scenarios
            num_workers_required = maximum(investor_scenarios)
        else
            num_workers_required = 0
        end
    end

    if num_workers_required > 0
        if hpc
            nodes = split(ENV["SLURM_NODELIST"], ",")
            num_procs = min(Int(ceil(num_workers_required / length(nodes))), 4)
            threads_per_worker = max(1, Sys.CPU_THREADS ÷ num_procs)
            node_pairs = [(n, num_procs) for n in nodes]
            Distributed.addprocs(node_pairs; exeflags = "--threads=$(threads_per_worker)")
        else
            num_workers = min(Int(num_workers_required), 4)
            threads_per_worker = max(1, Sys.CPU_THREADS ÷ num_workers)
            Distributed.addprocs(
                num_workers;
                lazy = false,
                exeflags = "--threads=$(threads_per_worker)",
            )
        end
    end

    return
end

"""
This function runs price prediction if investors are parallelized but scenarios are sequential.
"""
function parallelize_only_investors(investor::Investor,
    sys_data_dir::String,
    expected_portfolio::Vector{<: Project{<: BuildPhase}},
    rps_target::String,
    reserve_penalty::String,
    resource_adequacy::ResourceAdequacy,
    irm_scalar::Float64,
    zones::Vector{String},
    lines::Vector{ZonalLine},
    peak_load::Float64,
    average_capital_cost_multiplier::Float64,
    iteration_year::Int64,
    yearly_horizon::Int64,
    solver::JuMP.MOI.OptimizerWithAttributes,
    sys_results_dir::String,
    investor_name::String)
    investor_name,
    investor_dir,
    market_names,
    carbon_tax,
    reserve_products,
    ordc_products,
    rep_period_interval,
    rep_hour_weight,
    avg_block_size,
    fixed_block_size,
    chron_weights,
    scenarios = gather_prediction_parameters(investor, sys_data_dir, iteration_year)

    for scenario in scenarios
        create_expected_marketdata(investor_dir,
            sys_data_dir,
            market_names,
            carbon_tax,
            reserve_products,
            ordc_products,
            rps_target,
            reserve_penalty,
            resource_adequacy,
            irm_scalar,
            expected_portfolio,
            zones,
            lines,
            peak_load,
            rep_period_interval,
            rep_hour_weight,
            avg_block_size,
            fixed_block_size,
            chron_weights,
            average_capital_cost_multiplier,
            scenario,
            iteration_year,
            yearly_horizon,
            solver,
            sys_results_dir,
            investor_name)
    end

    return
end

"""
This function runs the construct_ordc function in parallel for different scenarios.
"""
function parallelize_ordc_construction(args)
    scenario,
    sys_UC,
    data_dir,
    investors,
    representative_periods,
    rep_period_interval,
    case,
    iteration_year,
    rolling_horizon,
    simulation_years = args
    for sim_year in
        collect(iteration_year:min(iteration_year + rolling_horizon - 1, simulation_years))
        construct_ordc(
            sys_UC,
            data_dir,
            scenario,
            sim_year,
            investors,
            0,
            representative_periods[scenario][sim_year],
            rep_period_interval,
            get_ordc_curved(case),
            get_ordc_unavailability_method(case),
            get_reserve_penalty(case),
        )
    end
end

"""
This function runs the update_delta_irm! function in parallel for different scenarios.
"""
function parallelize_update_delta_irm!(args)
    scenario,
    sys_PRAS,
    active_projects,
    capacity_forward_years,
    resource_adequacy,
    peak_load,
    static_capacity_bool,
    iteration_year,
    simulation_years,
    data_dir,
    rt_resolution,
    results_dir,
    outage_dir = args
    resource_adequacy = update_delta_irm!(
        sys_PRAS[scenario],
        active_projects,
        capacity_forward_years,
        resource_adequacy[scenario],
        peak_load[scenario][min(
            iteration_year + capacity_forward_years - 1,
            simulation_years,
        )],
        static_capacity_bool,
        scenario,
        iteration_year,
        data_dir,
        rt_resolution,
        results_dir,
        outage_dir,
        simulation_years,
    )
    return (scenario, resource_adequacy)
end

"""
This function creates repeated arguments for parallel runs
"""
function repeat_arguments(num_scenarios::Int, args...)
    repeated_args = []
    for arg in args
        push!(repeated_args, fill(arg, num_scenarios))
    end
    return repeated_args
end

"""
This function runs the update_simulation_derating_data! function in parallel for different scenarios.
"""
function parallelize_update_derating_data(args)
    scenario, simulation, iteration_year, methodology, ra_metric, marginal_cc = args
    update_simulation_derating_data!(
        simulation,
        scenario,
        iteration_year;
        methodology = methodology,
        ra_metric = ra_metric,
        marginal_cc = marginal_cc,
    )
    return
end
