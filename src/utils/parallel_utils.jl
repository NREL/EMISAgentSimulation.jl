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
                file_name = joinpath(investor_dir, "markets_data", "scenario_multiplier_data.csv")
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

          node_pairs = [(n, num_procs) for n in  nodes]
          Distributed.addprocs(node_pairs)
        else
          Distributed.addprocs(min(Int(num_workers_required), 4), lazy=false)
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
    rep_hour_weight,
    scenarios = gather_prediction_parameters(investor, sys_data_dir, iteration_year)

    for scenario in scenarios
        create_expected_marketdata(investor_dir,
                                market_names,
                                carbon_tax,
                                reserve_products,
                                ordc_products,
                                rps_target,
                                reserve_penalty,
                                expected_portfolio,
                                zones,
                                lines,
                                peak_load,
                                rep_hour_weight,
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
