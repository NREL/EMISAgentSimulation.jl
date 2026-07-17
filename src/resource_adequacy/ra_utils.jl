function add_outage_info!(
    PSY_gen::T,
    tech::Union{ThermalTech, RenewableTech, HydroTech, BatteryTech},
) where {T <: Union{PSY.Generator, PSY.Storage}}
    (λ, μ) = outage_to_rate((get_FOR(tech), get_MTTR(tech)))
    ext = PSY.get_ext(PSY_gen)
    ext["outage_probability"] = λ
    ext["recovery_probability"] = μ

    return
end

function get_generators(generator_type::String, projects)
    return sort!(
        filter(
            project -> occursin(generator_type, string(PSY.get_prime_mover_type(project))),
            projects,
        );
        by = x -> get_device_size(x),
    )
end

function get_availability_df(
    timeseries_data_dir::String,
    scenario::String,
    simulation_years::Int64,
    availability_type::String,
)
    availability_df = DataFrames.DataFrame()
    for sim_year in 1:simulation_years
        availability_df = vcat(
            availability_df,
            read_data(
                joinpath(
                    timeseries_data_dir,
                    scenario,
                    "sim_year_$(sim_year)",
                    "Availability",
                    "$(availability_type)_availability.csv",
                ),
            ),
        )
    end
    return availability_df
end

function slice_year(df::DataFrames.DataFrame, year::Int, simulation_years::Int)
    rows_per_year = DataFrames.nrow(df) ÷ simulation_years
    return df[((year - 1) * rows_per_year + 1):(year * rows_per_year), :]
end

function calculate_RA_metrics(sys::PSY.System,
    exportoutage::Bool,
    base_dir::String,
    outage_dir::String,
    iteration_year::Int64;
    samples::Int64 = PRAS_N_SAMPLES,
    seed::Int64 = 42,
    simulation_years::Int64 = 15)
    system_period_of_interest = range(1; length = DEFAULT_HOURS_PER_YEAR * simulation_years);
    # correlated_outage_csv_location = joinpath(outage_dir, "ThermalFOR_scenario_1_new.csv")

    @info "Calculating RA metrics for iteration year: $(iteration_year) with system_period_of_interest: $(system_period_of_interest)"

    # Build PRAS.SystemModel on the main process (needs PSY.System's live SQLite connection).
    # PRAS.SystemModel is plain arrays — safe to serialize and send to a remote worker.
    # generate_pras_system is in SiennaPRASInterface (SPI), not in PRASCore (PRAS).
    @timeit EMIS_TIMER "attach_outage_data" attach_outage_data_from_ext!(sys)
    pras_system =
        @timeit EMIS_TIMER "generate_pras_system" SPI.generate_pras_system(sys, PSY.Area)

    resultspec = Dict{String, Any}("shortfall" => PRAS.Shortfall())
    if exportoutage == true
        resultspec["gens_avail"] = PRAS.GeneratorAvailability()
    end

    ra_metrics = Dict{String, Float64}()
    @info "Running PRAS.assess ($(samples) samples, seed=$(seed))..."
    if !isnothing(PRAS_WORKER[])
        t_pras = @elapsed @timeit EMIS_TIMER "PRAS.assess (remote)" results_tuple =
            Distributed.remotecall_fetch(
                PRAS_WORKER[], pras_system, samples, seed,
            ) do pras_system, samples, seed
                PRAS.assess(pras_system,
                    PRAS.SequentialMonteCarlo(samples = samples, seed = seed),
                    values(resultspec)...)
            end
    else
        t_pras = @elapsed @timeit EMIS_TIMER "PRAS.assess (local)" results_tuple =
            PRAS.assess(pras_system,
                PRAS.SequentialMonteCarlo(samples = samples, seed = seed),
                values(resultspec)...)
    end

    results = Dict{String, Any}(zip(keys(resultspec), results_tuple))
    shortfall = results["shortfall"]
    eue_overall = PRAS.EUE(shortfall)
    lole_overall = PRAS.LOLE(shortfall)
    neue_overall = PRAS.NEUE(shortfall)

    ra_metrics["LOLE"] = val(lole_overall) / simulation_years
    ra_metrics["NEUE"] = val(neue_overall)
    @info "Finished PRAS simulation (elapsed: $(round(t_pras; digits=2))s)"
    @info "LOLE: $(ra_metrics["LOLE"])"
    @info "NEUE: $(ra_metrics["NEUE"])"

    if exportoutage
        gens_avail = results["gens_avail"]
        @info "Export outage profile from PRAS simulation... "
        scenarionum = 1
        df_outage = DataFrames.DataFrame()
        for (j, asset_name) in enumerate(gens_avail.generators)
            df_outage[!, asset_name] = Int.(gens_avail.available[j, :, scenarionum])
        end
        outage_csv_location=joinpath(base_dir, "GeneratorOutage")
        CSV.write(
            joinpath(outage_csv_location, "1/Generator_year$(iteration_year+1).csv"),
            df_outage;
            writeheader = true,
        )
    end

    PSY.set_units_base_system!(sys, PSY.IS.UnitSystem.DEVICE_BASE)
    return ra_metrics, shortfall
end

"""
This function does nothing if Device is not of RenewableGen type.
"""
function add_capacity_market_device_forecast!(sys_PRAS::PSY.System,
    device_PRAS::D,
    availability_raw::Vector{Float64},
    rt_resolution::Int64) where {D <: Union{PSY.ThermalGen, PSY.HydroGen, PSY.Storage}}
    return
end

"""
This function adds forecast timeseries to the future capacity market system if Device is of RenewableGen type.
"""
function add_capacity_market_device_forecast!(sys_PRAS::PSY.System,
    device_PRAS::D,
    availability_raw::Vector{Float64},
    rt_resolution::Int64) where {D <: PSY.RenewableGen}

    ######### Adding to PRAS System##########
    # time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
    #                                                 PSY.SingleTimeSeries,
    #                                                 first(PSY.get_components(PSY.ElectricLoad, sys_PRAS)),
    #                                                 "max_active_power"
    #                                                 )))
    sys_interval = PSY.get_forecast_interval(sys_PRAS)
    sys_horizon = PSY.get_forecast_horizon(sys_PRAS)
    sys_horizon = Dates.Hour(sys_horizon).value
    forecast_count = PSY.get_forecast_window_count(sys_PRAS)
    sys_resolution = PSY.get_time_series_resolutions(sys_PRAS)[1]
    start_datetime = PSY.get_forecast_initial_timestamp(sys_PRAS)
    finish_datetime =
        start_datetime + Dates.Hour((
            forecast_count * sys_interval/sys_resolution +
            (sys_horizon - sys_interval/sys_resolution) - 1
        ))
    time_stamps = StepRange(start_datetime, Dates.Hour(1), finish_datetime);
    additional_timestep = length(time_stamps) - DEFAULT_HOURS_PER_YEAR

    append!(
        availability_raw,
        availability_raw[(length(availability_raw) - additional_timestep + 1):end],
    )
    data = Dict(
        time_stamps[i] => availability_raw[i:(i + sys_horizon - 1)] for i in
        1:Int(sys_interval / sys_resolution):(length(time_stamps) - sys_horizon + 1)
    )
    forecast = PSY.Deterministic("max_active_power", data, Dates.Minute(rt_resolution))
    PSY.add_time_series!(sys_PRAS, device_PRAS, forecast)
    return
end

function add_capacity_market_project!(capacity_market_system::PSY.System,
    project::Project,
    simulation_dir::String,
    scenario::String,
    target_year::Int64,
    rt_resolution::Int64,
    simulation_years::Int64,
    timeseries_data_dir::String,
    availability_df_rt::DataFrame)

    # @info "Adding project $(get_name(project)) to capacity market system - scenario $(scenario) for year $(target_year)"

    PSY_project = create_PSY_generator(project, capacity_market_system)
    PSY.add_component!(capacity_market_system, PSY_project)
    add_nominal_outage_to_component!(capacity_market_system, PSY_project)

    for product in get_products(project)
        add_device_services!(capacity_market_system, PSY_project, product)
    end

    type = get_type(get_tech(project))
    zone = get_zone(get_tech(project))
    project_name = "$(type)_$(zone)"

    availability_raw_rt = ones(size(availability_df_rt, 1))
    if in(get_name(project), names(availability_df_rt))
        availability_raw_rt = availability_df_rt[:, Symbol(get_name(project))]
    elseif in(project_name, names(availability_df_rt))
        availability_raw_rt = availability_df_rt[:, Symbol(project_name)]
    else
        @warn "No availability data found for $(project_name), using default availability of 1.0"
    end

    add_device_forecast_PRAS!(
        capacity_market_system,
        PSY_project,
        availability_raw_rt,
        rt_resolution,
        simulation_years,
    )
    return
end

function create_capacity_mkt_system(initial_system::PSY.System,
    active_projects::Vector{Project},
    capacity_forward_years::Int64,
    scenario::String,
    iteration_year::Int64,
    simulation_dir::String,
    rt_resolution::Int64,
    simulation_years::Int64,
    timeseries_data_dir::String,
    availability_df_rt::DataFrame)
    @info "Creating Forward Capacity Market System"
    capacity_market_system = deepcopy(initial_system)

    capacity_market_year = iteration_year + capacity_forward_years - 1
    capacity_market_projects = Project[]
    option_leaftypes = leaftypes(Project{Option})
    non_option_projects =
        filter(project -> !in(typeof(project), option_leaftypes), active_projects)

    for project in non_option_projects
        # @info "Evaluating project $(get_name(project)) for inclusion in capacity market system for year $(capacity_market_year)"
        end_life_year = get_end_life_year(project)
        construction_year = get_construction_year(project)

        if end_life_year >= capacity_market_year &&
           construction_year <= capacity_market_year
            push!(capacity_market_projects, project)
            if !(get_name(project) in PSY.get_name.(get_all_techs(capacity_market_system)))
                add_capacity_market_project!(capacity_market_system, project,
                    simulation_dir, scenario,
                    capacity_market_year, rt_resolution, simulation_years,
                    timeseries_data_dir, availability_df_rt)
            end
        end
    end

    for device in get_all_techs(capacity_market_system)
        if !(PSY.get_name(device) in get_name.(capacity_market_projects))
            PSY.remove_component!(capacity_market_system, device)
        end
    end

    #=
    nodal_loads = PSY.get_components(PSY.StandardLoad, capacity_market_system)
    for load in nodal_loads
        zone = "zone_$(PSY.get_name(PSY.get_area(PSY.get_bus(load))))"
        scaled_active_power = deepcopy(PSY.get_max_active_power(load)) * (1 + load_growth["load_$(zone)"]) ^ (capacity_forward_years)
        PSY.set_max_active_power!(load, scaled_active_power)

    end
    =#
    return capacity_market_system
end

function check_ra_conditions(
    ra_targets::Dict{String, Float64},
    ra_metrics::Dict{String, Float64},
)
    metrics = keys(ra_targets)

    adequacy_conditions = falses(length(metrics))
    scarcity_conditions = falses(length(metrics))
    for (idx, metric) in enumerate(metrics)
        if ra_metrics[metric] <= ra_targets[metric]
            adequacy_conditions[idx] = true
        else
            scarcity_conditions[idx] = true
        end
    end

    adequacy_conditions_met = prod(adequacy_conditions)
    scarcity_conditions_met = true
    if sum(scarcity_conditions) < 1
        scarcity_conditions_met = false
    end

    return adequacy_conditions_met, scarcity_conditions_met
end

function update_delta_irm!(initial_system::PSY.System,
    active_projects::Vector{Project},
    capacity_forward_years::Int64,
    resource_adequacy::ResourceAdequacy,
    forward_peak_load::Float64,
    static_capacity_market::Bool,
    scenario::String,
    iteration_year::Int64,
    simulation_dir::String,
    rt_resolution::Int64,
    results_dir::String,
    outage_dir::String,
    simulation_years::Int64)
    timeseries_data_dir = joinpath(results_dir, "timeseries_data_files")
    availability_df_rt =
        get_availability_df(timeseries_data_dir, scenario, simulation_years, "REAL_TIME")

    if !(static_capacity_market)
        capacity_market_year = iteration_year + capacity_forward_years - 1

        capacity_market_system = create_capacity_mkt_system(initial_system,
            active_projects,
            capacity_forward_years,
            scenario,
            iteration_year,
            simulation_dir,
            rt_resolution,
            simulation_years,
            timeseries_data_dir,
            availability_df_rt)

        ra_targets = get_targets(resource_adequacy)
        delta_irm = 0.0
        all_capacity_market_projects = get_all_techs(capacity_market_system)
        removeable_projects = PSY.Generator[]

        CT_generators = get_generators("CT", all_capacity_market_projects)
        append!(removeable_projects, CT_generators)
        CC_generators = get_generators("CC", all_capacity_market_projects)
        append!(removeable_projects, CC_generators)

        @info "Updating delta IRM for scenario: $(scenario) - Year: $(iteration_year)"
        @time begin
            if !isempty(ra_targets)
                ra_metrics, shortfall = calculate_RA_metrics(capacity_market_system,
                    false, results_dir, outage_dir,
                    iteration_year;
                    samples = PRAS_N_SAMPLES,
                    simulation_years = simulation_years,
                )

                @info "RA metrics: $(ra_metrics)"

                adequacy_conditions_met, scarcity_conditions_met =
                    check_ra_conditions(ra_targets, ra_metrics)
                @info "Adequacy conditions met: $(adequacy_conditions_met), Scarcity conditions met: $(scarcity_conditions_met)"
                count = 1
                total_added_capacity = 0.0
                total_removed_capacity = 0.0
                removed_capacity = 0.0
                ct_project_template =
                    first(filter(p -> occursin("new_CT", get_name(p)), active_projects))

                if !(adequacy_conditions_met)
                    while !(adequacy_conditions_met)
                        @info "RA metrics: $(ra_metrics)"
                        @info "Updating delta IRM for scenario: $(scenario) - Year: $(iteration_year)"
                        scalar = 2
                        ratio = 0
                        for metric in keys(ra_targets)
                            @info "Updating $(metric)"
                            ratio +=
                                (ra_metrics[metric] - ra_targets[metric]) /
                                ra_targets[metric]
                        end
                        ratio = max(1, scalar * ratio / length(keys(ra_targets)))

                        for i in 1:ceil(ratio)
                            # @info i
                            @info "Updating delta IRM for scenario: $(scenario) - Year: $(iteration_year)"
                            incremental_project = deepcopy(ct_project_template)
                            set_name!(incremental_project, "addition_CT_project_$(count)")
                            total_added_capacity += get_maxcap(incremental_project)
                            add_capacity_market_project!(
                                capacity_market_system,
                                incremental_project,
                                simulation_dir,
                                scenario,
                                capacity_market_year,
                                rt_resolution,
                                simulation_years,
                                timeseries_data_dir,
                                availability_df_rt,
                            )
                            count += 1
                        end

                        ra_metrics, shortfall = calculate_RA_metrics(
                            capacity_market_system,
                            false,
                            results_dir,
                            outage_dir,
                            iteration_year;
                            samples = PRAS_N_SAMPLES,
                            simulation_years = simulation_years,
                        )
                        adequacy_conditions_met, scarcity_conditions_met =
                            check_ra_conditions(ra_targets, ra_metrics)
                    end

                elseif !(scarcity_conditions_met)
                    while !(scarcity_conditions_met) && (total_removed_capacity <= 400)
                        if !(isempty(removeable_projects))
                            removed_project = removeable_projects[1]
                            popfirst!(removeable_projects)
                            removed_capacity =
                                get_device_size(removed_project) *
                                PSY.get_base_power(removed_project)
                            total_removed_capacity += removed_capacity
                            PSY.remove_component!(capacity_market_system, removed_project)
                            ra_metrics, shortfall = calculate_RA_metrics(
                                capacity_market_system,
                                false,
                                results_dir,
                                outage_dir,
                                iteration_year;
                                samples = PRAS_N_SAMPLES,
                                simulation_years = simulation_years,
                            )
                            #println(ra_metrics)
                            adequacy_conditions_met, scarcity_conditions_met =
                                check_ra_conditions(ra_targets, ra_metrics)
                            count += 1
                        end
                    end
                    total_removed_capacity -= removed_capacity
                end

                delta_irm =
                    (total_added_capacity - total_removed_capacity) / forward_peak_load
            end
        end

        @info "delta_irm is $(delta_irm)."
        set_delta_irm!(resource_adequacy, iteration_year, delta_irm)
    end

    return resource_adequacy
end

function create_base_system(initial_system::PSY.System,
    active_projects::Vector{Project},
    capacity_forward_years::Int64,
    scenario::String,
    resource_adequacy::ResourceAdequacy,
    iteration_year::Int64,
    simulation_dir::String,
    outage_dir::String,
    rt_resolution::Int64,
    simulation::Union{AgentSimulation, AgentSimulationData},
    availability_df_rt::DataFrame)
    timeseries_data_dir = joinpath(get_results_dir(simulation), "timeseries_data_files")
    simulation_years = get_simulation_years(get_case(simulation))
    capacity_market_year = iteration_year + capacity_forward_years - 1

    case = get_case(simulation)
    simulation_years = get_simulation_years(case)

    capacity_market_system = create_capacity_mkt_system(initial_system,
        active_projects,
        capacity_forward_years,
        scenario,
        iteration_year,
        simulation_dir,
        rt_resolution,
        get_total_horizon(get_case(simulation)),
        timeseries_data_dir,
        availability_df_rt,
    )

    ra_targets = get_targets(resource_adequacy)

    @info "RA targets: $(ra_targets)"

    all_capacity_market_projects = get_all_techs(capacity_market_system)
    removeable_projects = PSY.Generator[]

    CT_generators = get_generators("CT", all_capacity_market_projects)
    append!(removeable_projects, CT_generators)
    CC_generators = get_generators("CC", all_capacity_market_projects)
    append!(removeable_projects, CC_generators)

    @time begin
        if !isempty(ra_targets)
            ra_metrics, shortfall = calculate_RA_metrics(
                capacity_market_system,
                false,
                get_results_dir(simulation),
                outage_dir,
                iteration_year;
                samples = PRAS_N_SAMPLES,
                simulation_years = simulation_years,
            )

            @info "RA metrics: $(ra_metrics)"
            adequacy_conditions_met, scarcity_conditions_met =
                check_ra_conditions(ra_targets, ra_metrics)

            count = 1
            total_added_capacity = 0.0
            total_removed_capacity = 0.0
            removed_capacity = 0.0
            ct_project_template =
                first(filter(p -> occursin("new_CT", get_name(p)), active_projects))

            if !(adequacy_conditions_met)
                while !(adequacy_conditions_met)
                    scalar = 2
                    ratio = 0
                    for metric in keys(ra_targets)
                        ratio +=
                            (ra_metrics[metric] - ra_targets[metric]) / ra_targets[metric]
                    end
                    ratio = max(1, scalar * ratio / length(keys(ra_targets)))
                    @info "adding $ratio new CTs to system to meet adequacy targets"
                    for i in 1:ceil(ratio)
                        incremental_project = deepcopy(ct_project_template)
                        set_name!(incremental_project, "addition_CT_project_$(count)")
                        total_added_capacity += get_maxcap(incremental_project)
                        add_capacity_market_project!(
                            capacity_market_system,
                            incremental_project,
                            simulation_dir,
                            scenario,
                            capacity_market_year,
                            rt_resolution,
                            get_total_horizon(get_case(simulation)),
                            timeseries_data_dir,
                            availability_df_rt,
                        )
                        count += 1
                    end

                    ra_metrics, shortfall = calculate_RA_metrics(
                        capacity_market_system,
                        false,
                        get_results_dir(simulation),
                        outage_dir,
                        iteration_year;
                        samples = PRAS_N_SAMPLES,
                        simulation_years = simulation_years,
                    )
                    @info "Added Capacity"
                    @info "RA metrics after adding capacity: $(ra_metrics)"
                    adequacy_conditions_met, scarcity_conditions_met =
                        check_ra_conditions(ra_targets, ra_metrics)
                end

            elseif !(scarcity_conditions_met)
                while !(scarcity_conditions_met)
                    if !(isempty(removeable_projects))
                        removed_project = removeable_projects[1]
                        popfirst!(removeable_projects)
                        removed_capacity =
                            get_device_size(removed_project) *
                            PSY.get_base_power(removed_project)
                        total_removed_capacity += removed_capacity
                        @info "Removing $(get_name(removed_project)) with capacity $(removed_capacity) MW to meet scarcity targets"
                        PSY.remove_component!(capacity_market_system, removed_project)
                        ra_metrics, shortfall = calculate_RA_metrics(
                            capacity_market_system,
                            false,
                            get_results_dir(simulation),
                            outage_dir,
                            iteration_year;
                            samples = PRAS_N_SAMPLES,
                            simulation_years = simulation_years,
                        )
                        @info "Base System: removed capacity: $(removed_capacity) MW from project $(get_name(removed_project))"
                        @info "RA metrics after removal: $(ra_metrics)"
                        adequacy_conditions_met, scarcity_conditions_met =
                            check_ra_conditions(ra_targets, ra_metrics)
                        count += 1
                    end
                end
                total_removed_capacity -= removed_capacity
            end
        end
    end

    return capacity_market_system
end
