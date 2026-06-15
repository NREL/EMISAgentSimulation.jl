function add_outage_info!(
                          PSY_gen::T,
                          tech::Union{ThermalTech, RenewableTech, HydroTech, BatteryTech}
                          ) where T <: Union{PSY.Generator, PSY.Storage}
    (λ, μ) = outage_to_rate((get_FOR(tech), get_MTTR(tech)))
    ext = PSY.get_ext(PSY_gen)
    ext["outage_probability"] = λ
    ext["recovery_probability"] = μ

    return
end


function calculate_RA_metrics(sys::PSY.System,
                              exportoutage::Bool,
                              base_dir::String,
                              outage_dir::String,
                              iteration_year::Int64;
                              samples::Int64 = 100,
                              seed::Int64 = 42)

    system_period_of_interest = range(1, length = 8760 * 15);
    correlated_outage_csv_location = joinpath(outage_dir, "ThermalFOR_scenario_1_new.csv")
    pras_system = make_pras_system(sys,
                                    system_model="Single-Node",
                                    aggregation="Area",
                                    period_of_interest = system_period_of_interest,
                                    outage_flag=false,
                                    lump_pv_wind_gens=false,
                                    availability_flag=true,
                                    outage_csv_location = correlated_outage_csv_location);

    total_load = calculate_total_load(sys, 60)

    ra_metrics = Dict{String, Float64}()
    if !isnothing(PRAS_WORKER[])
        shortfall, gens_avail = Distributed.remotecall_fetch(
            PRAS_WORKER[], pras_system, samples, seed
        ) do pras_system, samples, seed
            PRAS.assess(pras_system,
                PRAS.SequentialMonteCarlo(samples = samples, seed = seed),
                PRAS.Shortfall(), PRAS.GeneratorAvailability())
        end
    else
        shortfall, gens_avail = @time PRAS.assess(pras_system,
            PRAS.SequentialMonteCarlo(samples = samples, seed = seed),
            PRAS.Shortfall(), PRAS.GeneratorAvailability())
    end
    @info "Finished PRAS simulation... "
    eue_overall = PRAS.EUE(shortfall)
    lole_overall = PRAS.LOLE(shortfall)

    if exportoutage == true
        @info "Export outage profile from PRAS simulation... "
        scenarionum = 1
        df_outage = DataFrames.DataFrame()
        for (j,asset_name) in enumerate(gens_avail.generators)
            df_outage[!,asset_name] = Int.(gens_avail.available[j,:,scenarionum])
        end
        outage_csv_location=joinpath(base_dir,"GeneratorOutage") #get_base_dir(case)
        CSV.write(joinpath(outage_csv_location,"1/Generator_year$(iteration_year+1).csv"), df_outage,writeheader = true)
    end

    ra_metrics["LOLE"] = val(lole_overall) / 15
    ra_metrics["NEUE"] = val(eue_overall) * 1e6 / total_load

    PSY.set_units_base_system!(sys, PSY.IS.UnitSystem. DEVICE_BASE)

    return ra_metrics, shortfall

end

"""
This function does nothing if Device is not of RenewableGen type.
"""
function add_capacity_market_device_forecast!(sys_PRAS::PSY.System,
                                                device_PRAS::D,
                                                availability_raw::Vector{Float64},
                                                rt_resolution::Int64) where D <: Union{PSY.ThermalGen, PSY.HydroGen, PSY.Storage}

    return
end

"""
This function adds forecast timeseries to the future capacity market system if Device is of RenewableGen type.
"""
function add_capacity_market_device_forecast!(sys_PRAS::PSY.System,
                                                device_PRAS::D,
                                                availability_raw::Vector{Float64},
                                                rt_resolution::Int64) where D <: PSY.RenewableGen

    ######### Adding to PRAS System##########
    # time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
    #                                                 PSY.SingleTimeSeries,
    #                                                 first(PSY.get_components(PSY.ElectricLoad, sys_PRAS)),
    #                                                 "max_active_power"
    #                                                 )))
    sys_interval = sys_PRAS.data.time_series_params.forecast_params.interval
    sys_horizon = sys_PRAS.data.time_series_params.forecast_params.horizon
    forecast_count = sys_PRAS.data.time_series_params.forecast_params.count
    sys_resolution = sys_PRAS.data.time_series_params.resolution
    start_datetime = sys_PRAS.data.time_series_params.forecast_params.initial_timestamp
    finish_datetime = start_datetime + Dates.Hour((forecast_count * sys_interval/sys_resolution + (sys_horizon - sys_interval/sys_resolution) - 1))
    time_stamps = StepRange(start_datetime, Dates.Hour(1), finish_datetime);

    additional_timestep = length(time_stamps) - 8760

    # intervals = Int(36 * 60 / rt_resolution)
    append!(availability_raw, availability_raw[(length(availability_raw) - additional_timestep + 1):end])
    data = Dict(time_stamps[i] => availability_raw[i:(i + sys_horizon - 1)] for i in 1:Int(sys_interval/sys_resolution):(length(time_stamps)-sys_horizon + 1))
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
                                    simulation_years::Int64)
    PSY_project = create_PSY_generator(project, capacity_market_system)

    PSY.add_component!(capacity_market_system, PSY_project)

    for product in get_products(project)
        add_device_services!(capacity_market_system, PSY_project, product)
    end

    type = get_type(get_tech(project))
    zone = get_zone(get_tech(project))

    #max_year = maximum(map(s -> parse(Int, filter(x -> !isempty(x) && all(isdigit, x), split(s, "_"))[end]), readdir(joinpath(simulation_dir, "timeseries_data_files", scenario))))

    availability_df_rt = DataFrames.DataFrame()

    for sim_year in 1:simulation_years
        availability_df_rt = vcat(availability_df_rt, read_data(joinpath(simulation_dir, "timeseries_data_files", scenario, "sim_year_$(sim_year)", "Availability", "REAL_TIME_availability.csv")))
    end

    if in(get_name(project), names(availability_df_rt))
        availability_raw_rt = availability_df_rt[:, Symbol(get_name(project))]
    elseif in("$(type)_$(zone)", names(availability_df_rt))
        availability_raw_rt = availability_df_rt[:, Symbol("$(type)_$(zone)")]
    end

    add_device_forecast_PRAS!(capacity_market_system, PSY_project, availability_raw_rt, rt_resolution)

    return
end

function create_capacity_mkt_system(initial_system::PSY.System,
                                    active_projects::Vector{Project},
                                    capacity_forward_years::Int64,
                                    scenario::String,
                                    iteration_year::Int64,
                                    simulation_dir::String,
                                    rt_resolution::Int64,
                                    simulation_years::Int64)

    println("Creating Forward Capacity Market System")
    capacity_market_system = deepcopy(initial_system)

    capacity_market_year = iteration_year + capacity_forward_years - 1
    capacity_market_projects = Project[]
    option_leaftypes = leaftypes(Project{Option})
    non_option_projects = filter(project -> !in(typeof(project), option_leaftypes), active_projects)

    for project in non_option_projects
        end_life_year = get_end_life_year(project)
        construction_year = get_construction_year(project)
        if end_life_year >= capacity_market_year && construction_year <= capacity_market_year
            push!(capacity_market_projects, project)
            if !(get_name(project) in PSY.get_name.(get_all_techs(capacity_market_system)))
                add_capacity_market_project!(capacity_market_system, project, simulation_dir, scenario, capacity_market_year, rt_resolution, simulation_years)
            end
        end

    end

    for device in get_all_techs(capacity_market_system)
        if !(PSY.get_name(device) in get_name.(capacity_market_projects))
            PSY.remove_component!(capacity_market_system, device)
        end
    end

    nodal_loads = PSY.get_components(PSY.PowerLoad, capacity_market_system)

    @warn "UPDATE PSY TIMESERIES!"
    #= for load in nodal_loads
        zone = "zone_$(PSY.get_name(PSY.get_area(PSY.get_bus(load))))"
        scaled_active_power = deepcopy(PSY.get_max_active_power(load)) * (1 + load_growth["load_$(zone)"]) ^ (capacity_forward_years)
        PSY.set_max_active_power!(load, scaled_active_power)

    end
 =#
    return capacity_market_system

end

function check_ra_conditions(ra_targets::Dict{String, Float64}, ra_metrics::Dict{String, Float64})

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
    

    if !(static_capacity_market)

        capacity_market_year = iteration_year + capacity_forward_years - 1

        capacity_market_system = create_capacity_mkt_system(initial_system,
                                                            active_projects,
                                                            capacity_forward_years,
                                                            scenario,
                                                            iteration_year,
                                                            simulation_dir,
                                                            rt_resolution,
                                                            simulation_years)

        ra_targets = get_targets(resource_adequacy)
        delta_irm = 0.0

        all_capacity_market_projects = get_all_techs(capacity_market_system)
        removeable_projects = PSY.Generator[]

        CT_generators = sort!(filter(project -> occursin("CT", string(PSY.get_prime_mover_type(project))), all_capacity_market_projects), by = x -> get_device_size(x))
        append!(removeable_projects, CT_generators)
        CC_generators = sort!(filter(project -> occursin("CC", string(PSY.get_prime_mover_type(project))), all_capacity_market_projects), by = x -> get_device_size(x))
        append!(removeable_projects, CC_generators)

        @time begin
        if !isempty(ra_targets)
            ra_metrics, shortfall = calculate_RA_metrics(capacity_market_system, false, results_dir, outage_dir, iteration_year)
            #println(ra_metrics)
            adequacy_conditions_met, scarcity_conditions_met = check_ra_conditions(ra_targets, ra_metrics)

            count = 1
            total_added_capacity = 0.0
            total_removed_capacity = 0.0
            removed_capacity = 0.0

            if !(adequacy_conditions_met)
                while !(adequacy_conditions_met)

                    scalar = 2
                    ratio = 0
                    for metric in keys(ra_targets)
                        ratio += (ra_metrics[metric] - ra_targets[metric]) / ra_targets[metric] 
                    end
                    ratio = max(1, scalar * ratio /  length(keys(ra_targets)))

                    for i in 1:ceil(ratio)
                        incremental_project = deepcopy(first(filter(p -> occursin("new_CT", get_name(p)), active_projects)))
                        set_name!(incremental_project, "addition_CT_project_$(count)")
                        total_added_capacity += get_maxcap(incremental_project)
                        add_capacity_market_project!(capacity_market_system, incremental_project, simulation_dir, scenario, capacity_market_year, rt_resolution, simulation_years)
                        count += 1
                    end
                    
                    ra_metrics, shortfall = calculate_RA_metrics(capacity_market_system, false, results_dir, outage_dir, iteration_year)
                    #println(ra_metrics)
                    adequacy_conditions_met, scarcity_conditions_met = check_ra_conditions(ra_targets, ra_metrics)
                end

            elseif !(scarcity_conditions_met)
                while !(scarcity_conditions_met) && (total_removed_capacity <= 400)
                    if !(isempty(removeable_projects))
                        removed_project = removeable_projects[1]
                        popfirst!(removeable_projects)
                        removed_capacity = get_device_size(removed_project) * PSY.get_base_power(removed_project)
                        total_removed_capacity += removed_capacity
                        PSY.remove_component!(capacity_market_system, removed_project)
                        ra_metrics, shortfall = calculate_RA_metrics(capacity_market_system, false, results_dir, outage_dir, iteration_year)
                        #println(ra_metrics)
                        adequacy_conditions_met, scarcity_conditions_met = check_ra_conditions(ra_targets, ra_metrics)
                        count += 1
                    end
                end
                total_removed_capacity -= removed_capacity
            end

            delta_irm = (total_added_capacity - total_removed_capacity) / forward_peak_load
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
    simulation::Union{AgentSimulation,AgentSimulationData})

    capacity_market_year = iteration_year + capacity_forward_years - 1

    capacity_market_system = create_capacity_mkt_system(initial_system,
                                                        active_projects,
                                                        capacity_forward_years,
                                                        scenario,
                                                        iteration_year,
                                                        simulation_dir,
                                                        rt_resolution,
                                                        get_total_horizon(get_case(simulation)))

    ra_targets = get_targets(resource_adequacy)
    println(ra_targets)

    all_capacity_market_projects = get_all_techs(capacity_market_system)
    removeable_projects = PSY.Generator[]

    CT_generators = sort!(filter(project -> occursin("CT", string(PSY.get_prime_mover_type(project))), all_capacity_market_projects), by = x -> get_device_size(x))
    append!(removeable_projects, CT_generators)
    CC_generators = sort!(filter(project -> occursin("CC", string(PSY.get_prime_mover_type(project))), all_capacity_market_projects), by = x -> get_device_size(x))
    append!(removeable_projects, CC_generators)

    @time begin
    if !isempty(ra_targets)
        ra_metrics, shortfall = calculate_RA_metrics(capacity_market_system, false, get_results_dir(simulation), outage_dir, iteration_year;samples = 100)

        println(ra_metrics)
        adequacy_conditions_met, scarcity_conditions_met = check_ra_conditions(ra_targets, ra_metrics)

        count = 1
        total_added_capacity = 0.0
        total_removed_capacity = 0.0
        removed_capacity = 0.0

        if !(adequacy_conditions_met)
            while !(adequacy_conditions_met)
                scalar = 2
                ratio = 0
                for metric in keys(ra_targets)
                    ratio += (ra_metrics[metric] - ra_targets[metric]) / ra_targets[metric] 
                end
                ratio = max(1, scalar * ratio /  length(keys(ra_targets)))

                for i in 1:ceil(ratio)
                    incremental_project = deepcopy(first(filter(p -> occursin("new_CT", get_name(p)), active_projects)))
                    set_name!(incremental_project, "addition_CT_project_$(count)")
                    total_added_capacity += get_maxcap(incremental_project)
                    add_capacity_market_project!(capacity_market_system, incremental_project, simulation_dir, scenario, capacity_market_year, rt_resolution, get_total_horizon(get_case(simulation)))
                    count += 1
                end
                
                ra_metrics, shortfall = calculate_RA_metrics(capacity_market_system, false,get_results_dir(simulation), outage_dir, iteration_year;samples = 100)
                println("Added Capacity")
                println(ra_metrics)
                adequacy_conditions_met, scarcity_conditions_met = check_ra_conditions(ra_targets, ra_metrics)
                
            end

        elseif !(scarcity_conditions_met)
            while !(scarcity_conditions_met)
                if !(isempty(removeable_projects))
                    removed_project = removeable_projects[1]
                    popfirst!(removeable_projects)
                    removed_capacity = get_device_size(removed_project) * PSY.get_base_power(removed_project)
                    total_removed_capacity += removed_capacity
                    PSY.remove_component!(capacity_market_system, removed_project)
                    ra_metrics, shortfall = calculate_RA_metrics(capacity_market_system, false,get_results_dir(simulation), outage_dir, iteration_year;samples = 100)
                    println("Removed Capacity")
                    println(ra_metrics)
                    adequacy_conditions_met, scarcity_conditions_met = check_ra_conditions(ra_targets, ra_metrics)
                    count += 1
                end
            end
            total_removed_capacity -= removed_capacity
        end
    end
    end

    return capacity_market_system
end
