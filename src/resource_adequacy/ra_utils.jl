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


function calculate_RA_metrics(sys::PSY.System)

    system_period_of_interest = range(1, length = 8760);

    pras_system = make_pras_system(sys,
                                   system_model = "Single-Node",
                                   aggregation = "Area",
                                   period_of_interest = system_period_of_interest,
                                   outage_flag = false);

    total_load = calculate_total_load(sys, 60)

    ra_metrics = Dict{String, Float64}()
    seed = 3
    shortfall, = @time PRAS.assess(pras_system,  PRAS.SequentialMonteCarlo(samples = 100, seed=seed),  PRAS.Shortfall())
    @info "Finished PRAS simulation... "
    eue_overall = PRAS.EUE(shortfall)
    lole_overall = PRAS.LOLE(shortfall)

    ra_metrics["LOLE"] = val(lole_overall)
    ra_metrics["NEUE"] = val(eue_overall) * 1e6 / total_load

    PSY.set_units_base_system!(sys, PSY.IS.UnitSystem. DEVICE_BASE)

    return ra_metrics

end

"""
This function does nothing if Device is not of RenewableGen type.
"""
function add_capacity_market_device_forecast!(sys_UC::PSY.System,
                                                device_UC::D,
                                                availability_raw::Vector{Float64},
                                                da_resolution::Int64) where D <: Union{PSY.ThermalGen, PSY.HydroGen, PSY.Storage}

    return
end

"""
This function adds forecast timeseries to the future capacity market system if Device is of RenewableGen type.
"""
function add_capacity_market_device_forecast!(sys_UC::PSY.System,
                                                device_UC::D,
                                                availability_raw::Vector{Float64},
                                                da_resolution::Int64) where D <: PSY.RenewableGen

    ######### Adding to UC##########
    time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
                                                    PSY.SingleTimeSeries,
                                                    first(PSY.get_components(PSY.ElectricLoad, sys_UC)),
                                                    "max_active_power"
                                                    )))

    intervals = Int(24 * 60 / da_resolution)
    append!(availability_raw, availability_raw[(length(availability_raw) - intervals + 1):end])
    data = Dict(time_stamps[i] => availability_raw[i:(i + intervals - 1)] for i in 1:intervals:length(time_stamps))
    forecast = PSY.Deterministic("max_active_power", data, Dates.Minute(da_resolution))
    PSY.add_time_series!(sys_UC, device_UC, forecast)

    return
end

function add_capacity_market_project!(capacity_market_system::PSY.System,
                                    project::Project,
                                    simulation_dir::String,
                                    da_resolution::Int64)
    PSY_project = create_PSY_generator(project, capacity_market_system)

    PSY.add_component!(capacity_market_system, PSY_project)

    for product in get_products(project)
        add_device_services!(capacity_market_system, PSY_project, product)
    end

    type = get_type(get_tech(project))
    zone = get_zone(get_tech(project))

    availability_df = read_data(joinpath(simulation_dir, "timeseries_data_files", "Availability", "DAY_AHEAD_availability.csv"))

    if in(get_name(project), names(availability_df))
        availability_raw = availability_df[:, Symbol(get_name(project))]
    elseif in("$(type)_$(zone)", names(availability_df))
            availability_raw = availability_df[:, Symbol("$(type)_$(zone)")]
    end

    add_capacity_market_device_forecast!(capacity_market_system, PSY_project, availability_raw, da_resolution)
end

function create_capacity_mkt_system(initial_system::PSY.System,
                                    active_projects::Vector{Project},
                                    capacity_forward_years::Int64,
                                    iteration_year::Int64,
                                    load_growth::AxisArrays.AxisArray{Float64, 1},
                                    simulation_dir::String,
                                    da_resolution::Int64)

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
                add_capacity_market_project!(capacity_market_system, project, simulation_dir, da_resolution)
            end
        end

    end

    for device in get_all_techs(capacity_market_system)
        if !(PSY.get_name(device) in get_name.(capacity_market_projects))
            PSY.remove_component!(capacity_market_system, device)
        end
    end

    nodal_loads = PSY.get_components(PSY.ElectricLoad, capacity_market_system)

    for load in nodal_loads
        zone = "zone_$(PSY.get_name(PSY.get_area(PSY.get_bus(load))))"
        scaled_active_power = deepcopy(PSY.get_max_active_power(load)) * (1 + load_growth["load_$(zone)"]) ^ (capacity_forward_years)
        PSY.set_max_active_power!(load, scaled_active_power)

    end

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
                            peak_load::Float64,
                            static_capacity_market::Bool,
                            iteration_year::Int64,
                            load_growth::AxisArrays.AxisArray{Float64, 1},
                            simulation_dir::String,
                            da_resolution::Int64)

    if !(static_capacity_market)
        forward_peak_load = peak_load * (1 + Statistics.mean(load_growth)) ^ (capacity_forward_years)
        capacity_market_system = create_capacity_mkt_system(initial_system,
                                                            active_projects,
                                                            capacity_forward_years,
                                                            iteration_year,
                                                            load_growth,
                                                            simulation_dir,
                                                            da_resolution)

        ra_targets = get_targets(resource_adequacy)
        delta_irm = 0.0

        all_capacity_market_projects = get_all_techs(capacity_market_system)
        removeable_projects = PSY.Generator[]

        CT_generators = sort!(filter(project -> occursin("CT", string(PSY.get_prime_mover(project))), all_capacity_market_projects), by = x -> get_device_size(x))
        append!(removeable_projects, CT_generators)
        CC_generators = sort!(filter(project -> occursin("CC", string(PSY.get_prime_mover(project))), all_capacity_market_projects), by = x -> get_device_size(x))
        append!(removeable_projects, CC_generators)

        @time begin
        if !isempty(ra_targets)
            ra_metrics = calculate_RA_metrics(capacity_market_system)
            println(ra_metrics)
            adequacy_conditions_met, scarcity_conditions_met = check_ra_conditions(ra_targets, ra_metrics)

            count = 1
            total_added_capacity = 0.0
            total_removed_capacity = 0.0
            removed_capacity = 0.0

            if !(adequacy_conditions_met)
                while !(adequacy_conditions_met)
                    incremental_project = deepcopy(first(filter(p -> occursin("new_CT", get_name(p)), active_projects)))
                    set_name!(incremental_project, "addition_CT_project_$(count)")
                    total_added_capacity += get_maxcap(incremental_project)
                    add_capacity_market_project!(capacity_market_system, incremental_project, simulation_dir, da_resolution)
                    ra_metrics = calculate_RA_metrics(capacity_market_system)
                    println(ra_metrics)
                    adequacy_conditions_met, scarcity_conditions_met = check_ra_conditions(ra_targets, ra_metrics)
                    count += 1
                end

            elseif !(scarcity_conditions_met)
                while !(scarcity_conditions_met) && (total_removed_capacity <= 400)
                    if !(isempty(removeable_projects))
                        removed_project = removeable_projects[1]
                        popfirst!(removeable_projects)
                        removed_capacity = get_device_size(removed_project) * PSY.get_base_power(removed_project)
                        total_removed_capacity += removed_capacity
                        PSY.remove_component!(capacity_market_system, removed_project)
                        ra_metrics = calculate_RA_metrics(capacity_market_system)
                        println(ra_metrics)
                        adequacy_conditions_met, scarcity_conditions_met = check_ra_conditions(ra_targets, ra_metrics)
                        count += 1
                    end
                end
                total_removed_capacity -= removed_capacity
            end

            delta_irm = (total_added_capacity - total_removed_capacity) / forward_peak_load
        end
        end

        set_delta_irm!(resource_adequacy, iteration_year, delta_irm)
    end

    return
end
