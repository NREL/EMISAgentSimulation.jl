"""
function for ifelse elementwise treatment
"""
function elementwise_ifelse(x, y)
    if x <= 0
        z = x
    else
        z = y
    end
    return z
end

"""
This function calculates the derating data for existing and new renewable generation
based on top 100 net-load hour methodology.
"""
function calculate_derating_data(simulation_dir::String,
                                active_projects::Vector{Project},
                                derating_scale::Float64)

    cap_mkt_params = read_data(joinpath(simulation_dir, "markets_data", "Capacity.csv"))

    renewable_existing = filter(p -> typeof(p) == RenewableGenEMIS{Existing}, active_projects)
    renewable_options = filter(p -> typeof(p) == RenewableGenEMIS{Option}, active_projects)

    zones = unique(get_zone.(get_tech.(renewable_existing)))
    types = unique(get_type.(get_tech.(renewable_existing)))

    load_n_vg_data =  read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))
    availability_data = read_data(joinpath(simulation_dir, "timeseries_data_files", "Availability", "DAY_AHEAD_availability.csv"))

    num_hours = DataFrames.nrow(load_n_vg_data)
    num_top_hours = cap_mkt_params.num_top_hours[1]

    existing_vg_power = zeros(num_hours)

    load = vec(sum(Matrix(load_n_vg_data[:, r"load"]), dims=2))

    for g in renewable_existing
        existing_vg_power += load_n_vg_data[!,get_name(g)]
    end

    net_load_df = load_n_vg_data[:, 1:4]
    net_load_df[:, "net_load"] = load - existing_vg_power

    net_load_sorted_df = deepcopy(DataFrames.sort(net_load_df, "net_load", rev = true))

    derating_factors = read_data(joinpath(simulation_dir, "markets_data", "derating_dict.csv"))

    type_zone_max_cap = Dict{String, Float64}()
    for zone in zones
        for type in types
            type_zone_id = "$(type)_$(zone)"
            type_zone_max_cap[type_zone_id] = 0.0
            net_load_df[:, "net_load_w/o_existing_$(type_zone_id)"] = deepcopy(net_load_df[:, "net_load"])
            for g in renewable_existing
                gen_name = get_name(g)
                tech = get_tech(g)
                if "$(get_type(tech))_$(get_zone(tech))" == type_zone_id
                    net_load_df[:, "net_load_w/o_existing_$(type_zone_id)"] += load_n_vg_data[:, gen_name]
                    type_zone_max_cap[type_zone_id] += get_maxcap(g)
                end
            end
        end
    end

    for zone in zones
        for type in types
            type_zone_id = "$(type)_$(zone)"
            gen_sorted_df = deepcopy(DataFrames.sort(net_load_df, "net_load_w/o_existing_$(type_zone_id)", rev = true))

            load_reduction = gen_sorted_df[1:num_top_hours, "net_load_w/o_existing_$(type_zone_id)"] - gen_sorted_df[1:num_top_hours, "net_load"]
            derating_factors[:, "existing_$(type_zone_id)"] .= min(sum(load_reduction) * derating_scale / type_zone_max_cap[type_zone_id] / num_top_hours, 1.0)
        end
    end

    for g in renewable_options
        gen_name = get_name(g)
        tech = get_tech(g)
        type_zone_id = "$(get_type(tech))_$(get_zone(tech))"
        gen_cap = get_maxcap(g)

        net_load_df[:, "net_load_with_$(gen_name)"] =  deepcopy(net_load_df[:, "net_load"] - availability_data[:, "$(type_zone_id)"] * gen_cap)
        gen_sorted_df = deepcopy(DataFrames.sort(net_load_df, "net_load_with_$(gen_name)", rev = true))

        load_reduction = net_load_sorted_df[1:num_top_hours, "net_load"] - gen_sorted_df[1:num_top_hours, "net_load_with_$(gen_name)"]

        derating_factors[:, "new_$(type_zone_id)"] .= min(sum(load_reduction) * derating_scale / gen_cap / num_top_hours, 1.0)
    end

    # Storage CC script
    battery_existing = filter(p -> typeof(p) == BatteryEMIS{Existing}, active_projects)
    peak_reductions = [sum(get_maxcap.(battery_existing))]
    if length(battery_existing) > 0
        eff_charge = Statistics.mean(get_efficiency(get_tech(battery_existing[1])))
        stor_duration = Int(round(Statistics.mean(get_storage_capacity(get_tech(battery_existing[1]))[:max] ./ get_maxcap(battery_existing[1]))))
    else
        eff_charge = 0.85
        stor_duration = 4
    end

    stor_buffer_minutes = cap_mkt_params.stor_buffer_minutes[1]

    inc = 1 #set to 1 for system-wide, but will need to be replaced with zonal level array if/when convert to zonal

    max_demands = repeat((maximum(net_load_df[:,"net_load"]) .- peak_reductions), num_hours)

    batt_powers = repeat(peak_reductions, num_hours)

    poss_charges = min.(batt_powers .* eff_charge, (max_demands - net_load_df[:,"net_load"]) .* eff_charge)

    necessary_discharges = (max_demands - net_load_df[:,"net_load"])

    poss_batt_changes = zeros(size(necessary_discharges)[1])
    for n in collect(1:1:size(necessary_discharges)[1])
        poss_batt_changes[n] = elementwise_ifelse(necessary_discharges[n], poss_charges[n])
    end

    batt_e_level = zeros((inc, num_hours))
    batt_e_level[1] = min(poss_batt_changes[1], 0)
    for n in collect(2:1:num_hours)
        batt_e_level[n] = batt_e_level[n - 1] + poss_batt_changes[n]
        batt_e_level[n] = min(batt_e_level[n], 0.0)
    end

    required_MWhs = -minimum(batt_e_level)

    # This line of code will implement a buffer on all storage duration
    # requirements, i.e. if the stor_buffer_minutes is set to 60 minutes
    # then a 2-hour peak would be served by a 3-hour device, a 3-hour peak
    # by a 4-hour device, etc.
    stor_buffer_hrs = stor_buffer_minutes / 60
    required_MWhs = required_MWhs + (batt_powers[1] * stor_buffer_hrs)[1]

    if length(battery_existing) > 0
        stor_CC = peak_reductions[1] * stor_duration[1] / required_MWhs
        stor_CC = min(stor_CC, 1.0)
    else
        stor_CC = 1.0
    end

    derating_factors[:, "BA_$(stor_duration)"] .= stor_CC
    write_data(joinpath(simulation_dir, "markets_data"), "derating_dict.csv", derating_factors)
    return
end

"""
This function calculates the derating data for existing and new renewable generation
and storage based on PRAS outcomes.
"""

function calculate_derating_factors(
    simulation::Union{AgentSimulation,AgentSimulationData},
    iteration_year::Int64,
    derating_scale::Float64,
    methodology::String,
    ra_matric::String)

    if methodology == "ELCC"
        methodology = PRAS.ELCC
    elseif methodology == "EFC"
        methodology = PRAS.EFC
    else
        @error "Capacity Accreditation methodology should be either ELCC, EFC or TopNetLoad"
    end

    if ra_matric == "LOLE"
        ra_matric = PRAS.LOLE
    elseif ra_matric == "EUE"
        ra_matric = PRAS.EUE
    else
        @error "Resource Adequacy metric should be either LOLE or EUE"
    end


    simulation_dir = get_data_dir(get_case(simulation))
    da_resolution = get_da_resolution(get_case(simulation))
    zones = get_zones(simulation)

    derating_factors = read_data(joinpath(simulation_dir, "markets_data", "derating_dict.csv"))

    load_growth = get_annual_growth(simulation)[:, iteration_year]

    active_projects = get_activeprojects(simulation)
    
    existing = filter(p -> typeof(p) == RenewableGenEMIS{Existing}, active_projects)
    options = filter(p -> typeof(p) == RenewableGenEMIS{Option}, active_projects)

    existing_types = unique(get_type.(get_tech.(existing)))
    new_types = unique(get_type.(get_tech.(options)))

    capacity_forward_years = 3

    resource_adequacy = get_resource_adequacy(simulation)

    sys_UC = get_system_UC(simulation)

    base_sys = deepcopy(sys_UC)

    # create adjusted base system (by iteratively adding or removing generators) such that it meets the RA targets
    adjusted_base_system = create_base_system(base_sys,
        active_projects,
        capacity_forward_years,
        resource_adequacy,
        iteration_year,
        load_growth,
        simulation_dir,
        da_resolution,
        simulation)

    system_period_of_interest = range(1, length = 8760)
    correlated_outage_csv_location = "/kfs2/projects/gmlcmarkets/Phase2_EMIS_Analysis/Correlated Outages/ThermalFOR_2011.csv"

    # create "Base" PRAS system to be used for calculation of ELCC or EFC.
    base_pras_system = make_pras_system(adjusted_base_system,
                                system_model="Single-Node",
                                aggregation="Area",
                                period_of_interest = system_period_of_interest,
                                outage_flag=false,
                                lump_pv_wind_gens=false,
                                availability_flag=true,
                                outage_csv_location = correlated_outage_csv_location)

    for zone in zones
        
        for type in new_types
            idx = findfirst(x -> ((get_type(get_tech(x)) == type) && (get_zone(get_tech(x)) == zone)), options)
            if !isnothing(idx)
                
                augmented_sys = deepcopy(adjusted_base_system)
                build_size = 4 # set to 4 considering that there are 4 investors, so if a project is viable, there could be 4 such units coming online together.
                max_cap = get_maxcap(options[idx]) * build_size
                for i in 1:build_size
                    new_project = deepcopy(options[idx])
                    set_name!(new_project, "$(get_name(new_project))_$i")
                    add_capacity_market_project!(augmented_sys, new_project, simulation_dir, da_resolution)
                end
                
                augmented_pras_system = make_pras_system(augmented_sys,
                                    system_model="Single-Node",
                                    aggregation="Area",
                                    period_of_interest = system_period_of_interest,
                                    outage_flag=false,
                                    lump_pv_wind_gens=false,
                                    availability_flag=true,
                                    outage_csv_location = correlated_outage_csv_location)

                # Call PRAS accreditation methodology. Adjust sample size, seed, etc. here.
                cc_result  =  PRAS.assess(base_pras_system,  augmented_pras_system,  methodology{ra_matric}(Int(ceil(max_cap)), "Region"), PRAS.SequentialMonteCarlo(samples = 10, seed = 42))
                cc_lower,  cc_upper  =  extrema(cc_result) 
                cc_final = (cc_lower + cc_upper) * derating_scale / (2 * max_cap)
                derating_factors[!, "new_$(type)_$(zone)"] .= cc_final
            end
        end
    end
    
    # For average ELCC/EFC, existing units are removed. The new system with reduced units now becomes the base PRAS system.
    augmented_sys = deepcopy(adjusted_base_system)
    augmented_pras_system = make_pras_system(augmented_sys,
                                    system_model="Single-Node",
                                    aggregation="Area",
                                    period_of_interest = system_period_of_interest,
                                    outage_flag=false,
                                    lump_pv_wind_gens=false,
                                    availability_flag=true,
                                    outage_csv_location = correlated_outage_csv_location)

    for zone in zones
        for type in existing_types
            pruned_based_sys = deepcopy(adjusted_base_system)
            total_capacity = 0.0
            zone_tech_units = existing[findall(x -> ((get_type(get_tech(x)) == type) && (get_zone(get_tech(x)) == zone)), existing)]
            if !isempty(zone_tech_units)
                for project in zone_tech_units
                    remove_system_component!(pruned_based_sys, project)
                    total_capacity += get_maxcap(project)
                end

                @assert total_capacity > 0
                pruned_base_pras_system = make_pras_system(pruned_based_sys,
                                        system_model="Single-Node",
                                        aggregation="Area",
                                        period_of_interest = system_period_of_interest,
                                        outage_flag=false,
                                        lump_pv_wind_gens=false,
                                        availability_flag=true,
                                        outage_csv_location = correlated_outage_csv_location)
                #  Call PRAS accreditation methodology. Adjust sample size, seed, etc. here.
                cc_result  =  PRAS.assess(pruned_base_pras_system, augmented_pras_system, PRAS.ELCC{ra_matric}(Int(ceil(total_capacity)), "Region"), PRAS.SequentialMonteCarlo(samples = 10, seed = 42))
                cc_lower,  cc_upper  =  extrema(cc_result) 
                cc_final = (cc_lower + cc_upper) * derating_scale / (2 * total_capacity)

                derating_factors[!, "existing_$(type)_$(zone)"] .= cc_final
            end
        end
    end
       
    pruned_based_sys = deepcopy(adjusted_base_system)
    battery_projects = filter(p -> typeof(p) == BatteryEMIS{Existing}, active_projects)

    # Battery capacuty credits are based on average ELCC/EFC.
    if isempty(battery_projects)
        battery_projects = filter(p -> typeof(p) == BatteryEMIS{Option}, active_projects)
    end

    if !isempty(battery_projects)
        stor_duration = Int(round(Statistics.mean(get_storage_capacity(get_tech(battery_projects[1]))[:max] ./ get_maxcap(battery_projects[1]))))
        total_capacity = 0.0
        for project in battery_projects
            remove_system_component!(pruned_based_sys, project)
            total_capacity += get_maxcap(project)
        end
        
        pruned_base_pras_system = make_pras_system(pruned_based_sys,
            system_model="Single-Node",
            aggregation="Area",
            period_of_interest = system_period_of_interest,
            outage_flag=false,
            lump_pv_wind_gens=false,
            availability_flag=true,
            outage_csv_location = correlated_outage_csv_location)
        
            # Call PRAS accreditation methodology. Adjust sample size, seed, etc. here.
        cc_result  =  PRAS.assess(pruned_base_pras_system, augmented_pras_system, PRAS.ELCC{ra_matric}(Int(ceil(total_capacity)), "Region"), PRAS.SequentialMonteCarlo(samples = 10, seed = 42))
        cc_lower,  cc_upper  =  extrema(cc_result) 
        cc_final = (cc_lower + cc_upper) * derating_scale / (2 * total_capacity)
        derating_factors[!, "BA_$(stor_duration)"] .= cc_final
    end

    # Overwrite file with new derating factors.
    write_data(joinpath(simulation_dir, "markets_data"), "derating_dict.csv", derating_factors)
end


"""
This function does nothing is project is not of ThermalGenEMIS, HydroGenEMIS, RenewableGenEMIS or BatteryEMIS type.
"""
function update_derating_factor!(project::P,
                               simulation_dir::String,
                               derating_scale::Float64
                               ) where P <: Project{<:BuildPhase}
    return
end

"""
This function updates the derating factors of ThermalGenEMIS and HydroGenEMIS projects.
"""
function update_derating_factor!(project::Union{ThermalGenEMIS{<:BuildPhase}, HydroGenEMIS{<:BuildPhase}},
                               simulation_dir::String,
                               derating_scale::Float64
                               )

    derating_data = read_data(joinpath(simulation_dir, "markets_data", "derating_dict.csv"))
    derating_factor = derating_data[1, get_type(get_tech(project))]
    for product in get_products(project)
        set_derating!(product, derating_factor)
    end
    return
end

"""
This function updates the derating factors of existing RenewableGenEMIS projects.
"""
function update_derating_factor!(project::RenewableGenEMIS{Existing},
                               simulation_dir::String,
                               derating_scale::Float64
                               )

    derating_data = read_data(joinpath(simulation_dir, "markets_data", "derating_dict.csv"))
    name = get_name(project)
    tech = get_tech(project)
    type_zone_id = "$(get_type(tech))_$(get_zone(tech))"

    if  in("existing_$(type_zone_id)", names(derating_data))
        derating_factor = derating_data[1, "existing_$(type_zone_id)"]
    else
        error("Derating data not found")
    end

    for product in get_products(project)
        set_derating!(product, derating_factor)
    end

    return
end

"""
This function updates the derating factors of new RenewableGenEMIS projects.
"""
function update_derating_factor!(project::RenewableGenEMIS{<:BuildPhase},
                               simulation_dir::String,
                               derating_scale::Float64
                               )

    derating_data = read_data(joinpath(simulation_dir, "markets_data", "derating_dict.csv"))
    name = get_name(project)
    tech = get_tech(project)
    type_zone_id = "$(get_type(tech))_$(get_zone(tech))"

    if in("new_$(type_zone_id)", names(derating_data))
        derating_factor = derating_data[1, "new_$(type_zone_id)"]
    else
        error("Derating data not found")
    end

    for product in get_products(project)
        set_derating!(product, derating_factor)
    end

    return
end

"""
This function updates the derating factors of BatteryEMIS projects.
"""
function update_derating_factor!(project::BatteryEMIS{<:BuildPhase},
                               simulation_dir::String,
                               derating_scale::Float64
                               )
    tech = get_tech(project)
    duration = Int(round(get_storage_capacity(tech)[:max] / get_maxcap(project)))
    project_type = "$(get_type(tech))_$(duration)"

    derating_data = read_data(joinpath(simulation_dir, "markets_data", "derating_dict.csv"))
    derating_factor = derating_data[1, project_type]
    derating_factor = min(derating_factor * derating_scale, 1.0)
    for product in get_products(project)
        set_derating!(product, derating_factor)
    end
    return
end

"""
This function updates the derating factors of all active projects in the simulation.
"""
function update_simulation_derating_data!(
    simulation::Union{AgentSimulation,AgentSimulationData},
    iteration_year::Int64,
    derating_scale::Float64;
    methodology::String = "ELCC",
    ra_metric::String = "LOLE")

    data_dir = get_data_dir(get_case(simulation))
    active_projects = get_activeprojects(simulation)

    if methodology == "TopNetLoad"
        calculate_derating_data(data_dir, active_projects, derating_scale)
    else
        calculate_derating_factors(simulation, iteration_year, derating_scale, methodology, ra_metric)
    end
    for project in active_projects
        update_derating_factor!(project, data_dir, derating_scale)
    end
    return
end


