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
    write_data(joinpath(simulation_dir, "markets_data"), "derating_dict.csv", derating_factors)
    return
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
    duration = Int(get_storage_capacity(tech)[:max] / get_maxcap(project))
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
function update_simulation_derating_data!(simulation::Union{AgentSimulation, AgentSimulationData}, derating_scale::Float64)
    data_dir = get_data_dir(get_case(simulation))
    active_projects = get_activeprojects(simulation)

    calculate_derating_data(data_dir, active_projects, derating_scale)
    for project in active_projects
        update_derating_factor!(project, data_dir, derating_scale)
    end
    return
end
