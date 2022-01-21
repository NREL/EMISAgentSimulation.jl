"""
This function does nothing if project is not in queue.
Returns nothing.
"""
function start_construction!(projects::Vector{<: Project{<: BuildPhase}},
                             index::Int64,
                             project::P,
                             iteration_year::Int64) where P <: Project{<: BuildPhase}
    return
end

"""
This function converts a project from queue to planned if
queue time has been completed.
Returns nothing.
"""
function start_construction!(projects::Vector{<: Project{<: BuildPhase}},
                             index::Int64,
                             project::P,
                             iteration_year::Int64) where P <: Project{Queue}

    queue_time = length(get_queue_cost(get_finance_data(project)))

     if get_decision_year(project) + queue_time == iteration_year
        println("CONSTRUCTING:")
        println(get_name(project))
        projects[index] = convert(Project{Planned}, project)
     end
end

"""
This function does nothing if project is not a planned project.
Returns nothing.
"""
function finish_construction!(projects::Vector{<: Project{<: BuildPhase}},
                             index::Int64,
                             project::P,
                             sys_UC::Union{Nothing, PSY.System},
                             sys_ED::Union{Nothing, PSY.System},
                             simulation_dir::String,
                             iteration_year::Int64,
                             da_resolution::Int64,
                             rt_resolution::Int64) where P <: Project{<: BuildPhase}
    return
end

"""
This function converts a project from planned to existing if
construction time has been completed.
Returns nothing.
"""
function finish_construction!(projects::Vector{<: Project{<: BuildPhase}},
                             index::Int64,
                             project::P,
                             sys_UC::Nothing,
                             sys_ED::Nothing,
                             simulation_dir::String,
                             iteration_year::Int64,
                             da_resolution::Int64,
                             rt_resolution::Int64) where P <: Project{Planned}

     if get_construction_year(project) == iteration_year
        projects[index] = convert(Project{Existing}, project)

        if typeof(project) == RenewableGenEMIS{Planned}
            type = get_type(get_tech(project))
            zone = get_zone(get_tech(project))

            availability_df = read_data(joinpath(simulation_dir, "timeseries_data_files", "Availability", "DAY_AHEAD_availability.csv"))
            availability_df_rt = read_data(joinpath(simulation_dir, "timeseries_data_files", "Availability", "REAL_TIME_availability.csv"))

            if in(get_name(project), names(availability_df))
                availability_raw = availability_df[:, Symbol(get_name(project))]
                availability_raw_rt = availability_df_rt[:, Symbol(get_name(project))]
            elseif in("$(type)_$(zone)", names(availability_df))
                availability_raw = availability_df[:, Symbol("$(type)_$(zone)")]
                availability_raw_rt = availability_df[:, Symbol("$(type)_$(zone)")]
            end

            load_n_vg_df =  read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))
            load_n_vg_df[:, get_name(project)] = availability_raw * get_maxcap(project)

            load_n_vg_df_rt =  read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data_rt.csv"))
            load_n_vg_df_rt[:, get_name(project)] = availability_raw_rt * get_maxcap(project)

            write_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data"), "load_n_vg_data.csv", load_n_vg_df)
            write_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data"), "load_n_vg_data_rt.csv", load_n_vg_df_rt)
        end
     end

    return

end

"""
This function converts a project from planned to existing if
construction time has been completed.
Returns nothing.
"""
function finish_construction!(projects::Vector{<: Project{<: BuildPhase}},
                             index::Int64,
                             project::P,
                             sys_UC:: PSY.System,
                             sys_ED:: PSY.System,
                             simulation_dir::String,
                             iteration_year::Int64,
                             da_resolution::Int64,
                             rt_resolution::Int64) where P <: Project{Planned}
     if get_construction_year(project) == iteration_year
        projects[index] = convert(Project{Existing}, project)
        PSY_project_UC = create_PSY_generator(project, sys_UC)
        PSY_project_ED = create_PSY_generator(project, sys_ED)

        PSY.add_component!(sys_UC, PSY_project_UC)
        PSY.add_component!(sys_ED, PSY_project_ED)

        for product in get_products(project)
            add_device_services!(sys_UC, PSY_project_UC, product)
            add_device_services!(sys_ED, PSY_project_ED, product)
        end

        type = get_type(get_tech(project))
        zone = get_zone(get_tech(project))

        availability_df = read_data(joinpath(simulation_dir, "timeseries_data_files", "Availability", "DAY_AHEAD_availability.csv"))
        availability_df_rt = read_data(joinpath(simulation_dir, "timeseries_data_files", "Availability", "REAL_TIME_availability.csv"))

        if in(get_name(project), names(availability_df))
            availability_raw = availability_df[:, Symbol(get_name(project))]
            availability_raw_rt = availability_df_rt[:, Symbol(get_name(project))]
        elseif in("$(type)_$(zone)", names(availability_df))
            availability_raw = availability_df[:, Symbol("$(type)_$(zone)")]
            availability_raw_rt = availability_df_rt[:, Symbol("$(type)_$(zone)")]
        end

        add_device_forecast!(simulation_dir, sys_UC, sys_ED, PSY_project_UC, PSY_project_ED, availability_raw, availability_raw_rt, da_resolution, rt_resolution)

        if type == "NU_ST" || type == "RE_CT"
            convert_to_thermal_clean_energy!(PSY_project_UC, sys_UC)
            convert_to_thermal_clean_energy!(PSY_project_ED, sys_ED)
        elseif type == "CT"
            convert_to_thermal_fast_start!(PSY_project_UC, sys_UC)
            convert_to_thermal_fast_start!(PSY_project_ED, sys_ED)
        end


        add_clean_energy_contribution!(sys_UC, PSY_project_UC)
     end
     println("FINISHED CONSTRUCTING: $(get_name(project))")

     return
end

"""
This function does nothing if project is not an existing project.
Returns nothing.
"""
function retire_old!(projects::Vector{<: Project{<: BuildPhase}},
                     index::Int64,
                     project::P,
                     sys_UC::Union{Nothing, PSY. System},
                     sys_ED::Union{Nothing, PSY. System},
                     simulation_dir::String,
                     iteration_year::Int64) where P <: Project{<: BuildPhase}
    return false
end

"""
This function converts a project from existing to retired if
lifetime has ended.
Returns nothing.
"""
function retire_old!(projects::Vector{<: Project{<: BuildPhase}},
                     index::Int64,
                     project::P,
                     sys_UC::Union{Nothing, PSY. System},
                     sys_ED::Union{Nothing, PSY. System},
                     simulation_dir::String,
                     iteration_year::Int64) where P <: Project{Existing}

    if get_end_life_year(project) == iteration_year
        set_retirement_year!(projects[index], iteration_year + 1)
        projects[index] = convert(Project{Retired}, project)
        remove_system_component!(sys_UC, project)
        remove_system_component!(sys_ED, project)
        remove_renewable_gen_data!(project, simulation_dir)
    end
end

