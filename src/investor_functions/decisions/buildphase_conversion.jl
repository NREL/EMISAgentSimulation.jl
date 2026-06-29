"""
This function does nothing if project is not in queue.
Returns nothing.
"""
function start_construction!(projects::Vector{<: Project{<: BuildPhase}},
                             index::Int64,
                             project::P,
                             iteration_year::Int64,
                             step_size::Int64) where P <: Project{<: BuildPhase}
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
                             iteration_year::Int64,
                             step_size::Int64) where P <: Project{Queue}

    queue_time = length(get_queue_cost(get_finance_data(project)))

    # check if project construction start year is within this iteration step
    if iteration_year <= get_decision_year(project) + queue_time <= iteration_year + step_size - 1
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
                             sys_MDs::Union{Nothing, Vector{PSY.System}},
                             sys_UCs::Union{Nothing, Vector{PSY.System}},
                             sys_EDs::Union{Nothing, Vector{PSY.System}},
                             sys_PRAS::Dict{String, PSY.System},
                             simulation_dir::String,
                             iteration_year::Int64,
                             step_size::Int64,
                             pcm_scenario::String,
                             simulation_years::Int64,
                             scenario_names::Vector{String},
                             da_resolution::Int64,
                             rt_resolution::Int64,
                             timeseries_data_dir::String) where P <: Project{<: BuildPhase}
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
                             sys_MDs::Nothing,
                             sys_UCs::Nothing,
                             sys_EDs::Nothing,
                             sys_PRAS::Dict{String, PSY.System},
                             simulation_dir::String,
                             iteration_year::Int64,
                             step_size::Int64,
                             pcm_scenario::String,
                             simulation_years::Int64,
                             scenario_names::Vector{String},
                             da_resolution::Int64,
                             rt_resolution::Int64,
                             timeseries_data_dir::String) where P <: Project{Planned}

    # check if project construction end year is within this iteration step
     if iteration_year <= get_construction_year(project) <= iteration_year + step_size - 1
        projects[index] = convert(Project{Existing}, project)

        if typeof(project) == RenewableGenEMIS{Planned}
            type = get_type(get_tech(project))
            zone = get_zone(get_tech(project))

            for scenario in scenario_names
                
                availability_df = read_data(joinpath(timeseries_data_dir, scenario, "sim_year_$(iteration_year)", "Availability", "DAY_AHEAD_availability.csv"))
                availability_df_rt = read_data(joinpath(timeseries_data_dir, scenario, "sim_year_$(iteration_year)", "Availability", "REAL_TIME_availability.csv"))

                if in(get_name(project), names(availability_df))
                    availability_raw = availability_df[:, Symbol(get_name(project))]
                    availability_raw_rt = availability_df_rt[:, Symbol(get_name(project))]
                elseif in("$(type)_$(zone)", names(availability_df))
                    availability_raw = availability_df[:, Symbol("$(type)_$(zone)")]
                    availability_raw_rt = availability_df[:, Symbol("$(type)_$(zone)")]
                end

                for year in 1:simulation_years
                    load_n_vg_df =  read_data(joinpath(timeseries_data_dir, scenario, "sim_year_$(year)", "Net Load Data", "load_n_vg_data.csv"))
                    load_n_vg_df[:, get_name(project)] = availability_raw * get_maxcap(project)

                    load_n_vg_df_rt =  read_data(joinpath(timeseries_data_dir, scenario, "sim_year_$(year)", "Net Load Data", "load_n_vg_data_rt.csv"))
                    load_n_vg_df_rt[:, get_name(project)] = availability_raw_rt * get_maxcap(project)

                    write_data(joinpath(timeseries_data_dir, scenario, "sim_year_$(year)", "Net Load Data"), "load_n_vg_data.csv", load_n_vg_df)
                    write_data(joinpath(timeseries_data_dir, scenario, "sim_year_$(year)", "Net Load Data"), "load_n_vg_data_rt.csv", load_n_vg_df_rt)
                end
            end
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
                             sys_MDs::Vector{PSY.System},
                             sys_UCs::Vector{PSY.System},
                             sys_EDs::Vector{PSY.System},
                             sys_PRAS::Dict{String, PSY.System},
                             simulation_dir::String,
                             iteration_year::Int64,
                             step_size::Int64,
                             pcm_scenario::String,
                             simulation_years::Int64,
                             scenario_names::Vector{String},
                             da_resolution::Int64,
                             rt_resolution::Int64,
                             timeseries_data_dir::String) where P <: Project{Planned}
                             
    # check if project construction end year is within this iteration step
     if iteration_year <= get_construction_year(project) <= iteration_year + step_size - 1
        projects[index] = convert(Project{Existing}, project)
        @info "Finishing construction of project $(get_name(project)) in year $(iteration_year)."
        PSY_project_MD = create_PSY_generator(project, sys_MDs[iteration_year])
        PSY_project_UC = create_PSY_generator(project, sys_UCs[iteration_year])
        PSY_project_ED = create_PSY_generator(project, sys_EDs[iteration_year])

        for y in iteration_year:simulation_years
            PSY_project_MD_iteration_year = create_PSY_generator(project, sys_MDs[y])
            PSY_project_UC_iteration_year = create_PSY_generator(project, sys_UCs[y])
            PSY_project_ED_iteration_year = create_PSY_generator(project, sys_EDs[y])

            PSY.add_component!(sys_MDs[y], PSY_project_MD_iteration_year)
            PSY.add_component!(sys_UCs[y], PSY_project_UC_iteration_year)
            PSY.add_component!(sys_EDs[y], PSY_project_ED_iteration_year)
            add_nominal_outage_to_component!(sys_MDs[y], PSY_project_MD_iteration_year)
            add_nominal_outage_to_component!(sys_UCs[y], PSY_project_UC_iteration_year)
            add_nominal_outage_to_component!(sys_EDs[y], PSY_project_ED_iteration_year)
        end

        for y in iteration_year:simulation_years
            for product in get_products(project)
                add_device_services!(sys_MDs[y], PSY.get_component(typeof(PSY_project_MD), sys_MDs[y], PSY.get_name(PSY_project_MD)), product)
                add_device_services!(sys_UCs[y], PSY.get_component(typeof(PSY_project_UC), sys_UCs[y], PSY.get_name(PSY_project_UC)), product)
                add_device_services!(sys_EDs[y], PSY.get_component(typeof(PSY_project_ED), sys_EDs[y], PSY.get_name(PSY_project_ED)), product)
            end
        end

        for scenario in keys(sys_PRAS)
            PSY_project_PRAS = create_PSY_generator(project, sys_PRAS[scenario])
            bus_name = PSY.get_name(PSY.get_bus(PSY_project_PRAS))
            @info "Adding project $(get_name(project)) to PRAS system for scenario $(scenario). - project bus name is $(bus_name)"
            PSY.add_component!(sys_PRAS[scenario], PSY_project_PRAS)
            add_nominal_outage_to_component!(sys_PRAS[scenario], PSY_project_PRAS)
        end

            for product in get_products(project)
                @info "Adding product $(product) to project $(get_name(project)) in PRAS system for scenario $(scenario)."
                add_device_services!(sys_PRAS[scenario], PSY.get_component(typeof(PSY_project_PRAS), sys_PRAS[scenario], PSY.get_name(PSY_project_PRAS)), product)
            end            

        end

        project_type = get_type(get_tech(project))
        project_zone = get_zone(get_tech(project))

        for y in 1:simulation_years            
            availability_df = read_data(joinpath(timeseries_data_dir, pcm_scenario, "sim_year_$(y)", "Availability", "DAY_AHEAD_availability.csv"))
            availability_df_rt = read_data(joinpath(timeseries_data_dir, pcm_scenario, "sim_year_$(y)", "Availability", "REAL_TIME_availability.csv"))

            project_name = get_name(project)
            availability_df_names = names(availability_df)
            project_type_zone = "$(project_type)_$(project_zone)"
            @info "Finding availability data for project $(project_name) in year $(y) for project type $(project_type) and zone $(project_zone)."
            if in(project_name, availability_df_names)
                @info "Found availability data for project $(project_name) in year $(y)."
                availability_raw = availability_df[:, Symbol(project_name)]
                availability_raw_rt = availability_df_rt[:, Symbol(project_name)]
            elseif in(project_type_zone, availability_df_names)
                @info "Found availability data for project type $(project_type) and zone $(project_zone) in year $(y)."
                availability_raw = availability_df[:, Symbol(project_type_zone)]
                availability_raw_rt = availability_df_rt[:, Symbol(project_type_zone)]
            else
                @error "Could not find availability data for project $(project_name) or project type and zone combination $(project_type_zone) in year $(y)."
            end

            if y >= iteration_year
                add_device_forecast!(sys_MDs[y], sys_UCs[y], sys_EDs[y],
                                     PSY.get_component(typeof(PSY_project_MD), sys_MDs[y],
                                     PSY.get_name(PSY_project_MD)), PSY.get_component(typeof(PSY_project_UC),
                                     sys_UCs[y], PSY.get_name(PSY_project_UC)), PSY.get_component(typeof(PSY_project_ED),
                                     sys_EDs[y], PSY.get_name(PSY_project_ED)), availability_raw, availability_raw_rt,
                                     da_resolution, rt_resolution)
                
                if project_type == "NU_ST" || project_type == "RE_CT" || typeof(PSY_project_UC) == PSY.RenewableDispatch || typeof(PSY_project_UC) == PSY.HydroTurbine || typeof(PSY_project_UC) == PSY.HydroDispatch
                    # convert_to_thermal_clean_energy!(PSY_project_UC, sys_UCs[y])
                    # convert_to_thermal_clean_energy!(PSY_project_ED, sys_EDs[y])
                    add_clean_energy_contribution!(sys_MDs[y], PSY.get_component(typeof(PSY_project_MD), sys_MDs[y], PSY.get_name(PSY_project_MD)))
                    add_clean_energy_contribution!(sys_UCs[y], PSY.get_component(typeof(PSY_project_UC), sys_UCs[y], PSY.get_name(PSY_project_UC)))
                    add_clean_energy_contribution!(sys_EDs[y], PSY.get_component(typeof(PSY_project_ED), sys_EDs[y], PSY.get_name(PSY_project_ED)))
                    if project_type == "RE_CT"
                        convert_to_thermal_fast_start!(PSY.get_component(typeof(PSY_project_MD), sys_MDs[y], PSY.get_name(PSY_project_MD)), sys_MDs[y])
                        convert_to_thermal_fast_start!(PSY.get_component(typeof(PSY_project_UC), sys_UCs[y], PSY.get_name(PSY_project_UC)), sys_UCs[y])
                        convert_to_thermal_fast_start!(PSY.get_component(typeof(PSY_project_ED), sys_EDs[y], PSY.get_name(PSY_project_ED)), sys_EDs[y])
                    end
                elseif project_type == "CT"
                    convert_to_thermal_fast_start!(PSY.get_component(typeof(PSY_project_MD), sys_MDs[y], PSY.get_name(PSY_project_MD)), sys_MDs[y])
                    convert_to_thermal_fast_start!(PSY.get_component(typeof(PSY_project_UC), sys_UCs[y], PSY.get_name(PSY_project_UC)), sys_UCs[y])
                    convert_to_thermal_fast_start!(PSY.get_component(typeof(PSY_project_ED), sys_EDs[y], PSY.get_name(PSY_project_ED)), sys_EDs[y])
                end   
            end     

            write_vg_data(sys_UCs[y], sys_EDs[y], PSY_project_UC, PSY_project_ED, simulation_dir, availability_raw, availability_raw_rt, scenario_names, y, timeseries_data_dir)

        end

        for scenario in keys(sys_PRAS)
            PSY_project_PRAS = create_PSY_generator(project, sys_PRAS[scenario])
            availability_df_rt = DataFrames.DataFrame()

            for sim_year in 1:simulation_years
                availability_df_rt = vcat(availability_df_rt, read_data(joinpath(timeseries_data_dir, scenario, "sim_year_$(sim_year)", "Availability", "REAL_TIME_availability.csv")))
            end

            if in(get_name(project), names(availability_df_rt))
                availability_raw_rt = availability_df_rt[:, Symbol(get_name(project))]
            elseif in("$(project_type)_$(project_zone)", names(availability_df_rt))
                availability_raw_rt = availability_df_rt[:, Symbol("$(project_type)_$(project_zone)")]
            end

            add_device_forecast_PRAS!(sys_PRAS[scenario], 
                                    PSY.get_component(typeof(PSY_project_PRAS), sys_PRAS[scenario], PSY.get_name(PSY_project_PRAS)),
                                    availability_raw_rt,
                                    rt_resolution,
                                    simulation_years)

            if project_type == "NU_ST" || project_type == "RE_CT" || typeof(PSY_project_PRAS) == PSY.RenewableDispatch || typeof(PSY_project_PRAS) == PSY.HydroTurbine || typeof(PSY_project_PRAS) == PSY.HydroDispatch
                if project_type == "RE_CT"
                    convert_to_thermal_fast_start!(PSY.get_component(typeof(PSY_project_PRAS), sys_PRAS[scenario], PSY.get_name(PSY_project_PRAS)), sys_PRAS[scenario])
                end
            elseif project_type == "CT"
                convert_to_thermal_fast_start!(PSY.get_component(typeof(PSY_project_PRAS), sys_PRAS[scenario], PSY.get_name(PSY_project_PRAS)), sys_PRAS[scenario])
            end
        end

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
                     sys_MDs::Union{Nothing, Vector{PSY.System}},
                     sys_UCs::Union{Nothing, Vector{PSY.System}},
                     sys_EDs::Union{Nothing, Vector{PSY.System}},
                     sys_PRAS::Dict{String, PSY.System},
                     simulation_dir::String,
                     iteration_year::Int64,
                     step_size::Int64,
                     scenario_names::Vector{String},
                     total_horizon::Int64,
                     timeseries_data_dir::String) where P <: Project{<: BuildPhase}
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
                     sys_MDs::Union{Nothing, Vector{PSY.System}},
                     sys_UCs::Union{Nothing, Vector{PSY.System}},
                     sys_EDs::Union{Nothing, Vector{PSY.System}},
                     sys_PRAS::Dict{String, PSY.System},
                     simulation_dir::String,
                     iteration_year::Int64,
                     step_size::Int64,
                     scenario_names::Vector{String},
                     total_horizon::Int64,
                     timeseries_data_dir::String
                     ) where P <: Project{Existing}

    # check if project end of life is within this iteration step
    if iteration_year <= get_end_life_year(project) <= iteration_year + step_size - 1
        set_retirement_year!(projects[index], iteration_year + step_size)
        projects[index] = convert(Project{Retired}, project)
        for y in iteration_year:total_horizon
            remove_system_component!(sys_MDs[y], project)
            remove_system_component!(sys_UCs[y], project)
            remove_system_component!(sys_EDs[y], project)
        end

        for scenario in keys(sys_PRAS)
            remove_system_component!(sys_PRAS[scenario], project)
        end

        remove_renewable_gen_data!(project, simulation_dir, iteration_year, total_horizon, scenario_names, timeseries_data_dir)
    end
end

