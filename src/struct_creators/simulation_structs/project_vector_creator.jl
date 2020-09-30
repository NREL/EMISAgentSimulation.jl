"""
This function returns an empty vector of existing projects,
if data for exsting projects does not exist.
"""
function create_project_existing(projectdata::Nothing,
                                simulation_data::AgentSimulationData,
                                sys_UC::Union{Nothing, PSY.System},
                                investor_name::String,
                                investor_dir::String,
                                scenario_names::Vector{String})
    return Project{Existing}[]
end

"""
This function returns a vector of existing projects for the investor.
"""
function create_project_existing(projectdata::DataFrames.DataFrame,
                                simulation_data::AgentSimulationData,
                                sys_UC::Nothing,
                                investor_name::String,
                                investor_dir::String,
                                scenario_names::Vector{String})

    project_existing = Vector{Project{Existing}}(undef, size(projectdata,1))

    for i in 1:size(projectdata,1)
        project_existing[i] = create_project(projectdata[i,:],
                                           simulation_data,
                                           investor_name,
                                           investor_dir,
                                           scenario_names,
                                           false)
        println("Created project $(get_name(project_existing[i])) for $(investor_name)")
    end

    return project_existing
end

"""
This function returns a vector of existing projects for the investor if PSY data exists.
"""
function create_project_existing(projectdata::DataFrames.DataFrame,
                                simulation_data::AgentSimulationData,
                                sys_UC::PSY.System,
                                investor_name::String,
                                investor_dir::String,
                                scenario_names::Vector{String})

    existing_gens = PSY.get_components(PSY.Generator, sys_UC)
    existing_storage = PSY.get_components(PSY.Storage, sys_UC)

    existing_techs = union(existing_gens, existing_storage)

    project_existing = Vector{Project{Existing}}(undef, DataFrames.nrow(projectdata))

    for i in 1:DataFrames.nrow(projectdata)

        device = PSY.get_components_by_name(PSY.Device, sys_UC, projectdata[i, "GEN_UID"])

        project_existing[i] = create_project(projectdata[i, :],
                                           device[1],
                                           simulation_data,
                                           investor_name,
                                           investor_dir,
                                           scenario_names,
                                           false)

        println("Created project $(get_name(project_existing[i])) for $(investor_name)")
    end

    return project_existing
end

"""
This function returns an empty vector of option projects,
if data for option projects does not exist.
"""
function create_project_options(projectdata::Nothing,
                                   simulation_data::AgentSimulationData,
                                   investor_name::String,
                                   investor_dir::String,
                                   scenario_names::Vector{String})
    return Project{Option}[]
end

"""
This function returns a vector of project options for the investor.
"""
function create_project_options(projectdata::DataFrames.DataFrame,
                                   simulation_data::AgentSimulationData,
                                   investor_name::String,
                                   investor_dir::String,
                                   scenario_names::Vector{String})
    project_option = Vector{Project{Option}}(undef, size(projectdata,1))
    for i in 1:size(projectdata,1)
        project_option[i] = create_project(projectdata[i,:],
                                         simulation_data,
                                         investor_name,
                                         investor_dir,
                                         scenario_names,
                                         true)

        println("Created project $(get_name(project_option[i])) for $(investor_name)")
    end
        return project_option
end
