"""
This struct contains all the data for the simulation to be run.
    case: CaseDefinition for creating and running the simulation.
    results_dir: Path to the directory where results are stored.
    iteration_year: Index of current iteration year.
    system_UC: SIIP PSY System for Unit Commitment.
    system_ED: SIIP PSY System for Economic Dispatch.
    zones: Names of modeled zones.
    lines: Inter-zonal ransmission lines modeled.
    hour_weight: Weight allocated to each hour of actual market clearing.
    peak_load: Annual system-wide peak load for capacity market clearing.
    markets: Boolean for selecting which markets are modeled.
    investors: Vectors of all investors created in the simulation.
    derating_data: Derating data foe each technology type.
    annual_growth: Annual growth rate data of different parameters.
"""

mutable struct AgentSimulation
    case::CaseDefinition
    results_dir::String
    iteration_year::Int64
    system_UC::Union{Nothing, PSY.System}
    system_ED::Union{Nothing, PSY.System}
    zones::Vector{String}
    lines::Vector{ZonalLine}
    hour_weight::Vector{Float64}
    peak_load::Float64
    markets::Dict{Symbol, Bool}
    investors::Vector{Investor}
    derating_data::DataFrames.DataFrame
    annual_growth::AxisArrays.AxisArray{Float64, 2}
end

get_case(sim::AgentSimulation) = sim.case
get_results_dir(sim::AgentSimulation) = sim.results_dir
get_iteration_year(sim::AgentSimulation) = sim.iteration_year
get_system_UC(sim::AgentSimulation) = sim.system_UC
get_system_ED(sim::AgentSimulation) = sim.system_ED
get_zones(sim::AgentSimulation) = sim.zones
get_lines(sim::AgentSimulation) = sim.lines
get_hour_weight(sim::AgentSimulation) = sim.hour_weight
get_peak_load(sim::AgentSimulation) = sim.peak_load
get_markets(sim::AgentSimulation) = sim.markets
get_investors(sim::AgentSimulation) = sim.investors
get_derating_data(sim::AgentSimulation) = sim.derating_data
get_annual_growth(sim::AgentSimulation) = sim.annual_growth

function set_iteration_year!(simulation::AgentSimulation, iteration_year::Int64)
    simulation.iteration_year = iteration_year
end

function set_peak_load!(simulation::AgentSimulation, peak_load::Float64)
    simulation.peak_load = peak_load
end


"""
This struct temporarily holds the data for creating a simulation.
"""
mutable struct AgentSimulationData
    case::CaseDefinition
    iteration_year::Int64
    system_UC::Union{Nothing, PSY.System}
    system_ED::Union{Nothing, PSY.System}
    zones::Vector{String}
    lines::Vector{ZonalLine}
    hour_weight::Vector{Float64}
    rep_hour_weight::Vector{Float64}
    peak_load::Float64
    markets::Dict{Symbol, Bool}
    investors::Union{Nothing, Vector{Investor}}
    queue_cost_data::DataFrames.DataFrame
    derating_data::DataFrames.DataFrame
    annual_growth::AxisArrays.AxisArray{Float64, 2}
end

function AgentSimulationData(case::CaseDefinition,
                        system_UC::Union{Nothing, PSY.System},
                        system_ED::Union{Nothing, PSY.System},
                        zones::Vector{String},
                        lines::Vector{ZonalLine},
                        hour_weight::Vector{Float64},
                        rep_hour_weight::Vector{Float64},
                        peak_load::Float64,
                        markets::Dict{Symbol, Bool},
                        queue_cost_data::DataFrames.DataFrame,
                        derating_data::DataFrames.DataFrame,
                        annual_growth::AxisArrays.AxisArray{Float64, 2})
    return AgentSimulationData(case,
                          1,
                          system_UC,
                          system_ED,
                          zones,
                          lines,
                          hour_weight,
                          rep_hour_weight,
                          peak_load,
                          markets,
                          nothing,
                          queue_cost_data,
                          derating_data,
                          annual_growth)
end

get_case(sim::AgentSimulationData) = sim.case
get_iteration_year(sim::AgentSimulationData) = sim.iteration_year
get_system_UC(sim::AgentSimulationData) = sim.system_UC
get_system_ED(sim::AgentSimulationData) = sim.system_ED
get_zones(sim::AgentSimulationData) = sim.zones
get_lines(sim::AgentSimulationData) = sim.lines
get_hour_weight(sim::AgentSimulationData) = sim.hour_weight
get_rep_hour_weight(sim::AgentSimulationData) = sim.rep_hour_weight
get_peak_load(sim::AgentSimulationData) = sim.peak_load
get_markets(sim::AgentSimulationData) = sim.markets
get_investors(sim::AgentSimulationData) = sim.investors
get_queue_cost_data(sim::AgentSimulationData) = sim.queue_cost_data
get_derating_data(sim::AgentSimulationData) = sim.derating_data
get_annual_growth(sim::AgentSimulationData) = sim.annual_growth

function set_investors!(simulation_data::AgentSimulationData,
                        investors::Vector{Investor})
    simulation_data.investors = investors
    return
end

function get_allprojects(simulation::Union{AgentSimulation, AgentSimulationData})
    projects = Project[]
    append!(projects, vcat(get_projects.(get_investors(simulation))...))
    return projects
end

# Get all investors' projects except those which have been retired.
function get_activeprojects(simulation::Union{AgentSimulation, AgentSimulationData})
    projects = get_allprojects(simulation)
    active_projects = filter(x -> !isa(x, Project{Retired}) , projects)
    return active_projects
end

# Get all investors' projects except options.
function get_decidedprojects(simulation::Union{AgentSimulation, AgentSimulationData})
    projects = get_allprojects(simulation)
    decided_projects = filter(x -> !isa(x, Project{Option}) , projects)
    return decided_projects
end
