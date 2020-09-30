abstract type DeviceEMIS end

abstract type Project{T <: BuildPhase} <: DeviceEMIS end    #Projects, parameterized by buildphase
abstract type GeneratorEMIS{T <: BuildPhase} <: Project{T} end  #Generators - subtype of Projects
abstract type StorageEMIS{T <: BuildPhase} <: Project{T} end  #Storage - subtype of Projects

get_name(project::P) where P <: Project{<: BuildPhase} = project.name
get_tech(project::P) where P <: Project{<: BuildPhase}  = project.tech
get_decision_year(project::P) where P <: Project{<: BuildPhase}  = project.decision_year
get_construction_year(project::P) where P <: Project{<: BuildPhase}  = project.construction_year
get_retirement_year(project::P) where P <: Project{<: BuildPhase}  = project.retirement_year
get_end_life_year(project::P) where P <: Project{<: BuildPhase}  = project.end_life_year
get_finance_data(project::P) where P <: Project{<: BuildPhase}  = project.finance_data
get_products(project::P) where P <: Project{<: BuildPhase} = project.products
get_mincap(project::P) where P <: GeneratorEMIS{<: BuildPhase}  = get_active_power_limits(get_tech(project))[:min]
get_maxcap(project::P) where P <: GeneratorEMIS{<: BuildPhase}  = get_active_power_limits(get_tech(project))[:max]
get_mincap(project::P) where P <: StorageEMIS{<: BuildPhase}  = get_output_active_power_limits(get_tech(project))[:min]
get_maxcap(project::P) where P <: StorageEMIS{<: BuildPhase}  = get_output_active_power_limits(get_tech(project))[:max]

function set_name!(project::P,
                   name::String) where P <: Project{<: BuildPhase}
        project.name = name
    return
end

function set_decision_year!(project::P, year::Int64) where P <: Project{<: BuildPhase}
    project.decision_year = year
    return
end

function set_construction_year!(project::P, year::Int64) where P <: Project{<: BuildPhase}
    project.construction_year = year
    return
end

function set_retirement_year!(project::P, year) where P <: Project{<: BuildPhase}
    project.retirement_year = year
end

function set_end_life_year!(project::P, year::Int64) where P <: Project{<: BuildPhase}
    project.end_life_year = year
    return
end

function set_investment_cost!(project::P,
                            investment_cost::Vector{Float64}) where P <: Project{<: BuildPhase}
    return
end

function set_investment_cost!(project::P,
                             investment_cost::Vector{Float64}) where P <: Project{Option}
   project.finance_data.investment_cost = investment_cost
end

function set_effective_investment_cost!(project::P,
                                      investment_cost::Float64) where P <: Project{<: BuildPhase}
    project.finance_data.effective_investment_cost = investment_cost
end
