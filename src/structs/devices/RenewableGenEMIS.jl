"""
This struct contains the technical data of a renewable generator.
    type: Generator technology type.
    active_power_limits: Named tuple of minimum and maximum generation limits.
    ramp_limits: Ramp limits (MW/hr).
    operation_cost: Operation cost based on PowerSystems.jl definitions.
    bus: Name of the bus on which the generator is located.
    zone: Name of zone in which the generator is located.
"""
struct RenewableTech <: Tech
    type::String
    active_power_limits::NamedTuple{(:min, :max),Tuple{Float64,Float64}}
    ramp_limits::Union{Nothing, NamedTuple{(:up, :down), Tuple{Float64, Float64}}}
    operation_cost::Union{Nothing, PSY.TwoPartCost}
    bus::Int64
    zone::String
    FOR::Float64
end

get_type(tech::RenewableTech) = tech.type
get_active_power_limits(tech::RenewableTech) = tech.active_power_limits
get_ramp_limits(tech::RenewableTech) = tech.ramp_limits
get_operation_cost(tech::RenewableTech) = tech.operation_cost
get_bus(tech::RenewableTech) = tech.bus
get_zone(tech::RenewableTech) = tech.zone
get_FOR(tech::RenewableTech) = tech.FOR

"""
This struct contains all the data for a renewable generator.
    name: Generator's name.
    tech: Struct containing technical data.
    decision_year: The year in which the investment decision is made.
    construction_year: The year in which construction is completed.
    retirement_year: The year in which the generator is retired.
    end_life_year: Age-based retirement year.
    products: Vector containing the data of all products provided by the generator.
    finance_data: Struct containing financial data.
"""
mutable struct RenewableGenEMIS{T <: BuildPhase} <: GeneratorEMIS{T}
    name::String
    tech::RenewableTech
    decision_year::Int64
    construction_year::Int64
    retirement_year::Int64
    end_life_year::Int64
    products::Vector{Product}
    finance_data::Finance
end
