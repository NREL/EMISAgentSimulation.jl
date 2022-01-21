"""
This struct contains the technical data of a hydropower generator.
    type: Generator technology type.
    active_power_limits: Named tuple of minimum and maximum generation limits.
    ramp_limits: Ramp limits (MW/hr).
    time_limits: Minimum up and Minimum down time limits (hr).
    operation_cost: Operation cost based on PowerSystems.jl definitions.
    bus: Name of the bus on which the generator is located.
    zone: Name of zone in which the generator is located.
"""
struct HydroTech <: Tech
    type::String
    active_power_limits::NamedTuple{(:min, :max),Tuple{Float64,Float64}}
    ramp_limits::Union{Nothing, NamedTuple{(:up, :down), Tuple{Float64, Float64}}}
    time_limits::Union{Nothing, NamedTuple{(:up, :down), Tuple{Float64, Float64}}}
    operation_cost::Union{Nothing,PSY.TwoPartCost}
    bus::Int64
    zone::String
    FOR::Float64
    MTTR::Int64
end

get_type(tech::HydroTech) = tech.type
get_active_power_limits(tech::HydroTech) = tech.active_power_limits
get_time_limits(tech::HydroTech) = tech.time_limits
get_ramp_limits(tech::HydroTech) = tech.ramp_limits
get_operation_cost(tech::HydroTech) = tech.operation_cost
get_bus(tech::HydroTech) = tech.bus
get_zone(tech::HydroTech) = tech.zone
get_FOR(tech::HydroTech) = tech.FOR
get_MTTR(tech::HydroTech) = tech.MTTR

"""
This struct contains all the data for a hydropower generator.
    name: Generator's name.
    tech: Struct containing technical data.
    decision_year: The year in which the investment decision is made.
    construction_year: The year in which construction is completed.
    retirement_year: The year in which the generator is retired.
    end_life_year: Age-based retirement year.
    products: Vector containing the data of all products provided by the generator.
    finance_data: Struct containing financial data.
"""
mutable struct HydroGenEMIS{T <: BuildPhase} <: GeneratorEMIS{T}
    name::String
    tech::HydroTech
    decision_year::Int64
    construction_year::Int64
    retirement_year::Int64
    end_life_year::Int64
    products::Vector{Product}
    finance_data::Finance
end
