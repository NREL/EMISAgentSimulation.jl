"""
This struct contains the technical data of a thermal generator.
    type: Generator technology type.
    fuel: Fuel type.
    active_power_limits: Named tuple of minimum and maximum generation limits.
    ramp_limits: Ramp limits (MW/hr).
    time_limits: Minimum up and Minimum down time limits (hr).
    operation_cost: Operation cost based on PowerSystems.jl definitions.
    fuel_cost: Fuel cost
    heat_rate_curve: Heat Rate curve
    bus: Name of the bus on which the generator is located.
    zone: Name of zone in which the generator is located.
"""
struct ThermalTech <: Tech
    type::String
    fuel::String
    active_power_limits::NamedTuple{(:min, :max),Tuple{Float64,Float64}}
    ramp_limits::Union{Nothing, NamedTuple{(:up, :down), Tuple{Float64, Float64}}}
    time_limits::Union{Nothing, NamedTuple{(:up, :down), Tuple{Float64, Float64}}}
    operation_cost::Union{Nothing, PSY.ThreePartCost}
    fuel_cost::Float64
    heat_rate_curve::Vector{Tuple{Float64, Float64}}
    bus::Int64
    zone::String
    FOR::Float64
end

get_type(tech::ThermalTech) = tech.type
get_fuel(tech::ThermalTech) = tech.fuel
get_active_power_limits(tech::ThermalTech) = tech.active_power_limits
get_ramp_limits(tech::ThermalTech) = tech.ramp_limits
get_time_limits(tech::ThermalTech) = tech.time_limits
get_operation_cost(tech::ThermalTech) = tech.operation_cost
get_fuel_cost(tech::ThermalTech) = tech.fuel_cost
get_bus(tech::ThermalTech) = tech.bus
get_zone(tech::ThermalTech) = tech.zone
get_FOR(tech::ThermalTech) = tech.FOR
get_heat_rate_curve(tech::ThermalTech) = tech.heat_rate_curve

"""
This struct contains all the data for a thermal generator.
    name: Generator's name.
    tech: Struct containing technical data.
    decision_year: The year in which the investment decision is made.
    construction_year: The year in which construction is completed.
    retirement_year: The year in which the generator is retired.
    end_life_year: Age-based retirement year.
    products: Vector containing the data of all products provided by the generator.
    finance_data: Struct containing financial data.
"""
mutable struct ThermalGenEMIS{T <: BuildPhase} <: GeneratorEMIS{T}
    name::String
    tech::ThermalTech
    decision_year::Int64
    construction_year::Int64
    retirement_year::Int64
    end_life_year::Int64
    products::Vector{Product}
    finance_data::Finance
end
