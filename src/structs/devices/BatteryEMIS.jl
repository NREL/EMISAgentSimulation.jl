"""
This struct contains the technical data of a battery.
    type: Battery technology type.
    input_active_power_limits: Named tuple of minimum and maximum input power limits.
    output_active_power_limits: Named tuple of minimum and maximum output power limits.
    ramp_limits: Ramp limits (MW/hr).
    storage_capacity: Named tuple of minimum and maximum storage capacity.
    soc: State of charge - Important for first time-step of the iteration.
    efficiency: Charging and discharging efficiencies.
    bus: Name of bus at which the battery is located.
    zone: Name of zone in which the battery is located.
"""

struct BatteryTech
    type::String
    input_active_power_limits::NamedTuple{(:min, :max),Tuple{Float64,Float64}}
    output_active_power_limits::NamedTuple{(:min, :max),Tuple{Float64,Float64}}
    ramp_limits::Union{Nothing, NamedTuple{(:up, :down), Tuple{Float64, Float64}}}
    storage_capacity::NamedTuple{(:min, :max),Tuple{Float64,Float64}}
    soc::Float64
    efficiency::NamedTuple{(:in, :out), Tuple{Float64, Float64}}
    bus::Int64
    zone::String
end

get_type(tech::BatteryTech) = tech.type
get_input_active_power_limits(tech::BatteryTech) = tech.input_active_power_limits
get_output_active_power_limits(tech::BatteryTech) = tech.output_active_power_limits
get_ramp_limits(tech::BatteryTech) = tech.ramp_limits
get_storage_capacity(tech::BatteryTech) = tech.storage_capacity
get_soc(tech::BatteryTech) = tech.soc
get_efficiency(tech::BatteryTech) = tech.efficiency
get_bus(tech::BatteryTech) = tech.bus
get_zone(tech::BatteryTech) = tech.zone

"""
This struct contains all the data for battery storage.
    name: Battery's name.
    tech: Struct containing technical data.
    decision_year: The year in which the investment decision is made.
    construction_year: The year in which construction is completed.
    retirement_year: The year in which the generator is retired.
    end_life_year: Age-based retirement year.
    products: Vector containing the data of all products provided by the generator.
    finance_data: Struct containing financial data.
"""
mutable struct BatteryEMIS{T <: BuildPhase} <: StorageEMIS{T}
    name::String
    tech::BatteryTech
    decision_year::Int64
    construction_year::Int64
    retirement_year::Int64
    end_life_year::Int64
    products::Vector{Product}
    finance_data::Finance
end
