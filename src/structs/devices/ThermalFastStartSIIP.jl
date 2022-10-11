"""
    mutable struct ThermalFastStartSIIP <: ThermalGen
        name::String
        available::Bool
        status::Bool
        bus::Bus
        active_power::Float64
        reactive_power::Float64
        rating::Float64
        active_power_limits::NamedTuple{(:min, :max), Tuple{Float64, Float64}}
        reactive_power_limits::Union{Nothing, Min_Max}
        ramp_limits::Union{Nothing, NamedTuple{(:up, :down), Tuple{Float64, Float64}}}
        operation_cost::OperationalCost
        base_power::Float64
        time_limits::Union{Nothing, NamedTuple{(:up, :down), Tuple{Float64, Float64}}}
        prime_mover::PrimeMovers
        fuel::ThermalFuels
        services::Vector{Service}
        time_at_status::Float64
        dynamic_injector::Union{Nothing, DynamicInjection}
        ext::Dict{String, Any}
        time_series_container::InfrastructureSystems.TimeSeriesContainer
        internal::InfrastructureSystemsInternal
    end
Data Structure for thermal generation technologies.
# Arguments
- `name::String`
- `available::Bool`
- `status::Bool`
- `bus::Bus`
- `active_power::Float64`, validation range: `active_power_limits`, action if invalid: `warn`
- `reactive_power::Float64`, validation range: `reactive_power_limits`, action if invalid: `warn`
- `rating::Float64`: Thermal limited MVA Power Output of the unit. <= Capacity, validation range: `(0, nothing)`, action if invalid: `error`
- `active_power_limits::NamedTuple{(:min, :max), Tuple{Float64, Float64}}`, validation range: `(0, nothing)`, action if invalid: `warn`
- `reactive_power_limits::Union{Nothing, Min_Max}`
- `ramp_limits::Union{Nothing, NamedTuple{(:up, :down), Tuple{Float64, Float64}}}`: ramp up and ramp down limits in MW (in component base per unit) per minute, validation range: `(0, nothing)`, action if invalid: `error`
- `operation_cost::OperationalCost`
- `base_power::Float64`: Base power of the unit in MVA, validation range: `(0, nothing)`, action if invalid: `warn`
- `time_limits::Union{Nothing, NamedTuple{(:up, :down), Tuple{Float64, Float64}}}`: Minimum up and Minimum down time limits in hours, validation range: `(0, nothing)`, action if invalid: `error`
- `prime_mover::PrimeMovers`: Prime mover technology according to EIA 923
- `fuel::ThermalFuels`: Prime mover fuel according to EIA 923
- `services::Vector{Service}`: Services that this device contributes to
- `time_at_status::Float64`
- `dynamic_injector::Union{Nothing, DynamicInjection}`: corresponding dynamic injection device
- `ext::Dict{String, Any}`
- `time_series_container::InfrastructureSystems.TimeSeriesContainer`: internal time_series storage
- `internal::InfrastructureSystemsInternal`: power system internal reference, do not modify
"""
mutable struct ThermalFastStartSIIP <: PSY.ThermalGen
    name::String
    available::Bool
    status::Bool
    bus::PSY.Bus
    active_power::Float64
    reactive_power::Float64
    "Thermal limited MVA Power Output of the unit. <= Capacity"
    rating::Float64
    active_power_limits::NamedTuple{(:min, :max), Tuple{Float64, Float64}}
    reactive_power_limits::Union{Nothing, PSY.Min_Max}
    "ramp up and ramp down limits in MW (in component base per unit) per minute"
    ramp_limits::Union{Nothing, NamedTuple{(:up, :down), Tuple{Float64, Float64}}}
    operation_cost::PSY.OperationalCost
    "Base power of the unit in MVA"
    base_power::Float64
    "Minimum up and Minimum down time limits in hours"
    time_limits::Union{Nothing, NamedTuple{(:up, :down), Tuple{Float64, Float64}}}
    "Prime mover technology according to EIA 923"
    prime_mover::PSY.PrimeMovers
    "Prime mover fuel according to EIA 923"
    fuel::PSY.ThermalFuels
    "Services that this device contributes to"
    services::Vector{PSY.Service}
    time_at_status::Float64
    "corresponding dynamic injection device"
    dynamic_injector::Union{Nothing, PSY.DynamicInjection}
    ext::Dict{String, Any}
    "internal time_series storage"
    time_series_container::InfrastructureSystems.TimeSeriesContainer
    "power system internal reference, do not modify"
    internal::InfrastructureSystems.InfrastructureSystemsInternal
end

function ThermalFastStartSIIP(name, available, status, bus, active_power, reactive_power, rating, active_power_limits, reactive_power_limits, ramp_limits, operation_cost, base_power, time_limits=nothing, prime_mover=PSY.PrimeMovers.OT, fuel=PSY.ThermalFuels.OTHER, services=PSY.Device[], time_at_status=PSY.INFINITE_TIME, dynamic_injector=nothing, ext=Dict{String, Any}(), time_series_container=InfrastructureSystems.TimeSeriesContainer(), )
    ThermalFastStartSIIP(name, available, status, bus, active_power, reactive_power, rating, active_power_limits, reactive_power_limits, ramp_limits, operation_cost, base_power, time_limits, prime_mover, fuel, services, time_at_status, dynamic_injector, ext, time_series_container, InfrastructureSystemsInternal(), )
end

function ThermalFastStartSIIP(; name, available, status, bus, active_power, reactive_power, rating, active_power_limits, reactive_power_limits, ramp_limits, operation_cost, base_power, time_limits=nothing, prime_mover=PSY.PrimeMovers.OT, fuel=PSY.ThermalFuels.OTHER, services=PSY.Device[], time_at_status=PSY.INFINITE_TIME, dynamic_injector=nothing, ext=Dict{String, Any}(), time_series_container=InfrastructureSystems.TimeSeriesContainer(), internal=InfrastructureSystems.InfrastructureSystemsInternal(), )
    ThermalFastStartSIIP(name, available, status, bus, active_power, reactive_power, rating, active_power_limits, reactive_power_limits, ramp_limits, operation_cost, base_power, time_limits, prime_mover, fuel, services, time_at_status, dynamic_injector, ext, time_series_container, internal, )
end

# Constructor for demo purposes; non-functional.
function ThermalFastStartSIIP(::Nothing)
    ThermalFastStartSIIP(;
        name="init",
        available=false,
        status=false,
        bus=Bus(nothing),
        active_power=0.0,
        reactive_power=0.0,
        rating=0.0,
        active_power_limits=(min=0.0, max=0.0),
        reactive_power_limits=nothing,
        ramp_limits=nothing,
        operation_cost=ThreePartCost(nothing),
        base_power=0.0,
        time_limits=nothing,
        prime_mover=PrimeMovers.OT,
        fuel=ThermalFuels.OTHER,
        services=Device[],
        time_at_status=INFINITE_TIME,
        dynamic_injector=nothing,
        ext=Dict{String, Any}(),
        time_series_container=InfrastructureSystems.TimeSeriesContainer(),
    )
end

"""Get [`ThermalFastStartSIIP`](@ref) `name`."""
PSY.get_name(value::ThermalFastStartSIIP) = value.name
"""Get [`ThermalFastStartSIIP`](@ref) `available`."""
PSY.get_available(value::ThermalFastStartSIIP) = value.available
"""Get [`ThermalFastStartSIIP`](@ref) `status`."""
PSY.get_status(value::ThermalFastStartSIIP) = value.status
"""Get [`ThermalFastStartSIIP`](@ref) `bus`."""
PSY.get_bus(value::ThermalFastStartSIIP) = value.bus
"""Get [`ThermalFastStartSIIP`](@ref) `active_power`."""
PSY.get_active_power(value::ThermalFastStartSIIP) = get_value(value, value.active_power)
"""Get [`ThermalFastStartSIIP`](@ref) `reactive_power`."""
PSY.get_reactive_power(value::ThermalFastStartSIIP) = get_value(value, value.reactive_power)
"""Get [`ThermalFastStartSIIP`](@ref) `rating`."""
PSY.get_rating(value::ThermalFastStartSIIP) = get_value(value, value.rating)
"""Get [`ThermalFastStartSIIP`](@ref) `active_power_limits`."""
PSY.get_active_power_limits(value::ThermalFastStartSIIP) = get_value(value, value.active_power_limits)
"""Get [`ThermalFastStartSIIP`](@ref) `reactive_power_limits`."""
PSY.get_reactive_power_limits(value::ThermalFastStartSIIP) = get_value(value, value.reactive_power_limits)
"""Get [`ThermalFastStartSIIP`](@ref) `ramp_limits`."""
PSY.get_ramp_limits(value::ThermalFastStartSIIP) = get_value(value, value.ramp_limits)
"""Get [`ThermalFastStartSIIP`](@ref) `operation_cost`."""
PSY.get_operation_cost(value::ThermalFastStartSIIP) = value.operation_cost
"""Get [`ThermalFastStartSIIP`](@ref) `base_power`."""
PSY.get_base_power(value::ThermalFastStartSIIP) = value.base_power
"""Get [`ThermalFastStartSIIP`](@ref) `time_limits`."""
PSY.get_time_limits(value::ThermalFastStartSIIP) = value.time_limits
"""Get [`ThermalFastStartSIIP`](@ref) `prime_mover`."""
PSY.get_prime_mover(value::ThermalFastStartSIIP) = value.prime_mover
"""Get [`ThermalFastStartSIIP`](@ref) `fuel`."""
PSY.get_fuel(value::ThermalFastStartSIIP) = value.fuel
"""Get [`ThermalFastStartSIIP`](@ref) `services`."""
PSY.get_services(value::ThermalFastStartSIIP) = value.services
"""Get [`ThermalFastStartSIIP`](@ref) `time_at_status`."""
PSY.get_time_at_status(value::ThermalFastStartSIIP) = value.time_at_status
"""Get [`ThermalFastStartSIIP`](@ref) `dynamic_injector`."""
PSY.get_dynamic_injector(value::ThermalFastStartSIIP) = value.dynamic_injector
"""Get [`ThermalFastStartSIIP`](@ref) `ext`."""
PSY.get_ext(value::ThermalFastStartSIIP) = value.ext
"""Get [`ThermalFastStartSIIP`](@ref) `time_series_container`."""
PSY.get_time_series_container(value::ThermalFastStartSIIP) = value.time_series_container
"""Get [`ThermalFastStartSIIP`](@ref) `internal`."""
PSY.get_internal(value::ThermalFastStartSIIP) = value.internal

"""Set [`ThermalFastStartSIIP`](@ref) `name`."""
PSY.set_name!(value::ThermalFastStartSIIP, val) = value.name = val
"""Set [`ThermalFastStartSIIP`](@ref) `available`."""
PSY.set_available!(value::ThermalFastStartSIIP, val) = value.available = val
"""Set [`ThermalFastStartSIIP`](@ref) `status`."""
PSY.set_status!(value::ThermalFastStartSIIP, val) = value.status = val
"""Set [`ThermalFastStartSIIP`](@ref) `bus`."""
PSY.set_bus!(value::ThermalFastStartSIIP, val) = value.bus = val
"""Set [`ThermalFastStartSIIP`](@ref) `active_power`."""
PSY.set_active_power!(value::ThermalFastStartSIIP, val) = value.active_power = set_value(value, val)
"""Set [`ThermalFastStartSIIP`](@ref) `reactive_power`."""
PSY.set_reactive_power!(value::ThermalFastStartSIIP, val) = value.reactive_power = set_value(value, val)
"""Set [`ThermalFastStartSIIP`](@ref) `rating`."""
PSY.set_rating!(value::ThermalFastStartSIIP, val) = value.rating = set_value(value, val)
"""Set [`ThermalFastStartSIIP`](@ref) `active_power_limits`."""
PSY.set_active_power_limits!(value::ThermalFastStartSIIP, val) = value.active_power_limits = set_value(value, val)
"""Set [`ThermalFastStartSIIP`](@ref) `reactive_power_limits`."""
PSY.set_reactive_power_limits!(value::ThermalFastStartSIIP, val) = value.reactive_power_limits = set_value(value, val)
"""Set [`ThermalFastStartSIIP`](@ref) `ramp_limits`."""
PSY.set_ramp_limits!(value::ThermalFastStartSIIP, val) = value.ramp_limits = set_value(value, val)
"""Set [`ThermalFastStartSIIP`](@ref) `operation_cost`."""
PSY.set_operation_cost!(value::ThermalFastStartSIIP, val) = value.operation_cost = val
"""Set [`ThermalFastStartSIIP`](@ref) `base_power`."""
PSY.set_base_power!(value::ThermalFastStartSIIP, val) = value.base_power = val
"""Set [`ThermalFastStartSIIP`](@ref) `time_limits`."""
PSY.set_time_limits!(value::ThermalFastStartSIIP, val) = value.time_limits = val
"""Set [`ThermalFastStartSIIP`](@ref) `prime_mover`."""
PSY.set_prime_mover!(value::ThermalFastStartSIIP, val) = value.prime_mover = val
"""Set [`ThermalFastStartSIIP`](@ref) `fuel`."""
PSY.set_fuel!(value::ThermalFastStartSIIP, val) = value.fuel = val
"""Set [`ThermalFastStartSIIP`](@ref) `services`."""
PSY.set_services!(value::ThermalFastStartSIIP, val) = value.services = val
"""Set [`ThermalFastStartSIIP`](@ref) `time_at_status`."""
PSY.set_time_at_status!(value::ThermalFastStartSIIP, val) = value.time_at_status = val
"""Set [`ThermalFastStartSIIP`](@ref) `ext`."""
PSY.set_ext!(value::ThermalFastStartSIIP, val) = value.ext = val
"""Set [`ThermalFastStartSIIP`](@ref) `time_series_container`."""
PSY.set_time_series_container!(value::ThermalFastStartSIIP, val) = value.time_series_container = val
