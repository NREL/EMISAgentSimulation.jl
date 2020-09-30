"""
This struct contains all the data for the transmission lines between zones.
    name: Path to the directory where input data is stored.
    from_zone: Name of originating zone.
    to_zone: Name of destination zone.
    active_power_limit: Maximum flow capacity (MW).
"""
mutable struct ZonalLine
    name::String
    from_zone::String
    to_zone::String
    active_power_limit::Float64
end

get_name(line::ZonalLine) = line.name
get_from_zone(line::ZonalLine) = line.from_zone
get_to_zone(line::ZonalLine) = line.to_zone
get_active_power_limit(line::ZonalLine) = line.active_power_limit

function  set_active_power_limit!(line::ZonalLine, active_power_limit::Float64)
    line.active_power_limit = active_power_limit
    return
end
