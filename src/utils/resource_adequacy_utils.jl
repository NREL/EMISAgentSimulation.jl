# Helper functions for resource adequacy workflows using PRAS and SiennaPRASInterface (SPI).

"""
    get_regional_load_shares(system::PRAS.SystemModel) -> Dict{String, Float64}

Returns a dictionary mapping each region name to its share of total system load,
computed as each region's cumulative load across all timestamps divided by the
grand total — equivalent to the time-averaged load proportion.
"""
function get_regional_load_shares(system::PRAS.SystemModel)
    region_names = system.regions.names
    regional_load = system.regions.load  # (num_regions, num_timestamps)
    regional_load_totals = sum(regional_load, dims=2)[:, 1]
    total_load = sum(regional_load_totals)
    return Dict(region_names[i] => regional_load_totals[i] / total_load for i in 1:length(region_names))
end

