"""
This struct contains all data for running a CEM for price forecasts or for endoogeneous Economic Dispatch.
    zones: Name of zones included in the market clearing problem.
    lines: Vector of inter-zonal transmission lines included in the market clearing problem.
    capital_cost_multiplier: Average capital cost multiplier for new builds.
    inv_periods: Vector of all markets modeled for the entire CEM planning horizon.
    carbon_tax: Vector of carbon tax for the planning horizon
    projects: Vector of all projects included in the market clearing problem.
    rep_period_interval: Total number of hours in representative period selection e.g., 24 for representative days, 168 for representative weeks, etc. 
    rep_hour_weight: Weight of representative hours in the market clearing problem.
    avg_block_size: Number of hours in CEM aggregated time blocks. Default = 4. Set to 1 if hourly granularity is needed.
    fixed_block_size: Whether the block size in CEM time block aggregation is fixed, i.e., all blocks are of the same length = avg_block_size. If set to FALSE, the model will select variable size time blocks based on chronological clustering.
    chron_weights: Weights for representative hours for modeling chronological constraints
"""

struct MarketClearingProblem{Z, T}

    zones::Vector{String}
    lines::Vector{ZonalLine}
    capital_cost_multiplier::Float64
    inv_periods::Vector{MarketCollection{Z, T}}
    carbon_tax::Vector{Float64}
    projects::Vector{MarketProject}
    rep_period_interval::Int64
    rep_hour_weight::Vector{Float64}
    avg_block_size::Int64
    fixed_block_size::Bool
    chron_weights::Matrix{Int64}

    function MarketClearingProblem(
        zones::Vector{String},
        lines::Vector{ZonalLine},
        capital_cost_multiplier::Float64,
        inv_periods::Vector{MarketCollection{Z, T}},
        carbon_tax::Vector{Float64},
        projects::Vector{MarketProject},
        rep_period_interval::Int64,
        rep_hour_weight::Vector{Float64},
        avg_block_size::Int64,
        fixed_block_size::Bool,
        chron_weights::Matrix{Int64}
    ) where {Z, T}
        @assert T == length(rep_hour_weight)
        @assert T == size(chron_weights, 2)
        @assert length(zones) ==  Z
        @assert length(inv_periods) == length(carbon_tax)
        @assert allunique(getproperty.(projects, :name))
        new{Z, T}(zones, lines, capital_cost_multiplier, inv_periods, carbon_tax, projects, rep_period_interval, rep_hour_weight, avg_block_size, fixed_block_size, chron_weights)
    end

end
