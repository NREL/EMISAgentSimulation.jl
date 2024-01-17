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
        chron_weights::Matrix{Int64}
    ) where {Z, T}
        @assert T == length(rep_hour_weight)
        @assert T == size(chron_weights, 2)
        @assert length(zones) ==  Z
        @assert length(inv_periods) == length(carbon_tax)
        @assert allunique(getproperty.(projects, :name))
        new{Z, T}(zones, lines, capital_cost_multiplier, inv_periods, carbon_tax, projects, rep_period_interval, rep_hour_weight, chron_weights)
    end

end
