"""
This struct contains the reserve down market parameters
for price forcasts using CEM and endogeneous Economic Dispatch.
The reserve down market assumes inelastic demand.
"""

struct ReserveDownMarket{T} #

    demand::Vector{Float64}         # MW
    price_cap::Float64              # $/MWh
    zones::Vector{String}
    eligible_projects::Vector{String}

    function ReserveDownMarket{}(demand::Vector{Float64}, price_cap::Float64, zones::Vector{String}, eligible_projects::Vector{String})
        @assert all(demand .>= 0)
        @assert price_cap >= 0
        T = length(demand)
        new{T}(demand, price_cap, zones, eligible_projects)
    end

end
