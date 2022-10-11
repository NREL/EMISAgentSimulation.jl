"""
This struct contains the reserve up market parameters
for price forcasts using CEM and endogeneous Economic Dispatch.
The Reserve Up Market is implemented using a piece-wise linear Operating Reserve Demand Curve (ORDC),
with break-points and price-points parameterized within this struct.
"""

struct ReserveUpMarket{T} #

    demand::Vector{Float64}         # MW
    price_cap::Float64              # $/MWh
    zones::Vector{String}
    eligible_projects::Vector{String}

    function ReserveUpMarket{}(demand::Vector{Float64}, price_cap::Float64, zones::Vector{String}, eligible_projects::Vector{String})
        @assert all(demand .>= 0)
        @assert price_cap >= 0
        T = length(demand)
        new{T}(demand, price_cap, zones, eligible_projects)
    end

end
