"""
This struct contains the reserve down market parameters
for price forcasts using CEM and endogeneous Economic Dispatch.
The reserve down market assumes inelastic demand.
"""

struct ReserveDownMarket{Z, T} #

    demand::AxisArrays.AxisArray{Float64, 2}         # MW
    price_cap::AxisArrays.AxisArray{Float64, 1}       # $/MWh

    function ReserveDownMarket{}(demand::AxisArrays.AxisArray{Float64, 2}, price_cap::AxisArrays.AxisArray{Float64, 1})
        @assert all(demand .>= 0)
        @assert all(price_cap .>= 0)
        Z = size(demand, 1)
        T = size(demand, 2)
        new{Z, T}(demand, price_cap)
    end

end
