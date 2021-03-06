"""
This struct contains the reserve up market parameters
for price forcasts using CEM and endogeneous Economic Dispatch.
The Reserve Up Market is implemented using a piece-wise linear Operating Reserve Demand Curve (ORDC),
with break-points and price-points parameterized within this struct.
"""

struct ReserveUpMarket

    break_points::AxisArrays.AxisArray{Vector{Float64}}   # ORDC break points
    price_points::AxisArrays.AxisArray{Vector{Float64}}    # $/MWh
    demand::AxisArrays.AxisArray{Float64, 2}

    function ReserveUpMarket(b, p, d)
        @assert length(p) == length(b)
        @assert all(d .>= 0)
        for i = 1:length(b)
            @assert all(b[i] .>= 0)
            @assert all(p[i] .>= 0)

            @assert length(b[i]) > 1
            @assert length(b[i]) == length(p[i])
        end
        new(b, p, d)
    end

end
