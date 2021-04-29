"""
This struct contains the reserve up market parameters
for price forcasts using CEM and endogeneous Economic Dispatch.
The Reserve ORDC Market is implemented using a piece-wise linear Operating Reserve Demand Curve (ORDC),
with break-points and price-points parameterized within this struct.
"""

struct ReserveORDCMarket{T}

    break_points::AxisArrays.AxisArray{Vector{Float64}}   # ORDC break points
    price_points::AxisArrays.AxisArray{Vector{Float64}}    # $/MWh
    zones::Vector{String}
    eligible_projects::Vector{String}
    stepped::Bool

    function ReserveORDCMarket(b, p, z, e, s)
        T = length(b)
        @assert length(p) == length(b)
        for i = 1:length(b)
            @assert all(b[i] .>= 0)
            @assert all(p[i] .>= 0)

            @assert length(b[i]) > 1
            @assert length(b[i]) == length(p[i])
        end
        new{T}(b, p, z, e, s)
    end

end