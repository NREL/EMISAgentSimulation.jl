"""
This struct contains the capacity market parameters
for price forcasts using CEM and endogeneous Economic Dispatch.
The Capacity market is implemented using a piece-wise Capacity Demand Curve,
with break-points and price-points parameterized within this struct.
"""
struct CapacityMarket # Assumes a piece-wise linear demand curve

    break_points::Vector{Float64}   # demand curve break points
    price_points::Vector{Float64}    # $/MW-year

    function CapacityMarket(b, p)
        @assert all(b .>= 0)
        @assert all(p .>= 0)
        @assert length(b) > 1
        @assert length(p) == length(b)
        new(b, p)
    end

end
