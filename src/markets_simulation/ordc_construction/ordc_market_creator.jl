"""
This function creates the ReserveORDCMarket struct
based on the vector of point tuples
"""
function create_ordc_market(points::Union{Vector{String}, PooledArrays.PooledArray{String,UInt32,1,Array{UInt32,1}}}, parameters::DataFrames.DataFrame, eligible_products::Vector{String})
    T = length(points)
    break_points = AxisArrays.AxisArray([Vector{Float64}() for t in 1:T], (1:T))
    price_points = AxisArrays.AxisArray([Vector{Float64}() for t in 1:T], (1:T))

    for t = 1:T
        tuples = split.(chop.(split(chop(points[t], head = 1, tail = 2), "), "), head = 1, tail = 0), ", ")
        break_points[t] = [parse.(Float64, tuple)[1] for tuple in tuples]
        price_points[t] = [parse.(Float64, tuple)[2] for tuple in tuples]
    end

    zones = ["zone_$(n)" for n in split(parameters[1, "eligible_zones"], ";")]

    stepped = parameters[1, "stepped"]

    ordc_market = ReserveORDCMarket(break_points, price_points, zones, eligible_products, stepped)

    return ordc_market
end
