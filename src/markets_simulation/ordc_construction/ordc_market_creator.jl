"""
This function creates the ReserveORDCMarket struct
based on the vector of point tuples
"""
function create_ordc_market(points::Vector{String})
    T = length(points)
    break_points = AxisArrays.AxisArray([Vector{Float64}() for t in 1:T], (1:T))
    price_points = AxisArrays.AxisArray([Vector{Float64}() for t in 1:T], (1:T))

    for t = 1:T
        tuples = split.(chop.(split(chop(points[t], head = 1, tail = 2), "), "), head = 1, tail = 0), ", ")
        break_points[t] = [parse.(Float64, tuple)[1] for tuple in tuples]
        price_points[t] = [parse.(Float64, tuple)[2] for tuple in tuples]
    end

    ordc_market = ReserveORDCMarket(break_points, price_points)

    return ordc_market
end
