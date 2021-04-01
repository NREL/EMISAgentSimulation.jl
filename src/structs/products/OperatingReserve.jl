"""
This struct contains data for the operating reserve product
parameterized by the reserve direction (up or down).
    name: ReserveUp or ReserveDown.
    max_limit: Maximum percentage of reserve that can be provided by the project.
    marginal_cost: Marginal cost per unit of reserve provided.
"""
struct OperatingReserve{T <: ReserveDirection} <: ReserveEMIS{T}
    name::Symbol
    max_limit::Float64
    marginal_cost::Float64
    zonal::Bool
end

get_max_limit(prod::OperatingReserve) = prod.max_limit
