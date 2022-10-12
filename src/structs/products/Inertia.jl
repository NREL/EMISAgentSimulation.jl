"""
This struct contains data for the capacity market product.
    name: Inertia.
    synchronous: Whether unit provides synchronous inertia.
    h_constant: H-constant of the unit.
"""
mutable struct Inertia <: OperatingProduct
    name::Symbol
    synchronous::Bool
    h_constant::Float64
    marginal_cost::Float64
end

# Derating factors, accepted percentage and capacity bids only returned when product is of type Capacity
get_synchronous(prod::Product) = nothing
get_synchronous(prod::Inertia) = prod.synchronous
get_h_constant(prod::Product) = nothing
get_h_constant(prod::Inertia) = prod.h_constant

