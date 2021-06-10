"""
This struct contains data for the energy market product.
    name: Energy.
    capacity_factors: Dictionary of expected capacity factors for each scenario.
    marginal_cost: Marginal cost per unit of energy provided.
"""
mutable struct Energy <: OperatingProduct
    name::Symbol
    capacity_factors::Dict{String, Array{Float64, 2}}
    marginal_cost::Float64
    expected_production::Float64
end

# Capacity factors only returned when product is of type Energy
function get_capacity_factors(prod::T) where T<:Product
    return
end

get_capacity_factors(prod::Energy) = prod.capacity_factors


# Capacity factors only set when product is of type Energy
function set_capacity_factors!(product::T,
                              scenario_name::String,
                              capacity_factors::Array{Float64, 2}) where T <: Product
    return
end

function set_capacity_factors!(product::Energy,
                              scenario_name::String,
                              capacity_factors::Array{Float64, 2})
    product.capacity_factors[scenario_name] = capacity_factors
    return
end


# Expected energy production only returned when product is of type Energy
function get_expected_production(prod::T) where T<:Product
    return
end

get_expected_production(prod::Energy) = prod.expected_production

# Expected energy production only set when product is of type Energy
function set_expected_production!(product::T,
                                expected_production::Float64) where T <: Product
    return
end

function set_expected_production!(product::Energy,
                                expected_production::Float64)
    product.expected_production = expected_production
    return
end
