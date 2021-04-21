"""
This struct contains data for the energy market product.
    name: CarbonTax.
    emission: Emission Intensity (ton/MMBTU).
    avg_heat_rate: Heat Rate (MMBTU/MWh)
    fuel_cost: Fuel Cost (\$/MMBTU)
"""
mutable struct CarbonTax <: Product
    name::Symbol
    emission_intensity::Float64
    avg_heat_rate::Float64
    fuel_cost::Float64
    total_emission::Vector{Float64}
end

# Emission only returned when product is of type CarbonTax
function get_emission_intensity(prod::T) where T<:Product
    return
end

get_emission_intensity(prod::CarbonTax) = prod.emission_intensity

# Heat Rate only returned when product is of type CarbonTax
function get_avg_heat_rate(prod::T) where T<:Product
    return
end

get_avg_heat_rate(prod::CarbonTax) = prod.avg_heat_rate

# Fuel Cost only returned when product is of type CarbonTax
function get_fuel_cost(prod::T) where T<:Product
    return
end

get_fuel_cost(prod::CarbonTax) = prod.fuel_cost

# Total emissionsy only returned when product is of type CarbonTax
function get_total_emission(prod::T) where T<:Product
    return
end

get_total_emission(prod::CarbonTax) = prod.total_emission


