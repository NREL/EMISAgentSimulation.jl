"""
This struct contains data for the capacity market product.
    name: Capacity.
    derating: The derating factor for capacity market participation based on technology type.
    accepted_perc: The percentage of capacity cleared in capacity markets for each scenario.
    capacity_bid: Capacity market bid placed by the project.
"""
mutable struct Capacity <: Product
    name::Symbol
    derating::Float64
    accepted_perc::Dict{String, Vector{Float64}}
    capacity_bid::Float64
end

# Derating factors, accepted percentage and capacity bids only returned when product is of type Capacity
get_derating(prod::Product) = nothing
get_derating(prod::Capacity) = prod.derating
get_accepted_perc(prod::Product) = nothing
get_accepted_perc(prod::Capacity) = prod.accepted_perc
get_capacity_bid(prod::Product) = nothing
get_capacity_bid(prod::Capacity) = prod.capacity_bid

# Derating factors only set when product is of type Capacity
function set_derating!(prod::T, derating_factor) where T <: Product
    return
end

function set_derating!(prod::Capacity, derating_factor)
    prod.derating = derating_factor
    return
end

# Capacity bids only set when product is of type Capacity
function set_capacity_bid!(prod::T, capacity_bid) where T <: Product
    return
end

function set_capacity_bid!(prod::Capacity, capacity_bid)
    prod.capacity_bid = capacity_bid
    return
end

# Capacity accepted percentage only set when product is of type Capacity
function set_accepted_perc!(product::T,
                            scenario_name::String,
                            capacity_accepted_perc::Array{Float64, 1}) where T <: Product
    return
end

function set_accepted_perc!(product::Capacity,
                            scenario_name::String,
                            capacity_accepted_perc::Array{Float64, 1})
    product.accepted_perc[scenario_name] = capacity_accepted_perc
    return
end
