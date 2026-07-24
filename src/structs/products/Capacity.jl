"""
This struct contains data for the capacity market product.
    name: Capacity.
    derating: The derating factor for capacity market participation based on technology type for each scenario and season. Annual (non-seasonal) cases use a single "annual" season key.
    accepted_perc: The percentage of capacity cleared in capacity markets for each scenario and season. Annual (non-seasonal) cases use a single "annual" season key.
    capacity_bid: Capacity market bid placed by the project.
"""
mutable struct Capacity <: Product
    name::Symbol
    derating::Dict{String, Dict{String, Float64}}
    accepted_perc::Dict{String, Dict{String, Vector{Float64}}}
    capacity_bid::Float64
end

# Derating factors, accepted percentage and capacity bids only returned when product is of type Capacity
get_derating(prod::Product) = nothing
get_derating(prod::Product, scenario::String, season::String) = nothing
get_derating(prod::Capacity) = prod.derating
get_derating(prod::Capacity, scenario::String, season::String) = prod.derating[scenario][season]
get_derating(prod::Capacity, scenario::String) = prod.derating[scenario]["annual"]
get_accepted_perc(prod::Product) = nothing
get_accepted_perc(prod::Product, scenario::String, season::String) = nothing
get_accepted_perc(prod::Capacity) = prod.accepted_perc
get_accepted_perc(prod::Capacity, scenario::String, season::String) = prod.accepted_perc[scenario][season]
get_accepted_perc(prod::Capacity, scenario::String) = prod.accepted_perc[scenario]["annual"]
get_capacity_bid(prod::Product) = nothing
get_capacity_bid(prod::Capacity) = prod.capacity_bid

# Derating factors only set when product is of type Capacity
function set_derating!(prod::T, scenario::String, season::String, derating_factor) where T <: Product
    return
end

function set_derating!(prod::Capacity, scenario::String, season::String, derating_factor)
    if !haskey(prod.derating, scenario)
        prod.derating[scenario] = Dict{String, Float64}()
    end
    prod.derating[scenario][season] = derating_factor
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
                            season::String,
                            capacity_accepted_perc::Array{Float64, 1}) where T <: Product
    return
end

function set_accepted_perc!(product::Capacity,
                            scenario_name::String,
                            season::String,
                            capacity_accepted_perc::Array{Float64, 1})
    if !haskey(product.accepted_perc, scenario_name)
        product.accepted_perc[scenario_name] = Dict{String, Vector{Float64}}()
    end
    product.accepted_perc[scenario_name][season] = capacity_accepted_perc
    return
end
