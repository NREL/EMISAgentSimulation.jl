"""
This struct holds the market price data for each product.
    energy_price: Dictionary of zonal energy prices with hourly resolution for each scenario.
    reserve_price: Dictionary of reserve prices with hourly resolution for each scenario.
    capacity_price: Dictionary of annual capacity prices for each scenario.
    rec_price: Dictionary of annual REC prices for each scenario.
"""

mutable struct MarketPrices
    energy_price::Union{Nothing, Dict{String, AxisArrays.AxisArray{Float64, 3}}}
    reserve_price::Union{Nothing, Dict{String, Dict{String, Array{Float64, 2}}}}
    capacity_price::Union{Nothing, Dict{String, AxisArrays.AxisArray{Float64, 1}}}
    rec_price::Union{Nothing, Dict{String, AxisArrays.AxisArray{Float64, 1}}}
end

function MarketPrices()
    return MarketPrices(nothing, nothing, nothing, nothing)
end

get_energy_price(prices::MarketPrices) = prices.energy_price
get_reserve_price(prices::MarketPrices) = prices.reserve_price
get_capacity_price(prices::MarketPrices) = prices.capacity_price
get_rec_price(prices::MarketPrices) = prices.rec_price

get_prices(prices::MarketPrices, prod::T) where T<: Product = nothing
get_prices(prices::MarketPrices, prod::Energy) = get_energy_price(prices)

function get_prices(prices::MarketPrices, prod::Union{OperatingReserve{ReserveUpEMIS}, OperatingReserve{ReserveDownEMIS}})
    product_name = get_name(prod)
        price_data = get_reserve_price(prices)[String(product_name)]
    return price_data
end

get_prices(prices::MarketPrices, prod::Capacity) = get_capacity_price(prices)
get_prices(prices::MarketPrices, prod::REC) = get_rec_price(prices)

function set_energy_price!(prices::MarketPrices, scenario_name::String, energy_price::AxisArrays.AxisArray{Float64, 3})
    if !isnothing(prices.energy_price)
        prices.energy_price[scenario_name] = energy_price
    else
        prices.energy_price = Dict(scenario_name => energy_price)
    end
    return
end

function set_reserve_price!(prices::MarketPrices, scenario_name::String, reserve_price::Dict{String, Array{Float64, 2}})

    reserve_price_dict = Dict(product => Dict(scenario_name => reserve_price[product]) for product in keys(reserve_price))
    prices.reserve_price = reserve_price_dict

    return
end

function set_capacity_price!(prices::MarketPrices, scenario_name::String, capacity_price::AxisArrays.AxisArray{Float64, 1})
    if !isnothing(prices.capacity_price)
        prices.capacity_price[scenario_name] = capacity_price
    else
        prices.capacity_price = Dict(scenario_name => capacity_price)
    end
    return
end

function set_rec_price!(prices::MarketPrices, scenario_name::String, rec_price::AxisArrays.AxisArray{Float64, 1})
    if !isnothing(prices.rec_price)
        prices.rec_price[scenario_name] = rec_price
    else
        prices.rec_price = Dict(scenario_name => rec_price)
    end
    return
end
