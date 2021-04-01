"""
This struct holds the market price data for each product.
    energy_price: Dictionary of zonal energy prices with hourly resolution for each scenario.
    reserveup_price: Dictionary of zonal reserve up prices with hourly resolution for each scenario.
    reservedown_price: Dictionary of zonal reserve down prices with hourly resolution for each scenario.
    synchronous_reserve_price: Dictionary synchronous reserve prices with hourly resolution for each scenario.
    capacity_price: Dictionary of annual capacity prices for each scenario.
    rec_price: Dictionary of annual REC prices for each scenario.
"""

mutable struct MarketPrices
    energy_price::Union{Nothing, Dict{String, AxisArrays.AxisArray{Float64, 3}}}
    reserveup_price::Union{Nothing, Dict{String, AxisArrays.AxisArray{Float64, 3}}}
    reservedown_price::Union{Nothing, Dict{String, AxisArrays.AxisArray{Float64, 3}}}
    synchronous_reserve_price::Union{Nothing, Dict{String, AxisArrays.AxisArray{Float64, 2}}}
    capacity_price::Union{Nothing, Dict{String, AxisArrays.AxisArray{Float64, 1}}}
    rec_price::Union{Nothing, Dict{String, AxisArrays.AxisArray{Float64, 1}}}
end

function MarketPrices()
    return MarketPrices(nothing, nothing, nothing, nothing, nothing, nothing)
end

get_energy_price(prices::MarketPrices) = prices.energy_price
get_reserveup_price(prices::MarketPrices) = prices.reserveup_price
get_reservedown_price(prices::MarketPrices) = prices.reservedown_price
get_synchronous_reserve_price(prices::MarketPrices) = prices.synchronous_reserve_price
get_capacity_price(prices::MarketPrices) = prices.capacity_price
get_rec_price(prices::MarketPrices) = prices.rec_price

get_prices(prices::MarketPrices, prod::T) where T<: Product = nothing
get_prices(prices::MarketPrices, prod::Energy) = get_energy_price(prices)
get_prices(prices::MarketPrices, prod::OperatingReserve{ReserveDownEMIS}) = get_reservedown_price(prices)

function get_prices(prices::MarketPrices, prod::OperatingReserve{ReserveUpEMIS})
    product_name = get_name(prod)
    if product_name == :ReserveUp
        price_data = get_reserveup_price(prices)
    elseif product_name == :SynchronousReserve
        price_data = get_synchronous_reserve_price(prices)
    end
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

function set_reserveup_price!(prices::MarketPrices, scenario_name::String, reserveup_price::AxisArrays.AxisArray{Float64, 3})
    if !isnothing(prices.reserveup_price)
        prices.reserveup_price[scenario_name] = reserveup_price
    else
        prices.reserveup_price = Dict(scenario_name => reserveup_price)
    end
    return
end

function set_reservedown_price!(prices::MarketPrices, scenario_name::String, reservedown_price::AxisArrays.AxisArray{Float64, 3})
    if !isnothing(prices.reservedown_price)
        prices.reservedown_price[scenario_name] = reservedown_price
    else
        prices.reservedown_price = Dict(scenario_name => reservedown_price)
    end
    return
end

function set_synchronous_reserve_price!(prices::MarketPrices, scenario_name::String, synchronous_reserve_price::AxisArrays.AxisArray{Float64, 2})
    if !isnothing(prices.synchronous_reserve_price)
        prices.synchronous_reserve_price[scenario_name] = synchronous_reserve_price
    else
        prices.synchronous_reserve_price = Dict(scenario_name => synchronous_reserve_price)
    end
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
