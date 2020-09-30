"""
This function calculates the operating market profits.
"""
function calculate_realized_operating_profit(prices::AxisArrays.AxisArray{Float64, 2},
                                            marginal_cost::Float64,
                                            output::Array{Float64, 2},
                                            realized_hour_weight::Vector{Float64})

    profit = sum(realized_hour_weight[i] *
                ((prices[1, i] - marginal_cost) * output[1, i])
                 for i in 1:length(realized_hour_weight))

    return profit

end

"""
This function does nothing if the method is not specified for a product.
"""
function calculate_realized_profit(project_name::String,
                                   size::Float64,
                                   zone::String,
                                   product::T,
                                   market_prices::MarketPrices,
                                   capacity_factors::Dict{String, Array{Float64, 2}},
                                   reserve_up_perc::Dict{String, Array{Float64, 2}},
                                   reserve_down_perc::Dict{String, Array{Float64, 2}},
                                   capacity_accepted_bids::Dict{String, Float64},
                                   rec_accepted_bids::Dict{String, Float64},
                                   realized_hour_weight::Vector{Float64},
                                   iteration_year::Int64,
                                   capacity_forward_years::Int64) where T <: Product
    return nothing
end

"""
This function calculates the realized energy market profits.
"""
function calculate_realized_profit(project_name::String,
                                   size::Float64,
                                   zone::String,
                                   product::Energy,
                                   market_prices::MarketPrices,
                                   capacity_factors::Dict{String, Array{Float64, 2}},
                                   reserve_up_perc::Dict{String, Array{Float64, 2}},
                                   reserve_down_perc::Dict{String, Array{Float64, 2}},
                                   capacity_accepted_bids::Dict{String, Float64},
                                   rec_accepted_bids::Dict{String, Float64},
                                   realized_hour_weight::Vector{Float64},
                                   iteration_year::Int64,
                                   capacity_forward_years::Int64)

    update_year = iteration_year

    if in(project_name, keys(capacity_factors))
        output = size * capacity_factors[project_name]

        profit = calculate_realized_operating_profit(get_prices(market_prices, product)["realized"][zone, :, :],
                                                get_marginal_cost(product),
                                                output,
                                                realized_hour_weight)

        return profit, update_year
    else
        return nothing, update_year
    end
end

"""
This function calculates the realized reserve up market profits.
"""
function calculate_realized_profit(project_name::String,
                                  size::Float64,
                                  zone::String,
                                  product::OperatingReserve{ReserveUpEMIS},
                                  market_prices::MarketPrices,
                                  capacity_factors::Dict{String, Array{Float64, 2}},
                                  reserve_up_perc::Dict{String, Array{Float64, 2}},
                                  reserve_down_perc::Dict{String, Array{Float64, 2}},
                                  capacity_accepted_bids::Dict{String, Float64},
                                  rec_accepted_bids::Dict{String, Float64},
                                  realized_hour_weight::Vector{Float64},
                                  iteration_year::Int64,
                                  capacity_forward_years::Int64)

    update_year = iteration_year

    if in(project_name, keys(reserve_up_perc))
        output = size * reserve_up_perc[project_name]

        profit =calculate_realized_operating_profit(get_prices(market_prices, product)["realized"][zone, :, :],
                                                get_marginal_cost(product),
                                                output,
                                                realized_hour_weight)

        return profit, update_year
    else
        return nothing, update_year
    end
end

"""
This function calculates the realized reserve down market profits.
"""
function calculate_realized_profit(project_name::String,
                                  size::Float64,
                                  zone::String,
                                  product::OperatingReserve{ReserveDownEMIS},
                                  market_prices::MarketPrices,
                                  capacity_factors::Dict{String, Array{Float64, 2}},
                                  reserve_up_perc::Dict{String, Array{Float64, 2}},
                                  reserve_down_perc::Dict{String, Array{Float64, 2}},
                                  capacity_accepted_bids::Dict{String, Float64},
                                  rec_accepted_bids::Dict{String, Float64},
                                  realized_hour_weight::Vector{Float64},
                                  iteration_year::Int64,
                                  capacity_forward_years::Int64)

    update_year = iteration_year

    if in(project_name, keys(reserve_down_perc))
        output = size * reserve_down_perc[project_name]

        profit = calculate_realized_operating_profit(get_prices(market_prices, product)["realized"][zone, :, :],
                                                get_marginal_cost(product),
                                                output,
                                                realized_hour_weight)

        return profit, update_year
    else
        return nothing, update_year
    end
end

"""
This function calculates the realized capacity market profits.
"""
function calculate_realized_profit(project_name::String,
                                  size::Float64,
                                  zone::String,
                                  product::Capacity,
                                  market_prices::MarketPrices,
                                  capacity_factors::Dict{String, Array{Float64, 2}},
                                  reserve_up_perc::Dict{String, Array{Float64, 2}},
                                  reserve_down_perc::Dict{String, Array{Float64, 2}},
                                  capacity_accepted_bids::Dict{String, Float64},
                                  rec_accepted_bids::Dict{String, Float64},
                                  realized_hour_weight::Vector{Float64},
                                  iteration_year::Int64,
                                  capacity_forward_years::Int64)

    update_year = iteration_year + capacity_forward_years - 1

    if in(project_name, keys(capacity_accepted_bids))
        profit = size *
                get_derating(product) *
                get_prices(market_prices, product)["realized"][1] *
                capacity_accepted_bids[project_name]

    return profit, update_year
    else
        return nothing, update_year
    end
end

"""
This function calculates the realized REC market profits.
"""
function calculate_realized_profit(project_name::String,
                                  size::Float64,
                                  zone::String,
                                  product::REC,
                                  market_prices::MarketPrices,
                                  capacity_factors::Dict{String, Array{Float64, 2}},
                                  reserve_up_perc::Dict{String, Array{Float64, 2}},
                                  reserve_down_perc::Dict{String, Array{Float64, 2}},
                                  capacity_accepted_bids::Dict{String, Float64},
                                  rec_accepted_bids::Dict{String, Float64},
                                  realized_hour_weight::Vector{Float64},
                                  iteration_year::Int64,
                                  capacity_forward_years::Int64)

    update_year = iteration_year

    if in(project_name, keys(rec_accepted_bids))
        profit = get_prices(market_prices, product)["realized"][1] *
                rec_accepted_bids[project_name]

    return profit, update_year
    else
        return nothing, update_year
    end
end
