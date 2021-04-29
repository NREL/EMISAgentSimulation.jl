"""
This function calculates the operating market profits.
"""
function calculate_realized_operating_profit(prices::AxisArrays.AxisArray{Float64, 2},
                                            marginal_cost::Float64,
                                            output::Array{Float64, 2},
                                            realized_hour_weight::Vector{Float64})

    profit = sum(realized_hour_weight[i] *
                ((prices[1, i] - marginal_cost) * output[1, i])
                 for i in 1:size(output, 2))

    return profit

end

"""
This function calculates the operating market profits.
"""
function calculate_realized_operating_profit(prices::AxisArrays.AxisArray{Float64, 2},
                                            marginal_cost::Vector{Float64},
                                            output::Array{Float64, 2},
                                            carbon_emissions::Vector{Float64},
                                            carbon_tax::Float64,
                                            realized_hour_weight::Vector{Float64},
                                            resolution::Int64)

    profit = sum(realized_hour_weight[i] *
                ((prices[1, i] - marginal_cost[i]) * output[1, i] - carbon_emissions[i] * carbon_tax)
                 for i in 1:size(output, 2))


    return profit

end

"""
This function calculates the operating market profits.
"""
function calculate_realized_operating_profit(prices::Array{Float64, 2},
                                            marginal_cost::Float64,
                                            output::Array{Float64, 2},
                                            realized_hour_weight::Vector{Float64})

    profit = sum(realized_hour_weight[i] *
                ((prices[1, i] - marginal_cost) * output[1, i])
                 for i in 1:size(output, 2))
    return profit

end

function calculate_carbon_emissions(emission_intensity::Float64,
                                    heat_rate_curve::Nothing,
                                    output::Array{Float64, 2})
    emissions = zeros(size(output, 2))
    return emissions
end

function calculate_carbon_emissions(emission_intensity::Float64,
                                    heat_rate_curve::Vector{Tuple{Float64, Float64}},
                                    output::Array{Float64, 2})

    emissions = zeros(size(output, 2))

    for t in 1:length(emissions)
        for s in 2:length(heat_rate_curve)
            if (output[1, t] > heat_rate_curve[s - 1][2]) &&  (output[1, t] <= heat_rate_curve[s][2])
                slope = (heat_rate_curve[s][1] - heat_rate_curve[s - 1][1])/(heat_rate_curve[s][2] - heat_rate_curve[s - 1][2])
                ihr = heat_rate_curve[s - 1][1] + slope * (output[1, t] - heat_rate_curve[s - 1][2])
                emissions[t] = ihr * emission_intensity * output[1, t]
            end
        end
    end
    return emissions
end


function calculate_marginal_cost(heat_rate_curve::Nothing,
                                  output::Array{Float64, 2},
                                  fuel_cost::Float64)
    marginal_cost = zeros(size(output, 2))
    return marginal_cost
end

function calculate_marginal_cost(heat_rate_curve::Vector{Tuple{Float64, Float64}},
                                output::Array{Float64, 2},
                                fuel_cost::Float64)

    marginal_cost = zeros(size(output, 2))

    for t in 1:length(marginal_cost)
        for s in 2:length(heat_rate_curve)
            if (output[1, t] > heat_rate_curve[s - 1][2]) &&  (output[1, t] <= heat_rate_curve[s][2])
                slope = (heat_rate_curve[s][1] - heat_rate_curve[s - 1][1])/(heat_rate_curve[s][2] - heat_rate_curve[s - 1][2])
                ihr = heat_rate_curve[s - 1][1] + slope * (output[1, t] - heat_rate_curve[s - 1][2])
                marginal_cost[t] = ihr * fuel_cost
            end
        end
    end
    return marginal_cost
end

"""
This function does nothing if the method is not specified for a product.
"""
function calculate_realized_profit(project::Project,
                                   product::T,
                                   market_prices::MarketPrices,
                                   capacity_factors::Dict{String, Array{Float64, 2}},
                                   reserve_perc::Dict{String, Dict{String, Array{Float64, 2}}},
                                   capacity_accepted_bids::Dict{String, Float64},
                                   rec_accepted_bids::Dict{String, Float64},
                                   realized_hour_weight::Vector{Float64},
                                   iteration_year::Int64,
                                   capacity_forward_years::Int64,
                                   carbon_tax::Float64,
                                   da_resolution::Int64,
                                   rt_resolution::Int64,
                                   rt_products::Vector{String}) where T <: Product

    return nothing, 0
end

"""
This function calculates the realized energy market profits.
"""
function calculate_realized_profit(project::Project,
                                   product::Energy,
                                   market_prices::MarketPrices,
                                   capacity_factors::Dict{String, Array{Float64, 2}},
                                   reserve_perc::Dict{String, Dict{String, Array{Float64, 2}}},
                                   capacity_accepted_bids::Dict{String, Float64},
                                   rec_accepted_bids::Dict{String, Float64},
                                   realized_hour_weight::Vector{Float64},
                                   iteration_year::Int64,
                                   capacity_forward_years::Int64,
                                   carbon_tax::Float64,
                                   da_resolution::Int64,
                                   rt_resolution::Int64,
                                   rt_products::Vector{String})

    project_name = get_name(project)
    size = get_maxcap(project)
    zone = get_zone(get_tech(project))

    update_year = iteration_year
    if in(project_name, keys(capacity_factors))
        output = size * capacity_factors[project_name]
        emission_intensity = get_emission_intensity(project)
        heat_rate_curve = get_heat_rate_curve(get_tech(project))

        carbon_emissions = calculate_carbon_emissions(emission_intensity, heat_rate_curve, output) .* (rt_resolution / 60)
        marginal_cost = calculate_marginal_cost(heat_rate_curve, output, get_fuel_cost(get_tech(project))) .* (rt_resolution / 60)

        profit = calculate_realized_operating_profit(get_prices(market_prices, product)["realized"][zone, :, :],
                                                marginal_cost,
                                                output,
                                                carbon_emissions,
                                                carbon_tax,
                                                realized_hour_weight,
                                                rt_resolution)

        for product in get_products(project)
            set_total_emission!(product, carbon_emissions)
        end
        return profit, update_year
    else
        return nothing, update_year
    end
end

"""
This function calculates the realized reserve up market profits.
"""
function calculate_realized_profit(project::Project,
                                  product::Union{OperatingReserve{ReserveUpEMIS}, OperatingReserve{ReserveDownEMIS}},
                                  market_prices::MarketPrices,
                                  capacity_factors::Dict{String, Array{Float64, 2}},
                                  reserve_perc::Dict{String, Dict{String, Array{Float64, 2}}},
                                  capacity_accepted_bids::Dict{String, Float64},
                                  rec_accepted_bids::Dict{String, Float64},
                                  realized_hour_weight::Vector{Float64},
                                  iteration_year::Int64,
                                  capacity_forward_years::Int64,
                                  carbon_tax::Float64,
                                  da_resolution::Int64,
                                  rt_resolution::Int64,
                                  rt_products::Vector{String})

    project_name = get_name(project)
    size = get_maxcap(project)

    if get_name(product) in rt_products
        scale = rt_resolution / 60
    else
        scale = da_resolution / 60
    end

    update_year = iteration_year
    if in(project_name, keys(reserve_perc))
        output = size * reserve_perc[project_name][String(get_name(product))]

        profit = calculate_realized_operating_profit(get_prices(market_prices, product)["realized"],
                                                get_marginal_cost(product) * scale,
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
function calculate_realized_profit(project::Project,
                                  product::Capacity,
                                  market_prices::MarketPrices,
                                  capacity_factors::Dict{String, Array{Float64, 2}},
                                  reserve_perc::Dict{String, Dict{String, Array{Float64, 2}}},
                                  capacity_accepted_bids::Dict{String, Float64},
                                  rec_accepted_bids::Dict{String, Float64},
                                  realized_hour_weight::Vector{Float64},
                                  iteration_year::Int64,
                                  capacity_forward_years::Int64,
                                  carbon_tax::Float64,
                                  da_resolution::Int64,
                                  rt_resolution::Int64,
                                  rt_products::Vector{String})

    project_name = get_name(project)
    size = get_maxcap(project)

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
function calculate_realized_profit(project::Project,
                                  product::REC,
                                  market_prices::MarketPrices,
                                  capacity_factors::Dict{String, Array{Float64, 2}},
                                  reserve_perc::Dict{String, Dict{String, Array{Float64, 2}}},
                                  capacity_accepted_bids::Dict{String, Float64},
                                  rec_accepted_bids::Dict{String, Float64},
                                  realized_hour_weight::Vector{Float64},
                                  iteration_year::Int64,
                                  capacity_forward_years::Int64,
                                  carbon_tax::Float64,
                                  da_resolution::Int64,
                                  rt_resolution::Int64,
                                  rt_products::Vector{String})

    project_name = get_name(project)
    update_year = iteration_year
    if in(project_name, keys(rec_accepted_bids))
        profit = get_prices(market_prices, product)["realized"][1] *
                rec_accepted_bids[project_name]

    return profit, update_year
    else
        return nothing, update_year
    end
end
