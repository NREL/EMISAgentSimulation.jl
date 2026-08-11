"""
This function calculates the operating market profits.
"""
function calculate_realized_operating_profit(
        prices_md::Array{Float64, 2},
        prices_uc::Array{Float64, 2},
        prices_ed::Array{Float64, 2},
        marginal_cost::Float64,
        output_md::Array{Float64, 2},
        output_uc::Array{Float64, 2},
        output_ed::Array{Float64, 2},
        realized_hour_weight::Dict{String, Dict{Int64, Vector{Float64}}})

    replace!(prices_md, NaN => 0.0)
    replace!(prices_uc, NaN => 0.0)
    replace!(prices_ed, NaN => 0.0)
    replace!(output_md, NaN => 0.0)
    replace!(output_uc, NaN => 0.0)
    replace!(output_ed, NaN => 0.0)

    profit = sum(prices_md[1, i] * output_md[1, i] + prices_uc[1, i] * (output_uc[1, i] - output_md[1, i]) + prices_ed[1, i] * (output_ed[1, i] - output_uc[1, i]) 
                - marginal_cost * output_ed[1, i] for i in 1:size(output_ed, 2))

    return profit
end



function calculate_realized_operating_profit(prices::AxisArrays.AxisArray{Float64, 2},
                                            marginal_cost::Float64,
                                            output::Array{Float64, 2},
                                            realized_hour_weight::Dict{String, Dict{Int64, Vector{Float64}}})

    replace!(prices, NaN => 0.0)
    replace!(output, NaN => 0.0)

    profit = sum(((prices[1, i] - marginal_cost) * output[1, i])
                 for i in 1:size(output, 2))


    return profit
end

"""
This function calculates the operating market profits.
"""
function calculate_realized_operating_profit(
                                            prices_md::AxisArrays.AxisArray{Float64, 2},
                                            prices_uc::AxisArrays.AxisArray{Float64, 2},
                                            prices_ed::AxisArrays.AxisArray{Float64, 2},
                                            marginal_cost::Vector{Float64},
                                            output_md::Array{Float64, 2},
                                            output_uc::Array{Float64, 2},
                                            output_ed::Array{Float64, 2},
                                            carbon_emissions::Vector{Float64},
                                            carbon_tax::Float64,
                                            realized_hour_weight::Dict{String, Dict{Int64, Vector{Float64}}},
                                            resolution::Int64)

    replace!(prices_md, NaN => 0.0)
    replace!(prices_uc, NaN => 0.0)
    replace!(prices_ed, NaN => 0.0)
    replace!(output_md, NaN => 0.0)
    replace!(output_uc, NaN => 0.0)
    replace!(output_ed, NaN => 0.0)
    replace!(marginal_cost, NaN => 0.0)

    # profit = sum((max(0.0, (prices[1, i] - marginal_cost[i]) * output[1, i] - carbon_emissions[i] * carbon_tax))
    #              for i in 1:size(output_ed, 2))

    profit = sum(prices_md[1, i] * output_md[1, i] + prices_uc[1, i] * (output_uc[1, i] - output_md[1, i]) + prices_ed[1, i] * (output_ed[1, i] - output_uc[1, i]) 
                - marginal_cost[i] * output_ed[1, i] - carbon_emissions[i] * carbon_tax for i in 1:size(output_ed, 2))


    return profit

end

"""
This function calculates the operating market profits.
"""
function calculate_realized_operating_profit(prices::Array{Float64, 2},
                                            marginal_cost::Float64,
                                            output::Array{Float64, 2},
                                            realized_hour_weight::Dict{String, Dict{Int64, Vector{Float64}}})


    replace!(prices, NaN => 0.0)
    replace!(output, NaN => 0.0)

    profit = sum((max(0.0, (prices[1, i] - marginal_cost) * output[1, i]))
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
                                   capacity_factors_md::Dict{String, Array{Float64, 2}},
                                   capacity_factors_uc::Dict{String, Array{Float64, 2}},
                                   capacity_factors_ed::Dict{String, Array{Float64, 2}},
                                   reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
                                   reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
                                   reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
                                   inertia_perc::Dict{String, Array{Float64, 2}},
                                   capacity_accepted_bids::Dict{String, Dict{String, Float64}},
                                   rec_accepted_bids::Dict{String, Float64},
                                   realized_hour_weight::Dict{String, Dict{Int64, Vector{Float64}}},
                                   iteration_year::Int64,
                                   capacity_forward_years::Int64,
                                   carbon_tax::Float64,
                                   da_resolution::Int64,
                                   rt_resolution::Int64,
                                   rt_products::Vector{String},
                                   pcm_scenario::String) where T <: Product

    return nothing, 0
end

"""
This function calculates the realized energy market profits.
"""
function calculate_realized_profit(project::Project,
                                   product::Energy,
                                   market_prices::MarketPrices,
                                   capacity_factors_md::Dict{String, Array{Float64, 2}},
                                   capacity_factors_uc::Dict{String, Array{Float64, 2}},
                                   capacity_factors_ed::Dict{String, Array{Float64, 2}},
                                   reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
                                   reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
                                   reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
                                   inertia_perc::Dict{String, Array{Float64, 2}},
                                   capacity_accepted_bids::Dict{String, Dict{String, Float64}},
                                   rec_accepted_bids::Dict{String, Float64},
                                   realized_hour_weight::Dict{String, Dict{Int64, Vector{Float64}}},
                                   iteration_year::Int64,
                                   capacity_forward_years::Int64,
                                   carbon_tax::Float64,
                                   da_resolution::Int64,
                                   rt_resolution::Int64,
                                   rt_products::Vector{String},
                                   pcm_scenario::String)

    project_name = get_name(project)
    size = get_maxcap(project)
    zone = get_zone(get_tech(project))

    update_year = iteration_year
    if in(project_name, union(keys(capacity_factors_md), keys(capacity_factors_uc), keys(capacity_factors_ed)))
        output_md = size * capacity_factors_md[project_name]
        output_uc = size * capacity_factors_uc[project_name]
        output_ed = size * capacity_factors_ed[project_name]
        emission_intensity = get_emission_intensity(project)
        heat_rate_curve = get_heat_rate_curve(get_tech(project))

        carbon_emissions = calculate_carbon_emissions(emission_intensity, heat_rate_curve, output_ed) .* (rt_resolution / 60)
        marginal_cost = calculate_marginal_cost(heat_rate_curve, output_ed, get_fuel_cost(get_tech(project))) .* (rt_resolution / 60)

        profit = calculate_realized_operating_profit(
                                                get_prices(market_prices, product)["realized-md"][zone, :, :],
                                                get_prices(market_prices, product)["realized-uc"][zone, :, :],
                                                get_prices(market_prices, product)["realized-ed"][zone, :, :],
                                                marginal_cost,
                                                output_md,
                                                output_uc,
                                                output_ed,
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
                                  capacity_factors_md::Dict{String, Array{Float64, 2}},
                                  capacity_factors_uc::Dict{String, Array{Float64, 2}},
                                  capacity_factors_ed::Dict{String, Array{Float64, 2}},
                                  reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
                                  reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
                                  reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
                                  inertia_perc::Dict{String, Array{Float64, 2}},
                                  capacity_accepted_bids::Dict{String, Dict{String, Float64}},
                                  rec_accepted_bids::Dict{String, Float64},
                                  realized_hour_weight::Dict{String, Dict{Int64, Vector{Float64}}},
                                  iteration_year::Int64,
                                  capacity_forward_years::Int64,
                                  carbon_tax::Float64,
                                  da_resolution::Int64,
                                  rt_resolution::Int64,
                                  rt_products::Vector{String},
                                  pcm_scenario::String)

    project_name = get_name(project)
    size = get_maxcap(project)

    if get_name(product) in rt_products
        scale = rt_resolution / 60
    else
        scale = da_resolution / 60
    end

    update_year = iteration_year
    if in(project_name, union(keys(reserve_perc_md), keys(reserve_perc_uc), keys(reserve_perc_ed)))
        output_md = size * reserve_perc_md[project_name][String(get_name(product))]
        output_uc = size * reserve_perc_uc[project_name][String(get_name(product))]
        output_ed = size * reserve_perc_ed[project_name][String(get_name(product))]

        profit = calculate_realized_operating_profit(
                                                get_prices(market_prices, product)["realized-md"],
                                                get_prices(market_prices, product)["realized-uc"],
                                                get_prices(market_prices, product)["realized-ed"],
                                                get_marginal_cost(product) * scale,
                                                output_md,
                                                output_uc,
                                                output_ed,
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
                                  capacity_factors_md::Dict{String, Array{Float64, 2}},
                                  capacity_factors_uc::Dict{String, Array{Float64, 2}},
                                  capacity_factors_ed::Dict{String, Array{Float64, 2}},
                                  reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
                                  reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
                                  reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
                                  inertia_perc::Dict{String, Array{Float64, 2}},
                                  capacity_accepted_bids::Dict{String, Dict{String, Float64}},
                                  rec_accepted_bids::Dict{String, Float64},
                                  realized_hour_weight::Dict{String, Dict{Int64, Vector{Float64}}},
                                  iteration_year::Int64,
                                  capacity_forward_years::Int64,
                                  carbon_tax::Float64,
                                  da_resolution::Int64,
                                  rt_resolution::Int64,
                                  rt_products::Vector{String},
                                  pcm_scenario::String)

    project_name = get_name(project)
    size = get_maxcap(project)

    update_year = iteration_year + capacity_forward_years - 1

    # Sum realized capacity revenue across seasons. `capacity_prices` is season-keyed
    # (Dict{season => AxisArray}); `capacity_accepted_bids` is Dict{season => Dict{name => fraction}}.
    # In annual mode there is a single "annual" season, reproducing the prior scalar result.
    capacity_prices = get_prices(market_prices, product)["realized"]
    profit = 0.0
    cleared_any = false
    for (season, season_price) in capacity_prices
        season_bids = capacity_accepted_bids[season]
        if in(project_name, keys(season_bids))
            cleared_any = true
            profit += size *
                      get_derating(product, pcm_scenario, season) *
                      season_price[1] *
                      season_bids[project_name]
        end
    end

    if cleared_any
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
                                  capacity_factors_md::Dict{String, Array{Float64, 2}},
                                  capacity_factors_uc::Dict{String, Array{Float64, 2}},
                                  capacity_factors_ed::Dict{String, Array{Float64, 2}},
                                  reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
                                  reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
                                  reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
                                  inertia_perc::Dict{String, Array{Float64, 2}},
                                  capacity_accepted_bids::Dict{String, Dict{String, Float64}},
                                  rec_accepted_bids::Dict{String, Float64},
                                  realized_hour_weight::Dict{String, Dict{Int64, Vector{Float64}}},
                                  iteration_year::Int64,
                                  capacity_forward_years::Int64,
                                  carbon_tax::Float64,
                                  da_resolution::Int64,
                                  rt_resolution::Int64,
                                  rt_products::Vector{String},
                                  pcm_scenario::String)

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

"""
This function calculates the realized Inertia market profits.
"""
function calculate_realized_profit(project::Project,
                                  product::Inertia,
                                  market_prices::MarketPrices,
                                  capacity_factors_md::Dict{String, Array{Float64, 2}},
                                  capacity_factors_uc::Dict{String, Array{Float64, 2}},
                                  capacity_factors_ed::Dict{String, Array{Float64, 2}},
                                  reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
                                  reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
                                  reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
                                  inertia_perc::Dict{String, Array{Float64, 2}},
                                  capacity_accepted_bids::Dict{String, Dict{String, Float64}},
                                  rec_accepted_bids::Dict{String, Float64},
                                  realized_hour_weight::Dict{String, Dict{Int64, Vector{Float64}}},
                                  iteration_year::Int64,
                                  capacity_forward_years::Int64,
                                  carbon_tax::Float64,
                                  da_resolution::Int64,
                                  rt_resolution::Int64,
                                  rt_products::Vector{String},
                                  pcm_scenario::String)

    project_name = get_name(project)
    size = get_maxcap(project)

    if get_name(product) in rt_products
        scale = rt_resolution / 60
    else
        scale = da_resolution / 60
    end

    update_year = iteration_year
    if in(project_name, keys(inertia_perc))
        output = size * inertia_perc[project_name]

        profit = calculate_realized_operating_profit(get_prices(market_prices, product)["realized"],
                                                get_marginal_cost(product) * scale,
                                                output,
                                                realized_hour_weight)


        return profit, update_year
    else
        return nothing, update_year
    end
end