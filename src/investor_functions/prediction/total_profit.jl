"""
This function updates the expected profit of the project from all products.
Returns nothing.
"""
function update_expected_profit!(project::P,
                                 scenario_data::Vector{Scenario},
                                 market_prices::MarketPrices,
                                 carbon_tax_data::Vector{Float64},
                                 price_years::NamedTuple{(:start_year, :end_year),
                                              Tuple{Int64, Int64}},
                                 update_years::NamedTuple{(:start_year, :end_year),
                                              Tuple{Int64, Int64}},
                                 rep_hour_weight::Vector{Float64},
                                 queue_cost::Vector{Float64},
                                 capacity_forward_years::Int64,
                                 iteration_year::Int64,
                                 solver::JuMP.MOI.OptimizerWithAttributes) where P <: Project{<: BuildPhase}

    extend_profit_array(project, get_name.(scenario_data), iteration_year)

    carbon_tax_vector = carbon_tax_data[iteration_year:end]

    profit_array_length = price_years[:end_year] - price_years[:start_year] + 1

    operating_profit_array = AxisArrays.AxisArray{Float64, 2}[] # Operating profit array for each scenario
    capacity_revenue_array = Vector{Float64}[]  # Capacity revenue array for each scenario
    rec_revenue_array = Vector{Float64}[]  # REC revenue array for each scenario

    expected_energy_production = zeros(profit_array_length)

    capacity_prices = get_capacity_price(market_prices)
    rec_prices = get_rec_price(market_prices)

    sz = price_years[:end_year] - price_years[:start_year] + 1

    total_annual_cost = calculate_annualized_costs(project,
                                                   sz,
                                                   queue_cost)

    for scenario in scenario_data
        operating_profit,
        energy_production = calculate_operating_profit(project,
                                                      get_name(scenario),
                                                      market_prices,
                                                      carbon_tax_vector,
                                                      price_years[:start_year],
                                                      price_years[:end_year],
                                                      rep_hour_weight,
                                                      solver)
        push!(operating_profit_array, operating_profit)

        # Calculate expected energy production (NOTE: For storage, we save expected energy consumption (instead of production) here, for adding to clean energy (REC) requirement constraint)
        expected_energy_production += get_probability(scenario) * energy_production

        capacity_revenues = zeros(profit_array_length)
        rec_revenues = zeros(profit_array_length)


        for product in get_products(project)
            if !isnothing(capacity_prices)
                capacity_revenues = update_capacity_revenues!(
                                                              product,
                                                              get_name(scenario),
                                                              price_years,
                                                              capacity_revenues,
                                                              capacity_prices,
                                                              get_maxcap(project)
                                                            )
            end

            if !isnothing(rec_prices)
                rec_revenues = update_rec_revenues!(
                                                    product,
                                                    price_years,
                                                    rec_revenues,
                                                    rec_prices[get_name(scenario)],
                                                    energy_production
                                                    )
            end
        end

        annual_operating_profit =  [sum(operating_profit[:, y])
                               for y in price_years[:start_year]:price_years[:end_year]]

        capacity_going_forward_cost = total_annual_cost - annual_operating_profit - rec_revenues

        for product in get_products(project)
            update_profit!(project,
                           get_name(scenario),
                           product,
                           operating_profit,
                           capacity_revenues,
                           capacity_forward_years,
                           rec_revenues,
                           price_years,
                           update_years,
                           iteration_year)
        end

        push!(capacity_revenue_array, capacity_revenues)
        push!(rec_revenue_array, rec_revenues)
    end

    profit_axis_values = AxisArrays.axisvalues(operating_profit_array[1])

    expected_operating_profit = operating_profit_array[1]

    for x in profit_axis_values[1]
        for y in profit_axis_values[2]
            expected_operating_profit[x, y] = sum(get_probability(scenario_data[i]) * operating_profit_array[i][x, y] for i in 1:length(operating_profit_array))
        end
    end

    # Calculate expected profits
    expected_annual_operating_profit = [sum(expected_operating_profit[x, y] for x in profit_axis_values[1]) for y in profit_axis_values[2]]

    expected_capacity_revenues = [sum(get_probability(scenario_data[i]) * capacity_revenue_array[i][y] for i in 1:length(capacity_revenue_array))
                                 for y in 1:profit_array_length]

    expected_rec_revenues = [sum(get_probability(scenario_data[i]) * rec_revenue_array[i][y] for i in 1:length(rec_revenue_array))
                                 for y in 1:profit_array_length]

    #println(get_name(project))
    #println("expected_annual_operating_profit is $(expected_annual_operating_profit) & expected_capacity_revenues is $(expected_capacity_revenues) & expected_rec_revenues is $(expected_rec_revenues)")

    #Calculate capacity market going forward cost.
    expected_capacity_going_forward_cost = total_annual_cost - expected_annual_operating_profit - expected_rec_revenues
    capacity_market_bid = 0.0

    capacity_revenue_start_year = max(1, capacity_forward_years - price_years[:start_year] + 1)

    if capacity_revenue_start_year <= length(expected_capacity_going_forward_cost)
        capacity_market_bid = round(max(0.0, expected_capacity_going_forward_cost[capacity_revenue_start_year] / get_maxcap(project)), digits = 4)
    end

    #Calculate REC market bid.
    if length(expected_annual_operating_profit) >= 1
        if expected_energy_production[1] <= 1e-3
            rec_market_bid = 0.0
        else
            rec_market_bid = max(0.0, (total_annual_cost[1] - expected_annual_operating_profit[1] - expected_capacity_revenues[1]) / max(1e-3, expected_energy_production[1]))
        end
        for product in get_products(project)
            update_bid!(product,
                        capacity_market_bid,
                        rec_market_bid,
                        expected_energy_production[1],
                        iteration_year
                        )
            set_expected_production!(product, expected_energy_production[1])
        end
    else
        rec_market_bid = 0.0
        for product in get_products(project)
            update_bid!(product,
                        capacity_market_bid,
                        rec_market_bid,
                        0.0,
                        iteration_year
                        )
        end
    end

    # Sets the project profit values beyond the horizon equal to the last value.
    if get_end_life_year(project) > update_years[:end_year]
        for scenario_name in get_name.(scenario_data)
            for product_name in get_name.(get_products(project))
                for y in (update_years[:end_year] + 1):get_end_life_year(project)
                    profit = get_scenario_profit(get_finance_data(project))[scenario_name][iteration_year][Symbol(product_name),
                                                                                               update_years[:end_year]]
                    set_scenario_profit!(get_finance_data(project), scenario_name, product_name, iteration_year, y, profit)
                end
            end
        end
    end

    return
end

"""
This function does nothing if the product is neither
an operating market or capacity or REC product.
The update_profit functions for each product is defined in files pertaning to those products.
"""
function update_profit!(project::P,
                        scenario_name::String,
                        product::T,
                        operating_profit::AxisArrays.AxisArray{Float64, 2},
                        capacity_revenues::Vector{Float64},
                        capacity_forward_years::Int64,
                        rec_revenues::Vector{Float64},
                        price_years::NamedTuple{(:start_year, :end_year),
                                    Tuple{Int64, Int64}},
                        update_years::NamedTuple{(:start_year, :end_year),
                                     Tuple{Int64, Int64}},
                        iteration_year::Int64) where {P <: Project{<: BuildPhase}, T <: Product}
    return
end

"""
This function does nothing if project is not an option type.
Returns nothing
"""
function extend_profit_array(project::P, scenario_names::Vector{String}, iteration_year::Int64) where P <: Project{<: BuildPhase}
    return
end

"""
This function extends the profit array of option projects
if their end life year exceeds existing length of the profit array.
"""
function extend_profit_array(project::P, scenario_names::Vector{String}, iteration_year::Int64) where P <: Project{Option}

    if get_end_life_year(project) > size(get_scenario_profit(get_finance_data(project))[scenario_names[1]][iteration_year], 2)
        for scenario_name in scenario_names
            profit = AxisArrays.AxisArray(zeros(length(get_products(project)), get_end_life_year(project)),
                                AxisArrays.Axis{:prod}(get_name.(get_products(project))),
                                AxisArrays.Axis{:year}(1:1:get_end_life_year(project)))
            set_scenario_profit!(get_finance_data(project), scenario_name, iteration_year, profit)
        end
    end
    return
end
