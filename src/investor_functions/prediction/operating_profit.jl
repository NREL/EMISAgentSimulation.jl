"""
This function returns the product of marginal cost and quantity provided.
"""
function calculate_cost(marginal_cost::Float64, Q::JuMP.VariableRef)
    return marginal_cost * Q
end

"""
This function returns the product of market price and quantity provided.
"""
function calculate_revenue(price::Float64, Q::JuMP.VariableRef)
    return price * Q
end

"""
This function returns the product of marginal cost and quantity provided.
"""
function calculate_cost(marginal_cost::Float64, Q::JuMP.GenericAffExpr{Float64, JuMP.VariableRef})
    return marginal_cost * Q
end

"""
This function returns the product of market price and quantity provided.
"""
function calculate_revenue(price::Float64, Q::JuMP.GenericAffExpr{Float64, JuMP.VariableRef})
    return price * Q
end

"""
This function does nothing if the product is not either energy, reserve up or reserve down.
"""
function add_to_generationlimits!(min_gen::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                 max_gen::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                 Q::JuMP.VariableRef,
                                 product::T) where T <: Product
return
end

"""
This function adds to the maximum generation expression if the
product is Reserve Up.
Returns nothing.
"""
function add_to_generationlimits!(min_gen::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                 max_gen::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                 Q::JuMP.VariableRef,
                                 product::OperatingReserve{ReserveUpEMIS})
    JuMP.add_to_expression!(max_gen, Q)
    return
end

"""
This function adds to the minimum generation expression if the
product is Reserve Down.
Returns nothing.
"""
function add_to_generationlimits!(min_gen::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                 max_gen::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                 Q::JuMP.VariableRef,
                                 product::OperatingReserve{ReserveDownEMIS})
    JuMP.add_to_expression!(min_gen, -Q)
    return
end

"""
This function adds to both minimum and maximum generation expressions if the
product is Energy.
Returns nothing.
"""
function add_to_generationlimits!(min_gen::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                 max_gen::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                 Q::JuMP.VariableRef,
                                 product::Energy)
    JuMP.add_to_expression!(max_gen, Q)
    JuMP.add_to_expression!(min_gen, Q)
    return
end

"""
This function calculates operating market profits for a generator
based on expected market prices stored in the investor struct.
Returns an axis array of profits for each operating market product.
"""
function calculate_operating_profit(gen::P,
                                    scenario_name::String,
                                    all_prices::MarketPrices,
                                    start_year::Int64,
                                    end_year::Int64,
                                    rep_hour_weight::Vector{Float64},
                                    solver::JuMP.MOI.OptimizerWithAttributes) where P <: GeneratorEMIS{<: BuildPhase}

    num_hours = length(rep_hour_weight)

    products = find_operating_products(get_products(gen))

    market_prices = Dict(get_name(prod) => get_prices(all_prices, prod)[scenario_name][get_zone(get_tech(gen)), :, :] for prod in products)

    max_perc = Dict(zip(get_name.(products), get_maxperc.(products, scenario_name, end_year, num_hours)))


    m = JuMP.Model(solver)

    JuMP.@variable(m, Q[i in get_name.(products),
                   y in start_year:end_year,
                   h in 1:num_hours] >= 0)

    JuMP.@variable(m, temp == 0)


    JuMP.@expression(m, profit[i in 1:length(products), y in start_year:end_year], 0 + temp)

    JuMP.@expression(m, min_gen[y in start_year:end_year, h in 1:num_hours], 0 + temp)

    JuMP.@expression(m, max_gen[y in start_year:end_year, h in 1:num_hours], 0 + temp)


    for y in start_year:end_year
        for h in 1:num_hours
            for (i, product) in enumerate(products)

                hourly_profit = rep_hour_weight[h] *
                         (calculate_revenue(market_prices[get_name(product)][y, h],
                                       Q[get_name(product), y, h]) -
                          calculate_cost(get_marginal_cost(product),
                                   Q[get_name(product), y, h]))

                JuMP.add_to_expression!(profit[i, y], hourly_profit)

                add_to_generationlimits!(min_gen[y, h],
                                        max_gen[y, h],
                                        Q[get_name(product), y, h],
                                        product)
                JuMP.@constraint(m, Q[get_name(product), y, h] <=
                                max_perc[get_name(product)][y, h] * get_maxcap(gen))
            end
            JuMP.@constraint(m, min_gen[y,h] >= get_mincap(gen) * 0)
            JuMP.@constraint(m, max_gen[y,h] <= get_maxcap(gen))
        end
    end


    JuMP.@objective(m, Max, sum(profit))

    JuMP.optimize!(m)

    profit_array = AxisArrays.AxisArray(value.(profit),
                   AxisArrays.Axis{:prod}(get_name.(products)),
                   AxisArrays.Axis{:year}(start_year:end_year))

    energy_production = AxisArrays.AxisArray([sum(value(Q[:Energy, y, h]) * rep_hour_weight[h] for h in 1:num_hours) for y in start_year:end_year],
                        AxisArrays.Axis{:year}(start_year:end_year))

    return profit_array, energy_production

end

"""
This function does nothing if the product is not of Energy, ReserveUp or ReserveDown type
"""
function storage_product_constraints(product::T,
                                    Q::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                    p_out::JuMP.VariableRef,
                                    p_in::JuMP.VariableRef,
                                    p_ru::JuMP.VariableRef,
                                    p_rd::JuMP.VariableRef
                                    ) where T <: OperatingProduct
    return
end

"""
This function calculates the net energy provided by the storage device.
"""
function storage_product_constraints(product::Energy,
                                    Q::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                    p_out::JuMP.VariableRef,
                                    p_in::JuMP.VariableRef,
                                    p_ru::JuMP.VariableRef,
                                    p_rd::JuMP.VariableRef
                                    )
    JuMP.add_to_expression!(Q, p_out - p_in)
    return
end

"""
This function calculates the net reserve up provided by the storage device.
"""
function storage_product_constraints(product::ReserveUpEMIS,
                                    Q::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                    p_out::JuMP.VariableRef,
                                    p_in::JuMP.VariableRef,
                                    p_ru::JuMP.VariableRef,
                                    p_rd::JuMP.VariableRef
                                    )
    JuMP.add_to_expression!(Q, p_ru)
    return
end

"""
This function calculates the net reserve down provided by the storage device.
"""
function storage_product_constraints(product::ReserveDownEMIS,
                                    Q::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                    p_out::JuMP.VariableRef,
                                    p_in::JuMP.VariableRef,
                                    p_ru::JuMP.VariableRef,
                                    p_rd::JuMP.VariableRef
                                    )
    JuMP.add_to_expression!(Q, p_rd)
    return
end


"""
This function calculates operating market profits for a storage device
based on expected market prices stored in the investor struct.
Returns an axis array of profits for each operating market product.
"""
function calculate_operating_profit(storage::P,
                                    scenario_name::String,
                                    all_prices::MarketPrices,
                                    start_year::Int64,
                                    end_year::Int64,
                                    rep_hour_weight::Vector{Float64},
                                    solver::JuMP.MOI.OptimizerWithAttributes ) where P <: StorageEMIS{<: BuildPhase}

    num_hours = length(rep_hour_weight)
    tech = get_tech(storage)

    products = find_operating_products(get_products(storage))
    energy_product = find_energy_product(products)

    capacity_factors = get_capacity_factors(energy_product[1])[scenario_name]

    market_prices = Dict(get_name(prod) => get_prices(all_prices, prod)[scenario_name][get_zone(tech), :, :] for prod in products)

    max_gen = get_maxcap(storage)                                   # maximum generation limit
    min_gen = get_mincap(storage)                                    # minimum generation limit

    input_power_limits = get_input_active_power_limits(tech)
    min_input = input_power_limits[:min]                            # minimum charging limit
    max_input = input_power_limits[:max]                            # maximum charging limit

    max_reserve_up = get_project_max_reserveup(storage)              # maximum reserve up limit
    max_reserve_down = get_project_max_reservedown(storage)           # maximum reserve down limit

    storage_capacity = get_storage_capacity(tech)
    min_storage = storage_capacity[:min]                            # minimum storage capacity
    max_storage = storage_capacity[:max]                            # maximum storage capacity

    efficiency = get_efficiency(tech)
    efficiency_in = efficiency[:in]                                 # charging efficiency
    efficiency_out = efficiency[:out]                               # discharging efficiency

    initstorage = get_soc(tech)
    m = JuMP.Model(solver)

    JuMP.@variable(m, p_e[p in start_year:end_year, t in 1:num_hours] >= 0) # Unit energy production [MW]
    JuMP.@variable(m, p_in[p in start_year:end_year, t in 1:num_hours] >= 0) # Storage charging [MW]

    JuMP.@variable(m, p_ru[p in start_year:end_year, t in 1:num_hours] >= 0) # Unit energy production [MW]
    JuMP.@variable(m, p_rd[p in start_year:end_year, t in 1:num_hours] >= 0) # Storage charging [MW]

    JuMP.@variable(m, p_in_ru[p in start_year:end_year, t in 1:num_hours] >= 0) # reserve up provided by storage charging [MW]
    JuMP.@variable(m, p_in_rd[p in start_year:end_year, t in 1:num_hours] >= 0) # reserve down provided by storage charging [MW]

    JuMP.@variable(m, p_out_ru[p in start_year:end_year, t in 1:num_hours] >= 0) # reserve up provided by storage discharging [MW]
    JuMP.@variable(m, p_out_rd[p in start_year:end_year, t in 1:num_hours] >= 0) # reserve down provided by storage discharging [MW]

    JuMP.@variable(m, temp == 0)

    # Quantity of products provided
    JuMP.@expression(m, Q[i in get_name.(products), p in start_year:end_year, t in 1:num_hours], 0 + temp)

    # Profit for each product
    JuMP.@expression(m, profit[i in 1:length(products), p in start_year:end_year], 0 + temp)

    # Storage level evolution
    JuMP.@expression(m, storage_level[p in start_year:end_year, t in 1:num_hours],
                    p_in[p, t] * efficiency_in - p_e[p, t] / efficiency_out)

    for p in start_year:end_year
        for t in 1:num_hours

            for (i, product) in enumerate(products)

                storage_product_constraints(product,
                                            Q[get_name(product), p, t],
                                            p_e[p, t],
                                            p_in[p, t],
                                            p_ru[p, t],
                                            p_rd[p, t])

                hourly_profit = rep_hour_weight[t] *
                                     (calculate_revenue(market_prices[get_name(product)][p, t],
                                                   Q[get_name(product), p, t]) -
                                      calculate_cost(get_marginal_cost(product),
                                               Q[get_name(product), p, t]))

                            JuMP.add_to_expression!(profit[i, p], hourly_profit)

            end

            if t == 1
                JuMP.add_to_expression!(storage_level[p, t], initstorage)
            else
                JuMP.add_to_expression!(storage_level[p, t], storage_level[p, t - 1])
            end

            # Storage technical constraints
            JuMP.@constraint(m, p_e[p, t] <= capacity_factors[p, t] * max_gen) # Capacity factor constraint

            JuMP.@constraint(m, p_e[p,t] - p_out_rd[p,t] >= min_gen * 0) # Minimum output dispatch
            JuMP.@constraint(m, p_e[p,t] + p_out_ru[p,t] <= max_gen) # Maximum output dispatch

            JuMP.@constraint(m, p_in[p,t] - p_in_ru[p,t] >=  min_input * 0) # Minimum input dispatch
            JuMP.@constraint(m, p_in[p,t] + p_in_rd[p,t] <= max_input) # Maximum input dispatch

            JuMP.@constraint(m, p_ru[p, t] <= p_out_ru[p, t] + p_in_ru[p, t]) # Total storage reserve up
            JuMP.@constraint(m, p_rd[p, t] <= p_out_rd[p, t] + p_in_rd[p, t]) # Total storage reserve down

            JuMP.@constraint(m, storage_level[p, t] >= min_storage) # Minimum storage level
            JuMP.@constraint(m, storage_level[p, t] <= max_storage) # Maximum storage level

        end
    end

    JuMP.@objective(m, Max, sum(profit))

    JuMP.optimize!(m)

    profit_array = AxisArrays.AxisArray(value.(profit),
                   AxisArrays.Axis{:prod}(get_name.(products)),
                   AxisArrays.Axis{:year}(start_year:end_year))

    energy_production = AxisArrays.AxisArray([sum(value(Q[:Energy, p, t]) * rep_hour_weight[t] for t in 1:num_hours) for p in start_year:end_year],
                        AxisArrays.Axis{:year}(start_year:end_year))

    return profit_array, energy_production

end

"""
This function updates the expected profit array of the project
if the product is an operating market product.
Returns nothing.
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
                        iteration_year::Int64) where {P <: Project{<: BuildPhase}, T <: OperatingProduct}

    profit = [operating_profit[get_name(product), y]
            for y in price_years[:start_year]:price_years[:end_year]]

    set_scenario_profit!(get_finance_data(project),
                scenario_name,
                get_name(product),
                iteration_year,
                update_years[:start_year],
                update_years[:end_year],
                profit)

    return
end
