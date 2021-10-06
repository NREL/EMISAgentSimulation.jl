"""
This function returns the total cost of generation.
"""
function calculate_cost(product::T,
                        marginal_cost::Float64,
                        carbon_tax::Float64,
                        emission_intensity::Float64,
                        inertia_constant::Float64,
                        Q::JuMP.VariableRef) where T <: OperatingProduct
    return marginal_cost * Q
end

"""
This function returns the total cost of generation.
"""
function calculate_cost(product::Energy,
                        marginal_cost::Float64,
                        carbon_tax::Float64,
                        emission_intensity::Float64,
                        inertia_constant::Float64,
                        Q::JuMP.VariableRef)
    return (marginal_cost + (carbon_tax * emission_intensity)) * Q
end

function calculate_cost(product::Inertia,
                        marginal_cost::Float64,
                        carbon_tax::Float64,
                        emission_intensity::Float64,
                        inertia_constant::Float64,
                        Q::JuMP.VariableRef)
    return marginal_cost * Q * inertia_constant
end

"""
This function returns the product of market price and quantity provided.
"""
function calculate_revenue(product::T, price::Float64, inertia_constant::Float64, Q::JuMP.VariableRef) where T <: OperatingProduct
    return price * Q
end

"""
This function returns the product of market price, quantity and inertia constant for Inertia constant.
"""
function calculate_revenue(product::Inertia, price::Float64, inertia_constant::Float64, Q::JuMP.VariableRef)
    return price * Q * inertia_constant
end

"""
This function returns the product of marginal cost and quantity provided.
"""
function calculate_cost(product::T, marginal_cost::Float64, inertia_constant::Float64, Q::JuMP.GenericAffExpr{Float64, JuMP.VariableRef}) where T <: OperatingProduct
    return marginal_cost * Q
end

"""
This function returns the product of marginal cost, quantity and inertia constant.
"""
function calculate_cost(product::Inertia, marginal_cost::Float64, inertia_constant::Float64, Q::JuMP.GenericAffExpr{Float64, JuMP.VariableRef})
    return marginal_cost * Q * inertia_constant
end

"""
This function returns the product of market price and quantity provided.
"""
function calculate_revenue(product::T,
                           price::Float64,
                           inertia_constant::Float64,
                           Q::JuMP.GenericAffExpr{Float64, JuMP.VariableRef}) where T <: OperatingProduct
    return price * Q
end

function calculate_revenue(product::Inertia,
                           price::Float64,
                           inertia_constant::Float64,
                           Q::JuMP.GenericAffExpr{Float64, JuMP.VariableRef})
    return price * Q * inertia_constant
end

"""
This function does nothing if the product is not either energy, reserve up or reserve down.
"""
function add_to_generationlimits!(min_gen::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                 max_gen::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                 Q::JuMP.VariableRef,
                                 product::T,
                                 synchronous_inertia::Bool) where T <: Product
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
                                 product::OperatingReserve{ReserveUpEMIS},
                                 synchronous_inertia::Bool)
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
                                 product::OperatingReserve{ReserveDownEMIS},
                                 synchronous_inertia::Bool)
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
                                 product::Energy,
                                 synchronous_inertia::Bool)
    JuMP.add_to_expression!(max_gen, Q)
    JuMP.add_to_expression!(min_gen, Q)
    return
end

"""
This function adds to the maximum generation expressions if the
product is Inertia and unit doesn't provide synchronous inertia.
Returns nothing.
"""
function add_to_generationlimits!(min_gen::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                 max_gen::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                 Q::JuMP.VariableRef,
                                 product::Inertia,
                                 synchronous_inertia::Bool)
    if !synchronous_inertia
        JuMP.add_to_expression!(max_gen, Q)
    end
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
                                    carbon_tax::Vector{Float64},
                                    start_year::Int64,
                                    end_year::Int64,
                                    rep_hour_weight::Vector{Float64},
                                    solver::JuMP.MOI.OptimizerWithAttributes) where P <: GeneratorEMIS{<: BuildPhase}

    num_hours = length(rep_hour_weight)

    all_products = get_products(gen)

    emission_intensity = 0.0

    for product in all_products
        emission_intensity = calculate_carbon_emissions(emission_intensity, product)
    end

    products = find_operating_products(all_products)

    market_prices = Dict{Symbol, Union{Array{Float64, 2}, AxisArrays.AxisArray{Float64, 2}}}()

    for prod in products
        if get_name(prod) == :Energy
            market_prices[get_name(prod)] = get_prices(all_prices, prod)[scenario_name][get_zone(get_tech(gen)), :, :]
        else
            market_prices[get_name(prod)] = get_prices(all_prices, prod)[scenario_name]
        end
    end

    max_perc = Dict(zip(get_name.(products), get_maxperc.(products, scenario_name, end_year, num_hours)))

    inertia_constant = get_inertia_constant(gen)
    synchronous_inertia = get_synchronous_inertia(gen)

    total_utilization = get_scenario_total_utilization(get_finance_data(gen))[scenario_name]

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
                         (calculate_revenue(product,market_prices[get_name(product)][y, h], inertia_constant, Q[get_name(product), y, h]) -
                          calculate_cost(product, get_marginal_cost(product), carbon_tax[y], emission_intensity, inertia_constant, Q[get_name(product), y, h]))

                JuMP.add_to_expression!(profit[i, y], hourly_profit)

                add_to_generationlimits!(min_gen[y, h],
                                        max_gen[y, h],
                                        Q[get_name(product), y, h],
                                        product,
                                        synchronous_inertia)

                JuMP.@constraint(m, Q[get_name(product), y, h] <=
                                max_perc[get_name(product)][y, h] * get_maxcap(gen))

            end

            JuMP.@constraint(m, Q[:Energy, y, h] == max_perc[:Energy][y, h] * get_maxcap(gen))

            if in(:Inertia, products) && synchronous_inertia
                JuMP.@constraint(m, Q[:Inertia, y, h] == Q[:Energy, y, h])
            end

            JuMP.@constraint(m, min_gen[y,h] >= get_mincap(gen) * 0)
            JuMP.@constraint(m, max_gen[y,h] <= get_maxcap(gen) * total_utilization[y, h])
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
                                    p_r::JuMP.VariableRef,
                                    p_inertia::JuMP.VariableRef
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
                                    p_r::JuMP.VariableRef,
                                    p_inertia::JuMP.VariableRef
                                    )
    JuMP.add_to_expression!(Q, p_out - p_in)
    return
end

"""
This function calculates the net reserve up provided by the storage device.
"""
function storage_product_constraints(product::OperatingReserve{ReserveUpEMIS},
                                    Q::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                    p_out::JuMP.VariableRef,
                                    p_in::JuMP.VariableRef,
                                    p_r::JuMP.VariableRef,
                                    p_inertia::JuMP.VariableRef
                                    )
        JuMP.add_to_expression!(Q, p_r)

    return
end

"""
This function calculates the net reserve down provided by the storage device.
"""
function storage_product_constraints(product::OperatingReserve{ReserveDownEMIS},
                                    Q::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                    p_out::JuMP.VariableRef,
                                    p_in::JuMP.VariableRef,
                                    p_r::JuMP.VariableRef,
                                    p_inertia::JuMP.VariableRef
                                    )
    JuMP.add_to_expression!(Q, p_r)
    return
end

"""
This function calculates the inertia provided by the storage device.
"""
function storage_product_constraints(product::Inertia,
                                    Q::JuMP.GenericAffExpr{Float64, JuMP.VariableRef},
                                    p_out::JuMP.VariableRef,
                                    p_in::JuMP.VariableRef,
                                    p_r::JuMP.VariableRef,
                                    p_inertia::JuMP.VariableRef
                                    )
    JuMP.add_to_expression!(Q, p_inertia)
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
                                    carbon_tax::Vector{Float64},
                                    start_year::Int64,
                                    end_year::Int64,
                                    rep_hour_weight::Vector{Float64},
                                    solver::JuMP.MOI.OptimizerWithAttributes ) where P <: StorageEMIS{<: BuildPhase}

    num_hours = length(rep_hour_weight)
    tech = get_tech(storage)

    products = find_operating_products(get_products(storage))
    energy_product = find_energy_product(products)

    capacity_factors = get_capacity_factors(energy_product[1])[scenario_name]

    market_prices = Dict{Symbol, Union{Array{Float64, 2}, AxisArrays.AxisArray{Float64, 2}}}()

    for prod in products
        if get_name(prod) == :Energy
            market_prices[get_name(prod)] = get_prices(all_prices, prod)[scenario_name][get_zone(get_tech(storage)), :, :]
        else
            market_prices[get_name(prod)] = get_prices(all_prices, prod)[scenario_name]
        end
    end

    max_gen = get_maxcap(storage)                                   # maximum generation limit
    min_gen = get_mincap(storage)                                    # minimum generation limit

    total_utilization = get_scenario_total_utilization(get_finance_data(storage))[scenario_name]

    input_power_limits = get_input_active_power_limits(tech)
    min_input = input_power_limits[:min]                            # minimum charging limit
    max_input = input_power_limits[:max]                            # maximum charging limit

    reserve_up_products = filter(p -> typeof(p) == OperatingReserve{ReserveUpEMIS}, products)
    reserve_down_products = filter(p -> typeof(p) == OperatingReserve{ReserveDownEMIS}, products)

    storage_capacity = get_storage_capacity(tech)
    min_storage = storage_capacity[:min]                            # minimum storage capacity
    max_storage = storage_capacity[:max]                            # maximum storage capacity

    efficiency = get_efficiency(tech)
    efficiency_in = efficiency[:in]                                 # charging efficiency
    efficiency_out = efficiency[:out]                               # discharging efficiency

    initstorage = get_soc(tech)

    inertia_constant = get_inertia_constant(storage)
    synchronous_inertia = get_synchronous_inertia(storage)

    m = JuMP.Model(solver)

    JuMP.@variable(m, p_e[p in start_year:end_year, t in 1:num_hours] >= 0) # Unit energy production [MW]
    JuMP.@variable(m, p_in[p in start_year:end_year, t in 1:num_hours] >= 0) # Storage charging [MW]

    JuMP.@variable(m, p_ru[rp in get_name.(reserve_up_products), p in start_year:end_year, t in 1:num_hours] >= 0) # Unit energy production [MW]
    JuMP.@variable(m, p_rd[rp in get_name.(reserve_down_products), p in start_year:end_year, t in 1:num_hours] >= 0) # Storage charging [MW]

    JuMP.@variable(m, p_inertia[p in start_year:end_year, t in 1:num_hours] >= 0) # Unit Inertia Provision [MW-s]

    JuMP.@variable(m, p_in_ru[rp in get_name.(reserve_up_products), p in start_year:end_year, t in 1:num_hours] >= 0) # reserve up provided by storage charging [MW]
    JuMP.@variable(m, p_in_rd[rp in get_name.(reserve_down_products), p in start_year:end_year, t in 1:num_hours] >= 0) # reserve down provided by storage charging [MW]

    JuMP.@variable(m, p_out_ru[rp in get_name.(reserve_up_products), p in start_year:end_year, t in 1:num_hours] >= 0) # reserve up provided by storage discharging [MW]
    JuMP.@variable(m, p_out_rd[rp in get_name.(reserve_down_products), p in start_year:end_year, t in 1:num_hours] >= 0) # reserve down provided by storage discharging [MW]

    JuMP.@variable(m, p_out_inertia[p in start_year:end_year, t in 1:num_hours] >= 0) # inertia provided by storage discharging [MW-s]
    JuMP.@variable(m, p_in_inertia[p in start_year:end_year, t in 1:num_hours] >= 0) # inertia provided by storage charging [MW-s]

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

                if product in reserve_up_products

                    storage_product_constraints(product,
                                                Q[get_name(product), p, t],
                                                p_e[p, t],
                                                p_in[p, t],
                                                p_ru[get_name(product), p, t],
                                                p_inertia[p, t])

                elseif product in reserve_down_products

                    storage_product_constraints(product,
                                                Q[get_name(product), p, t],
                                                p_e[p, t],
                                                p_in[p, t],
                                                p_rd[get_name(product), p, t],
                                                p_inertia[p, t])

                else
                    storage_product_constraints(product,
                                                Q[get_name(product), p, t],
                                                p_e[p, t],
                                                p_in[p, t],
                                                p_e[p, t],
                                                p_inertia[p, t])
                end

                hourly_profit = rep_hour_weight[t] *
                                     (calculate_revenue(product, market_prices[get_name(product)][p, t], inertia_constant, Q[get_name(product), p, t]))

                            JuMP.add_to_expression!(profit[i, p], hourly_profit)
                #=
                JuMP.@constraint(m, Q[get_name(product), y, h] <=
                                    max_perc[get_name(product)][y, h] * get_maxcap(gen))
                =#
            end

            if t == 1
                JuMP.add_to_expression!(storage_level[p, t], initstorage)
            else
                JuMP.add_to_expression!(storage_level[p, t], storage_level[p, t - 1])
            end

            # Storage technical constraints
            #JuMP.@constraint(m, p_e[p, t] <= capacity_factors[p, t] * max_gen) # Capacity factor constraint

            JuMP.@constraint(m, p_e[p,t] - sum(p_out_rd[rp, p, t]  for rp in get_name.(reserve_down_products)) >= min_gen * 0) # Minimum output dispatch
            JuMP.@constraint(m, p_e[p,t] + sum(p_out_ru[rp, p, t]  for rp in get_name.(reserve_up_products)) + p_out_inertia[p, t] == max_gen * total_utilization[p, t]) # Maximum output dispatch

            JuMP.@constraint(m, p_in[p,t] - sum(p_in_ru[rp, p, t]  for rp in get_name.(reserve_up_products)) - p_in_inertia[p, t] >=  min_input * 0) # Minimum input dispatch
            JuMP.@constraint(m, p_in[p,t] + sum(p_in_rd[rp, p, t]  for rp in get_name.(reserve_down_products)) <= max_input) # Maximum input dispatch

            for rp in get_name.(reserve_up_products)
                JuMP.@constraint(m, p_ru[rp, p, t] <= p_out_ru[rp, p, t]  + p_in_ru[rp, p, t]) # Total storage reserve up
            end

            for rp in get_name.(reserve_down_products)
                JuMP.@constraint(m, p_rd[rp, p, t] <= p_out_rd[rp, p, t] + p_in_rd[rp, p, t]) # Total storage reserve down
            end

            JuMP.@constraint(m, p_inertia[p, t] <= p_out_inertia[p, t] + p_in_inertia[p, t])

            JuMP.@constraint(m, storage_level[p, t] >= min_storage * 0) # Minimum storage level
            JuMP.@constraint(m, storage_level[p, t] <= max_storage) # Maximum storage level

        end
    end

    JuMP.@objective(m, Max, sum(profit))

    JuMP.optimize!(m)

    profit_array = AxisArrays.AxisArray(value.(profit),
                   AxisArrays.Axis{:prod}(get_name.(products)),
                   AxisArrays.Axis{:year}(start_year:end_year))


    energy_consumption = AxisArrays.AxisArray(zeros(end_year - start_year + 1), AxisArrays.Axis{:year}(start_year:end_year))

    energy_consumption = AxisArrays.AxisArray([sum(value(p_in[y, h]) * rep_hour_weight[h] for h in 1:num_hours) for y in start_year:end_year],
                        AxisArrays.Axis{:year}(start_year:end_year))


    return profit_array, energy_consumption

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
