"""
This function does nothing if product is not of Capacity type.
"""
function update_forward_profit!(product::T,
                               finance_data::Finance,
                               scenario_name::String,
                               update_year::Int64,
                               profit::Float64) where T <: Product
    return
end

"""
This function updates the capacity market forward reveneues.
"""
function update_forward_profit!(product::Capacity,
                               finance_data::Finance,
                               scenario_name::String,
                               update_year::Int64,
                               profit::Float64)

    for i = 1:length(get_scenario_profit(finance_data)[scenario_name])
        set_scenario_profit!(finance_data,
                        scenario_name,
                        get_name(product),
                        i,
                        update_year,
                        profit)
    end
    return
end


"""
This function does nothing if the project is retired.
"""
#=
function update_realized_profits!(project::P,
                                  market_prices::MarketPrices,
                                  capacity_factors::Dict{String, Array{Float64, 2}},
                                  reserve_perc::Dict{String, Dict{String, Array{Float64, 2}}},
                                  inertia_perc::Dict{String, Array{Float64, 2}},
                                  capacity_accepted_bids::Dict{String, Float64},
                                  rec_accepted_bids::Dict{String, Float64},
                                  rep_hour_weight::Vector{Float64},
                                  iteration_year::Int64,
                                  capacity_forward_years::Int64,
                                  carbon_tax::Float64,
                                  da_resolution::Int64,
                                  rt_resolution::Int64,
                                  rt_products::Vector{String}) where P <: Project{Retired}
    
end
=#
"""
This function updates the annual realized profit for active projects.
Returns nothing.
"""
function update_realized_profits!(project::P,
                                  market_prices::MarketPrices,
                                  capacity_factors::Dict{String, Array{Float64, 2}},
                                  reserve_perc::Dict{String, Dict{String, Array{Float64, 2}}},
                                  inertia_perc::Dict{String, Array{Float64, 2}},
                                  capacity_accepted_bids::Dict{String, Float64},
                                  rec_accepted_bids::Dict{String, Float64},
                                  realized_hour_weight::Vector{Float64},
                                  iteration_year::Int64,
                                  capacity_forward_years::Int64,
                                  carbon_tax::Float64,
                                  da_resolution::Int64,
                                  rt_resolution::Int64,
                                  rt_products::Vector{String}) where P <: Project{<: BuildPhase}

    for product in get_products(project)
        profit, update_year = calculate_realized_profit(project,
                                           product,
                                           market_prices,
                                           capacity_factors,
                                           reserve_perc,
                                           inertia_perc,
                                           capacity_accepted_bids,
                                           rec_accepted_bids,
                                           realized_hour_weight,
                                           iteration_year,
                                           capacity_forward_years,
                                           carbon_tax,
                                           da_resolution,
                                           rt_resolution,
                                           rt_products)

        finance_data =  get_finance_data(project)
        profit_array_length =  size(get_realized_profit(finance_data), 2)
        if !isnothing(profit) && update_year <= profit_array_length && update_year > 0
            
            set_realized_profit!(finance_data,
                    get_name(product),
                    update_year,
                    profit)

            for scenario_name in keys(get_scenario_profit(finance_data))
                update_forward_profit!(product, finance_data, scenario_name, update_year, profit)
            end
        end

    end

    return
end

"""
THis function updates the annual cash flow of Existing projects.
"""
function update_annual_cashflow!(project::Union{Project{Retired}, Project{Existing}}, iteration_year::Int64)

    finance_data = get_finance_data(project)
    annual_revenue = sum(get_realized_profit(finance_data)[:, iteration_year])
    if isnan(annual_revenue)
       println(get_name(project)) 
    end
    annual_cashflow = annual_revenue - get_fixed_OM_cost(finance_data)
    set_annual_cashflow!(finance_data, iteration_year, get_annual_cashflow(finance_data)[iteration_year] + annual_cashflow)
    return
end

"""
This function updates the annual cash flow of Queue projects
"""
function update_annual_cashflow!(project::Project{Queue}, iteration_year::Int64)
    finance_data = get_finance_data(project)
    queue_cost = get_queue_cost(finance_data)
    decision_year = get_decision_year(project)

    queue_year =  iteration_year - decision_year + 1

    annual_queuecost = queue_cost[queue_year]
    set_annual_cashflow!(finance_data, iteration_year, get_annual_cashflow(finance_data)[iteration_year] - annual_queuecost)
    return
end

"""
This function updates the annual cash flow of Planned projects
"""
function update_annual_cashflow!(project::Project{Planned}, iteration_year::Int64)
    finance_data = get_finance_data(project)

    queue_cost = get_queue_cost(finance_data)
    decision_year = get_decision_year(project)

    if decision_year + length(queue_cost) == iteration_year
        investment_cost = get_effective_investment_cost(finance_data)
        set_annual_cashflow!(finance_data, iteration_year, get_annual_cashflow(finance_data)[iteration_year] - investment_cost)
    end

    return
end

"""
This function does nothing if project is in Option or Retired phase.
"""
function update_annual_cashflow!(project:: Project{Option}, iteration_year::Int64)
    return
end

"""
This function keeps forecasts the same for all scenarios if Kalman Filter based updates are deactivated.
"""
function update_scenario_data!(scenario_data::Vector{Scenario}, kf::Nothing, iteration_year::Int64)
    for scenario in scenario_data
        parameter_values = deepcopy(get_parameter_values(scenario)[iteration_year])
        set_parameter_values!(scenario, iteration_year + 1, parameter_values)
    end

    return
end

"""
This function updates the forecasts for all scenarios using Kalman Filters and scaling factors.
"""
function update_scenario_data!(scenario_data::Vector{Scenario}, kf::KalmanFilter, iteration_year::Int64)
    state_estimate_data = get_state_estimate(kf)
    for scenario in scenario_data
        parameter_values = deepcopy(get_parameter_values(scenario)[iteration_year])

        axis_val = AxisArrays.axisvalues(parameter_values)
        parameter_multipliers = get_parameter_multipliers(scenario)
        for parameter in axis_val[1]
            for year in axis_val[2]
                parameter_values[parameter, year] = state_estimate_data[parameter] * parameter_multipliers[parameter]
            end
        end
        set_parameter_values!(scenario, iteration_year + 1, parameter_values)
    end

    return
end

"""
This function does nothing if forecast is Perfect.
"""
function update_forecast!(forecast::Perfect, measurement::AxisArrays.AxisArray{Float64, 1}, iteration_year::Int64)
    return
end

"""
This function updates the forecasts using Kalman Filters if they are Imperfect.
"""
function update_forecast!(forecast::Imperfect, measurement::AxisArrays.AxisArray{Float64, 1}, iteration_year::Int64)
    update_belief!(get_kalman_filter(forecast), measurement)
    update_scenario_data!(get_scenario_data(forecast), get_kalman_filter(forecast), iteration_year)
    return
end
