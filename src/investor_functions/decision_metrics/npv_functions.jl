"""
This function calculates the annualized capital, queue and fixed O&M costs
"""
function calculate_annualized_costs(project::P,
                                    sz::Int64,
                                    queue_cost::Vector{Float64}) where P <: Project{<: BuildPhase}

    finance = get_finance_data(project)
    r = get_discount_rate(finance)

    adjusted_investment_cost = get_effective_investment_cost(finance) * (1 + r) ^ (get_lag_time(finance) - 1)

    adjusted_queue_cost = 0

    for q in 1:length(queue_cost)
        adjusted_queue_cost += queue_cost[q] * (1 + r) ^ (get_lag_time(finance) + (length(queue_cost) - q))
    end

    α =  r / (1 - (1 + r) ^ (-get_capex_years(finance)))

    annualized_cap_cost = fill(α * (adjusted_investment_cost + adjusted_queue_cost), sz)

    annual_OM_cost = fill(get_fixed_OM_cost(finance), sz)

    total_annual_cost = annualized_cap_cost + annual_OM_cost

    return total_annual_cost

end

"""
This function returns the net present value of the existing and planned projects:
NPV for these buildphases doesn't include investment and queue costs.
"""
function calculate_npv(project::P,
                       scenario_name::String,
                       iteration_year::Int64,
                       yearly_horizon::Int64,
                       queue_time::Int64,
                       queue_cost::Vector{Float64},
                       update_start_year::Int64) where P <: Union{Project{Existing}, Project{Planned}}

    project_life_end = get_end_life_year(project)
    finance = get_finance_data(project)
    npv_years = project_life_end - iteration_year + 1
    discount_rate = get_discount_rate(finance)

    OM_cost_array = create_OMcost_array(finance, scenario_name, project_life_end, update_start_year, iteration_year)

    npv_cost_array = [OM_cost_array[y]
                     for y in iteration_year:project_life_end]

    npv_profit_array = [sum(get_scenario_profit(finance)[scenario_name][iteration_year][:, y])
                       for y in iteration_year:project_life_end]



    npv = sum(((1 / (1 + discount_rate)) ^ y) * (npv_profit_array[y] - npv_cost_array[y])
              for y in 1:npv_years)

    return npv
end


"""
This function returns the net present value of the queue and option projects:
NPV for these buildphases considers annualized capital and queue costs.
"""
function calculate_npv(project::P,
                       scenario_name::String,
                       iteration_year::Int64,
                       yearly_horizon::Int64,
                       queue_time::Int64,
                       queue_cost::Vector{Float64},
                       update_start_year::Int64) where P <: Union{Project{Option}, Project{Queue}}

    project_life_end = get_end_life_year(project)
    finance = get_finance_data(project)
    npv_years = min(project_life_end, yearly_horizon + iteration_year - 1)
    discount_rate = get_discount_rate(finance)

    queue_cost_array = create_queue_cost_array(size(get_scenario_profit(finance)[scenario_name][iteration_year], 2),
                                              get_decision_year(project),
                                              queue_time,
                                              queue_cost)

    npv_profit_array = [sum(get_scenario_profit(finance)[scenario_name][iteration_year][:, y])
                        for y in get_construction_year(project):npv_years]

    npv_cost_array = calculate_annualized_costs(project, length(npv_profit_array), queue_cost)

    npv = sum(((1 / (1 + discount_rate)) ^ (y + get_construction_year(project) - iteration_year )) * (npv_profit_array[y] - npv_cost_array[y])
              for y in 1:length(npv_profit_array))

    return npv
end

"""
This function updates the NPV of the project.
"""
function update_project_npv!(project::P,
                               scenario_data::Vector{Scenario},
                               market_prices::MarketPrices,
                               carbon_tax_data::Vector{Float64},
                               rep_hour_weight::Vector{Float64},
                               iteration_year::Int64,
                               yearly_horizon::Int64,
                               queue_cost::Vector{Float64},
                               capacity_forward_years::Int64,
                               solver::JuMP.MOI.OptimizerWithAttributes) where P <: Project{<: BuildPhase}

    price_years = find_price_years(get_construction_year(project),
                                   get_end_life_year(project),
                                   iteration_year,
                                   yearly_horizon)

    update_years = find_update_years(get_construction_year(project),
                                     get_end_life_year(project),
                                     iteration_year,
                                     yearly_horizon)

    update_expected_profit!(project,
                            scenario_data,
                            market_prices,
                            carbon_tax_data,
                            price_years,
                            update_years,
                            rep_hour_weight,
                            queue_cost,
                            capacity_forward_years,
                            iteration_year,
                            solver)

    expected_npv = 0.0

    finance_data = get_finance_data(project)

    for scenario in scenario_data
        scenario_name = get_name(scenario)
        scenario_npv = calculate_npv(project,
                            scenario_name,
                            iteration_year,
                            yearly_horizon,
                            length(queue_cost),
                            queue_cost,
                            update_years[:start_year])

        if -5 <= scenario_npv / get_maxcap(project) < 0   # Converting values from very small negative values to 0 to avoid rounding off issues
            scenario_npv = 0.0
        end
        set_scenario_npv!(finance_data, scenario_name, iteration_year, scenario_npv)
        expected_npv += scenario_npv * get_probability(scenario)
    end

    set_expected_npv!(finance_data, iteration_year, expected_npv)

    return
end

"""
This function does nothing if the project is retired.
Returns nothing.
"""
function update_project_npv!(project::P,
                               scenario_data::Vector{Scenario},
                               market_prices::MarketPrices,
                               carbon_tax_data::Vector{Float64},
                               rep_hour_weight::Vector{Float64},
                               iteration_year::Int64,
                               yearly_horizon::Int64,
                               queue_cost::Vector{Float64},
                               capacity_forward_years::Int64,
                               solver::JuMP.MOI.OptimizerWithAttributes ) where P <: Project{Retired}
    return
end
