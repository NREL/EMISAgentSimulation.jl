"""
This function updates the expected utility of a project for Risk Neutral investors.
"""
function update_project_utility!(project::P,
                               scenario_data::Vector{Scenario},
                               market_prices::MarketPrices,
                               carbon_tax_data::Vector{Float64},
                               risk_preference::RiskNeutral,
                               rep_hour_weight::Vector{Float64},
                               iteration_year::Int64,
                               yearly_horizon::Int64,
                               queue_cost::Vector{Float64},
                               capacity_forward_years::Int64,
                               solver::JuMP.MOI.OptimizerWithAttributes) where P <: Project{<: BuildPhase}

    update_project_npv!(project,
                        scenario_data,
                        market_prices,
                        carbon_tax_data,
                        rep_hour_weight,
                        iteration_year,
                        yearly_horizon,
                        queue_cost,
                        capacity_forward_years,
                        solver)

    finance_data = get_finance_data(project)

    for scenario in scenario_data
        scenario_name = get_name(scenario)
        set_scenario_utility!(finance_data, scenario_name, iteration_year, get_scenario_npv(finance_data)[scenario_name][iteration_year])
    end

    set_expected_utility!(finance_data, iteration_year, get_expected_npv(finance_data)[iteration_year])

    return

end

"""
This function updates the expected utility of a project for Risk Averse investors.
"""
function update_project_utility!(project::P,
                               scenario_data::Vector{Scenario},
                               market_prices::MarketPrices,
                               carbon_tax_data::Vector{Float64},
                               risk_preference::RiskAverse,
                               rep_hour_weight::Vector{Float64},
                               iteration_year::Int64,
                               yearly_horizon::Int64,
                               queue_cost::Vector{Float64},
                               capacity_forward_years::Int64,
                               solver::JuMP.MOI.OptimizerWithAttributes) where P <: Project{<: BuildPhase}

    update_project_npv!(project,
                        scenario_data,
                        market_prices,
                        carbon_tax_data,
                        rep_hour_weight,
                        iteration_year,
                        yearly_horizon,
                        queue_cost,
                        capacity_forward_years,
                        solver)

    expected_utility = 0.0
    finance_data = get_finance_data(project)

    for scenario in scenario_data
        scenario_name = get_name(scenario)
        scenario_npv = get_scenario_npv(finance_data)[scenario_name][iteration_year]
        scenario_utility = get_constant(risk_preference) -
                  (get_multiplier(risk_preference) * exp(-get_risk_coefficient(risk_preference) * scenario_npv))
        set_scenario_utility!(finance_data, scenario_name, iteration_year, scenario_utility)

        expected_utility += scenario_utility * get_probability(scenario)

    end

    set_expected_utility!(finance_data, iteration_year, expected_utility)

    return

end

"""
This function does nothing if a project is retired.
"""
function update_project_utility!(project::P,
                               scenario_data::Vector{Scenario},
                               market_prices::MarketPrices,
                               carbon_tax_data::Vector{Float64},
                               risk_preference::R,
                               rep_hour_weight::Vector{Float64},
                               iteration_year::Int64,
                               yearly_horizon::Int64,
                               queue_cost::Vector{Float64},
                               capacity_forward_years::Int64,
                               solver::JuMP.MOI.OptimizerWithAttributes) where {P <: Project{Retired}, R <: RiskPreference}

end
