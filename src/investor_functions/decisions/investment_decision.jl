"""
This function is a wrapper for running expected utility and first year profit calculations
for making investment decisions.
"""
function calculate_expected_utility(project::P,
                               scenario_data::Vector{Scenario},
                               market_prices::MarketPrices,
                               risk_preference::R,
                               rep_hour_weight::Vector{Float64},
                               iteration_year::Int64,
                               yearly_horizon::Int64,
                               queue_cost::Vector{Float64},
                               capacity_forward_years::Int64,
                               solver::JuMP.MOI.OptimizerWithAttributes) where {P <: Project{<: BuildPhase}, R <: RiskPreference}

    finance_data = get_finance_data(project)

    update_project_utility!(project,
                              scenario_data,
                              market_prices,
                              risk_preference,
                              rep_hour_weight,
                              iteration_year,
                              yearly_horizon,
                              queue_cost,
                              capacity_forward_years,
                              solver)

        annual_revenue = 0.0

        for scenario in scenario_data
            annual_revenue += get_probability(scenario) *
                             sum(get_scenario_profit(finance_data)[get_name(scenario)][iteration_year][:, get_construction_year(project)])
        end

        annual_profit = annual_revenue - get_fixed_OM_cost(finance_data)
        project_utility = get_expected_utility(finance_data)[iteration_year]

        return annual_profit, project_utility

end

"""
This function adds options of each technology type (also incorporating capital cost multiplier)
to the list of profitable options which can be sent to the queue.
Returns the updated list of profitable options.
"""
function add_profitable_option(projects::Vector{Project},
                                  max_new_options::Dict{String, Int64},
                                  profitable_options::Vector{Project},
                                  scenario_data::Vector{Scenario},
                                  market_prices::MarketPrices,
                                  risk_preference::R,
                                  capital_cost_multiplier::Float64,
                                  investor_name::String,
                                  iteration_year::Int64,
                                  yearly_horizon::Int64,
                                  simulation_years::Int64,
                                  rep_hour_weight::Vector{Float64},
                                  capacity_forward_years::Int64,
                                  solver::JuMP.MOI.OptimizerWithAttributes) where {R <: RiskPreference}

    if length(projects) >= 1
        if  get_construction_year(projects[1]) <= simulation_years

            total_counter = 1
            counter_by_zone = AxisArrays.AxisArray(ones(length(projects)), get_zone.(get_tech.(projects)))
            profitable_type_options = Project[]
            base_inv_cost = get_effective_investment_cost(get_finance_data(projects[1]))
            inv_cost_array_length = length(get_investment_cost(get_finance_data(projects[1])))
            queue_cost = get_queue_cost(get_finance_data(projects[1]))

            for project in projects

                annual_profit, project_utility = calculate_expected_utility(project,
                                                                    scenario_data,
                                                                    market_prices,
                                                                    risk_preference,
                                                                    rep_hour_weight,
                                                                    iteration_year,
                                                                    yearly_horizon,
                                                                    queue_cost,
                                                                    capacity_forward_years,
                                                                    solver)

                if project_utility >= 0 && annual_profit >= -1 && counter_by_zone[get_zone(get_tech(project))] <= max_new_options[get_name(project)]
                    new_option = deepcopy(project)
                    push!(profitable_type_options, new_option)
                end

            end

            while length(profitable_type_options) > 0
                sort!(profitable_type_options, by = x -> -get_expected_utility(get_finance_data(x))[iteration_year])
                most_profitable_option = deepcopy(profitable_type_options[1])
                zone = get_zone(get_tech(most_profitable_option))
                set_name!(most_profitable_option, "$(investor_name)_$(get_type(get_tech(most_profitable_option)))_$(zone)_year_$(iteration_year)_$(Int(counter_by_zone[zone]))")
                set_investment_cost!(most_profitable_option, fill(get_effective_investment_cost(get_finance_data(most_profitable_option)), inv_cost_array_length))
                push!(profitable_options, most_profitable_option)
                counter_by_zone[zone] += 1

                for option in profitable_type_options
                    set_effective_investment_cost!(option, base_inv_cost * (1 + capital_cost_multiplier * total_counter))
                end

                for (idx, option) in enumerate(profitable_type_options)
                    annual_profit, project_utility = calculate_expected_utility(option,
                                                                        scenario_data,
                                                                        market_prices,
                                                                        risk_preference,
                                                                        rep_hour_weight,
                                                                        iteration_year,
                                                                        yearly_horizon,
                                                                        queue_cost,
                                                                        capacity_forward_years,
                                                                        solver)
                    if project_utility < 0 || annual_profit < -1 || counter_by_zone[get_zone(get_tech(option))] > max_new_options[get_name(option)]
                        deleteat!(profitable_type_options, idx)
                    end
                end

                total_counter += 1

            end

        end
    end

    return profitable_options
end
"""
This function makes the investment decisions each year by sending a top N profitable options to the queue,
where N is the maximum number of new projects an investor can invest in within a year.
Returns nothing.
"""
function make_investments!(investor::Investor,
                           max_new_options::Dict{String, Int64},
                           iteration_year::Int64,
                           yearly_horizon::Int64,
                           simulation_years::Int64,
                           capacity_forward_years::Int64,
                           solver::JuMP.MOI.OptimizerWithAttributes)

    scenario_data = get_scenario_data(get_forecast(investor))
    market_prices = get_market_prices(investor)
    risk_preference = get_risk_preference(investor)
    projects = get_projects(investor)

    option_projects = get_options(investor)

    types = unique(get_type.(get_tech.(option_projects)))
    projects_by_type = Dict(type => Project[] for type in types)

    for type in types
        for project in option_projects
            if get_type(get_tech(project)) == type
                push!(projects_by_type[type], project)
            end
        end
    end

    profitable_options = Project[]
    for type in types
        profitable_options = add_profitable_option(projects_by_type[type],
                                                    max_new_options,
                                                    profitable_options,
                                                    scenario_data,
                                                    market_prices,
                                                    risk_preference,
                                                    get_cap_cost_multiplier(investor),
                                                    get_name(investor),
                                                    iteration_year,
                                                    yearly_horizon,
                                                    simulation_years,
                                                    get_rep_hour_weight(investor),
                                                    capacity_forward_years,
                                                    solver)
    end

        sort!(profitable_options, by = x -> -get_expected_utility(get_finance_data(x))[iteration_year])  # sort in descending order of project utility
        for i in 1:min(length(profitable_options), get_max_annual_projects(investor))
            push!(projects, convert(Project{Queue}, profitable_options[i]))
        end

    return
end
