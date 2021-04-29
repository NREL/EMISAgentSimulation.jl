"""
This function does nothing if the method is not specified for the project's buildphase.
"""
function retire_project!(project::P) where P <: Project{<: BuildPhase}
    return project
end

"""
This function converts Existing, Planned and Queued Projects to Retired buildphase.
"""
function retire_project!(project::P) where P <: Union{Project{Existing}, Project{Planned}, Project{Queue}}
    return convert(Project{Retired}, project)
end

"""
This function does nothing if PSY System is not defined.
"""
function remove_system_component!(sys::Nothing,
                                  project::P,
                                  ) where P <: Project{<: BuildPhase}
    return
end

"""
This function removes Existing projects from the list of Devices in PSY System.
"""
function remove_system_component!(sys::PSY.System,
                                   project::P,
                                  ) where P <: Project{Existing}
    component = PSY.get_components_by_name(PSY.Device, sys, get_name(project))[1]
    PSY.remove_component!(sys, component)
end

"""
This function does nothing if the method is not specified for the project's type and buildphase.
"""
function remove_renewable_gen_data!(project::P,
                                     simulation_dir::String,
                                     ) where P <: Project{<:BuildPhase}
    return
 end

"""
This function removes the PSY System timeseries data for Existing RenewableGen projects.
"""
  function remove_renewable_gen_data!(project::RenewableGenEMIS{Existing},
                                     simulation_dir::String,
                                     )
        load_n_vg_df =  read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))
        DataFrames.select!(load_n_vg_df, DataFrames.Not(get_name(project)))
        write_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data"), "load_n_vg_data.csv", load_n_vg_df)
    return
 end

 """
This function removes the future realized profits of projects which are being retired.
"""
function remove_future_profits!(project::Project, iteration_year::Int64)
    products = get_products(project)
    finance_data = get_finance_data(project)
    profit_array_length =  size(get_realized_profit(finance_data), 2)
    for year in iteration_year:profit_array_length
        for product in get_name.(products)
            set_realized_profit!(finance_data, product, year, 0.0)
        end
    end
    return
end


"""
This function does nothing if the method is not specified for the project's buildphase.
"""
function take_retirement_decision(project::P,
                                  index::Int64,
                                  projects::Vector{<: Project{<: BuildPhase}},
                                  scenario_data::Vector{Scenario},
                                  market_prices::MarketPrices,
                                  carbon_tax_data::Vector{Float64},
                                  risk_preference::R,
                                  sys_UC::Union{Nothing, PSY.System},
                                  sys_ED::Union{Nothing, PSY.System},
                                  simulation_dir::String,
                                  iteration_year::Int64,
                                  yearly_horizon::Int64,
                                  simulation_years::Int64,
                                  rep_hour_weight::Vector{Float64},
                                  capacity_forward_years::Int64,
                                  retirement_lookback::Int64,
                                  solver::JuMP.MOI.OptimizerWithAttributes) where {P <: Project{<: BuildPhase}, R <: RiskPreference}

    return
end

"""
This function decides whether or not to retire Existing projects.
"""
function take_retirement_decision(project::P,
                                  index::Int64,
                                  projects::Vector{<: Project{<: BuildPhase}},
                                  scenario_data::Vector{Scenario},
                                  market_prices::MarketPrices,
                                  carbon_tax_data::Vector{Float64},
                                  risk_preference::R,
                                  sys_UC::Union{Nothing, PSY.System},
                                  sys_ED::Union{Nothing, PSY.System},
                                  simulation_dir::String,
                                  iteration_year::Int64,
                                  yearly_horizon::Int64,
                                  simulation_years::Int64,
                                  rep_hour_weight::Vector{Float64},
                                  capacity_forward_years::Int64,
                                  retirement_lookback::Int64,
                                  solver::JuMP.MOI.OptimizerWithAttributes) where {P <: Project{Existing}, R <: RiskPreference}

    finance_data = get_finance_data(project)
    queue_cost = get_queue_cost(finance_data)

    if  (get_construction_year(project) <= simulation_years)

        update_project_utility!(project,
                              scenario_data,
                              market_prices,
                              carbon_tax_data,
                              risk_preference,
                              rep_hour_weight,
                              iteration_year,
                              yearly_horizon,
                              queue_cost,
                              capacity_forward_years,
                              solver)

        project_utility = get_expected_utility(finance_data)[iteration_year]

        annual_revenue = 0.0
        for scenario in scenario_data
            annual_revenue += get_probability(scenario) *
                             sum(get_scenario_profit(finance_data)[get_name(scenario)][iteration_year][:, iteration_year])
        end

        annual_profit = annual_revenue - get_fixed_OM_cost(finance_data)

        cash_flow = get_annual_cashflow(finance_data)
        if (iteration_year > retirement_lookback)
            previous_cashflow = cash_flow[(iteration_year - retirement_lookback):(iteration_year - 1)]
            if (annual_profit < -1) && (project_utility < 0) && (all(previous_cashflow .< -1))
                println("RETIRING:")
                println(get_name(project))
                set_retirement_year!(projects[index], iteration_year)
                projects[index] = retire_project!(project)
                remove_system_component!(sys_UC, project)
                remove_system_component!(sys_ED, project)
                remove_renewable_gen_data!(project, simulation_dir)
                remove_future_profits!(project, iteration_year)
            end
        end

    end

    return

end

"""
This function decides whether or not to retire Planned projects.
"""
function take_retirement_decision(project::P,
                                  index::Int64,
                                  projects::Vector{<: Project{<: BuildPhase}},
                                  scenario_data::Vector{Scenario},
                                  market_prices::MarketPrices,
                                  carbon_tax_data::Vector{Float64},
                                  risk_preference::R,
                                  sys_UC::Union{Nothing, PSY.System},
                                  sys_ED::Union{Nothing, PSY.System},
                                  simulation_dir::String,
                                  iteration_year::Int64,
                                  yearly_horizon::Int64,
                                  simulation_years::Int64,
                                  rep_hour_weight::Vector{Float64},
                                  capacity_forward_years::Int64,
                                  retirement_lookback::Int64,
                                  solver::JuMP.MOI.OptimizerWithAttributes) where {P <: Project{Planned}, R <: RiskPreference}

    finance_data = get_finance_data(project)
    queue_cost = get_queue_cost(finance_data)

    if  (get_construction_year(project) <= simulation_years)

        update_project_utility!(project,
                              scenario_data,
                              market_prices,
                              carbon_tax_data,
                              risk_preference,
                              rep_hour_weight,
                              iteration_year,
                              yearly_horizon,
                              queue_cost,
                              capacity_forward_years,
                              solver)
    end

    return

end

"""
This function decides whether or not to retire Queued projects.
"""
function take_retirement_decision(project::P,
                                  index::Int64,
                                  projects::Vector{<: Project{<: BuildPhase}},
                                  scenario_data::Vector{Scenario},
                                  market_prices::MarketPrices,
                                  carbon_tax_data::Vector{Float64},
                                  risk_preference::R,
                                  sys_UC::Union{Nothing, PSY.System},
                                  sys_ED::Union{Nothing, PSY.System},
                                  simulation_dir::String,
                                  iteration_year::Int64,
                                  yearly_horizon::Int64,
                                  simulation_years::Int64,
                                  rep_hour_weight::Vector{Float64},
                                  capacity_forward_years::Int64,
                                  retirement_lookback::Int64,
                                  solver::JuMP.MOI.OptimizerWithAttributes) where {P <: Project{Queue}, R <: RiskPreference}

    finance_data = get_finance_data(project)
    queue_cost = get_queue_cost(finance_data)

    if  get_construction_year(project) <= simulation_years

        update_project_utility!(project,
                              scenario_data,
                              market_prices,
                              carbon_tax_data,
                              risk_preference,
                              rep_hour_weight,
                              iteration_year,
                              yearly_horizon,
                              queue_cost,
                              capacity_forward_years,
                              solver)


        project_utility = get_expected_utility(finance_data)[iteration_year]

        if project_utility < 0
            set_retirement_year!(projects[index], iteration_year)
            projects[index] = retire_project!(project)
            remove_future_profits!(project, iteration_year)
        end
    end

    return

end

"""
This function makes the retirement decisions each year.
Returns nothing.
"""
function retire_unprofitable!(investor::Investor,
                              sys_UC::Union{Nothing, PSY.System},
                              sys_ED::Union{Nothing, PSY.System},
                              simulation_dir::String,
                              iteration_year::Int64,
                              yearly_horizon::Int64,
                              simulation_years::Int64,
                              capacity_forward_years::Int64,
                              solver::JuMP.MOI.OptimizerWithAttributes)

    scenario_data = get_scenario_data(get_forecast(investor))
    market_prices = get_market_prices(investor)
    carbon_tax_data = get_carbon_tax(investor)
    risk_preference = get_risk_preference(investor)
    projects = get_projects(investor)
        for i in 1:length(projects)
            take_retirement_decision(projects[i],
                                    i,
                                    projects,
                                    scenario_data,
                                    market_prices,
                                    carbon_tax_data,
                                    risk_preference,
                                    sys_UC,
                                    sys_ED,
                                    simulation_dir,
                                    iteration_year,
                                    yearly_horizon,
                                    simulation_years,
                                    get_rep_hour_weight(investor),
                                    capacity_forward_years,
                                    get_retirement_lookback(investor),
                                    solver)

        end
    return
end
