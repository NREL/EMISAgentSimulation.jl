function run_investor_iteration(investor::Investor,
                                 active_projects::Vector{Project},
                                 iteration_year::Int64,
                                 yearly_horizon::Int64,
                                 simulation_years::Int64,
                                 start_year::Int64,
                                 capacity_forward_years::Int64,
                                 sys_data_dir::String,
                                 sys_UC::Union{Nothing, PSY.System},
                                 services::Vector{PSY.Service},
                                 solver::JuMP.MOI.OptimizerWithAttributes
                            )

    investor_dir = get_data_dir(investor)
    projects = get_projects(investor)

    active_projects_copy = deepcopy(active_projects)

    market_names = get_markets(investor)

    # Create empty market prices struct
    market_prices = MarketPrices()

    option_projects = get_options(investor)
    max_new_options = Dict(get_name(project) => 0 for project in option_projects)

    scenarios = get_scenario_data(get_forecast(investor))

    for scenario in scenarios
        scenario_name = get_name(scenario)
        output_file = joinpath(investor_dir, "expected_market_data", "$(scenario_name)_year_$(iteration_year).jld2")
        expected_data = FileIO.load(output_file)

        set_energy_price!(market_prices, scenario_name, expected_data["energy_price"])

        if in(:ReserveUp, market_names)
            set_reserveup_price!(market_prices, scenario_name, expected_data["reserve_up_price"])
        end

        if in(:ReserveDown, market_names)
            set_reservedown_price!(market_prices, scenario_name, expected_data["reserve_down_price"])
        end

        if in(:Capacity, market_names)
            set_capacity_price!(market_prices, scenario_name, expected_data["capacity_price"])
        end

        if in(:REC, market_names)
            set_rec_price!(market_prices, scenario_name, expected_data["rec_price"])
        end

        for project in projects
            if in(get_name(project), get_name.(active_projects_copy))
                update_capacity_factors!(project, scenario_name, expected_data["capacity_factors"])
            end
            if in(get_name(project), get_name.(active_projects_copy))
                update_capacity_accepted_perc!(project, scenario_name, expected_data["capacity_accepted_perc"])
            end

        end

         max_new_options = update_max_new_options!(max_new_options, expected_data["new_options"], option_projects)
    end

    set_market_prices!(investor, market_prices)

    retire_unprofitable!(investor,
                         sys_UC,
                         sys_data_dir,
                         iteration_year,
                         yearly_horizon,
                         simulation_years,
                         capacity_forward_years,
                         solver)

    make_investments!(investor,
                      max_new_options,
                      iteration_year,
                      yearly_horizon,
                      simulation_years,
                      capacity_forward_years,
                      solver)

    for (i, project) in enumerate(projects)
        start_construction!(projects,
                            i,
                            project,
                            iteration_year)

        finish_construction!(projects,
                            i,
                            project,
                            sys_UC,
                            services,
                            sys_data_dir,
                            iteration_year,
                            start_year)

        update_lifecycle!(project,
                          iteration_year,
                          simulation_years)
    end

    return
end


