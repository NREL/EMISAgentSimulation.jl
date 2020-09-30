
function run_agent_simulation(simulation::AgentSimulation, simulation_years::Int64)

    total_horizon = get_total_horizon(get_case(simulation))
    rolling_horizon = get_rolling_horizon(get_case(simulation))

    installed_capacity =  zeros(simulation_years)

    capacity_forward_years = get_capacity_forward_years(simulation)

    # Set initial capacity market profits considering forward capacity auctions
    if get_markets(simulation)[:Capacity]
        initial_existing_projects = vcat(get_existing.(get_investors(simulation))...)
        initial_capacity_prices = [70000.0, 80000.0]              # initial capacity prices - arbitrarily selected here - #TODO: need some meachanism to generate these
        for y in 1:capacity_forward_years - 1
            for project in initial_existing_projects
                if get_end_life_year(project) >= y
                    for product in get_products(project)
                        update_initial_capacity_revenues!(project, product, initial_capacity_prices, y)
                    end
                end
            end
        end
    end

    sys_UC = get_system_UC(simulation)
    sys_ED = get_system_ED(simulation)

    investors = get_investors(simulation)

    average_capital_cost_multiplier = Statistics.mean(get_cap_cost_multiplier.(investors))

    services = get_system_services(sys_UC)

    for iteration_year = 1:simulation_years

        yearly_horizon = min(total_horizon - iteration_year + 1, rolling_horizon)

        println("Year $(iteration_year)")
        set_iteration_year!(simulation, iteration_year)

        active_projects = deepcopy(get_activeprojects(simulation))

        installed_capacity = update_installed_cap!(installed_capacity,
                                                   active_projects,
                                                   iteration_year,
                                                   simulation_years)

        create_investor_predictions(investors,
                                          active_projects,
                                          iteration_year,
                                          yearly_horizon,
                                          get_data_dir(get_case(simulation)),
                                          get_results_dir(simulation),
                                          average_capital_cost_multiplier,
                                          get_zones(simulation),
                                          get_lines(simulation),
                                          get_peak_load(simulation),
                                          get_solver(get_case(simulation)),
                                          get_parallel_investors(get_case(simulation)),
                                          get_parallel_scenarios(get_case(simulation))
                                          )

        for investor in investors
            run_investor_iteration(investor,
                                    active_projects,
                                    iteration_year,
                                    yearly_horizon,
                                    simulation_years,
                                    get_start_year(get_case(simulation)),
                                    capacity_forward_years,
                                    get_data_dir(get_case(simulation)),
                                    sys_UC,
                                    services,
                                    get_solver(get_case(simulation))
                                    )
        end

        update_simulation_derating_data!(simulation)

        #Get all existing projects to calculate realized profits for energy and REC markets.
        all_existing_projects = vcat(get_existing.(get_investors(simulation))...)

        # Get all projects which are expected to be online for the forward capacity market auction.
        capacity_market_year = iteration_year + capacity_forward_years - 1
        capacity_market_projects = Project[]

        for project in get_activeprojects(simulation)
            end_life_year = get_end_life_year(project)
            construction_year = get_construction_year(project)
            if end_life_year >= capacity_market_year && construction_year <= capacity_market_year
                push!(capacity_market_projects, project)
            end
        end

        installed_capacity = update_installed_cap!(installed_capacity,
                                                   all_existing_projects,
                                                   iteration_year,
                                                   simulation_years)

        println("Current Installed Capacity = $(round(installed_capacity[iteration_year])) MW")

        #Find which markets to simulate.
        markets = union(hcat(get_markets.(get_investors(simulation))...))

        #Create realzed market prices for existing projects.
        realized_market_prices,
        realized_capacity_factors,
        realized_reserve_up_perc,
        realized_reserve_down_perc,
        capacity_accepted_bids,
        rec_accepted_bids = create_realized_marketdata(simulation,
                                                             sys_UC,
                                                             sys_ED,
                                                             markets,
                                                             all_existing_projects,
                                                             capacity_market_projects,
                                                             capacity_forward_years,
                                                             iteration_year,
                                                             simulation_years,
                                                             get_solver(get_case(simulation)),
                                                             get_results_dir(simulation))

        existing_project_types = unique(get_type.(get_tech.(all_existing_projects)))

        #Update forecasts and realized profits of all existing projects for each investor.

        for investor in get_investors(simulation)
            if iteration_year < simulation_years
                update_forecast!(get_forecast(investor), get_annual_growth(simulation)[:, iteration_year], iteration_year)
            end

            projects = get_projects(investor)
            for (i, project) in enumerate(projects)

                update_realized_profits!(project,
                                         realized_market_prices,
                                         realized_capacity_factors,
                                         realized_reserve_up_perc,
                                         realized_reserve_down_perc,
                                         capacity_accepted_bids,
                                         rec_accepted_bids,
                                         get_hour_weight(simulation),
                                         iteration_year,
                                         capacity_forward_years)

                update_annual_cashflow!(project, iteration_year)

                retire_old!(projects,
                            i,
                            project,
                            sys_UC,
                            get_data_dir(get_case(simulation)),
                            iteration_year)

            end
        end

        println("COMPLETED YEAR $(iteration_year)")
    end

    final_portfolio = vcat(get_existing.(get_investors(simulation))...)

    for project in final_portfolio
        extrapolate_profits!(project, simulation_years)
    end

    FileIO.save(joinpath(get_results_dir(simulation), "simulation_data.jld2"), "simulation_data", simulation)

    return
end


