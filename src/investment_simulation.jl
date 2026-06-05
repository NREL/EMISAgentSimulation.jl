
function run_agent_simulation(simulation::AgentSimulation, simulation_years::Int64, current_siip_sim, siip_system, current_year::Int64)
    case = get_case(simulation)
    total_horizon = get_total_horizon(case)
    rolling_horizon = get_rolling_horizon(case)

    installed_capacity = zeros(simulation_years)
    capacity_forward_years = get_capacity_forward_years(simulation)

    # Set initial capacity market profits considering forward capacity auctions
    if current_year == 1
        if get_markets(simulation)[:Capacity]
            initial_existing_projects = vcat(get_existing.(get_investors(simulation))...)

            simulation_dir = get_data_dir(get_case(simulation))
            capacity_mkt_param_file = joinpath(simulation_dir, "markets_data", "Capacity.csv")
            capacity_mkt_params = read_data(capacity_mkt_param_file)[1, :]
            introduction_year = capacity_mkt_params["introduction_year"]

            if introduction_year >= 3
                initial_capacity_prices = [0.0, 0.0]              # initial capacity prices - arbitrarily selected here - #TODO: need some meachanism to generate these
            else
                initial_capacity_prices = [70000.0, 80000.0]
            end

            for y in 1:capacity_forward_years - 1
                for project in initial_existing_projects
                    if get_end_life_year(project) >= y
                        for product in get_products(project)
                            update_initial_capacity_revenues!(project, product, initial_capacity_prices, y, get_pcm_scenario(case))
                        end
                    end
                end
            end
        end
    end

    sys_MDs = get_system_MDs(simulation)
    sys_UCs = get_system_UCs(simulation)
    sys_EDs = get_system_EDs(simulation)
    sys_PRAS = get_system_PRAS(simulation)

    investors = get_investors(simulation)

    average_capital_cost_multiplier = Statistics.mean(get_cap_cost_multiplier.(investors))

    clean_energy_percentage_vector = zeros(simulation_years)

    for iteration_year = current_year:simulation_years

        yearly_horizon = min(total_horizon - iteration_year + 1, rolling_horizon)

        println("Year $(iteration_year)")
        set_iteration_year!(simulation, iteration_year)

        active_projects = deepcopy(get_activeprojects(simulation))

        installed_capacity = update_installed_cap!(installed_capacity,
                                                   active_projects,
                                                   iteration_year,
                                                   simulation_years)

        scenario_names = String.(get_all_scenario_names(get_data_dir(case)))

        simulation_dir = get_data_dir(get_case(simulation))

        # save existing net load csv file for potential checkpoint re-runs
        for scenario in scenario_names
            pre_update_da_net_load = joinpath(simulation_dir, "timeseries_data_files", scenario, "sim_year_$(iteration_year)", "Net Load Data", "load_n_vg_data_pre_update.csv")
            pre_update_rt_net_load = joinpath(simulation_dir, "timeseries_data_files", scenario, "sim_year_$(iteration_year)", "Net Load Data", "load_n_vg_data_rt_pre_update.csv")
            post_update_da_net_load = joinpath(simulation_dir, "timeseries_data_files", scenario, "sim_year_$(iteration_year)", "Net Load Data", "load_n_vg_data.csv")
            post_update_rt_net_load = joinpath(simulation_dir, "timeseries_data_files", scenario, "sim_year_$(iteration_year)", "Net Load Data", "load_n_vg_data_rt.csv")
            if isfile(pre_update_da_net_load)
                cp(pre_update_da_net_load, post_update_da_net_load, force=true)
                cp(pre_update_rt_net_load, post_update_rt_net_load, force=true)
            else
                cp(post_update_da_net_load, pre_update_da_net_load, force=true)
                cp(post_update_rt_net_load, pre_update_rt_net_load, force=true)
            end
        end

        if current_year == 1
            for scenario in scenario_names
                derating_factors = read_data(joinpath(get_data_dir(case), "markets_data", "derating_data", scenario, "derating_dict.csv"))

                output_file = joinpath(get_results_dir(simulation), "derating_data", scenario, "derating_data_year_$(iteration_year).h5")

                save_derating_factors(output_file, derating_factors)
            end
        end

        num_scenarios = length(scenario_names)
        sys_PRAS_list, active_projects_list, capacity_forward_years_list, resource_adequacies, peak_loads, static_capacity_bools, iteration_years, simulation_years_list, data_dirs, rt_resolutions, results_dirs, outage_dirs = repeat_arguments(num_scenarios, sys_PRAS, active_projects, capacity_forward_years, get_resource_adequacy(simulation), get_peak_load(simulation), get_static_capacity_market(case), iteration_year, simulation_years, get_data_dir(case), get_rt_resolution(case), get_results_dir(simulation), get_outage_dir(case))

        # Parallelize the processing of scenarios using Distributed.pmap
        @time resource_adequacy_tuples = Distributed.pmap(parallelize_update_delta_irm!, zip(scenario_names, sys_PRAS_list, active_projects_list, capacity_forward_years_list, resource_adequacies, peak_loads, static_capacity_bools, iteration_years, simulation_years_list, data_dirs, rt_resolutions, results_dirs, outage_dirs))
        set_resource_adequacy!(simulation, Dict(key => value for (key, value) in resource_adequacy_tuples))

        create_investor_predictions(investors,
                                    active_projects,
                                    iteration_year,
                                    yearly_horizon,
                                    get_data_dir(case),
                                    get_results_dir(simulation),
                                    average_capital_cost_multiplier,
                                    get_zones(simulation),
                                    get_lines(simulation),
                                    get_peak_load(simulation),
                                    get_rps_target(case),
                                    get_reserve_penalty(case),
                                    get_resource_adequacy(simulation),
                                    get_irm_scalar(case),
                                    get_solver(case),
                                    get_parallel_investors(case),
                                    get_parallel_scenarios(case)
                                    )

        for investor in investors
            run_investor_iteration(investor,
                                    active_projects,
                                    iteration_year,
                                    yearly_horizon,
                                    simulation_years,
                                    capacity_forward_years,
                                    sys_MDs,
                                    sys_UCs,
                                    sys_EDs,
                                    sys_PRAS,
                                    case,
                                    scenario_names
                                    )

        end

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

            # Update variable operation cost based on annual carbon tax for SIIP market clearing
            update_operation_cost!(project, sys_MDs[iteration_year], (get_carbon_tax(simulation)), iteration_year)
            update_operation_cost!(project, sys_UCs[iteration_year], (get_carbon_tax(simulation)), iteration_year)
            update_operation_cost!(project, sys_EDs[iteration_year], (get_carbon_tax(simulation)), iteration_year)      
            for scenario in keys(sys_PRAS)
                update_operation_cost!(project, sys_PRAS[scenario], (get_carbon_tax(simulation)), iteration_year)
            end

        end
        installed_capacity = update_installed_cap!(installed_capacity,
                                                   all_existing_projects,
                                                   iteration_year,
                                                   simulation_years)

        println("Current Installed Capacity = $(round(installed_capacity[iteration_year])) MW")

        #Find which markets to simulate.
        markets = union(hcat(get_markets.(get_investors(simulation))...))

        # for d in PSY.get_components(PSYE.ThermalCleanEnergy,sys_UCs[iteration_year])
        #     if "CleanEnergyConstraint" ∉ PSY.get_name.(d.services)
        #         println("$(PSY.get_name(d))")
        #         add_clean_energy_contribution!(sys_UCs[iteration_year], d)
        #     end
        # end

        #Create realzed market prices for existing projects.
        realized_market_prices,
        realized_capacity_factors_md,
        realized_capacity_factors_uc,
        realized_capacity_factors_ed,
        realized_reserve_perc_md,
        realized_reserve_perc_uc,
        realized_reserve_perc_ed,
        realized_inertia_perc,
        capacity_accepted_bids,
        rec_accepted_bids,
        clean_energy_percentage_vector[iteration_year],
        cet_achieved_ratio = create_realized_marketdata(simulation,
                            sys_MDs[iteration_year],
                            sys_UCs[iteration_year],
                            sys_EDs[iteration_year],
                            markets,
                            get_rps_target(case),
                            get_reserve_penalty(case),
                            get_ordc_curved(case),
                            all_existing_projects,
                            capacity_market_projects,
                            capacity_forward_years,
                            iteration_year,
                            simulation_years,
                            get_solver(case),
                            get_results_dir(simulation),
                            current_siip_sim,
                            siip_system)


        existing_project_types = unique(get_type.(get_tech.(all_existing_projects)))
        rt_products = String.(split(read_data(joinpath(get_data_dir(case), "markets_data", "reserve_products.csv"))[1,"rt_products"], "; "))

        if iteration_year < simulation_years

            update_rec_correction_factors!(get_activeprojects(simulation),
                                        realized_capacity_factors_ed,
                                        get_rt_resolution(case),
                                        iteration_year)

            if get_markets(simulation)[:CarbonTax]
                max_carbon_tax_increment = get_max_carbon_tax_increase(case)
                if cet_achieved_ratio == 0.0
                    delta_carbon_tax = max_carbon_tax_increment
                else
                    delta_carbon_tax = max_carbon_tax_increment * max(0.0, (1 - cet_achieved_ratio))
                end
                new_carbon_tax = max((get_carbon_tax(simulation)[iteration_year] + delta_carbon_tax), get_carbon_tax(simulation)[iteration_year + 1])
                simulation.carbon_tax[iteration_year + 1] = new_carbon_tax
            end
        end

        for scenario in keys(sys_PRAS)
            ra_metrics, shortfall = calculate_RA_metrics(deepcopy(sys_PRAS[scenario]),false,get_results_dir(simulation), get_outage_dir(case), iteration_year)
            save_shortfall_data(joinpath(get_results_dir(simulation), "shortfall_data_$(scenario)_year$(iteration_year).h5"), shortfall)
            println(ra_metrics)
            set_metrics!(get_resource_adequacy(simulation)[scenario], iteration_year, ra_metrics)
        end
       
        #Update forecasts and realized profits of all existing projects for each investor.

        for investor in get_investors(simulation)
            # DEPRECATED: Load growth forecast updates commented out due to changes in input timeseries structure.
            #= if iteration_year < simulation_years
                update_forecast!(get_forecast(investor), get_annual_growth(simulation)[:, iteration_year], iteration_year)
            end =#

            projects = get_projects(investor)
            for (i, project) in enumerate(projects)

                update_realized_profits!(project,
                                         realized_market_prices,
                                         realized_capacity_factors_md,
                                         realized_capacity_factors_uc,
                                         realized_capacity_factors_ed,
                                         realized_reserve_perc_md,
                                         realized_reserve_perc_uc,
                                         realized_reserve_perc_ed,
                                         realized_inertia_perc,
                                         capacity_accepted_bids,
                                         rec_accepted_bids,
                                         get_hour_weight(simulation),
                                         iteration_year,
                                         capacity_forward_years,
                                         get_carbon_tax(simulation)[iteration_year],
                                         get_da_resolution(case),
                                         get_rt_resolution(case),
                                         rt_products,
                                         get_pcm_scenario(case))

                update_annual_cashflow!(project, iteration_year)

                retire_old!(projects,
                            i,
                            project,
                            sys_MDs,
                            sys_UCs,
                            sys_EDs,
                            sys_PRAS,
                            get_data_dir(case),
                            iteration_year,
                            scenario_names,
                            total_horizon)

            end

            update_portfolio_preference_multipliers!(investor, iteration_year)

        end

        simulations, iteration_years, derating_scales, methodologies, ra_metric_list, marginal_cc_switches =  repeat_arguments(num_scenarios, simulation, iteration_year, get_derating_scale(case), get_accreditation_methodology(case), get_accreditation_metric(case), get_marginal_cc_switch(case))
        
        @time Distributed.pmap(parallelize_update_derating_data, zip(scenario_names, simulations, iteration_years, derating_scales, methodologies, ra_metric_list, marginal_cc_switches))

        for scenario in scenario_names
            derating_factors = read_data(joinpath(get_data_dir(case), "markets_data", "derating_data", scenario, "derating_dict.csv"))

            output_file = joinpath(get_results_dir(simulation), "derating_data", scenario, "derating_data_year_$(iteration_year+1).h5")

            save_derating_factors(output_file, derating_factors)
        end
    
        active_projects = get_activeprojects(simulation)

        for project in active_projects
            for scenario in scenario_names
                update_derating_factor!(project, get_data_dir(case), scenario, get_derating_scale(case), get_marginal_cc_switch(case))
            end
        end

        # reserve_ts_scaling_factor = calculate_reserve_scaling_factor(simulation)
        reserve_ts_scaling(simulation, iteration_year)

        println("COMPLETED YEAR $(iteration_year)")
        # FileIO.save(joinpath(get_results_dir(simulation), "shortfall_data_year$(iteration_year).jld2"), "shortfall_data", shortfall)
        save_simulation(joinpath(get_results_dir(simulation), "simulation_data_year$(iteration_year).h5"), simulation)
        save_clean_energy_percentage(joinpath(get_results_dir(simulation), "clean_energy_percentage_year$(iteration_year).h5"), clean_energy_percentage_vector)
    end

    final_portfolio = vcat(get_existing.(get_investors(simulation))...)

    for project in final_portfolio
        extrapolate_profits!(project, simulation_years)
    end

    save_clean_energy_percentage(joinpath(get_results_dir(simulation), "clean_energy_percentage.h5"), clean_energy_percentage_vector)
    save_simulation(joinpath(get_results_dir(simulation), "simulation_data.h5"), simulation)

    return
end
