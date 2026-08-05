
function run_agent_simulation(
    simulation::AgentSimulation,
    current_siip_sim,
    siip_system,
    current_year::Int64,
)
    case = get_case(simulation)
    total_horizon = get_total_horizon(case)
    rolling_horizon = get_rolling_horizon(case)
    step_size = get_step_size(case)
    simulation_years = get_simulation_years(case)

    installed_capacity = zeros(simulation_years)
    capacity_forward_years = get_capacity_forward_years(simulation)
    results_dir = get_results_dir(simulation)
    timeseries_data_dir = joinpath(results_dir, "timeseries_data_files")
    simulation_dir = get_data_dir(get_case(simulation))
    scenario_names = String.(get_all_scenario_names(get_data_dir(case)))
    total_sim_time = 0.0

    availability_rt_by_scenario = Dict(
        scenario => get_availability_df(
            timeseries_data_dir,
            scenario,
            simulation_years,
            "REAL_TIME",
        )
        for scenario in scenario_names
    )
    availability_by_scenario = Dict(
        scenario => get_availability_df(
            timeseries_data_dir,
            scenario,
            simulation_years,
            "DAY_AHEAD",
        )
        for scenario in scenario_names
    )

    # Set initial capacity market profits considering forward capacity auctions
    @info "Setting initial capacity market profits for existing projects based on forward capacity auctions"
    if current_year == 1
        if get_markets(simulation)[:Capacity]
            initial_existing_projects = vcat(get_existing.(get_investors(simulation))...)

            capacity_mkt_param_file =
                joinpath(simulation_dir, "markets_data", "Capacity.csv")
            capacity_mkt_params = read_data(capacity_mkt_param_file)[1, :]
            introduction_year = capacity_mkt_params["introduction_year"]

            if introduction_year >= 3
                initial_capacity_prices = [0.0, 0.0]              # initial capacity prices - arbitrarily selected here - #TODO: need some meachanism to generate these
            else
                initial_capacity_prices = [70000.0, 80000.0]
            end

            for y in 1:(capacity_forward_years - 1)
                for project in initial_existing_projects
                    if get_end_life_year(project) >= y
                        for product in get_products(project)
                            update_initial_capacity_revenues!(
                                project,
                                product,
                                initial_capacity_prices,
                                y,
                                get_pcm_scenario(case),
                            )
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
    carbon_tax = get_carbon_tax(simulation)

    # Update operation cost for all projects based on carbon tax in the first year
    for iteration_year in current_year:step_size:simulation_years
        t_start = time()
        yearly_horizon = min(total_horizon - iteration_year + 1, rolling_horizon)
        ts_now = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
        @info "Starting iteration year $(iteration_year) @ $(ts_now)"
        set_iteration_year!(simulation, iteration_year)

        active_projects = deepcopy(get_activeprojects(simulation))
        installed_capacity = update_installed_cap!(installed_capacity,
            active_projects,
            iteration_year,
            simulation_years)

        # save existing net load csv file for potential checkpoint re-runs
        for scenario in scenario_names
            pre_update_da_net_load = joinpath(
                timeseries_data_dir,
                scenario,
                "sim_year_$(iteration_year)",
                "Net Load Data",
                "load_n_vg_data_pre_update.csv",
            )
            pre_update_rt_net_load = joinpath(
                timeseries_data_dir,
                scenario,
                "sim_year_$(iteration_year)",
                "Net Load Data",
                "load_n_vg_data_rt_pre_update.csv",
            )
            post_update_da_net_load = joinpath(
                timeseries_data_dir,
                scenario,
                "sim_year_$(iteration_year)",
                "Net Load Data",
                "load_n_vg_data.csv",
            )
            post_update_rt_net_load = joinpath(
                timeseries_data_dir,
                scenario,
                "sim_year_$(iteration_year)",
                "Net Load Data",
                "load_n_vg_data_rt.csv",
            )
            # Restore from pre_update only when genuinely restarting at this exact year.
            # On a fresh run or for years beyond the restart year, always save so that
            # stale pre_update files copied from a prior run don't clobber retirements
            # and new builds applied in earlier years of this run.
            if isfile(pre_update_da_net_load) && iteration_year == current_year &&
               current_year > 1
                cp(pre_update_da_net_load, post_update_da_net_load; force = true)
                cp(pre_update_rt_net_load, post_update_rt_net_load; force = true)
            else
                cp(post_update_da_net_load, pre_update_da_net_load; force = true)
                cp(post_update_rt_net_load, pre_update_rt_net_load; force = true)
            end
        end

        if current_year == 1
            for scenario in scenario_names
                derating_factors = read_data(
                    joinpath(
                        get_data_dir(case),
                        "markets_data",
                        "derating_data",
                        scenario,
                        "derating_dict.csv",
                    ),
                )

                save_derating_factors(
                    joinpath(
                        results_dir,
                        "derating_data",
                        scenario,
                        "derating_data_year_$(iteration_year).h5",
                    ),
                    derating_factors,
                )
            end
        end

        # Kept for traceability with pg/grid-solutions (commented to avoid stale runtime setup).
        # num_scenarios = length(scenario_names)
        # sys_PRAS_list, active_projects_list, capacity_forward_years_list,
        # resource_adequacies, peak_loads, static_capacity_bools,
        # iteration_years, simulation_years_list, data_dirs,
        # rt_resolutions, results_dirs,
        # outage_dirs = repeat_arguments(num_scenarios, sys_PRAS, active_projects,
        #     capacity_forward_years, get_resource_adequacy(simulation),
        #     get_peak_load(simulation), get_static_capacity_market(case),
        #     iteration_year, simulation_years, get_data_dir(case),
        #     get_rt_resolution(case),
        #     get_results_dir(simulation), get_outage_dir(case))
        # @info "Resource adequacies: $(resource_adequacies)"

        # Parallelize the processing of scenarios using Distributed.pmap
        # NOTE: pmap risks OOM because each worker receives a full sys_PRAS dict copy
        # (all scenarios serialized to each worker). Use sequential loop instead and
        # rely on --threads auto for PRAS multi-threading in the main process.
        # @timeit EMIS_TIMER "update_delta_irm" resource_adequacy_tuples = Distributed.pmap(parallelize_update_delta_irm!,
        # zip(scenario_names, sys_PRAS_list, active_projects_list, capacity_forward_years_list,
        # resource_adequacies, peak_loads, static_capacity_bools, iteration_years, simulation_years_list,
        # data_dirs, rt_resolutions, results_dirs, outage_dirs))

        resource_adequacy_tuples = []
        for scenario in scenario_names
            @info "Updating resource adequacy for scenario: $scenario"
            resource_adequacy = @timeit EMIS_TIMER "update_delta_irm" update_delta_irm!(
                sys_PRAS[scenario],
                active_projects,
                capacity_forward_years,
                get_resource_adequacy(simulation)[scenario],
                get_peak_load(simulation)[scenario][min(
                    iteration_year + capacity_forward_years - 1,
                    simulation_years,
                )],
                get_static_capacity_market(case),
                scenario,
                iteration_year,
                get_data_dir(case),
                get_rt_resolution(case),
                get_results_dir(simulation),
                get_outage_dir(case),
                simulation_years,
            )
            @info "Resource adequacy for scenario $scenario: $resource_adequacy"
            push!(resource_adequacy_tuples, (scenario, resource_adequacy))
        end

        @info "Setting resource adequacy for all scenarios in the simulation"
        set_resource_adequacy!(
            simulation,
            Dict(key => value for (key, value) in resource_adequacy_tuples),
        )

        @info "Creating investor predictions for all investors based on updated resource adequacy and other market data"
        @timeit EMIS_TIMER "investor_predictions" create_investor_predictions(investors,
            active_projects,
            iteration_year,
            yearly_horizon,
            get_data_dir(case),
            get_results_dir(simulation),
            timeseries_data_dir,
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
            get_parallel_scenarios(case),
        )

        for investor in investors
            @timeit EMIS_TIMER "investor_iteration" run_investor_iteration(investor,
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
                scenario_names,
                timeseries_data_dir,
                availability_rt_by_scenario,
                availability_by_scenario,
            )
        end

        @info "Getting all existing projects to calculate realized profits for energy and REC markets."
        all_existing_projects = vcat(get_existing.(get_investors(simulation))...)

        # Get all projects which are expected to be online for the forward capacity market auction.
        capacity_market_year = iteration_year + capacity_forward_years - 1
        capacity_market_projects = Project[]

        for project in get_activeprojects(simulation)
            # @info "Year $(iteration_year): Updating operation costs for project $(get_name(project))"
            end_life_year = get_end_life_year(project)
            construction_year = get_construction_year(project)
            if end_life_year >= capacity_market_year &&
               construction_year <= capacity_market_year
                push!(capacity_market_projects, project)
            end

            # Update variable operation cost based on annual carbon tax for SIIP market clearing
            # @info "Updating variable operation"
            update_operation_cost!(
                project,
                sys_MDs[iteration_year],
                carbon_tax,
                iteration_year,
            )
            update_operation_cost!(
                project,
                sys_UCs[iteration_year],
                carbon_tax,
                iteration_year,
            )
            update_operation_cost!(
                project,
                sys_EDs[iteration_year],
                carbon_tax,
                iteration_year,
            )
            for scenario in keys(sys_PRAS)
                update_operation_cost!(
                    project,
                    sys_PRAS[scenario],
                    carbon_tax,
                    iteration_year,
                )
            end
        end
        installed_capacity = update_installed_cap!(installed_capacity,
            all_existing_projects,
            iteration_year,
            simulation_years)

        @info "Current Installed Capacity = $(round(installed_capacity[iteration_year])) MW"

        # Find which markets to simulate.
        markets = union(hcat(get_markets.(get_investors(simulation))...))

        # for d in PSY.get_components(PSYE.ThermalCleanEnergy,sys_UCs[iteration_year])
        #     if "CleanEnergyConstraint" ∉ PSY.get_name.(d.services)
        #         println("$(PSY.get_name(d))")
        #         add_clean_energy_contribution!(sys_UCs[iteration_year], d)
        #     end
        # end

        # Create realzed market prices for existing projects.
        @info "Creating realized market data for existing projects to calculate profits and update forecasts for next iteration"
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
        cet_achieved_ratio =
            @timeit EMIS_TIMER "realized_marketdata" create_realized_marketdata(simulation,
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
        rt_products = String.(
            split(
                read_data(
                    joinpath(get_data_dir(case), "markets_data", "reserve_products.csv"),
                )[
                    1,
                    "rt_products",
                ],
                "; ",
            ),
        )

        if iteration_year < simulation_years
            update_rec_correction_factors!(get_activeprojects(simulation),
                realized_capacity_factors_ed,
                get_rt_resolution(case),
                iteration_year,
                step_size)

            if get_markets(simulation)[:CarbonTax]
                max_carbon_tax_increment = get_max_carbon_tax_increase(case)
                if cet_achieved_ratio == 0.0
                    delta_carbon_tax = max_carbon_tax_increment
                else
                    delta_carbon_tax =
                        max_carbon_tax_increment * max(0.0, (1 - cet_achieved_ratio))
                end
                new_carbon_tax = max(
                    (get_carbon_tax(simulation)[iteration_year] + delta_carbon_tax),
                    get_carbon_tax(simulation)[iteration_year + step_size],
                )
                simulation.carbon_tax[iteration_year + step_size] = new_carbon_tax
            end
        end

        for scenario in keys(sys_PRAS)
            ra_metrics, shortfall = @timeit EMIS_TIMER "ra_metrics" calculate_RA_metrics(
                sys_PRAS[scenario],
                false,
                results_dir,
                get_outage_dir(case),
                iteration_year,
                simulation_years = simulation_years,
            )

            save_shortfall_data(
                joinpath(
                    results_dir,
                    "shortfall_data_$(scenario)_year$(iteration_year).h5",
                ),
                shortfall,
            )
            @info "RA Metrics for scenario $scenario in year $iteration_year: $ra_metrics"
            set_metrics!(
                get_resource_adequacy(simulation)[scenario],
                iteration_year,
                ra_metrics,
            )
        end

        #Update forecasts and realized profits of all existing projects for each investor.

        for investor in get_investors(simulation)
            projects = get_projects(investor)
            for (i, project) in enumerate(projects)
                # @info "$(i): Updating realized profits for $(get_name(project))"
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
                    step_size,
                    scenario_names,
                    total_horizon,
                    timeseries_data_dir)
            end

            update_portfolio_preference_multipliers!(investor, iteration_year)
        end

        @info "Updating derating data for all scenarios in the simulation based on updated resource adequacy and market conditions"
        # Parallelize the processing of scenarios using Distributed.pmap
        # NOTE: pmap can crash here because each worker receives serialized PSY.System
        # objects with process-local SQLite handles inside sys_PRAS.
        # simulations, iteration_years,
        # methodologies, ra_metric_list, marginal_cc_switches, timeseries_data_dir_list =
        #     repeat_arguments(num_scenarios,
        #         simulation, iteration_year,
        #         get_accreditation_methodology(case), get_accreditation_metric(case),
        #         get_marginal_cc_switch(case), timeseries_data_dir)
        # @timeit EMIS_TIMER "update_derating" Distributed.pmap(parallelize_update_derating_data,
        #     zip(scenario_names, simulations, iteration_years,
        #         methodologies, ra_metric_list, marginal_cc_switches,
        #         timeseries_data_dir_list))

        # Run sequentially to avoid Distributed serialization of PSY.System
        # objects that contain process-local SQLite handles.
        @timeit EMIS_TIMER "update_derating" for scenario in scenario_names
            update_simulation_derating_data!(
                simulation,
                scenario,
                iteration_year,
                timeseries_data_dir;
                methodology = get_accreditation_methodology(case),
                ra_metric = get_accreditation_metric(case),
                marginal_cc = get_marginal_cc_switch(case),
            )
        end

        for scenario in scenario_names
            derating_factors = read_data(
                joinpath(
                    get_data_dir(case),
                    "markets_data",
                    "derating_data",
                    scenario,
                    "derating_dict.csv",
                ),
            )

            save_derating_factors(
                joinpath(
                    results_dir,
                    "derating_data",
                    scenario,
                    "derating_data_year_$(iteration_year+step_size).h5",
                ),
                derating_factors,
            )
        end

        active_projects = get_activeprojects(simulation)

        for project in active_projects
            for scenario in scenario_names
                update_derating_factor!(
                    project,
                    get_data_dir(case),
                    scenario,
                    get_marginal_cc_switch(case),
                )
            end
        end

        # reserve_ts_scaling_factor = calculate_reserve_scaling_factor(simulation)
        @timeit EMIS_TIMER "reserve_ts_scaling" reserve_ts_scaling(
            simulation,
            iteration_year,
            step_size,
        )

        @timeit EMIS_TIMER "save_year_data" begin
            @info "COMPLETED ITERATION YEAR $(iteration_year)"

            save_clean_energy_percentage(
                joinpath(results_dir, "clean_energy_percentage_year$(iteration_year).h5"),
                clean_energy_percentage_vector,
            )

            save_simulation(simulation, results_dir, iteration_year)
        end
        t_end = time()
        iteration_time_hours = round((t_end - t_start) / 3600; digits = 2)
        total_sim_time += iteration_time_hours
        ts_now = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
        @info "Finished iteration year $(iteration_year) @ $(ts_now)"
        @info "Iteration year $(iteration_year) took $(iteration_time_hours) hours"
        @info "Total simulation time after completing iteration year $(iteration_year): $(round(total_sim_time, digits=2)) hours"
        print_timer(stderr, EMIS_TIMER)
    end

    final_portfolio = vcat(get_existing.(get_investors(simulation))...)

    @info "Extrapolating profits for all projects in the final portfolio"
    for project in final_portfolio
        extrapolate_profits!(project, simulation_years)
    end

    @info "Saving final simulation data and clean energy percentage vector"
    save_clean_energy_percentage(
        joinpath(results_dir, "clean_energy_percentage.h5"),
        clean_energy_percentage_vector,
    )

    save_simulation(simulation, results_dir, simulation_years)

    @info "EMIS SIMULATION COMPLETED!"
    print_timer(stderr, EMIS_TIMER)
    return
end
