""""
This function gathers data for making price and other market data predictions.
"""
function gather_prediction_parameters(investor::Investor,
                                       sys_data_dir::String,
                                       iteration_year::Int64)
    investor_name = get_name(investor)
    investor_dir = get_data_dir(investor)

    load_dir = joinpath(investor_dir, "timeseries_data_files", "Load")
    dir_exists(load_dir)
    cp(joinpath(sys_data_dir, "timeseries_data_files", "Load", "rep_load_$(iteration_year - 1).csv"),
    joinpath(load_dir, "load_$(iteration_year - 1).csv"), force = true)

    reserve_definition = read_data(joinpath(sys_data_dir, "markets_data", "reserve_products.csv"))

    reserve_products = String.(split(reserve_definition[1, "all_products"], "; "))
    ordc_products = String.(split(reserve_definition[1, "ordc_products"], "; "))

    reserve_dir = joinpath(investor_dir, "timeseries_data_files", "Reserves")
    dir_exists(reserve_dir)

    for reserve in reserve_products
        cp(joinpath(sys_data_dir, "timeseries_data_files", "Reserves", "rep_$(reserve)_$(iteration_year - 1).csv"),
        joinpath(reserve_dir, "$(reserve)_$(iteration_year - 1).csv"), force = true)
    end

    market_names = get_markets(investor)

    carbon_tax = get_carbon_tax(investor)

    rep_hour_weight = get_rep_hour_weight(investor)

    scenarios = get_scenario_data(get_forecast(investor))

    return investor_name, investor_dir, market_names, carbon_tax, reserve_products, ordc_products, rep_hour_weight, scenarios
end

"""
This function runs CEM for price predictions based on user-defined parallelization settings.
"""
function create_investor_predictions(investors::Vector{Investor},
                                          active_projects::Vector{Project},
                                          iteration_year::Int64,
                                          yearly_horizon::Int64,
                                          sys_data_dir::String,
                                          sys_results_dir::String,
                                          average_capital_cost_multiplier::Float64,
                                          zones::Vector{String},
                                          lines::Vector{ZonalLine},
                                          peak_load::Float64,
                                          rps_target::String,
                                          reserve_penalty::String,
                                          resource_adequacy::ResourceAdequacy,
                                          irm_scalar::Float64,
                                          solver::JuMP.MOI.OptimizerWithAttributes,
                                          parallelize_investors::Bool,
                                          parallelize_scenarios::Bool)

    if parallelize_investors

        if parallelize_scenarios

            scenarios_pmap = Scenario[]
            investor_name_pmap = String[]
            investor_dir_pmap = String[]
            market_names_pmap = Vector{Symbol}[]
            carbon_tax_pmap = Vector{Float64}[]
            reserve_products_pmap = Vector{String}[]
            ordc_products_pmap = Vector{String}[]
            rps_target_pmap = String[]
            reserve_penalty_pmap = String[]
            resource_adequacy_pmap = ResourceAdequacy[]
            irm_scalar_pmap = Float64[]
            rep_hour_weight_pmap = Vector{Float64}[]
            expected_portfolio_pmap = Vector{Project}[]

            for investor in investors

                investor_name,
                investor_dir,
                market_names,
                carbon_tax,
                reserve_products,
                ordc_products,
                rep_hour_weight,
                scenarios = gather_prediction_parameters(investor, sys_data_dir, iteration_year)

                for scenario in scenarios
                    push!(scenarios_pmap, scenario)
                    push!(investor_name_pmap, investor_name)
                    push!(investor_dir_pmap, investor_dir)
                    push!(market_names_pmap, market_names)
                    push!(carbon_tax_pmap, carbon_tax)
                    push!(reserve_products_pmap, reserve_products)
                    push!(ordc_products_pmap, ordc_products)
                    push!(rps_target_pmap, rps_target)
                    push!(reserve_penalty_pmap, reserve_penalty)
                    push!(resource_adequacy_pmap, resource_adequacy)
                    push!(irm_scalar_pmap, irm_scalar)
                    push!(rep_hour_weight_pmap, rep_hour_weight)
                    push!(expected_portfolio_pmap, active_projects)

                end

            end

            num_tasks = length(scenarios_pmap)
            Distributed.pmap(create_expected_marketdata,
                 investor_dir_pmap,
                 market_names_pmap,
                 carbon_tax_pmap,
                 reserve_products_pmap,
                 ordc_products_pmap,
                 rps_target_pmap,
                 reserve_penalty_pmap,
                 resource_adequacy_pmap,
                 irm_scalar_pmap,
                 expected_portfolio_pmap,
                 repeat([zones], num_tasks),
                 repeat([lines], num_tasks),
                 repeat([peak_load], num_tasks),
                 rep_hour_weight_pmap,
                 repeat([average_capital_cost_multiplier], num_tasks),
                 scenarios_pmap,
                 repeat([iteration_year], num_tasks),
                 repeat([yearly_horizon], num_tasks),
                 repeat([solver], num_tasks),
                 repeat([sys_results_dir], num_tasks),
                 investor_name_pmap)

        else

            num_tasks = length(investors)
            Distributed.pmap(parallelize_only_investors,
                            investors,
                            repeat([sys_data_dir], num_tasks),
                            repeat([active_projects], num_tasks),
                            repeat([rps_target], num_tasks),
                            repeat([reserve_penalty], num_tasks),
                            repeat([resource_adequacy], num_tasks),
                            repeat([irm_scalar], num_tasks),
                            repeat([zones], num_tasks),
                            repeat([lines], num_tasks),
                            repeat([peak_load], num_tasks),
                            repeat([average_capital_cost_multiplier], num_tasks),
                            repeat([iteration_year], num_tasks),
                            repeat([yearly_horizon], num_tasks),
                            repeat([solver], num_tasks),
                            repeat([sys_results_dir], num_tasks),
                            get_name.(investors))
        end

    else

        for investor in investors

            investor_name,
            investor_dir,
            market_names,
            carbon_tax,
            reserve_products,
            ordc_products,
            rep_hour_weight,
            scenarios = gather_prediction_parameters(investor, sys_data_dir, iteration_year)

            if parallelize_scenarios

                num_scenarios = length(scenarios)

                Distributed.pmap(create_expected_marketdata,
                    repeat([investor_dir], num_scenarios),
                    repeat([market_names], num_scenarios),
                    repeat([carbon_tax], num_scenarios),
                    repeat([reserve_products], num_scenarios),
                    repeat([ordc_products], num_scenarios),
                    repeat([rps_target], num_scenarios),
                    repeat([reserve_penalty], num_scenarios),
                    repeat([resource_adequacy], num_scenarios),
                    repeat([irm_scalar], num_scenarios),
                    repeat([active_projects], num_scenarios),
                    repeat([zones], num_scenarios),
                    repeat([lines], num_scenarios),
                    repeat([peak_load], num_scenarios),
                    repeat([rep_hour_weight], num_scenarios),
                    repeat([average_capital_cost_multiplier], num_scenarios),
                    scenarios,
                    repeat([iteration_year], num_scenarios),
                    repeat([yearly_horizon], num_scenarios),
                    repeat([solver], num_scenarios),
                    repeat([sys_results_dir], num_scenarios),
                    repeat([investor_name], num_scenarios))

            else

                for scenario in scenarios
                    create_expected_marketdata(investor_dir,
                                            market_names,
                                            carbon_tax,
                                            reserve_products,
                                            ordc_products,
                                            rps_target,
                                            reserve_penalty,
                                            resource_adequacy,
                                            irm_scalar,
                                            active_projects,
                                            zones,
                                            lines,
                                            peak_load,
                                            rep_hour_weight,
                                            average_capital_cost_multiplier,
                                            scenario,
                                            iteration_year,
                                            yearly_horizon,
                                            solver,
                                            sys_results_dir,
                                            investor_name)
                end

            end

        end

    end

    return
end
