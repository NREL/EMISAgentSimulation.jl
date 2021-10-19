"""
This function returns a vector of investors in the simulation.
"""
function create_investors(simulation_data::AgentSimulationData)
    horizon = get_total_horizon(get_case(simulation_data))
    simulation_data_dir = get_data_dir(get_case(simulation_data))
    dir_name = joinpath(get_data_dir(get_case(simulation_data)), "investors")
    investor_names = readdir(dir_name)
    investors = Vector{Investor}(undef, length(investor_names))
    for i = 1:length(investor_names)
        investor_dir = joinpath(dir_name, "$(investor_names[i])")

        rep_hour_weight = get_rep_hour_weight(simulation_data)

        # Read Investor Characteristics
        characteristics = read_data(joinpath(investor_dir, "characteristics.csv"))

        forecast_type = get_forecast_type(get_case(simulation_data))

        if forecast_type == "perfect"
            parameter_values = [get_annual_growth(simulation_data) for i in 1:horizon]
            scenario_data = [Scenario("perfect", 1.0, nothing, parameter_values)]
            forecast = Perfect(scenario_data)

        elseif forecast_type == "imperfect"

            # Read belief data file
            if get_info_symmetry(get_case(simulation_data))
                belief_filename = joinpath(simulation_data_dir, "markets_data", "symmetric_belief.csv")
            else
                belief_filename = joinpath(investor_dir, "markets_data", "investor_belief.csv")
            end

            @assert isfile(belief_filename)

            investor_belief_data = read_data(belief_filename)

            param_names= investor_belief_data.parameters

            @assert in("initial_estimate", names(investor_belief_data))
            initial_estimate = AxisArrays.AxisArray(investor_belief_data.initial_estimate, param_names)   # initial state initial_estimate

            # Create Kalman Filter for Belief Update
            if get_belief_update(get_case(simulation_data))
                @assert in("initial_error_cov", names(investor_belief_data))
                @assert in("process_cov", names(investor_belief_data))
                @assert in("measurement_cov", names(investor_belief_data))

                initial_error_covariance = AxisArrays.AxisArray(collect(LinearAlgebra.Diagonal(investor_belief_data.initial_error_cov)), param_names, param_names) # initial error covariance

                Q = AxisArrays.AxisArray(investor_belief_data.process_cov, param_names)        # process covariance
                R = AxisArrays.AxisArray(investor_belief_data.measurement_cov, param_names)    # measurement covariance

                belief = InvestorBelief(Q, R)
                kf = KalmanFilter(belief, initial_estimate, initial_error_covariance)

            else
                kf = nothing
            end

            # Populate scenario data
            scenario_data = Scenario[]
                # If user has provided parameter values for each scenario

            if get_uncertainty(get_case(simulation_data))
                # Assert that user has provided parameter multiplier values for each scenario
                if get_info_symmetry(get_case(simulation_data))
                    scenario_file_name = joinpath(simulation_data_dir, "markets_data", "symmetric_scenario_multiplier_data.csv")
                else
                    scenario_file_name = joinpath(investor_dir, "markets_data", "scenario_multiplier_data.csv")
                end

                @assert isfile(scenario_file_name)
                scenario_df = read_data(scenario_file_name)
                num_scenarios = DataFrames.nrow(scenario_df)
                for s in 1:num_scenarios
                    name = scenario_df[s, :scenario]
                    probability = scenario_df[s, :probability]
                    parameter_multipliers = Dict{String, Float64}()
                    parameter_values = [AxisArrays.AxisArray(zeros(length(param_names), horizon), param_names, collect(1:horizon)) for i in 1:horizon]
                    for param in param_names
                        if in(param, names(scenario_df))
                            parameter_multipliers[param] = scenario_df[s, Symbol(param)]
                            for i in 1:horizon
                                parameter_values[1][param, i] = parameter_multipliers[param] * initial_estimate[param]
                            end
                        else
                            parameter_multipliers[param] = 1.0
                            for i in 1:horizon
                                parameter_values[1][param, i] = initial_estimate[param]
                            end
                        end
                    end
                    push!(scenario_data, Scenario(name, probability, parameter_multipliers, parameter_values))
                end
            else
                name = "scenario_1"
                probability = 1.0
                parameter_multipliers = Dict(param => 1.0 for param in param_names)
                parameter_values = [AxisArrays.AxisArray(zeros(length(param_names), horizon), param_names, collect(1:horizon)) for i in 1:horizon]
                for param in param_names
                    for i in 1:horizon
                        parameter_values[1][param, i] = initial_estimate[param]
                    end
                end

                push!(scenario_data, Scenario(name, probability, parameter_multipliers, parameter_values))
            end

            forecast = Imperfect(kf, scenario_data)

        end
        #Empty vector of projects.
        projects = Project{<:BuildPhase}[]

        projectdata_existing = extract_projectdata(investor_dir, "projectexisting.csv")
        projectdata_options = extract_projectdata(investor_dir, "projectoptions.csv")

        sys_UC = get_system_UC(simulation_data)

        #Append existing and option projects.
        append!(projects, create_project_existing(projectdata_existing,
                                                  simulation_data,
                                                  sys_UC,
                                                  investor_names[i],
                                                  investor_dir,
                                                  get_name.(scenario_data)))

        append!(projects, create_project_options(projectdata_options,
                                                      simulation_data,
                                                      investor_names[i],
                                                      investor_dir,
                                                      get_name.(scenario_data)))

        add_investor_project_availability!(simulation_data_dir, projects)

        option_leaftypes = leaftypes(Project{Option})

        option_projects = filter(project -> in(typeof(project), option_leaftypes), projects)

        portfolio_preference_multipliers = Dict{Tuple{String, String}, Vector{Float64}}()

        for project in option_projects
            project_tech_specs = get_tech(project)
            tech = get_type(project_tech_specs)
            zone = get_zone(project_tech_specs)
            portfolio_preference_multipliers[(tech, zone)] = get_project_preference_multiplier(project)
        end

        preference_multiplier_range = (min = characteristics[1, "min_pref_multiplier"], max = characteristics[1, "max_pref_multiplier"])

        #Names of the markets in which the investor is participating.
        markets = Symbol[]
        simulation_markets = get_markets(simulation_data)

        for m in keys(simulation_markets)
            if simulation_markets[m]
                    push!(markets, m)
            end
        end

        #Carbon Tax Data
        simulation_years = get_total_horizon(get_case(simulation_data))
        start_year = get_start_year(get_case(simulation_data))

        carbon_tax = zeros(simulation_years)

        if in(:CarbonTax, markets)
            carbon_tax_data = read_data(joinpath(investor_dir, "markets_data", "CarbonTax.csv"))
            for y in 1:simulation_years
                carbon_tax[y] = carbon_tax_data[findfirst(x -> x == start_year + y - 1, carbon_tax_data[:, "Year"]), "\$/ton"]
            end
        end

        # Empty market prices struct
        market_prices = MarketPrices()

        capital_cost_multiplier = characteristics.capital_cost_multiplier[1]
        max_annual_projects = characteristics.max_annual_projects[1]

        # Risk Preference
        if get_risk_aversion(get_case(simulation_data))
            risk_preference_type = lowercase(characteristics.risk_preference[1])
            if risk_preference_type == "neutral"
                risk_preference = RiskNeutral()
            elseif risk_preference_type == "averse"
                risk_preference = RiskAverse(characteristics.uf_constant[1], characteristics.uf_multiplier[1], characteristics.uf_risk_coefficient[1])
            end
        else
            risk_preference = RiskNeutral()
        end

        retirement_lookback = characteristics.retire_lookback[1]

        investors[i] = Investor(investor_names[i],
                                investor_dir,
                                projects,
                                markets,
                                carbon_tax,
                                market_prices,
                                rep_hour_weight,
                                forecast,
                                capital_cost_multiplier,
                                preference_multiplier_range,
                                portfolio_preference_multipliers,
                                max_annual_projects,
                                risk_preference,
                                retirement_lookback)
    end

    return investors
end
