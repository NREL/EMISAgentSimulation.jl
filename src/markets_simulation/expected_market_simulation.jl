"""
This function creates the expected market data for each investor for each scenario using CEM.
"""
function create_expected_marketdata(investor_dir::String,
                                    market_names::Vector{Symbol},
                                    carbon_tax::Vector{Float64},
                                    reserve_products::Vector{String},
                                    ordc_products::Vector{String},
                                    expected_portfolio::Vector{<: Project{<: BuildPhase}},
                                    zones::Vector{String},
                                    lines::Vector{ZonalLine},
                                    peak_load::Float64,
                                    rep_hour_weight::Vector{Float64},
                                    average_capital_cost_multiplier::Float64,
                                    scenario::Scenario,
                                    iteration_year::Int64,
                                    yearly_horizon::Int64,
                                    solver::JuMP.MOI.OptimizerWithAttributes,
                                    sys_results_dir::String,
                                    investor_name::String)


    system = create_cem_mkt_clr_problem(investor_dir,
                                        market_names,
                                        carbon_tax,
                                        reserve_products,
                                        ordc_products,
                                        expected_portfolio,
                                        zones,
                                        lines,
                                        peak_load,
                                        rep_hour_weight,
                                        average_capital_cost_multiplier,
                                        scenario,
                                        iteration_year,
                                        yearly_horizon)

    capacity_price,
    energy_price,
    reserve_price,
    rec_price,
    capacity_factors,
    total_utilization,
    capacity_accepted_perc,
    new_options = cem(
                        system,
                        solver,
                        "C:/Users/manwar2/Documents/GitRepos/emt-tests/data/simulation_data/results/cem_results.txt"
                     )

    output_file = joinpath(investor_dir, "expected_market_data", "$(get_name(scenario))_year_$(iteration_year).jld2")

    FileIO.save(output_file,
                     "capacity_price", capacity_price,
                     "energy_price", energy_price,
                     "reserve_price", reserve_price,
                     "rec_price", rec_price,
                     "capacity_factors", capacity_factors,
                     "total_utilization", total_utilization,
                     "capacity_accepted_perc", capacity_accepted_perc,
                     "new_options", new_options
        )

    sys_results_file = joinpath(sys_results_dir, investor_name, "expected_market_data", "$(get_name(scenario))_year_$(iteration_year).jld2")

    FileIO.save(sys_results_file,
                     "capacity_price", capacity_price,
                     "energy_price", energy_price,
                     "reserve_price", reserve_price,
                     "rec_price", rec_price,
                     "capacity_factors", capacity_factors,
                     "total_utilization", total_utilization,
                     "capacity_accepted_perc", capacity_accepted_perc,
                     "new_options", new_options
        )

    return
end

"""
This function updates the maximum new options which can be built by the investor based on CEM predictions.
"""
function update_max_new_options!(max_new_options::Dict{String, Int64},
                                scenario_new_options::Dict{String, Int64},
                                option_projects::Vector{Project})
    for project in option_projects
        type = get_type(get_tech(project))
        zone = get_zone(get_tech(project))
        if scenario_new_options["option_$(type)_$(zone)"] > max_new_options[get_name(project)]
            max_new_options[get_name(project)] = scenario_new_options["option_$(type)_$(zone)"]
        end
    end
    return max_new_options
end
