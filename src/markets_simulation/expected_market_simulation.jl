"""
This function creates the expected market data for each investor for each scenario using CEM.
"""
function create_expected_marketdata(investor_dir::String,
                                    market_names::Vector{Symbol},
                                    carbon_tax::Vector{Float64},
                                    reserve_products::Vector{String},
                                    ordc_products::Vector{String},
                                    rps_target::String,
                                    reserve_penalty::String,
                                    resource_adequacy::ResourceAdequacy,
                                    irm_scalar::Float64,
                                    expected_portfolio::Vector{<: Project{<: BuildPhase}},
                                    zones::Vector{String},
                                    lines::Vector{ZonalLine},
                                    peak_load::Float64,
                                    rep_period_interval::Int64,
                                    rep_hour_weight::Vector{Float64},
                                    chron_weights::Matrix{Int64},
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
                                        rps_target,
                                        reserve_penalty,
                                        resource_adequacy,
                                        irm_scalar,
                                        expected_portfolio,
                                        zones,
                                        lines,
                                        peak_load,
                                        rep_period_interval,
                                        rep_hour_weight,
                                        chron_weights,
                                        average_capital_cost_multiplier,
                                        scenario,
                                        iteration_year,
                                        yearly_horizon)

    jump_model_dir = joinpath(sys_results_dir, investor_name, "expected_market_data")

    capacity_price,
    energy_price,
    reserve_price,
    rec_price,
    inertia_price,
    capacity_factors,
    total_utilization,
    capacity_accepted_perc,
    new_options,
    new_options_by_type,
    REC_slack,
    REC_supply,
    REC_demand,
    rps_compliant_projects,
    rec_requirement,
    investment,
    retirement,
    max_new_options,
    max_gen,
    demand_e,
    rep_hour_weight,
    p_e_print,
    in_flow_print,
    out_flow_print,
    v_e_print,
    p_in_print,
    linepowerlimit,
    rec_correction,
    p_e_storage_print,
    zone_storage,
    zone_projects,
    lines_to_zone,
    lines_from_zone,
    init_storage,
    p_e_detail,
    remaining_buildtime,
    projects,
    generator_projects,
    storage_projects,
    reserve_up_products,
    ordc_products,
    reserve_down_products,
    p_in_ru_detail,
    p_in_ordc_detail,
    p_in_inertia_detail,
    p_in_rd_detail,
    p_out_ru_detail,
    p_out_ordc_detail,
    p_out_inertia_detail,
    p_out_rd_detail = cem(
                        system,
                        solver,
                        jump_model_dir,
                        "C:/Users/manwar2/Documents/GitRepos/emt-tests/data/simulation_data/results/cem_results.txt"
                     )

    output_file = joinpath(investor_dir, "expected_market_data", "$(get_name(scenario))_year_$(iteration_year).jld2")

    FileIO.save(output_file,
                     "capacity_price", capacity_price,
                     "energy_price", energy_price,
                     "reserve_price", reserve_price,
                     "rec_price", rec_price,
                     "inertia_price", inertia_price,
                     "capacity_factors", capacity_factors,
                     "total_utilization", total_utilization,
                     "capacity_accepted_perc", capacity_accepted_perc,
                     "new_options", new_options,
                     "new_options_by_type", new_options_by_type,
                     "REC_slack", REC_slack,
                     "REC_supply", REC_supply,
                     "REC_demand", REC_demand,
                     "rps_compliant_projects", rps_compliant_projects,
                     "rec_requirement", rec_requirement,
                     "investment", investment,
                     "retirement", retirement,
                     "max_new_options", max_new_options,
                     "max_gen", max_gen,
                     "demand_e", demand_e,
                     "rep_hour_weight", rep_hour_weight,
                     "p_e_print", p_e_print,
                     "in_flow_print", in_flow_print,
                     "out_flow_print", out_flow_print,
                     "v_e_print", v_e_print,
                     "p_in_print", p_in_print,
                     "linepowerlimit", linepowerlimit,
                     "rec_correction", rec_correction,
                     "p_e_storage_print", p_e_storage_print,
                     "zone_storage", zone_storage,
                     "zone_projects", zone_projects,
                     "lines_to_zone", lines_to_zone,
                     "lines_from_zone", lines_from_zone,
                     "init_storage", init_storage,
                     "p_e_detail", p_e_detail,
                     "remaining_buildtime", remaining_buildtime,
                     "projects", projects,
                     "generator_projects", generator_projects,
                     "storage_projects", storage_projects,
                     "reserve_up_products", reserve_up_products,
                     "ordc_products", ordc_products,
                     "reserve_down_products", reserve_down_products,
                     "p_in_ru_detail", p_in_ru_detail,
                     "p_in_ordc_detail", p_in_ordc_detail,
                     "p_in_inertia_detail", p_in_inertia_detail,
                     "p_in_rd_detail", p_in_rd_detail,
                     "p_out_ru_detail", p_out_ru_detail,
                     "p_out_ordc_detail", p_out_ordc_detail,
                     "p_out_inertia_detail", p_out_inertia_detail,
                     "p_out_rd_detail", p_out_rd_detail,
        )

    sys_results_file = joinpath(sys_results_dir, investor_name, "expected_market_data", "$(get_name(scenario))_year_$(iteration_year).jld2")

    FileIO.save(sys_results_file,
                     "capacity_price", capacity_price,
                     "energy_price", energy_price,
                     "reserve_price", reserve_price,
                     "rec_price", rec_price,
                     "inertia_price", inertia_price,
                     "capacity_factors", capacity_factors,
                     "total_utilization", total_utilization,
                     "capacity_accepted_perc", capacity_accepted_perc,
                     "new_options", new_options,
                     "new_options_by_type", new_options_by_type,
                     "REC_slack", REC_slack,
                     "REC_supply", REC_supply,
                     "REC_demand", REC_demand,
                     "rps_compliant_projects", rps_compliant_projects,
                     "rec_requirement", rec_requirement,
                     "investment", investment,
                     "retirement", retirement,
                     "max_new_options", max_new_options,
                     "max_gen", max_gen,
                     "demand_e", demand_e,
                     "rep_hour_weight", rep_hour_weight,
                     "p_e_print", p_e_print,
                     "in_flow_print", in_flow_print,
                     "out_flow_print", out_flow_print,
                     "v_e_print", v_e_print,
                     "p_in_print", p_in_print,
                     "linepowerlimit", linepowerlimit,
                     "rec_correction", rec_correction,
                     "p_e_storage_print", p_e_storage_print,
                     "zone_storage", zone_storage,
                     "zone_projects", zone_projects,
                     "lines_to_zone", lines_to_zone,
                     "lines_from_zone", lines_from_zone,
                     "init_storage", init_storage,
                     "p_e_detail", p_e_detail,
                     "remaining_buildtime", remaining_buildtime,
                     "projects", projects,
                     "generator_projects", generator_projects,
                     "storage_projects", storage_projects,
                     "reserve_up_products", reserve_up_products,
                     "ordc_products", ordc_products,
                     "reserve_down_products", reserve_down_products,
                     "p_in_ru_detail", p_in_ru_detail,
                     "p_in_ordc_detail", p_in_ordc_detail,
                     "p_in_inertia_detail", p_in_inertia_detail,
                     "p_in_rd_detail", p_in_rd_detail,
                     "p_out_ru_detail", p_out_ru_detail,
                     "p_out_ordc_detail", p_out_ordc_detail,
                     "p_out_inertia_detail", p_out_inertia_detail,
                     "p_out_rd_detail", p_out_rd_detail,
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

"""
This function updates the maximum new options by technology type which can be built by the investor based on CEM predictions.
"""
function update_max_new_options_by_type!(max_new_options_by_type::Dict{String, Float64},
                                        scenario_new_options_by_type::Dict{String, Float64})
    for (type, value) in max_new_options_by_type
        if scenario_new_options_by_type[type] > value
            max_new_options_by_type[type] = scenario_new_options_by_type[type]
        end
    end
    return max_new_options_by_type
end
