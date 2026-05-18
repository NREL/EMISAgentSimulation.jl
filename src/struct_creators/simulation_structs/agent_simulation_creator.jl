"""
This function populates and returns the AgentSimulationData struct.
"""
function gather_data(case::CaseDefinition)

    data_dir = get_data_dir(case)
    test_system_dir = get_sys_dir(case)
    start_year = get_start_year(case)
    rep_period_interval = get_rep_period_interval(case)
    n_rep_periods = get_num_rep_periods(case)
    rep_checkpoint = get_rep_chronology_checkpoint(case)
    simulation_years = get_total_horizon(case)
    rolling_horizon = get_rolling_horizon(case)
    pcm_scenario = get_pcm_scenario(case)

    annual_growth_df = read_data(joinpath(data_dir, "markets_data", "annual_growth.csv"))
    annual_growth_df_simulation = filter(row -> row.year >= start_year, annual_growth_df)

    annual_growth_simulation = AxisArrays.AxisArray(collect(transpose(Matrix(annual_growth_df_simulation[:, 2:end]))),
                              names(annual_growth_df_simulation)[2:end],
                              1:DataFrames.nrow(annual_growth_df_simulation))

    zones = String[]
    zonal_lines = ZonalLine[]

    annual_growth_past_first = []

    scenarios = string.(get_all_scenario_names(data_dir))

    representative_periods = Dict(scenario => Dict{Int64, Union{Dict{Int64,Int64}, OrderedCollections.OrderedDict{Int64, Int64}}}() for scenario in scenarios)
    test_sys_hour_weight = Dict(scenario => Dict{Int64, Vector{Float64}}() for scenario in scenarios)
    rep_hour_weight = Dict(scenario => Dict{Int64, Vector{Float64}}() for scenario in scenarios)                    
    chron_weights = Dict(scenario => Dict{Int64, Matrix{Int64}}() for scenario in scenarios)
    system_peak_load = Dict(scenario => Dict{Int64, Float64}() for scenario in scenarios)
    
    for scenario in scenarios
        println(scenario)
        for sim_year in collect(1:simulation_years)
            println(sim_year)
            test_system_load_da = DataFrames.DataFrame(CSV.File(joinpath(test_system_dir, "RTS_Data", "timeseries_data_files", scenario, "sim_year_$(sim_year)", "Load", "DAY_AHEAD_regional_Load.csv")))
            test_system_load_rt = DataFrames.DataFrame(CSV.File(joinpath(test_system_dir, "RTS_Data", "timeseries_data_files", scenario, "sim_year_$(sim_year)", "Load", "REAL_TIME_regional_Load.csv")))

            base_year = test_system_load_da[1, "Year"]

            @assert base_year <= start_year

            annual_growth_df_past = filter(row -> row.year < start_year && row.year >= base_year, annual_growth_df)

            annual_growth_past = AxisArrays.AxisArray(collect(transpose(Matrix(annual_growth_df_past[:, 2:end]))),
                                names(annual_growth_df_past)[2:end],
                                1:DataFrames.nrow(annual_growth_df_past))

            if sim_year == 1
                annual_growth_past_first = annual_growth_past
            end

            zones,
            representative_periods[scenario][sim_year],
            rep_hour_weight[scenario][sim_year],
            chron_weights[scenario][sim_year], 
            system_peak_load[scenario][sim_year],
            test_sys_hour_weight[scenario][sim_year],
            zonal_lines = read_test_system(
                data_dir,
                test_system_dir,
                get_base_dir(case),
                scenario,
                test_system_load_da,
                test_system_load_rt,
                base_year,
                annual_growth_past,
                start_year,
                sim_year,
                rep_period_interval,
                n_rep_periods,
                rep_checkpoint)

            if isnothing(zones)
                zones = ["zone_1"]
            end

            if isnothing(zonal_lines)
                zonal_lines = [ZonalLine("line_1", zones[1], zones[1], 0.0)]
            end

        end
    end

    markets_dict = get_markets(case)
    
    if get_siip_market_clearing(case)
        base_power = 100.0
        sys_MDs, sys_UCs, sys_EDs, sys_PRAS, MD_horizon, MD_interval, UC_horizon, UC_interval, ED_horizon, ED_interval = 
        create_rts_sys(test_system_dir, base_power, data_dir, get_scratch_dir(case), scenarios, pcm_scenario, simulation_years, get_da_resolution(case), get_rt_resolution(case),
            get_md_horizon(case), get_md_interval(case), get_uc_horizon(case), get_uc_interval(case), get_ed_horizon(case), get_ed_interval(case),
            )
    else
        sys_MD = nothing
        sys_UC = nothing
        sys_ED = nothing
    end
    
    #updating past growth rate in PSY Systems
    # for sim_year in 1:simulation_years
    #     for y in 1:size(annual_growth_past_first)[2]
    #         apply_PSY_past_load_growth!(sys_MDs[sim_year], annual_growth_past_first[:, y], data_dir)
    #         apply_PSY_past_load_growth!(sys_UCs[sim_year], annual_growth_past_first[:, y], data_dir)
    #         apply_PSY_past_load_growth!(sys_EDs[sim_year], annual_growth_past_first[:, y], data_dir)
    #     end
    # end

    carbon_tax = zeros(simulation_years)

    if markets_dict[:CarbonTax]
        carbon_tax_data = read_data(joinpath(data_dir, "markets_data", "CarbonTax.csv"))
        for y in 1:simulation_years
            carbon_tax[y] = carbon_tax_data[findfirst(x -> x == start_year + y - 1, carbon_tax_data[:, "Year"]), "\$/ton"]
        end
    end

    rec_requirement = zeros(simulation_years)
    initial_rec_requirement = 0.0
    if markets_dict[:REC]
        rec_data = read_data(joinpath(data_dir, "markets_data", "REC_$(get_rps_target(case))_RPS.csv"))
        initial_rec_requirement = rec_data.rec_req[1]
        rec_increment = rec_data.annual_increment[1]
        rec_requirement = [initial_rec_requirement + y * rec_increment for y in 1:simulation_years]
    end

    queue_cost_df = read_data(joinpath(data_dir, "queue_cost_data.csv"))

    deratingdata = Dict(s => read_data(joinpath(data_dir, "markets_data", "derating_data", s, "derating_dict.csv")) for s in scenarios)

    ra_target_file = joinpath(data_dir, "markets_data", "resource_adequacy_targets.csv")
    ra_targets = Dict{String, Float64}()
    ra_metrics = Dict{String, Float64}()


    if isfile(ra_target_file)
        for row in eachrow(read_data(ra_target_file))
            ra_targets[row["Metric"]] = row["Target"]
        end
    end

    
    resource_adequacy = Dict(s => ResourceAdequacy(ra_targets, zeros(simulation_years), [ra_metrics for i in 1:simulation_years]) for s in scenarios)
    
    results_dir = make_results_dir(case)

    simulation_data = AgentSimulationData(case,
                                        results_dir,
                                        sys_MDs,
                                        sys_UCs,
                                        sys_EDs,
                                        sys_PRAS,
                                        zones,
                                        zonal_lines,
                                        representative_periods,
                                        rep_period_interval,
                                        test_sys_hour_weight,
                                        rep_hour_weight,
                                        chron_weights,
                                        system_peak_load,
                                        markets_dict,
                                        carbon_tax,
                                        rec_requirement,
                                        queue_cost_df,
                                        deratingdata,
                                        resource_adequacy)

    investors = create_investors(simulation_data)
    set_investors!(simulation_data, investors)
    
    iteration_year = 1 

    # Parallelize the processing of scenarios using Distributed.pmap
    num_scenarios = length(scenarios)
    sys_UC_list, data_dirs, investors_list, representative_periods_list, rep_period_intervals, cases, iteration_years, rolling_horizons, simulation_years_list = repeat_arguments(num_scenarios, deepcopy(sys_UCs[1]), data_dir, investors, representative_periods, rep_period_interval, case, iteration_year, rolling_horizon, simulation_years)
    @time Distributed.pmap(parallelize_ordc_construction, zip(scenarios, sys_UC_list, data_dirs, investors_list, representative_periods_list, rep_period_intervals, cases, iteration_years, rolling_horizons, simulation_years_list))
    

    for y in 1:simulation_years
        # convert_thermal_clean_energy!(sys_MDs[y])
        # convert_thermal_clean_energy!(sys_UCs[y])
        # convert_thermal_clean_energy!(sys_EDs[y])

        convert_thermal_fast_start!(sys_MDs[y])
        convert_thermal_fast_start!(sys_UCs[y])
        convert_thermal_fast_start!(sys_EDs[y])
        
        add_psy_ordc!(data_dir, markets_dict, sys_MDs[y], "MD", pcm_scenario, 1, get_da_resolution(case), get_rt_resolution(case), get_reserve_penalty(case), get_ordc_curved(case))
        add_psy_ordc!(data_dir, markets_dict, sys_UCs[y], "UC", pcm_scenario, 1, get_da_resolution(case), get_rt_resolution(case), get_reserve_penalty(case), get_ordc_curved(case))
        add_psy_ordc!(data_dir, markets_dict, sys_EDs[y], "ED", pcm_scenario, 1, get_da_resolution(case), get_rt_resolution(case), get_reserve_penalty(case), get_ordc_curved(case))

        if markets_dict[:Inertia]
            add_psy_inertia!(data_dir, sys_MDs[y], "MD", get_reserve_penalty(case), system_peak_load)
            add_psy_inertia!(data_dir, sys_UCs[y], "UC", get_reserve_penalty(case), system_peak_load)
            add_psy_inertia!(data_dir, sys_EDs[y], "ED", get_reserve_penalty(case), system_peak_load)
        end
        
        # TODO: need to update this for MD
        add_psy_clean_energy_constraint!(sys_UCs[y], initial_rec_requirement)

        # NG: this function works for ORDC because ORDC has SingleTimeSeries
        transform_psy_timeseries!(sys_MDs[y], sys_UCs[y], sys_EDs[y], get_da_resolution(case), get_rt_resolution(case), MD_horizon, UC_horizon, ED_horizon, MD_interval, UC_interval, ED_interval)    
    end
    
    for scenario in scenarios
        #convert_thermal_clean_energy!(sys_PRAS[scenario])
        convert_thermal_fast_start!(sys_PRAS[scenario])
        add_psy_ordc!(data_dir, markets_dict, sys_PRAS[scenario], "PRAS", scenario, 1, get_da_resolution(case), get_rt_resolution(case), get_reserve_penalty(case), get_ordc_curved(case))

        if markets_dict[:Inertia]
            add_psy_inertia!(data_dir, sys_PRAS[scenario], "PRAS", get_reserve_penalty(case), system_peak_load)
        end

        PSY.transform_single_time_series!(sys_PRAS[scenario], Int(ED_horizon * 60 / get_rt_resolution(case)), Dates.Hour(ED_interval))   
    end
 
    # Adding representative days availability data
    for scenario in scenarios
        for sim_year in collect(1:simulation_years)
            system_availability_data = DataFrames.DataFrame(CSV.File(joinpath(data_dir, "timeseries_data_files", scenario, "sim_year_$(sim_year)", "Availability", "DAY_AHEAD_availability.csv")))

            system_availability_data[!, "Period_Number"] = 1:size(system_availability_data, 1)
            system_availability_data[!, "Representative_Period"] = add_representative_period.(system_availability_data[:, "Period_Number"], rep_period_interval)

            rep_projects_availability = filter(row -> in(row["Representative_Period"], keys(representative_periods[scenario][sim_year])), system_availability_data)

            write_data(joinpath(data_dir, "timeseries_data_files", scenario, "sim_year_$(sim_year)","Availability"), "rep_DAY_AHEAD_availability.csv", rep_projects_availability)
        end
    end

    simulations, iteration_years, derating_scales, methodologies, ra_metric_list, marginal_cc_switches =  repeat_arguments(num_scenarios, simulation_data, iteration_year, get_derating_scale(case), get_accreditation_methodology(case), get_accreditation_metric(case), get_marginal_cc_switch(case))
   
    @time Distributed.pmap(parallelize_update_derating_data, zip(scenarios, simulations, iteration_years, derating_scales, methodologies, ra_metric_list, marginal_cc_switches))
    
    active_projects = get_activeprojects(simulation_data)

    for project in active_projects
        for scenario in scenarios
            update_derating_factor!(project, data_dir, scenario, get_derating_scale(case), get_marginal_cc_switch(case))
        end
    end

    return simulation_data
end

"""
This function creates the data directory for the simulated case.
"""
function make_case_data_dir(case::CaseDefinition)
    base_dir = get_base_dir(case)
    if get_heterogeneity(case)
        sys_data_dir = joinpath(base_dir, "Heterogeneous")

    else
        sys_data_dir = joinpath(base_dir, "Homogeneous")
    end

    case_dir = get_data_dir(case)
    dir_exists(case_dir)
    cp(sys_data_dir, case_dir, force=true, follow_symlinks=true)

    return
end

"""
This function creates the results directory for the simulated case.
"""
function make_results_dir(case::CaseDefinition)

    case_name = get_name(case)

    results_dir = joinpath(".", "Results", case_name)
    dir_exists(results_dir)

    return results_dir
end

"""
This function returns the AgentSimulation struct which contains all the required data for running the simulation.
"""
function create_agent_simulation(case::CaseDefinition)
    simulation_data = gather_data(case)
    simulation = AgentSimulation(case,
                            get_results_dir(simulation_data),
                            1,
                            get_system_MDs(simulation_data),
                            get_system_UCs(simulation_data),
                            get_system_EDs(simulation_data),
                            get_system_PRAS(simulation_data),
                            get_zones(simulation_data),
                            get_lines(simulation_data),
                            get_rep_periods(simulation_data),
                            get_rep_period_interval(simulation_data),
                            get_hour_weight(simulation_data),
                            get_peak_load(simulation_data),
                            get_markets(simulation_data),
                            get_carbon_tax(simulation_data),
                            get_rec_requirement(simulation_data),
                            get_investors(simulation_data),
                            get_derating_data(simulation_data),
                            get_resource_adequacy(simulation_data))

    return simulation
end
