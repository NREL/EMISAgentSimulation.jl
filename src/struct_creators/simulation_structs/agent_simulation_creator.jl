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

    annual_growth_df = read_data(joinpath(data_dir, "markets_data", "annual_growth.csv"))

    test_system_load_da = DataFrames.DataFrame(CSV.File(joinpath(test_system_dir, "RTS_Data", "timeseries_data_files", "Load", "DAY_AHEAD_regional_Load.csv")))
    test_system_load_rt = DataFrames.DataFrame(CSV.File(joinpath(test_system_dir, "RTS_Data", "timeseries_data_files", "Load", "REAL_TIME_regional_Load.csv")))

    base_year = test_system_load_da[1, "Year"]

    @assert base_year <= start_year

    annual_growth_df_past = filter(row -> row.year < start_year && row.year >= base_year, annual_growth_df)
    annual_growth_df_simulation = filter(row -> row.year >= start_year, annual_growth_df)

    annual_growth_past = AxisArrays.AxisArray(collect(transpose(Matrix(annual_growth_df_past[:, 2:end]))),
                              names(annual_growth_df_past)[2:end],
                              1:DataFrames.nrow(annual_growth_df_past))

    annual_growth_simulation = AxisArrays.AxisArray(collect(transpose(Matrix(annual_growth_df_simulation[:, 2:end]))),
                              names(annual_growth_df_simulation)[2:end],
                              1:DataFrames.nrow(annual_growth_df_simulation))

    zones,
    representative_periods,
    rep_hour_weight,
    chron_weights,
    system_peak_load,
    test_sys_hour_weight,
    zonal_lines = read_test_system(
        data_dir,
        test_system_dir,
        get_base_dir(case),
        test_system_load_da,
        test_system_load_rt,
        base_year,
        annual_growth_past,
        start_year,
        rep_period_interval,
        n_rep_periods,
        rep_checkpoint)

    if isnothing(zones)
        zones = ["zone_1"]
    end

    if isnothing(zonal_lines)
        zonal_lines = [ZonalLine("line_1", zones[1], zones[1], 0.0)]
    end

    markets_dict = get_markets(case)
    
    if get_siip_market_clearing(case)
        base_power = 100.0
        sys_UC, sys_ED = create_rts_sys(test_system_dir, base_power, data_dir, get_da_resolution(case), get_rt_resolution(case))
    else
        sys_UC = nothing
        sys_ED = nothing
    end

    #updating past growth rate in PSY Systems
    for y in 1:size(annual_growth_past)[2]
        apply_PSY_past_load_growth!(sys_UC, annual_growth_past[:, y], data_dir)
        apply_PSY_past_load_growth!(sys_ED, annual_growth_past[:, y], data_dir)
    end

    simulation_years = get_total_horizon(case)

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

    deratingdata = read_data(joinpath(data_dir, "markets_data", "derating_dict.csv"))

    ra_target_file = joinpath(data_dir, "markets_data", "resource_adequacy_targets.csv")
    ra_targets = Dict{String, Float64}()
    ra_metrics = Dict{String, Float64}()


    if isfile(ra_target_file)
        for row in eachrow(read_data(ra_target_file))
            ra_targets[row["Metric"]] = row["Target"]
        end
    end

    resource_adequacy = ResourceAdequacy(ra_targets, zeros(simulation_years), [ra_metrics for i in 1:simulation_years])
    
    results_dir = make_results_dir(case)

    simulation_data = AgentSimulationData(case,
                                        results_dir,
                                        sys_UC,
                                        sys_ED,
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
                                        annual_growth_simulation,
                                        resource_adequacy)

    investors = create_investors(simulation_data)
    set_investors!(simulation_data, investors)

    # convert_thermal_clean_energy!(sys_UC)
    # convert_thermal_clean_energy!(sys_ED)

    convert_thermal_fast_start!(sys_UC)
    convert_thermal_fast_start!(sys_ED)

    construct_ordc(deepcopy(sys_UC), data_dir, investors, 0, representative_periods, rep_period_interval, get_ordc_curved(case), get_ordc_unavailability_method(case), get_reserve_penalty(case))
    add_psy_ordc!(data_dir, markets_dict, sys_UC, "UC", 1, get_da_resolution(case), get_rt_resolution(case), get_reserve_penalty(case))
    add_psy_ordc!(data_dir, markets_dict, sys_ED, "ED", 1, get_da_resolution(case), get_rt_resolution(case), get_reserve_penalty(case))

    if markets_dict[:Inertia]
        add_psy_inertia!(data_dir, sys_UC, "UC", get_reserve_penalty(case), system_peak_load)
        add_psy_inertia!(data_dir, sys_ED, "ED", get_reserve_penalty(case), system_peak_load)
    end

    add_psy_clean_energy_constraint!(sys_UC, initial_rec_requirement)

    transform_psy_timeseries!(sys_UC, sys_ED, get_da_resolution(case), get_rt_resolution(case), 36, 2)

    # Adding representative days availability data to investor folders
    system_availability_data = DataFrames.DataFrame(CSV.File(joinpath(data_dir, "timeseries_data_files", "Availability", "DAY_AHEAD_availability.csv")))

    system_availability_data[!, "Period_Number"] = 1:size(system_availability_data, 1)
    system_availability_data[!, "Representative_Period"] = add_representative_period.(system_availability_data[:, "Period_Number"], rep_period_interval)

    rep_projects_availability = filter(row -> in(row["Representative_Period"], keys(representative_periods)), system_availability_data)

    for dir in get_data_dir.(investors)
        write_data(joinpath(dir, "timeseries_data_files", "Availability"), "DAY_AHEAD_availability.csv", rep_projects_availability)
    end

    update_simulation_derating_data!(
        simulation_data,
        1,
        get_derating_scale(case),
        methodology = get_accreditation_methodology(case),
        ra_metric = get_accreditation_metric(case))

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
    make_case_data_dir(case)
    simulation_data = gather_data(case)
    simulation = AgentSimulation(case,
                            get_results_dir(simulation_data),
                            1,
                            get_system_UC(simulation_data),
                            get_system_ED(simulation_data),
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
                            get_annual_growth(simulation_data),
                            get_resource_adequacy(simulation_data))

    return simulation
end
