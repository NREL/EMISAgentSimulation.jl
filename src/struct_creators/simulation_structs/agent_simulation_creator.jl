"""
This function populates and returns the AgentSimulationData struct.
"""
### NY_change
function gather_data(case::CaseDefinition; results_dir::Union{String, Nothing} = nothing)
    reset_timer!(EMIS_TIMER)
    data_dir = get_data_dir(case)
    test_system_dir = get_sys_dir(case)
    ntp_timeseries_data_dir = get_timeseries_data_dir(case)
    start_year = get_start_year(case)
    rep_period_interval = get_rep_period_interval(case)
    n_rep_periods = get_num_rep_periods(case)
    rep_checkpoint = get_rep_chronology_checkpoint(case)
    simulation_years = get_total_horizon(case)
    rolling_horizon = get_rolling_horizon(case)
    pcm_scenario = get_pcm_scenario(case)
    results_dir = results_dir === nothing ? make_results_dir(case) : results_dir
    da_resolution = get_da_resolution(case)
    rt_resolution = get_rt_resolution(case)
    md_horizon = get_md_horizon(case)
    md_interval = get_md_interval(case)
    uc_horizon = get_uc_horizon(case)
    uc_interval = get_uc_interval(case)
    ed_horizon = get_ed_horizon(case)
    ed_interval = get_ed_interval(case)
    outage_dir = get_outage_dir(case)
    reserve_penalty = get_reserve_penalty(case)
    rps_target = get_rps_target(case)
    base_dir = get_base_dir(case)
    siip_market_clearing = get_siip_market_clearing(case)
    scratch_dir = get_scratch_dir(case)
    derating_scale = get_derating_scale(case)
    accreditation_methodology = get_accreditation_methodology(case)
    accreditation_metric = get_accreditation_metric(case)
    marginal_cc_switch = get_marginal_cc_switch(case)

    timeseries_data_dir = joinpath(results_dir, "timeseries_data_files")
    annual_growth_df = read_data(joinpath(data_dir, "markets_data", "annual_growth.csv"))
    annual_growth_df_simulation = filter(row -> row.year >= start_year, annual_growth_df)

    annual_growth_simulation = AxisArrays.AxisArray(
        collect(transpose(Matrix(annual_growth_df_simulation[:, 2:end]))),
        names(annual_growth_df_simulation)[2:end],
        1:DataFrames.nrow(annual_growth_df_simulation))

    zones = String[]
    zonal_lines = ZonalLine[]
    annual_growth_past_first = []
    scenarios = string.(get_all_scenario_names(data_dir))

    representative_periods = Dict(
        scenario => Dict{
            Int64,
            Union{Dict{Int64, Int64}, OrderedCollections.OrderedDict{Int64, Int64}},
        }() for scenario in scenarios
    )
    test_sys_hour_weight =
        Dict(scenario => Dict{Int64, Vector{Float64}}() for scenario in scenarios)
    rep_hour_weight =
        Dict(scenario => Dict{Int64, Vector{Float64}}() for scenario in scenarios)
    chron_weights = Dict(scenario => Dict{Int64, Matrix{Int64}}() for scenario in scenarios)
    system_peak_load = Dict(scenario => Dict{Int64, Float64}() for scenario in scenarios)

    for scenario in scenarios
        @info "Processing scenario: $scenario"
        for sim_year in collect(1:simulation_years)
            @info "Processing simulation year: $sim_year"
            ### NY_change: copied load data over from /kfs2/projects/gmlcmarkets/Phase2_EMIS_Analysis/Feb2024_ERCOT_2011_MARKET_Test_NGUO_LDES/RTS-GMLC_NY/nys_psy/Data/zonal_load_profile.csv (zonal_model_tscost branch)
            test_system_load_da = DataFrames.DataFrame(
                CSV.File(
                    joinpath(
                        test_system_dir,
                        "RTS_Data",
                        "timeseries_data_files",
                        scenario,
                        "sim_year_$(sim_year)",
                        "Load",
                        "DAY_AHEAD_regional_Load.csv",
                    ),
                ),
            )
            test_system_load_rt = DataFrames.DataFrame(
                CSV.File(
                    joinpath(
                        test_system_dir,
                        "RTS_Data",
                        "timeseries_data_files",
                        scenario,
                        "sim_year_$(sim_year)",
                        "Load",
                        "REAL_TIME_regional_Load.csv",
                    ),
                ),
            )

            base_year = test_system_load_da[1, "Year"]
            @assert base_year <= start_year

            annual_growth_df_past = filter(
                row -> row.year < start_year && row.year >= base_year,
                annual_growth_df,
            )

            annual_growth_past = AxisArrays.AxisArray(
                collect(transpose(Matrix(annual_growth_df_past[:, 2:end]))),
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
            zonal_lines = @timeit EMIS_TIMER "setup/read_test_system" read_test_system(
                data_dir,
                test_system_dir,
                base_dir,
                scenario,
                test_system_load_da,
                test_system_load_rt,
                base_year,
                annual_growth_past,
                start_year,
                sim_year,
                rep_period_interval,
                n_rep_periods,
                rep_checkpoint,
                timeseries_data_dir)

            if isnothing(zones)
                zones = ["zone_1"]
            end

            if isnothing(zonal_lines)
                zonal_lines = [ZonalLine("line_1", zones[1], zones[1], 0.0)]
            end
        end
    end

    markets_dict = get_markets(case)
    sys_MDs = nothing
    sys_UCs = nothing
    sys_EDs = nothing
    sys_PRAS = Dict{String, PSY.System}()

    if siip_market_clearing
        base_power = BASE_POWER
        sys_MDs, sys_UCs, sys_EDs, sys_PRAS,
        MD_horizon, MD_interval, UC_horizon,
        UC_interval, ED_horizon, ED_interval =
            @timeit EMIS_TIMER "setup/create_rts_sys" create_rts_sys(test_system_dir,
                base_power, data_dir,
                scratch_dir, ntp_timeseries_data_dir, scenarios,
                pcm_scenario, simulation_years, da_resolution,
                rt_resolution, md_horizon,
                md_interval, uc_horizon,
                uc_interval, ed_horizon,
                ed_interval, outage_dir,
            )
    end

    carbon_tax = zeros(simulation_years)
    if markets_dict[:CarbonTax]
        carbon_tax_data = read_data(joinpath(data_dir, "markets_data", "CarbonTax.csv"))
        for y in 1:simulation_years
            carbon_tax[y] = carbon_tax_data[
                findfirst(x -> x == start_year + y - 1, carbon_tax_data[:, "Year"]),
                "\$/ton",
            ]
        end
    end

    rec_requirement = zeros(simulation_years)
    initial_rec_requirement = 0.0
    if markets_dict[:REC]
        rec_data = read_data(
            joinpath(data_dir, "markets_data", "REC_$(rps_target)_RPS.csv"),
        )
        initial_rec_requirement = rec_data.rec_req[1]
        rec_increment = rec_data.annual_increment[1]
        rec_requirement =
            [initial_rec_requirement + y * rec_increment for y in 1:simulation_years]
    end

    queue_cost_df = read_data(joinpath(data_dir, "queue_cost_data.csv"))
    deratingdata = Dict(
        s => read_data(
            joinpath(data_dir, "markets_data", "derating_data", s, "derating_dict.csv"),
        ) for s in scenarios
    )

    ra_target_file = joinpath(data_dir, "markets_data", "resource_adequacy_targets.csv")
    ra_targets = Dict{String, Float64}()
    ra_metrics = Dict{String, Float64}()

    if isfile(ra_target_file)
        for row in eachrow(read_data(ra_target_file))
            ra_targets[row["Metric"]] = row["Target"]
        end
    end

    resource_adequacy = Dict(
        s => ResourceAdequacy(
            ra_targets,
            zeros(simulation_years),
            [ra_metrics for i in 1:simulation_years],
        ) for s in scenarios
    )

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

    investors = @timeit EMIS_TIMER "setup/create_investors" create_investors(
        simulation_data,
        timeseries_data_dir,
    )
    set_investors!(simulation_data, investors)

    iteration_year = 1

    # Parallelize the processing of scenarios using Distributed.pmap
    num_scenarios = length(scenarios)
    ### NY_change: didn't work for only 1 scenario
    sys_UC_list,
    data_dirs,
    investors_list,
    representative_periods_list,
    rep_period_intervals,
    cases,
    iteration_years,
    rolling_horizons,
    simulation_years_list,
    timeseries_data_dir_list = repeat_arguments(
        num_scenarios,
        deepcopy(sys_UCs[1]),
        data_dir,
        investors,
        representative_periods,
        rep_period_interval,
        case,
        iteration_year,
        rolling_horizon,
        simulation_years,
        timeseries_data_dir,
    )
    @timeit EMIS_TIMER "setup/ordc_construction" Distributed.pmap(
        parallelize_ordc_construction,
        zip(
            scenarios,
            sys_UC_list,
            data_dirs,
            investors_list,
            representative_periods_list,
            rep_period_intervals,
            cases,
            iteration_years,
            rolling_horizons,
            simulation_years_list,
            timeseries_data_dir_list,
        ),
    )

    @timeit EMIS_TIMER "setup/transform_timeseries" for y in 1:simulation_years
        # convert_thermal_clean_energy!(sys_MDs[y])
        # convert_thermal_clean_energy!(sys_UCs[y])
        # convert_thermal_clean_energy!(sys_EDs[y])

        convert_thermal_fast_start!(sys_MDs[y])
        convert_thermal_fast_start!(sys_UCs[y])
        convert_thermal_fast_start!(sys_EDs[y])

        add_psy_ordc!(
            data_dir,
            markets_dict,
            sys_MDs[y],
            "MD",
            pcm_scenario,
            1,
            da_resolution,
            rt_resolution,
            reserve_penalty,
            timeseries_data_dir,
        )
        add_psy_ordc!(
            data_dir,
            markets_dict,
            sys_UCs[y],
            "UC",
            pcm_scenario,
            1,
            da_resolution,
            rt_resolution,
            reserve_penalty,
            timeseries_data_dir,
        )
        add_psy_ordc!(
            data_dir,
            markets_dict,
            sys_EDs[y],
            "ED",
            pcm_scenario,
            1,
            da_resolution,
            rt_resolution,
            reserve_penalty,
            timeseries_data_dir,
        )

        if markets_dict[:Inertia]
            add_psy_inertia!(
                data_dir,
                sys_MDs[y],
                "MD",
                reserve_penalty,
                system_peak_load,
            )
            add_psy_inertia!(
                data_dir,
                sys_UCs[y],
                "UC",
                reserve_penalty,
                system_peak_load,
            )
            add_psy_inertia!(
                data_dir,
                sys_EDs[y],
                "ED",
                reserve_penalty,
                system_peak_load,
            )
        end

        # TODO: need to update this for MD
        add_psy_clean_energy_constraint!(sys_UCs[y], initial_rec_requirement)

        # NG: this function works for ORDC because ORDC has SingleTimeSeries
        transform_psy_timeseries!(
            sys_MDs[y],
            sys_UCs[y],
            sys_EDs[y],
            da_resolution,
            rt_resolution,
            MD_horizon,
            UC_horizon,
            ED_horizon,
            MD_interval,
            UC_interval,
            ED_interval,
        )
    end

    @timeit EMIS_TIMER "setup/pras_transforms" for scenario in scenarios
        #convert_thermal_clean_energy!(sys_PRAS[scenario])
        PSY.transform_single_time_series!(
            sys_PRAS[scenario],
            Dates.Hour(Int(ED_horizon * 60 / rt_resolution)),
            Dates.Hour(ED_interval),
        )

        convert_thermal_fast_start!(sys_PRAS[scenario])

        add_psy_ordc!(data_dir, markets_dict, sys_PRAS[scenario],
            "PRAS", scenario, 1, da_resolution,
            rt_resolution, reserve_penalty, timeseries_data_dir)

        if markets_dict[:Inertia]
            add_psy_inertia!(
                data_dir,
                sys_PRAS[scenario],
                "PRAS",
                reserve_penalty,
                system_peak_load,
            )
        end
    end

    # Adding representative days availability data
    @timeit EMIS_TIMER "setup/availability_data" for scenario in scenarios
        for sim_year in collect(1:simulation_years)
            system_availability_data = DataFrames.DataFrame(
                CSV.File(
                    joinpath(
                        timeseries_data_dir,
                        scenario,
                        "sim_year_$(sim_year)",
                        "Availability",
                        "DAY_AHEAD_availability.csv",
                    ),
                ),
            )

            system_availability_data[!, "Period_Number"] =
                1:size(system_availability_data, 1)
            system_availability_data[!, "Representative_Period"] =
                add_representative_period.(
                    system_availability_data[:, "Period_Number"],
                    rep_period_interval,
                )

            rep_projects_availability = filter(
                row -> in(
                    row["Representative_Period"],
                    keys(representative_periods[scenario][sim_year]),
                ),
                system_availability_data,
            )

            write_data(
                joinpath(
                    timeseries_data_dir,
                    scenario,
                    "sim_year_$(sim_year)",
                    "Availability",
                ),
                "rep_DAY_AHEAD_availability.csv",
                rep_projects_availability,
            )
        end
    end

    simulations,
    iteration_years,
    derating_scales,
    methodologies,
    ra_metric_list,
    marginal_cc_switches = repeat_arguments(
        num_scenarios,
        simulation_data,
        iteration_year,
        derating_scale,
        accreditation_methodology,
        accreditation_metric,
        marginal_cc_switch,
    )

    @timeit EMIS_TIMER "setup/update_derating" Distributed.pmap(
        parallelize_update_derating_data,
        zip(
            scenarios,
            simulations,
            iteration_years,
            derating_scales,
            methodologies,
            ra_metric_list,
            marginal_cc_switches,
            timeseries_data_dir_list,
        ),
    )

    # update_simulation_derating_data!(
    #     simulation_data,
    #     scenarios[1],
    #     iteration_year,
    #     get_derating_scale(case),
    #     methodology = get_accreditation_methodology(case),
    #     ra_metric = get_accreditation_metric(case),
    #     marginal_cc = get_marginal_cc_switch(case)
    # )

    active_projects = get_activeprojects(simulation_data)

    for project in active_projects
        for scenario in scenarios
            update_derating_factor!(
                project,
                data_dir,
                scenario,
                derating_scale,
                marginal_cc_switch,
            )
        end
    end

    return simulation_data
end

"""
This function creates the data directory for the simulated case.
"""
function make_case_data_dir(case::CaseDefinition)
    case_dir = get_data_dir(case)
    dir_exists(case_dir)
    if !isdir(case_dir)
        projects_type = get_heterogeneity(case) ? "Heterogeneous" : "Homogeneous"
        @info "Copying system data for case $(get_name(case)) from $(get_base_dir(case)) to $case_dir: $projects_type"
        sys_data_dir = joinpath(get_base_dir(case), projects_type)
        cp(sys_data_dir, case_dir; force = true, follow_symlinks = true)
    end
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
function create_agent_simulation(
    case::CaseDefinition;
    results_dir::Union{String, Nothing} = nothing,
)
    simulation_data = gather_data(case; results_dir = results_dir)
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
