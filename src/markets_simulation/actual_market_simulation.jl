
"""
This function creates the realized market data based on actual market clearing.
"""
function create_realized_marketdata(simulation::AgentSimulation,
    sys_MD::Union{Nothing, PSY.System},
    sys_UC::Union{Nothing, PSY.System},
    sys_ED::Union{Nothing, PSY.System},
    market_names::Vector{Symbol},
    rps_target::String,
    reserve_penalty::String,
    ordc_curved::Bool,
    existing_projects::Vector{<: Project{<: BuildPhase}},
    capacity_market_projects::Vector{<: Project{<: BuildPhase}},
    capacity_forward_years::Int64,
    iteration_year::Int64,
    simulation_years::Int64,
    solver::JuMP.MOI.OptimizerWithAttributes,
    results_dir::String,
    current_siip_sim,
    siip_system)
    num_invperiods = 1

    simulation_dir = get_data_dir(get_case(simulation))
    zones = get_zones(simulation)
    lines = get_lines(simulation)
    case = get_case(simulation)
    pcm_scenario = get_pcm_scenario(case)
    all_scenarios = String.(get_all_scenario_names(simulation_dir))
    md_market_bool = get_md_market(case)
    single_stage_bool = get_single_stage(case)
    step_size = get_step_size(case)

    # DEPRECATED: Commenting out all code related to load growth.

    #annual_growth = get_annual_growth(simulation)[:, iteration_year]

    hour_weight = get_hour_weight(simulation)[pcm_scenario][iteration_year]

    # Get actual value of load growth
    #= load_growth = AxisArrays.AxisArray(zeros(length(zones)), zones)
    for (idx, zone) in enumerate(zones)
        load_growth[idx] = annual_growth["load_$(zone)"]
    end

    energy_annual_increment = load_growth
    reserveup_annual_increment = load_growth
    reservedown_annual_increment = load_growth
    =#
    ######## Energy and Ancillary Services market clearing ############################################

    zonal_load, system = create_economic_dispatch_problem(simulation,
        sys_UC,
        market_names,
        num_invperiods,
        existing_projects,
        pcm_scenario,
        iteration_year)

    num_hours = size(zonal_load)[2]
    num_days = Int(num_hours/24)

    rec_perc_requirement = get_rec_requirement(simulation)[iteration_year]

    energy_price_ed,
    energy_price_uc,
    energy_price_md,
    reserve_price_ed,
    reserve_price_uc,
    reserve_price_md,
    inertia_price,
    capacity_factors_md,
    capacity_factors_uc,
    capacity_factors_ed,
    reserve_perc_md,
    reserve_perc_uc,
    reserve_perc_ed,
    inertia_perc,
    start_up_costs,
    shut_down_costs,
    energy_voll,
    energy_voll_uc,
    energy_voll_md,
    reserve_voll,
    reserve_voll_uc,
    reserve_voll_md,
    inertia_voll = energy_mkt_clearing(
        sys_MD,
        sys_UC,
        sys_ED,
        system,
        simulation_dir,
        reserve_penalty,
        rec_perc_requirement,
        zones,
        num_days,
        pcm_scenario,
        iteration_year,
        get_da_resolution(case),
        get_rt_resolution(case),
        get_name(case),
        solver,
        get_base_dir(case),
        simulation,
        current_siip_sim,
        md_market_bool,
        single_stage_bool,
        siip_system,
        case,
    )
    @info "Finsihed Energy Market Clearing for year $(iteration_year)"
    @info "Clean energy requirement for this year is $(get_rec_requirement(simulation)[iteration_year] * 100) percent"
    total_production = 0.0
    total_cec_production = 0.0
    day = 0
    rt_resolution = get_rt_resolution(get_case(simulation))
    for time in 1:Int(24 * 60 / rt_resolution):(Int(24 * 60 / rt_resolution) * 365)
        day += 1
        daily_total_production = 0.0
        daily_cec_production = 0.0
        for gen in get_all_techs(sys_ED)
            name = PSY.get_name(gen)
            prime_mover = string(PSY.get_prime_mover_type(gen))

            clean_energy_project_check =
                any(occursin.(["WT", "PV", "HY"], prime_mover)) || (
                    occursin("ST", prime_mover) &&
                    occursin("NUCLEAR", string(PSY.get_fuel(gen)))
                )

            if !(occursin("BA", prime_mover))
                energy_production =
                    sum(
                        capacity_factors_ed[name][time:(time + Int(
                            24 * 60 / rt_resolution,
                        ) - 1)],
                    ) * get_device_size(gen)
                total_production += energy_production
                daily_total_production += energy_production

                if clean_energy_project_check
                    total_cec_production += energy_production
                    daily_cec_production += energy_production
                end
            end
        end
        #println("Clean energy contribution for day $(day) is $(round(daily_cec_production * 100.0 / daily_total_production, digits = 2)) percent")
    end

    @info "Total Annual clean energy contribution is $(round(total_cec_production * 100.0 / total_production, digits = 2)) percent"

    cet_achieved_ratio =
        round(total_cec_production / total_production; digits = 2) /
        get_rec_requirement(simulation)[iteration_year]

    # Replace energy_mkt_clearing(nothing, nothing, system, load_growth, zones, num_days, solver) with
    # energy_mkt_clearing(sys_UC, sys_ED, system, load_growth, zones, num_days, solver) to run SIIP production cost model

    # Create empty market prices struct
    market_prices = MarketPrices()

    set_energy_price!(market_prices, "realized-ed", energy_price_ed)
    set_energy_price!(market_prices, "realized-uc", energy_price_uc)
    set_energy_price!(market_prices, "realized-md", energy_price_md)

    set_reserve_price!(market_prices, "realized-ed", reserve_price_ed)
    set_reserve_price!(market_prices, "realized-uc", reserve_price_uc)
    set_reserve_price!(market_prices, "realized-md", reserve_price_md)

    if in(:Inertia, market_names)
        set_inertia_price!(market_prices, "realized", inertia_price)
    end

    ######## Capacity market clearing #####################################################################

    capacity_market_bool = false
    if in(:Capacity, market_names)
        capacity_market_bool = true
    end

    capacity_mkt_param_file = joinpath(simulation_dir, "markets_data", "Capacity.csv")
    peak_load = get_peak_load(simulation)[pcm_scenario]
    #capacity_annual_increment = load_growth

    capacity_mkt_params_all = read_data(capacity_mkt_param_file)
    capacity_mkt_params = capacity_mkt_params_all[1, :]

    # Load canonical season definitions (primary key: name -> months).
    # Absent file / disabled toggle both degrade to Dict("annual" => 1:12).
    seasons_file = get_capacity_seasons_file(simulation_dir)
    season_months = load_season_months(seasons_file)
    seasonal_capacity_market = get_seasonal_capacity_market(get_case(simulation))
    seasons = resolve_capacity_seasons(season_months, seasonal_capacity_market)

    # Index Capacity.csv rows by season (foreign key into capacity_seasons.csv).
    cap_params_by_season = index_capacity_params_by_season(capacity_mkt_params_all, seasons)

    introduction_year = capacity_mkt_params["introduction_year"]
    discontinuation_year = capacity_mkt_params["discontinuation_year"]
    capacity_year = iteration_year + capacity_forward_years - 1
    capacity_active =
        (capacity_year >= introduction_year && capacity_year < discontinuation_year) ? 1 : 0

    #average_load_growth = Statistics.mean(load_growth)

    delta_irm =
        get_delta_irm(get_resource_adequacy(simulation)[pcm_scenario], iteration_year)
    irm_scalar = get_irm_scalar(get_case(simulation))

    capacity_price_dict = Dict{String, AxisArrays.AxisArray{Float64, 1}}(
        season => AxisArrays.AxisArray(reshape([0.0], 1), [1]) for season in seasons
    )
    capacity_accepted_bids_dict =
        Dict{String, Dict{String, Float64}}(season => Dict("no_accepted_bids" => 0.0) for season in seasons)

    if in(:Capacity, market_names) &&
       iteration_year + capacity_forward_years - 1 <= simulation_years
        #system_peak_load = (1 + average_load_growth) ^ (capacity_forward_years) * peak_load
        system_peak_load = peak_load[iteration_year + capacity_forward_years - 1]
        capacity_active_bool = Bool(capacity_active * capacity_market_bool)

        for season in seasons
            # Build the season-specific supply curve.
            seasonal_supply_curve = Vector{Union{String, Float64}}[]
            for project in capacity_market_projects
                for product in get_products(project)
                    seasonal_supply_curve = update_capacity_supply_curve!(
                        seasonal_supply_curve,
                        product,
                        project,
                        pcm_scenario,
                        season,
                    )
                end
            end
            sort!(seasonal_supply_curve; by = x -> x[3])   # Sort by capacity bid

            seasonal_demand_curve = create_capacity_demand_curve(
                cap_params_by_season[season],
                system_peak_load,
                irm_scalar,
                delta_irm,
                capacity_active_bool,
            )

            capacity_price_dict[season], capacity_accepted_bids_dict[season] =
                capacity_market_clearing(seasonal_demand_curve, seasonal_supply_curve, solver)
        end
    end

    set_capacity_price!(market_prices, "realized", capacity_price_dict)

    # Phase 4d: realized-profit (`realized_profits_calculator.jl`) now consumes the full
    # season-keyed accepted-bids dict (season -> name -> fraction), so pass it through directly.
    capacity_accepted_bids = capacity_accepted_bids_dict

    # Phase 4e boundary: `save_realized_market_data` still writes single un-seasoned
    # `capacity_price::AxisArray` and `capacity_accepted_bids::Dict{String,Float64}`.
    # Bridge on the representative season ("annual" when present, else the first season)
    # until Phase 4e generalizes the HDF5 writer.
    representative_season = haskey(capacity_price_dict, "annual") ? "annual" : first(seasons)
    capacity_price = capacity_price_dict[representative_season]
    capacity_accepted_bids_flat = capacity_accepted_bids_dict[representative_season]

    ######### REC market clearing ############################################################################

    rec_market_bool = false
    if in(:REC, market_names)
        rec_market_bool = true
    end

    REC_mkt_params =
        read_data(joinpath(simulation_dir, "markets_data", "REC_$(rps_target)_RPS.csv"))
    pricecap_rec = REC_mkt_params.price_cap[1]
    rec_req = REC_mkt_params.rec_req[1] * rec_market_bool
    rec_annual_increment = REC_mkt_params.annual_increment[1] * rec_market_bool
    rec_non_binding_years = REC_mkt_params.non_binding_years[1] * rec_market_bool

    rec_price = AxisArrays.AxisArray(reshape([pricecap_rec], 1), [1])
    rec_accepted_bids = Dict{String, Float64}()

    total_demand = 0.0
    for z in zones
        for t in 1:num_hours
            #total_demand += zonal_load[z, t] * (1 + energy_annual_increment[z]) * hour_weight[t]
            total_demand += zonal_load[z, t] * hour_weight[t]
        end
    end

    total_clean_production = 0.0

    rec_supply_curve = Vector{Union{String, Float64}}[]

    total_storage_consumption = 0.0
    for project in existing_projects
        # Populate REC market supply curves
        clean_production = 0.0
        for product in get_products(project)
            rec_supply_curve = update_rec_supply_curve!(rec_supply_curve, product, project)
            clean_production += find_clean_energy_production(product, project)
            total_storage_consumption += find_storage_energy_consumption(product, project)
        end
        total_clean_production += clean_production
    end

    clean_energy_percentage = min(1.0, (total_clean_production / total_demand))
    #println(clean_energy_percentage)

    if in(:REC, market_names)
        if length(rec_supply_curve) >= 1
            rec_energy_requirment =
                total_demand * min(rec_req + (rec_annual_increment * iteration_year), 1)

            sort!(rec_supply_curve; by = x -> x[3])      # Sort REC supply curve by REC bid

            #rec_energy_requirment = min(total_clean_production, rec_energy_requirment)
            rec_price, rec_accepted_bids = rec_market_clearing_non_binding(
                rec_energy_requirment,
                pricecap_rec,
                rec_supply_curve,
                solver,
            )

            # if iteration_year <= rec_non_binding_years

            #     rec_energy_requirment = min(total_clean_production, rec_energy_requirment)
            #     #println(rec_energy_requirment)
            #     rec_price, rec_accepted_bids = rec_market_clearing_non_binding(rec_energy_requirment, pricecap_rec, rec_supply_curve, solver)
            # else
            #     total = 0
            #     for i in rec_supply_curve
            #         total += i[2]
            #     end
            #     #println(total)
            #     rec_energy_requirment = min(total_clean_production, rec_energy_requirment)
            #     #println(rec_energy_requirment)
            #     rec_price, rec_accepted_bids = rec_market_clearing_binding(rec_energy_requirment, pricecap_rec, rec_supply_curve, solver)
            # end

            set_rec_price!(market_prices, "realized", rec_price)
        end
    end

    ################# Write actual market clearing data ################################################
    @info "Saving expected market data for year $(iteration_year)"
    output_file =
        joinpath(results_dir, "realized_market_data", "year_$(iteration_year).h5")

    save_realized_market_data(output_file,
        capacity_price,
        energy_price_ed, energy_price_uc, energy_price_md,
        reserve_price_ed, reserve_price_uc, reserve_price_md,
        rec_price, inertia_price,
        capacity_factors_md, capacity_factors_uc, capacity_factors_ed,
        reserve_perc_md, reserve_perc_uc, reserve_perc_ed,
        capacity_accepted_bids_flat, rec_accepted_bids,
        inertia_perc, start_up_costs, shut_down_costs,
        energy_voll, energy_voll_uc, energy_voll_md,
        reserve_voll, reserve_voll_uc, reserve_voll_md,
        inertia_voll,
        rec_supply_curve, rec_energy_requirment, cet_achieved_ratio,
    )
    ################ Update realized load data and peak load #########################################
    num_scenarios = length(all_scenarios)

    sys_UC_list,
    data_dirs,
    investors_list,
    representative_periods_list,
    rep_period_intervals,
    cases,
    iteration_years,
    rolling_horizons,
    simulation_years_list, timeseries_data_dir_list = repeat_arguments(
        num_scenarios,
        deepcopy(sys_UC),
        simulation_dir,
        get_investors(simulation),
        get_rep_periods(simulation),
        get_rep_period_interval(simulation),
        case,
        iteration_year,
        get_rolling_horizon(case),
        get_total_horizon(case),
        joinpath(get_results_dir(simulation), "timeseries_data_files"),
    )
    @time Distributed.pmap(
        parallelize_ordc_construction,
        zip(
            all_scenarios,
            sys_UC_list,
            data_dirs,
            investors_list,
            representative_periods_list,
            rep_period_intervals,
            cases,
            iteration_years .+ step_size,
            rolling_horizons,
            simulation_years_list,
            timeseries_data_dir_list,
        ),
    )

    return market_prices,
    capacity_factors_md,
    capacity_factors_uc,
    capacity_factors_ed,
    reserve_perc_md,
    reserve_perc_uc,
    reserve_perc_ed,
    inertia_perc,
    capacity_accepted_bids,
    rec_accepted_bids,
    clean_energy_percentage,
    cet_achieved_ratio
end

function calculate_reserve_scaling_factor(simulation::AgentSimulation)
    # scaling factor based on PV+Wind capacity
    case = get_case(simulation)

    test_system_dir = get_sys_dir(case)
    existing_generator_data = DataFrames.DataFrame(
        CSV.File(joinpath(test_system_dir, "RTS_Data", "SourceData", "gen.csv")),
    )
    initial_pv_wind_capacity = sum(
        existing_generator_data[
            (existing_generator_data[
                :,
                "Unit Type",
            ] .== "WIND") .| (existing_generator_data[:, "Unit Type"] .== "PV"),
            "PMax MW",
        ],
    )

    current_portfolio = vcat(get_existing.(get_investors(simulation))...)

    current_pv_wind_capacity = 0.0

    for project in current_portfolio
        if get_type(get_tech(project)) == "WT" || get_type(get_tech(project)) == "PVe"
            current_pv_wind_capacity = current_pv_wind_capacity + get_maxcap(project)
        end
    end

    scaling_factor_non_ordc_reserves =
        (current_pv_wind_capacity - initial_pv_wind_capacity) / initial_pv_wind_capacity

    return scaling_factor_non_ordc_reserves
end

function reserve_ts_scaling(simulation::AgentSimulation,
    iteration_year::Int64, step_size::Int64)
    simulation_dir = get_data_dir(get_case(simulation))
    all_scenarios = String.(get_all_scenario_names(simulation_dir))
    results_dir = get_results_dir(simulation)
    timeseries_data_dir = joinpath(results_dir, "timeseries_data_files")

    reserve_products = split(
        read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[
            1,
            "all_products",
        ],
        "; ",
    )
    ordc_products = split(
        read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[
            1,
            "ordc_products",
        ],
        "; ",
    )
    non_ordc_products = filter(p -> !(p in ordc_products), reserve_products)

    if iteration_year <= get_total_horizon(get_case(simulation))-1
        for scenario in all_scenarios
            reserve_timeseries_data = Dict(
                r => read_data(
                    joinpath(
                        timeseries_data_dir,
                        scenario,
                        "sim_year_$(iteration_year + step_size)",
                        "Reserves",
                        "$(r).csv",
                    ),
                ) for r in non_ordc_products
            )
            rep_reserve_timeseries_data = Dict(
                r => read_data(
                    joinpath(
                        timeseries_data_dir,
                        scenario,
                        "sim_year_$(iteration_year + step_size)",
                        "Reserves",
                        "rep_$(r).csv",
                    ),
                ) for r in non_ordc_products
            )

            load_initial = read_data(
                joinpath(
                    timeseries_data_dir,
                    scenario,
                    "sim_year_1",
                    "Load",
                    "load.csv",
                ),
            )
            load_current = read_data(
                joinpath(
                    timeseries_data_dir,
                    scenario,
                    "sim_year_$(iteration_year + step_size)",
                    "Load",
                    "load.csv",
                ),
            )
            load_initial_total =
                sum(sum(eachcol(load_initial[:, Not(:Year, :Month, :Day, :Period)])))
            load_current_total =
                sum(sum(eachcol(load_current[:, Not(:Year, :Month, :Day, :Period)])))

            scaling_factor = (load_current_total - load_initial_total) / load_initial_total

            for product in non_ordc_products
                reserve_timeseries_data[product][:, product] =
                    reserve_timeseries_data[product][:, product] * (1 + scaling_factor)

                rep_reserve_timeseries_data[product][:, product] =
                    rep_reserve_timeseries_data[product][:, product] * (1 + scaling_factor)

                CSV.write(
                    joinpath(
                        timeseries_data_dir,
                        scenario,
                        "sim_year_$(iteration_year + step_size)",
                        "Reserves",
                        "$(product).csv",
                    ),
                    reserve_timeseries_data[product],
                )
                CSV.write(
                    joinpath(
                        timeseries_data_dir,
                        scenario,
                        "sim_year_$(iteration_year + step_size)",
                        "Reserves",
                        "rep_$(product).csv",
                    ),
                    rep_reserve_timeseries_data[product],
                )
            end
        end
    end
end
