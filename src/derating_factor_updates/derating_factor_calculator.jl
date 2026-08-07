"""
function for ifelse elementwise treatment
"""
function elementwise_ifelse(x, y)
    if x <= 0
        z = x
    else
        z = y
    end
    return z
end

"""
Returns the case-specific capacity season definition file path.
"""
function get_capacity_seasons_file(simulation_dir::String)
    primary = joinpath(simulation_dir, MARKETS_DATA_DIRNAME, CAPACITY_SEASONS_FILENAME)
    return primary
end

"""
Determines whether this run should execute in annual-only mode.
If `seasonal_capacity_market` is disabled in CaseDefinition, annual mode is forced
regardless of any season file content.
"""
function is_annual_only_seasons(
    season_months::Dict{String, Vector{Int64}},
    seasonal_capacity_market::Bool,
)
    if !seasonal_capacity_market
        return true
    end
    return length(season_months) == 1 && haskey(season_months, "annual")
end

"""
Initializes the output derating table for annual or seasonal mode.
In seasonal mode, creates one row per season and keeps existing column schema.
"""
function initialize_derating_output(
    derating_template::DataFrame,
    seasons::Vector{String},
    seasonal_mode::Bool,
)
    if !seasonal_mode
        return deepcopy(derating_template)
    end

    derating_factors = DataFrame()
    derating_factors[!, "season"] = seasons
    for col in names(derating_template)
        if col == "season"
            continue
        end
        col_type = Base.nonmissingtype(eltype(derating_template[!, col]))
        if col_type <: Number
            derating_factors[!, col] = zeros(Float64, length(seasons))
        else
            derating_factors[!, col] = fill(missing, length(seasons))
        end
    end
    return derating_factors
end

"""
Calculates raw (unscaled) capacity credit values for existing and new renewable generation
and storage using the top-N net-load-hour methodology, then writes them to
`derating_dict.csv`. Per-type scalars are applied later in `update_derating_factor!`.
"""
function calculate_derating_data(simulation::Union{AgentSimulation, AgentSimulationData},
    simulation_dir::String,
    scenario::String,
    iteration_year::Int64,
    active_projects::Vector{Project},
    marginal_cc::Bool,
    timeseries_data_dir::String)
    @info "Calculating derating data using top net load hour methodology - iteration year: $(iteration_year), scenario: $(scenario)"
    cap_mkt_params = read_data(joinpath(simulation_dir, "markets_data", "Capacity.csv"))

    seasons_file = get_capacity_seasons_file(simulation_dir)
    seasonal_capacity_market = get_seasonal_capacity_market(get_case(simulation))
    season_months = load_season_months(seasons_file)
    seasonal_mode = !is_annual_only_seasons(season_months, seasonal_capacity_market)
    if !seasonal_mode
        # Keep downstream loops uniform by using a single synthetic annual season.
        season_months = Dict{String, Vector{Int64}}("annual" => collect(1:12))
    end
    seasons = collect(keys(season_months))

    renewable_existing =
        filter(p -> typeof(p) == RenewableGenEMIS{Existing}, active_projects)
    renewable_options = filter(p -> typeof(p) == RenewableGenEMIS{Option}, active_projects)

    zones = unique(get_zone.(get_tech.(renewable_existing)))
    types = unique(get_type.(get_tech.(renewable_existing)))

    extract_year(str) = parse(Int, split(str, "_")[end])

    simulation_years = get_total_horizon(get_case(simulation))

    load_n_vg_data = vcat(
        [
            read_data(
                joinpath(
                    timeseries_data_dir,
                    scenario,
                    "sim_year_$(sim_year)",
                    "Net Load Data",
                    "load_n_vg_data_rt.csv",
                ),
            ) for sim_year in 1:simulation_years
        ]...,
    )
    availability_data =
        read_availability_df(timeseries_data_dir, scenario, simulation_years, "REAL_TIME")

    num_hours = DataFrames.nrow(load_n_vg_data)
    base_num_top_hours = cap_mkt_params.num_top_hours[1]

    month_vector = derive_month_vector(load_n_vg_data)

    derating_template = read_data(
        joinpath(
            simulation_dir,
            "markets_data",
            "derating_data",
            scenario,
            "derating_dict.csv",
        ),
    )
    derating_factors = initialize_derating_output(derating_template, seasons, seasonal_mode)
    season_row_index = Dict(season => idx for (idx, season) in enumerate(seasons))

    load_n_vg_cols = Set(names(load_n_vg_data))
    missing_cols =
        [get_name(g) for g in renewable_existing if !(get_name(g) in load_n_vg_cols)]
    if !isempty(missing_cols)
        @warn "calculate_derating_data (year=$iteration_year, scenario=$scenario): columns missing from net-load CSV: $missing_cols"
    end

    function calculate_average_storage_cc(
        stor_duration::Int64,
        peak_reductions_existing::Dict{Int64, Float64},
        init_CC::Dict{Int64, Float64},
        average_efficiency::Float64,
        net_load_df::DataFrame,
        num_hours::Int64,
        stor_buffer_minutes::Int64,
    )
        inc = 1 #set to 1 for system-wide, but will need to be replaced with zonal level array if/when convert to zonal

        peak_reduction = peak_reductions_existing[stor_duration]

        max_demands =
            repeat([maximum(net_load_df[:, "net_load"]) - peak_reduction], num_hours)

        batt_powers = repeat([peak_reduction], num_hours)

        poss_charges = min.(
            batt_powers .* average_efficiency,
            (max_demands - net_load_df[:, "net_load"]) .* average_efficiency,
        )

        necessary_discharges = (max_demands - net_load_df[:, "net_load"])

        poss_batt_changes = zeros(size(necessary_discharges)[1])
        for n in collect(1:1:size(necessary_discharges)[1])
            poss_batt_changes[n] =
                elementwise_ifelse(necessary_discharges[n], poss_charges[n])
        end

        batt_e_level = zeros((inc, num_hours))
        batt_e_level[1] = min(poss_batt_changes[1], 0)
        for n in collect(2:1:num_hours)
            batt_e_level[n] = batt_e_level[n - 1] + poss_batt_changes[n]
            batt_e_level[n] = min(batt_e_level[n], 0.0)
        end

        required_MWhs = -minimum(batt_e_level)

        # This line of code will implement a buffer on all storage duration
        # requirements, i.e. if the stor_buffer_minutes is set to 60 minutes
        # then a 2-hour peak would be served by a 3-hour device, a 3-hour peak
        # by a 4-hour device, etc.

        stor_buffer_hrs = stor_buffer_minutes / 60
        required_MWhs = required_MWhs + (batt_powers[1] * stor_buffer_hrs)[1]
        stor_CC = peak_reduction * stor_duration / required_MWhs
        stor_CC = min(stor_CC, 1.0)

        return stor_CC
    end

    function calculate_marginal_storage_cc(
        stor_duration::Int64,
        peak_reductions_existing::Dict{Int64, Float64},
        peak_reduction_new::Float64,
        init_CC::Union{Dict{Any, Any}, Dict{Int64, Float64}},
        average_efficiency::Float64,
        net_load_df::DataFrame,
        num_hours::Int64,
        stor_buffer_minutes::Int64,
    )
        inc = 1 #set to 1 for system-wide, but will need to be replaced with zonal level array if/when convert to zon

        existing_peak_reduction = sum(values(peak_reductions_existing))

        max_demands = repeat(
            [maximum(net_load_df[:, "net_load"]) - existing_peak_reduction],
            num_hours,
        )

        batt_powers = repeat([peak_reduction_new], num_hours)

        poss_charges = min.(
            batt_powers .* average_efficiency,
            (max_demands - net_load_df[:, "net_load"]) .* average_efficiency,
        )

        necessary_discharges = (max_demands - net_load_df[:, "net_load"])

        poss_batt_changes = zeros(size(necessary_discharges)[1])
        for n in collect(1:1:size(necessary_discharges)[1])
            poss_batt_changes[n] =
                elementwise_ifelse(necessary_discharges[n], poss_charges[n])
        end

        batt_e_level = zeros((inc, num_hours))
        batt_e_level[1] = min(poss_batt_changes[1], 0)
        for n in collect(2:1:num_hours)
            batt_e_level[n] = batt_e_level[n - 1] + poss_batt_changes[n]
            batt_e_level[n] = min(batt_e_level[n], 0.0)
        end

        required_MWhs = -minimum(batt_e_level)

        # This line of code will implement a buffer on all storage duration
        # requirements, i.e. if the stor_buffer_minutes is set to 60 minutes
        # then a 2-hour peak would be served by a 3-hour device, a 3-hour peak
        # by a 4-hour device, etc.

        stor_buffer_hrs = stor_buffer_minutes / 60
        required_MWhs = required_MWhs + (batt_powers[1] * stor_buffer_hrs)[1]
        stor_CC = peak_reduction_new * stor_duration / required_MWhs
        stor_CC = min(stor_CC, 1.0)

        return stor_CC
    end

    for season in seasons
        season_month_set = Set(season_months[season])
        # Filter full-year data to the subset of rows that belong to this season.
        season_mask = map(m -> m in season_month_set, month_vector)
        season_hours = count(season_mask)
        if season_hours == 0
            @warn "No rows found for season '$season' while calculating derating factors."
            continue
        end

        season_num_top_hours =
            seasonal_mode ? max(1, Int(round(base_num_top_hours * simulation_years * length(season_months[season]) / 12))) :
            Int(base_num_top_hours * simulation_years)

        season_load_n_vg_data = load_n_vg_data[season_mask, :]
        season_availability_data = availability_data[season_mask, :]

        existing_vg_power = zeros(season_hours)
        load = vec(sum(Matrix(season_load_n_vg_data[:, r"load"]); dims = 2))
        for g in renewable_existing
            existing_vg_power += season_load_n_vg_data[!, get_name(g)]
        end

        net_load_df = season_load_n_vg_data[:, 1:4]
        net_load_df[:, "net_load"] = load - existing_vg_power
        net_load_sorted_df = deepcopy(DataFrames.sort(net_load_df, "net_load"; rev = true))

        type_zone_max_cap = Dict{String, Float64}()
        for zone in zones
            for type in types
                type_zone_id = "$(type)_$(zone)"
                type_zone_max_cap[type_zone_id] = 0.0
                net_load_df[:, "net_load_w/o_existing_$(type_zone_id)"] =
                    deepcopy(net_load_df[:, "net_load"])
                for g in renewable_existing
                    gen_name = get_name(g)
                    tech = get_tech(g)
                    if "$(get_type(tech))_$(get_zone(tech))" == type_zone_id
                        net_load_df[:, "net_load_w/o_existing_$(type_zone_id)"] +=
                            season_load_n_vg_data[:, gen_name]
                        type_zone_max_cap[type_zone_id] += get_maxcap(g)
                    end
                end
            end
        end

        for zone in zones
            for type in types
                type_zone_id = "$(type)_$(zone)"
                if type_zone_max_cap[type_zone_id] <= 0
                    continue
                end

                gen_sorted_df = deepcopy(
                    DataFrames.sort(
                        net_load_df,
                        "net_load_w/o_existing_$(type_zone_id)";
                        rev = true,
                    ),
                )

                top_hours = min(season_num_top_hours, DataFrames.nrow(gen_sorted_df))
                if top_hours == 0
                    continue
                end

                load_reduction =
                    gen_sorted_df[1:top_hours, "net_load_w/o_existing_$(type_zone_id)"] -
                    gen_sorted_df[1:top_hours, "net_load"]
                derating_factors[season_row_index[season], "existing_$(type_zone_id)"] = min(
                    sum(load_reduction) / type_zone_max_cap[type_zone_id] /
                    top_hours,
                    1.0,
                )
            end
        end

        for g in renewable_options
            gen_name = get_name(g)
            tech = get_tech(g)
            type_zone_id = "$(get_type(tech))_$(get_zone(tech))"
            gen_cap = get_maxcap(g)
            if !(type_zone_id in names(season_availability_data))
                @warn "Availability data missing column $(type_zone_id) for new renewable option derating."
                continue
            end

            net_load_df[:, "net_load_with_$(gen_name)"] = deepcopy(
                net_load_df[:, "net_load"] -
                season_availability_data[:, "$(type_zone_id)"] * gen_cap,
            )
            gen_sorted_df =
                deepcopy(DataFrames.sort(net_load_df, "net_load_with_$(gen_name)"; rev = true))
            top_hours = min(season_num_top_hours, DataFrames.nrow(gen_sorted_df))
            if top_hours == 0
                continue
            end

            load_reduction =
                net_load_sorted_df[1:top_hours, "net_load"] -
                gen_sorted_df[1:top_hours, "net_load_with_$(gen_name)"]

            derating_factors[season_row_index[season], "new_$(type_zone_id)"] =
                min(sum(load_reduction) / gen_cap / top_hours, 1.0)
        end

        # Storage CC script
        stor_buffer_minutes = cap_mkt_params.stor_buffer_minutes[1]
        all_battery_existing = filter(p -> typeof(p) == BatteryEMIS{Existing}, active_projects)
        all_battery_options = filter(p -> typeof(p) == BatteryEMIS{Option}, active_projects)

        # Define a dictionary to store batteries with their corresponding storage durations
        existing_storage_duration_dict = Dict{Int, Vector{BatteryEMIS{Existing}}}()
        option_storage_duration_dict = Dict{Int, Vector{BatteryEMIS{Option}}}()

        # Iterate over each existing battery
        for battery in all_battery_existing
            stor_duration =
                Int(round(get_storage_capacity(get_tech(battery))[:max] / get_maxcap(battery)))
            # @info "Existing Battery $(get_name(battery)) has storage duration of $(stor_duration) hours: storage capacity $(get_storage_capacity(get_tech(battery))[:max]) MWh, max cap $(get_maxcap(battery)) MW"
            if haskey(existing_storage_duration_dict, stor_duration)
                push!(existing_storage_duration_dict[stor_duration], battery)
            else
                existing_storage_duration_dict[stor_duration] = [battery]
            end
        end

        # Iterate over each option battery
        for battery in all_battery_options
            stor_duration =
                Int(round(get_storage_capacity(get_tech(battery))[:max] / get_maxcap(battery)))
            # @info "Option Battery $(get_name(battery)) has storage duration of $(stor_duration) hours: storage capacity $(get_storage_capacity(get_tech(battery))[:max]) MWh, max cap $(get_maxcap(battery)) MW"
            if haskey(option_storage_duration_dict, stor_duration)
                push!(option_storage_duration_dict[stor_duration], battery)
            else
                option_storage_duration_dict[stor_duration] = [battery]
            end
        end

        peak_reductions_existing = Dict(
            sd => sum(get_maxcap.(existing_storage_duration_dict[sd])) for
            sd in keys(existing_storage_duration_dict)
        )
        init_CC = Dict(
            sd => ("STOR_$sd" in names(derating_factors) ? derating_factors[season_row_index[season], "STOR_$sd"] : 0.0) for
            sd in keys(existing_storage_duration_dict)
        )
        season_num_hours = DataFrames.nrow(net_load_df)

        for (stor_duration, battery_existing) in existing_storage_duration_dict
            efficiencies = get_efficiency.(get_tech.(battery_existing))
            average_efficiency =
                (
                    mean([eff.in for eff in efficiencies]) +
                    mean([eff.out for eff in efficiencies])
                ) / 2

            stor_CC = calculate_average_storage_cc(
                stor_duration,
                peak_reductions_existing,
                init_CC,
                average_efficiency,
                net_load_df,
                season_num_hours,
                stor_buffer_minutes,
            )

            derating_factors[season_row_index[season], "existing_STOR_$(stor_duration)"] =
                stor_CC
        end

        if marginal_cc
            for (stor_duration, battery_option) in option_storage_duration_dict
                efficiencies = get_efficiency.(get_tech.(battery_option))
                average_efficiency =
                    (
                        mean([eff.in for eff in efficiencies]) +
                        mean([eff.out for eff in efficiencies])
                    ) / 2
                peak_reduction_new = sum(get_maxcap.(battery_option))
                stor_CC = calculate_marginal_storage_cc(
                    stor_duration,
                    peak_reductions_existing,
                    peak_reduction_new,
                    init_CC,
                    average_efficiency,
                    net_load_df,
                    season_num_hours,
                    stor_buffer_minutes,
                )
                derating_factors[season_row_index[season], "new_STOR_$(stor_duration)"] =
                    stor_CC
            end

        else
            for (stor_duration, battery_option) in option_storage_duration_dict
                if "existing_STOR_$(stor_duration)" in names(derating_factors)
                    derating_factors[season_row_index[season], "new_STOR_$(stor_duration)"] =
                        derating_factors[season_row_index[season], "existing_STOR_$(stor_duration)"]
                else
                    derating_factors[season_row_index[season], "new_STOR_$(stor_duration)"] =
                        derating_factors[season_row_index[season], "STOR_$(stor_duration)"]
                end
            end
        end
    end

    write_data(
        joinpath(simulation_dir, "markets_data", "derating_data", scenario),
        "derating_dict.csv",
        derating_factors,
    )
    return
end

"""
Deepcopy base_system, remove all projects in projects_to_remove via
remove_system_component!, and return the resulting PRAS.SystemModel.
The intermediate PSY system is freed at function return, reducing peak memory
before any subsequent PRAS.assess call.
"""
function build_pruned_pras_system(
    base_system::PSY.System,
    projects_to_remove::AbstractVector{<:Project},
)::PRAS.SystemModel
    pruned_sys = deepcopy(base_system)
    for project in projects_to_remove
        remove_system_component!(pruned_sys, project)
    end
    return make_pras_system_spi(
        pruned_sys,
        PSY.Area,
        nothing;
        copper_plate = false,
        copy_system = false,
    )
end

"""
Deepcopy base_system, add all projects in projects_to_add via
add_capacity_market_project!, and return the resulting PRAS.SystemModel.
The caller is responsible for pre-configuring projects (deepcopy, set_name!, etc.)
before passing them. The intermediate PSY system is freed at function return,
reducing peak memory before any subsequent PRAS.assess call.
"""
function build_augmented_pras_system(
    base_system::PSY.System,
    projects_to_add::AbstractVector{<:Project},
    simulation_dir::String,
    scenario::String,
    capacity_market_year::Int64,
    rt_resolution,
    simulation_years,
    timeseries_data_dir::String,
    availability_df_rt::DataFrame,
)::PRAS.SystemModel
    augmented_sys = deepcopy(base_system)
    for project in projects_to_add
        add_capacity_market_project!(
            augmented_sys,
            project,
            simulation_dir,
            scenario,
            capacity_market_year,
            rt_resolution,
            simulation_years,
            timeseries_data_dir,
            availability_df_rt,
        )
    end
    return make_pras_system_spi(
        augmented_sys,
        PSY.Area,
        nothing;
        copper_plate = false,
        copy_system = false,
    )
end

"""
Calculates raw (unscaled) capacity credit values for existing and new renewable generation
and storage using PRAS (ELCC or EFC methodology), then writes them to `derating_dict.csv`.
Per-type scalars are applied later in `update_derating_factor!`.
"""

function calculate_derating_factors(
    simulation::Union{AgentSimulation, AgentSimulationData},
    scenario::String,
    iteration_year::Int64,
    methodology::String,
    ra_metric::String,
    marginal_cc::Bool,
    timeseries_data_dir::String)
    if methodology == "ELCC"
        methodology = PRAS.ELCC
    elseif methodology == "EFC"
        methodology = PRAS.EFC
    else
        @error "Capacity Accreditation methodology should be either ELCC, EFC or TopNetLoad"
    end

    if ra_metric == "LOLE"
        ra_metric = PRAS.LOLE
    elseif ra_metric == "EUE"
        ra_metric = PRAS.EUE
    else
        @error "Resource Adequacy metric should be either LOLE or EUE"
    end

    simulation_dir = get_data_dir(get_case(simulation))
    simulation_years = get_total_horizon(get_case(simulation))
    outage_dir = get_outage_dir(get_case(simulation))
    rt_resolution = get_rt_resolution(get_case(simulation))
    zones = get_zones(simulation)

    availability_df_rt =
        get_availability_df(timeseries_data_dir, scenario, simulation_years, "REAL_TIME")

    derating_factors = read_data(
        joinpath(
            simulation_dir,
            "markets_data",
            "derating_data",
            scenario,
            "derating_dict.csv",
        ),
    )

    active_projects = get_activeprojects(simulation)
    existing = filter(p -> typeof(p) == RenewableGenEMIS{Existing}, active_projects)
    options = filter(p -> typeof(p) == RenewableGenEMIS{Option}, active_projects)
    existing_types = unique(get_type.(get_tech.(existing)))
    new_types = unique(get_type.(get_tech.(options)))

    capacity_forward_years = get_capacity_forward_years(simulation)
    capacity_market_year = iteration_year + capacity_forward_years - 1
    resource_adequacy = get_resource_adequacy(simulation)
    sys_PRAS = get_system_PRAS(simulation)[scenario]

    # No deepcopy needed here: create_base_system -> create_capacity_mkt_system performs
    # deepcopy(initial_system) internally, so copying here would be redundant.
    base_sys = sys_PRAS

    # create adjusted base system (by iteratively adding or removing generators) such that it meets the RA targets
    adjusted_base_system = create_base_system(sys_PRAS,
        active_projects,
        capacity_forward_years,
        scenario,
        resource_adequacy[scenario],
        iteration_year,
        simulation_dir,
        outage_dir,
        rt_resolution,
        simulation,
        availability_df_rt,
    )

    # create "Base" PRAS system to be used for calculation of ELCC or EFC.
    base_pras_system = make_pras_system_spi(
        adjusted_base_system,
        PSY.Area,
        nothing;
        copper_plate = false,
        copy_system = false,
    )

    # Compute regional load shares once; reused in all PRAS assess calls below.
    regional_load_shares = collect(get_regional_load_shares(base_pras_system))

    if marginal_cc
        for zone in zones
            for type in new_types
                @info "Adding to capacity market: $(type)_$(zone)"
                idx = findfirst(
                    x -> (
                        (get_type(get_tech(x)) == type) && (get_zone(get_tech(x)) == zone)
                    ),
                    options,
                )
                if !isnothing(idx)
                    build_size = 4 # set to 4 considering that there are 4 investors, so if a project is viable, there could be 4 such units coming online together.
                    max_cap = get_maxcap(options[idx]) * build_size
                    new_projects = Project[]
                    for i in 1:build_size
                        p = deepcopy(options[idx])
                        set_name!(p, "$(get_name(p))_$i")
                        push!(new_projects, p)
                    end
                    augmented_pras_system = build_augmented_pras_system(
                        adjusted_base_system,
                        new_projects,
                        simulation_dir,
                        scenario,
                        capacity_market_year,
                        rt_resolution,
                        simulation_years,
                        timeseries_data_dir,
                        availability_df_rt,
                    )

                    # Call PRAS accreditation methodology. Adjust sample size, seed, etc. here.
                    cc_result = PRAS.assess(
                        base_pras_system,
                        augmented_pras_system,
                        methodology{ra_metric}(Int(ceil(max_cap)), regional_load_shares),
                        PRAS.SequentialMonteCarlo(;
                            samples = PRAS_N_SAMPLES,
                            seed = PRAS_MONTE_CARLO_SEED,
                        ),
                    )
                    cc_lower, cc_upper = extrema(cc_result)
                    cc_final = (cc_lower + cc_upper) / (2 * max_cap)
                    derating_factors[!, "new_$(type)_$(zone)"] .= cc_final
                end
            end
        end
    end

    # For average ELCC/EFC, existing units are removed. The new system with reduced units now becomes the base PRAS system.
    # No deepcopy needed here: SPI.generate_pras_system only reads the PSY system to build a
    # PRAS struct and does not mutate it. The resulting augmented_pras_system is a fresh object.
    augmented_pras_system = make_pras_system_spi(
        adjusted_base_system,
        PSY.Area,
        nothing;
        copper_plate = false,
        copy_system = false,
    )

    for zone in zones
        for type in existing_types
            zone_tech_units = existing[findall(
                x -> ((get_type(get_tech(x)) == type) && (get_zone(get_tech(x)) == zone)),
                existing,
            )]
            if !isempty(zone_tech_units)
                total_capacity = sum(get_maxcap.(zone_tech_units))
                @assert total_capacity > 0
                pruned_base_pras_system =
                    build_pruned_pras_system(adjusted_base_system, zone_tech_units)
                #  Call PRAS accreditation methodology. Adjust sample size, seed, etc. here.
                cc_result = PRAS.assess(
                    pruned_base_pras_system,
                    augmented_pras_system,
                    PRAS.ELCC{ra_metric}(Int(ceil(total_capacity)), regional_load_shares),
                    PRAS.SequentialMonteCarlo(;
                        samples = PRAS_N_SAMPLES,
                        seed = PRAS_MONTE_CARLO_SEED,
                    ),
                )
                cc_lower, cc_upper = extrema(cc_result)
                cc_final = (cc_lower + cc_upper) / (2 * total_capacity)

                derating_factors[!, "existing_$(type)_$(zone)"] .= cc_final
            end
        end
    end

    all_battery_existing = filter(p -> typeof(p) == BatteryEMIS{Existing}, active_projects)
    all_battery_options = filter(p -> typeof(p) == BatteryEMIS{Option}, active_projects)

    # Define a dictionary to store batteries with their corresponding storage durations
    existing_storage_duration_dict = Dict{Int, Vector{BatteryEMIS{Existing}}}()
    option_storage_duration_dict = Dict{Int, Vector{BatteryEMIS{Option}}}()

    # Iterate over each existing battery
    for battery in all_battery_existing
        stor_duration =
            Int(round(get_storage_capacity(get_tech(battery))[:max] / get_maxcap(battery)))
        @info "calculate_derating_factors - Existing Battery $(get_name(battery)) has storage duration of $(stor_duration) hours"
        if haskey(existing_storage_duration_dict, stor_duration)
            push!(existing_storage_duration_dict[stor_duration], battery)
        else
            existing_storage_duration_dict[stor_duration] = [battery]
        end
    end

    # Iterate over each option battery
    for battery in all_battery_options
        stor_duration =
            Int(round(get_storage_capacity(get_tech(battery))[:max] / get_maxcap(battery)))
        @info "calculate_derating_factors - Option Battery $(get_name(battery)) has storage duration of $(stor_duration) hours"
        if haskey(option_storage_duration_dict, stor_duration)
            push!(option_storage_duration_dict[stor_duration], battery)
        else
            option_storage_duration_dict[stor_duration] = [battery]
        end
    end

    for (stor_duration, battery_existing) in existing_storage_duration_dict
        total_capacity = sum(get_maxcap.(battery_existing))
        pruned_base_pras_system =
            build_pruned_pras_system(adjusted_base_system, battery_existing)

        # Call PRAS accreditation methodology. Adjust sample size, seed, etc. here.
        cc_result = PRAS.assess(
            pruned_base_pras_system,
            augmented_pras_system,
            PRAS.ELCC{ra_metric}(Int(ceil(total_capacity)), regional_load_shares),
            PRAS.SequentialMonteCarlo(;
                samples = PRAS_N_SAMPLES,
                seed = PRAS_MONTE_CARLO_SEED,
            ),
        )
        cc_lower, cc_upper = extrema(cc_result)
        cc_final = (cc_lower + cc_upper) / (2 * total_capacity)
        derating_factors[!, "existing_STOR_$(stor_duration)"] .= cc_final
    end

    # augmented_pras_system is no longer needed after the existing storage loop above.
    # Release it before the battery marginal CC block to reduce peak memory.
    augmented_pras_system = nothing

    if marginal_cc
        new_project_names = []
        max_cap = 0.0
        for (stor_duration, battery_options) in option_storage_duration_dict
            new_projects = Project[]
            for project in battery_options
                new_project = deepcopy(project)
                project_name = get_name(new_project)
                if !(project_name in new_project_names)
                    push!(new_project_names, project_name)
                    max_cap += get_maxcap(project)
                    push!(new_projects, new_project)
                end
            end
            augmented_pras_system = build_augmented_pras_system(
                adjusted_base_system,
                new_projects,
                simulation_dir,
                scenario,
                capacity_market_year,
                rt_resolution,
                simulation_years,
                timeseries_data_dir,
                availability_df_rt,
            )

            # Call PRAS accreditation methodology. Adjust sample size, seed, etc. here.
            cc_result = PRAS.assess(
                base_pras_system,
                augmented_pras_system,
                methodology{ra_metric}(Int(ceil(max_cap)), regional_load_shares),
                PRAS.SequentialMonteCarlo(;
                    samples = PRAS_N_SAMPLES,
                    seed = PRAS_MONTE_CARLO_SEED,
                ),
            )
            cc_lower, cc_upper = extrema(cc_result)
            cc_final = (cc_lower + cc_upper) / (2 * max_cap)
            derating_factors[!, "new_STOR_$(stor_duration)"] .= cc_final
        end

    else
        for (stor_duration, battery_options) in option_storage_duration_dict
            if "existing_STOR_$(stor_duration)" in names(derating_factors)
                derating_factors[:, "new_STOR_$(stor_duration)"] .=
                    derating_factors[:, "existing_STOR_$(stor_duration)"]
            else
                derating_factors[:, "new_STOR_$(stor_duration)"] .=
                    derating_factors[:, "STOR_$(stor_duration)"]
            end
        end
    end

    # Overwrite file with new derating factors.
    write_data(
        joinpath(simulation_dir, "markets_data", "derating_data", scenario),
        "derating_dict.csv",
        derating_factors,
    )
end

"""
Reads the per-type capacity credit scalar from `CC_SCALAR_FILENAME` for the given
scenario. Returns the value in column `type_key` (row 1) if the file and column exist,
otherwise returns `1.0` so that missing entries are a no-op.
"""
function read_cc_scalar(simulation_dir::String, scenario::String, type_key::String)::Float64
    filepath = joinpath(
        simulation_dir,
        "markets_data",
        "derating_data",
        scenario,
        CC_SCALAR_FILENAME,
    )
    if !isfile(filepath)
        return 1.0
    end
    df = read_data(filepath)
    if nrow(df) == 0 || !(type_key in names(df))
        return 1.0
    end
    val = df[1, type_key]
    return ismissing(val) ? 1.0 : Float64(val)
end

"""
This function does nothing is project is not of ThermalGenEMIS, HydroGenEMIS, RenewableGenEMIS or BatteryEMIS type.
"""
function update_derating_factor!(project::P,
    simulation_dir::String,
    scenario::String,
    marginal_cc::Bool,
) where {P <: Project{<:BuildPhase}}
    return
end

"""
Updates the derating factors of ThermalGenEMIS and HydroGenEMIS projects.
Reads the raw derating factor from `derating_dict.csv` and multiplies by the
per-type scalar from `cc_scalar.csv` (defaults to 1.0 if absent).
"""
function update_derating_factor!(
    project::Union{ThermalGenEMIS{<:BuildPhase}, HydroGenEMIS{<:BuildPhase}},
    simulation_dir::String,
    scenario::String,
    marginal_cc::Bool,
)
    derating_data = read_data(
        joinpath(
            simulation_dir,
            "markets_data",
            "derating_data",
            scenario,
            "derating_dict.csv",
        ),
    )
    seasonal = "season" in names(derating_data)
    type_key = get_type(get_tech(project))
    scalar = read_cc_scalar(simulation_dir, scenario, get_type(get_tech(project)))
    for row in eachrow(derating_data)
        season = seasonal ? String(row["season"]) : "annual"
        derating_factor = row[type_key]
        for product in get_products(project)
            set_derating!(product, scenario, season, derating_factor * scalar)
        end
    end
    return
end

"""
Updates the derating factors of existing RenewableGenEMIS projects.
Reads the raw CC from `derating_dict.csv` (written unscaled by the calculate step)
and multiplies by the per-type scalar from `cc_scalar.csv` (defaults to 1.0).
"""
function update_derating_factor!(project::RenewableGenEMIS{Existing},
    simulation_dir::String,
    scenario::String,
    marginal_cc::Bool,
)
    derating_data = read_data(
        joinpath(
            simulation_dir,
            "markets_data",
            "derating_data",
            scenario,
            "derating_dict.csv",
        ),
    )
    name = get_name(project)
    tech = get_tech(project)
    type_zone_id = "$(get_type(tech))_$(get_zone(tech))"

    if !in("existing_$(type_zone_id)", names(derating_data))
        error("Derating data not found")
    end

    seasonal = "season" in names(derating_data)
    scalar = read_cc_scalar(simulation_dir, scenario, get_type(tech))
    for row in eachrow(derating_data)
        season = seasonal ? String(row["season"]) : "annual"
        derating_factor = row["existing_$(type_zone_id)"]
        for product in get_products(project)
            set_derating!(product, scenario, season, derating_factor * scalar)
        end
    end

    return
end

"""
Updates the derating factors of new/option RenewableGenEMIS projects.
Reads the raw marginal or average CC from `derating_dict.csv` and multiplies by
the per-type scalar from `cc_scalar.csv` (defaults to 1.0).
"""
function update_derating_factor!(project::RenewableGenEMIS{<:BuildPhase},
    simulation_dir::String,
    scenario::String,
    marginal_cc::Bool,
)
    derating_data = read_data(
        joinpath(
            simulation_dir,
            "markets_data",
            "derating_data",
            scenario,
            "derating_dict.csv",
        ),
    )
    name = get_name(project)
    tech = get_tech(project)
    type_zone_id = "$(get_type(tech))_$(get_zone(tech))"

    col_name = marginal_cc ? "new_$(type_zone_id)" : "existing_$(type_zone_id)"
    if !in(col_name, names(derating_data))
        error("Derating data not found")
    end

    seasonal = "season" in names(derating_data)
    scalar = read_cc_scalar(simulation_dir, scenario, get_type(tech))
    for row in eachrow(derating_data)
        season = seasonal ? String(row["season"]) : "annual"
        derating_factor = row[col_name]
        for product in get_products(project)
            set_derating!(product, scenario, season, derating_factor * scalar)
        end
    end

    return
end

"""
Updates the derating factors of existing BatteryEMIS projects.
Reads the raw CC from `derating_dict.csv` (written unscaled) and multiplies by the
duration-based scalar `STOR_N` from `cc_scalar.csv` (defaults to 1.0). No cap
at 1.0 — scalars above 1.0 are supported.
"""
function update_derating_factor!(project::BatteryEMIS{Existing},
    simulation_dir::String,
    scenario::String,
    marginal_cc::Bool,
)
    tech = get_tech(project)
    duration = Int(round(get_storage_capacity(tech)[:max] / get_maxcap(project)))
    project_type = "existing_STOR_$(duration)"

    # @info "Existing battery derating factor update: $(get_name(project)), storage_capacity $(get_storage_capacity(tech)[:max]), max_cap $(get_maxcap(project)), duration of $(duration) hours"

    derating_data = read_data(
        joinpath(
            simulation_dir,
            "markets_data",
            "derating_data",
            scenario,
            "derating_dict.csv",
        ),
    )

    seasonal = "season" in names(derating_data)
    scalar = read_cc_scalar(simulation_dir, scenario, "STOR_$(duration)")
    for row in eachrow(derating_data)
        season = seasonal ? String(row["season"]) : "annual"
        derating_factor = row[project_type]
        for product in get_products(project)
            set_derating!(product, scenario, season, derating_factor * scalar)
        end
    end
    return
end

"""
Updates the derating factors of new/option BatteryEMIS projects.
Reads the raw CC from `derating_dict.csv` (written unscaled) and multiplies by the
duration-based scalar `STOR_N` from `cc_scalar.csv` (defaults to 1.0). No cap
at 1.0 — scalars above 1.0 are supported.
"""
function update_derating_factor!(project::BatteryEMIS{<:BuildPhase},
    simulation_dir::String,
    scenario::String,
    marginal_cc::Bool,
)
    tech = get_tech(project)
    duration = Int(round(get_storage_capacity(tech)[:max] / get_maxcap(project)))

    # @info "BuildPhase battery derating factor update: $(get_name(project)), storage_capacity $(get_storage_capacity(tech)[:max]), max_cap $(get_maxcap(project)), duration of $(duration) hours"

    if marginal_cc
        project_type = "new_STOR_$(duration)"
    else
        project_type = "existing_STOR_$(duration)"
    end

    derating_data = read_data(
        joinpath(
            simulation_dir,
            "markets_data",
            "derating_data",
            scenario,
            "derating_dict.csv",
        ),
    )
    seasonal = "season" in names(derating_data)
    scalar = read_cc_scalar(simulation_dir, scenario, "STOR_$(duration)")
    for row in eachrow(derating_data)
        season = seasonal ? String(row["season"]) : "annual"
        derating_factor = row[project_type]
        for product in get_products(project)
            set_derating!(product, scenario, season, derating_factor * scalar)
        end
    end
    return
end

"""
Writes raw (unscaled) capacity credit values to `derating_dict.csv` for a single
scenario and iteration year. Dispatches to `calculate_derating_data` (TopNetLoad)
or `calculate_derating_factors` (ELCC/EFC). Per-type scalars from `cc_scalar.csv`
are applied separately by the caller via `update_derating_factor!`.
"""
function update_simulation_derating_data!(
    simulation::Union{AgentSimulation, AgentSimulationData},
    scenario::String,
    iteration_year::Int64,
    timeseries_data_dir::String;
    methodology::String = "ELCC",
    ra_metric::String = "LOLE",
    marginal_cc::Bool = true)
    @info "Updating derating factors for scenario $(scenario) and iteration year $(iteration_year) using methodology $(methodology) and RA metric $(ra_metric). Marginal CC is set to $(marginal_cc)."
    data_dir = get_data_dir(get_case(simulation))
    active_projects = get_activeprojects(simulation)

    if methodology == "TopNetLoad"
        calculate_derating_data(
            simulation,
            data_dir,
            scenario,
            iteration_year,
            active_projects,
            marginal_cc,
            timeseries_data_dir,
        )
    else
        calculate_derating_factors(
            simulation,
            scenario,
            iteration_year,
            methodology,
            ra_metric,
            marginal_cc,
            timeseries_data_dir,
        )
    end

    return
end

