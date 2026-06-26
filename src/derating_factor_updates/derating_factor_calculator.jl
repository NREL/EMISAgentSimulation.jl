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
This function calculates the derating data for existing and new renewable generation
based on top 100 net-load hour methodology.
"""
function calculate_derating_data(simulation::Union{AgentSimulation, AgentSimulationData},
    simulation_dir::String,
    scenario::String,
    iteration_year::Int64,
    active_projects::Vector{Project},
    derating_scale::Float64,
    marginal_cc::Bool)
    @info "Calculating derating data using top net load hour methodology - iteration year: $(iteration_year), scenario: $(scenario)"
    cap_mkt_params = read_data(joinpath(simulation_dir, "markets_data", "Capacity.csv"))

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
                    simulation_dir,
                    "timeseries_data_files",
                    scenario,
                    "sim_year_$(sim_year)",
                    "Net Load Data",
                    "load_n_vg_data_rt.csv",
                ),
            ) for sim_year in simulation_years
        ]...,
    )
    availability_data = vcat(
        [
            read_data(
                joinpath(
                    simulation_dir,
                    "timeseries_data_files",
                    scenario,
                    "sim_year_$(sim_year)",
                    "Availability",
                    "REAL_TIME_availability.csv",
                ),
            ) for sim_year in simulation_years
        ]...,
    )

    num_hours = DataFrames.nrow(load_n_vg_data)
    num_top_hours = cap_mkt_params.num_top_hours[1] * simulation_years

    existing_vg_power = zeros(num_hours)

    load = vec(sum(Matrix(load_n_vg_data[:, r"load"]); dims = 2))

    for g in renewable_existing
        existing_vg_power += load_n_vg_data[!, get_name(g)]
    end

    net_load_df = load_n_vg_data[:, 1:4]
    net_load_df[:, "net_load"] = load - existing_vg_power

    net_load_sorted_df = deepcopy(DataFrames.sort(net_load_df, "net_load"; rev = true))

    derating_factors = read_data(
        joinpath(
            simulation_dir,
            "markets_data",
            "derating_data",
            scenario,
            "derating_dict.csv",
        ),
    )

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
                        load_n_vg_data[:, gen_name]
                    type_zone_max_cap[type_zone_id] += get_maxcap(g)
                end
            end
        end
    end

    for zone in zones
        for type in types
            type_zone_id = "$(type)_$(zone)"
            gen_sorted_df = deepcopy(
                DataFrames.sort(
                    net_load_df,
                    "net_load_w/o_existing_$(type_zone_id)";
                    rev = true,
                ),
            )

            load_reduction =
                gen_sorted_df[1:num_top_hours, "net_load_w/o_existing_$(type_zone_id)"] -
                gen_sorted_df[1:num_top_hours, "net_load"]
            derating_factors[:, "existing_$(type_zone_id)"] .= min(
                sum(load_reduction) * derating_scale / type_zone_max_cap[type_zone_id] /
                num_top_hours,
                1.0,
            )
        end
    end

    for g in renewable_options
        gen_name = get_name(g)
        tech = get_tech(g)
        type_zone_id = "$(get_type(tech))_$(get_zone(tech))"
        gen_cap = get_maxcap(g)

        net_load_df[:, "net_load_with_$(gen_name)"] = deepcopy(
            net_load_df[:, "net_load"] - availability_data[:, "$(type_zone_id)"] * gen_cap,
        )
        gen_sorted_df =
            deepcopy(DataFrames.sort(net_load_df, "net_load_with_$(gen_name)"; rev = true))

        load_reduction =
            net_load_sorted_df[1:num_top_hours, "net_load"] -
            gen_sorted_df[1:num_top_hours, "net_load_with_$(gen_name)"]

        derating_factors[:, "new_$(type_zone_id)"] .=
            min(sum(load_reduction) * derating_scale / gen_cap / num_top_hours, 1.0)
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
        sd => derating_factors[1, "STOR_$sd"] for
        sd in keys(existing_storage_duration_dict)
    )
    # @info "initial CC is $(init_CC)"
    # @info "peak reductions is $(peak_reductions_existing)"

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
            num_hours,
            stor_buffer_minutes,
        )

        derating_factors[:, "existing_STOR_$(stor_duration)"] .= stor_CC
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
                num_hours,
                stor_buffer_minutes,
            )
            derating_factors[:, "new_STOR_$(stor_duration)"] .= stor_CC
        end
    else
        for (stor_duration, battery_option) in option_storage_duration_dict
            if "existing_STOR_$(stor_duration)" in names(derating_factors)
                derating_factors[:, "new_STOR_$(stor_duration)"] .=
                    derating_factors[:, "existing_STOR_$(stor_duration)"]
            else
                derating_factors[:, "new_STOR_$(stor_duration)"] .=
                    derating_factors[:, "STOR_$(stor_duration)"]
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
    return SPI.generate_pras_system(pruned_sys, PSY.Area, false)
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
        )
    end
    return SPI.generate_pras_system(augmented_sys, PSY.Area, false)
end

"""
This function calculates the derating data for existing and new renewable generation
and storage based on PRAS outcomes.
"""

function calculate_derating_factors(
    simulation::Union{AgentSimulation, AgentSimulationData},
    scenario::String,
    iteration_year::Int64,
    derating_scale::Float64,
    methodology::String,
    ra_metric::String,
    marginal_cc::Bool)
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
    adjusted_base_system = create_base_system(base_sys,
        active_projects,
        capacity_forward_years,
        scenario,
        resource_adequacy[scenario],
        iteration_year,
        simulation_dir,
        outage_dir,
        rt_resolution,
        simulation)

    # create "Base" PRAS system to be used for calculation of ELCC or EFC.
    base_pras_system = SPI.generate_pras_system(adjusted_base_system,
        PSY.Area,
        false)

    # Compute regional load shares once; reused in all PRAS assess calls below.
    regional_load_shares = collect(get_regional_load_shares(base_pras_system))

    ##TODO: AA remove debug code after validation
    temp_dir = "/projects/gmlcmarkets/Phase2_EMIS_Analysis/GS_AAYAD/HPC_Analysis_Runs/20250310_no_sdes_High_RECT_Static_ORDC_RA_Cap_wo_md_storff_High_RPS/temp_data"
    @info "Debug: Saving PRAS system for scenario $(scenario) and iteration year $(iteration_year) to $(temp_dir) for debugging purposes."
    PSY.to_json(base_pras_system, joinpath(temp_dir, "base_pras_system_scenario_$(scenario)_year_$(iteration_year).json"))
    PSY.to_json(adjusted_base_system, joinpath(temp_dir, "adjusted_base_system_scenario_$(scenario)_year_$(iteration_year).json"))

    if marginal_cc
        for zone in zones
            for type in new_types
                @info "$(type)_$(zone)"
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
                    cc_final = (cc_lower + cc_upper) * derating_scale / (2 * max_cap)
                    derating_factors[!, "new_$(type)_$(zone)"] .= cc_final
                end
            end
        end
    end

    # For average ELCC/EFC, existing units are removed. The new system with reduced units now becomes the base PRAS system.
    # No deepcopy needed here: SPI.generate_pras_system only reads the PSY system to build a
    # PRAS struct and does not mutate it. The resulting augmented_pras_system is a fresh object.
    augmented_sys = adjusted_base_system
    augmented_pras_system = SPI.generate_pras_system(augmented_sys,
        PSY.Area,
        false)

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
                cc_final = (cc_lower + cc_upper) * derating_scale / (2 * total_capacity)

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
        cc_final = (cc_lower + cc_upper) * derating_scale / (2 * total_capacity)
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
            cc_final = (cc_lower + cc_upper) * derating_scale / (2 * max_cap)
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
This function does nothing is project is not of ThermalGenEMIS, HydroGenEMIS, RenewableGenEMIS or BatteryEMIS type.
"""
function update_derating_factor!(project::P,
    simulation_dir::String,
    scenario::String,
    derating_scale::Float64,
    marginal_cc::Bool,
) where {P <: Project{<:BuildPhase}}
    return
end

"""
This function updates the derating factors of ThermalGenEMIS and HydroGenEMIS projects.
"""
function update_derating_factor!(
    project::Union{ThermalGenEMIS{<:BuildPhase}, HydroGenEMIS{<:BuildPhase}},
    simulation_dir::String,
    scenario::String,
    derating_scale::Float64,
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
    derating_factor = derating_data[1, get_type(get_tech(project))]
    for product in get_products(project)
        set_derating!(product, scenario, derating_factor)
    end
    return
end

"""
This function updates the derating factors of existing RenewableGenEMIS projects.
"""
function update_derating_factor!(project::RenewableGenEMIS{Existing},
    simulation_dir::String,
    scenario::String,
    derating_scale::Float64,
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

    if in("existing_$(type_zone_id)", names(derating_data))
        derating_factor = derating_data[1, "existing_$(type_zone_id)"]
    else
        error("Derating data not found")
    end

    for product in get_products(project)
        set_derating!(product, scenario, derating_factor)
    end

    return
end

"""
This function updates the derating factors of new RenewableGenEMIS projects.
"""
function update_derating_factor!(project::RenewableGenEMIS{<:BuildPhase},
    simulation_dir::String,
    scenario::String,
    derating_scale::Float64,
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

    if marginal_cc
        if in("new_$(type_zone_id)", names(derating_data))
            derating_factor = derating_data[1, "new_$(type_zone_id)"]
        else
            error("Derating data not found")
        end
    else
        if in("existing_$(type_zone_id)", names(derating_data))
            derating_factor = derating_data[1, "existing_$(type_zone_id)"]
        else
            error("Derating data not found")
        end
    end

    for product in get_products(project)
        set_derating!(product, scenario, derating_factor)
    end

    return
end

"""
This function updates the derating factors of Existing BatteryEMIS projects.
"""
function update_derating_factor!(project::BatteryEMIS{Existing},
    simulation_dir::String,
    scenario::String,
    derating_scale::Float64,
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

    derating_factor = derating_data[1, project_type]
    derating_factor = min(derating_factor * derating_scale, 1.0)
    for product in get_products(project)
        set_derating!(product, scenario, derating_factor)
    end
    return
end

"""
This function updates the derating factors of BatteryEMIS projects.
"""
function update_derating_factor!(project::BatteryEMIS{<:BuildPhase},
    simulation_dir::String,
    scenario::String,
    derating_scale::Float64,
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
    derating_factor = derating_data[1, project_type]
    derating_factor = min(derating_factor * derating_scale, 1.0)
    for product in get_products(project)
        set_derating!(product, scenario, derating_factor)
    end
    return
end

"""
This function updates the derating factors of all active projects in the simulation.
"""
function update_simulation_derating_data!(
    simulation::Union{AgentSimulation, AgentSimulationData},
    scenario::String,
    iteration_year::Int64,
    derating_scale::Float64;
    methodology::String = "ELCC",
    ra_metric::String = "LOLE",
    marginal_cc::Bool = true)
    @info "Updating derating factors for scenario $(scenario) and iteration year $(iteration_year) using methodology $(methodology) and RA metric $(ra_metric). Marginal CC is set to $(marginal_cc). Derating scale is set to $(derating_scale)."
    data_dir = get_data_dir(get_case(simulation))
    active_projects = get_activeprojects(simulation)

    if methodology == "TopNetLoad"
        calculate_derating_data(
            simulation,
            data_dir,
            scenario,
            iteration_year,
            active_projects,
            derating_scale,
            marginal_cc,
        )
    else
        calculate_derating_factors(
            simulation,
            scenario,
            iteration_year,
            derating_scale,
            methodology,
            ra_metric,
            marginal_cc,
        )
    end

    return
end

