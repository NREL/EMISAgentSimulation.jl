function prune_system_devices!(
    sys::PSY.System,
    prune_dict::Dict{Type{<:PSY.Component}, Array{AbstractString}},
)
    for (type, device_names) in prune_dict
        for name in device_names
            device = PSY.get_component(type, sys, name)
            PSY.remove_component!(sys, device)
        end
    end
    return
end

"""
Use this function to specify which devices to remove from the RTS at the start of the simulation.
"""
function specify_pruned_units()
    pruned_unit = Dict{Type{<:PSY.Component}, Array{AbstractString}}()
    pruned_unit[PSY.ThermalStandard] =
        ["115_STEAM_1", "115_STEAM_2", "315_STEAM_1", "315_STEAM_2", "315_STEAM_3",
            "315_STEAM_4", "315_STEAM_5",
            "101_CT_1", "101_CT_2", "102_CT_1", "102_CT_2",
            "201_CT_1", "201_CT_2", "202_CT_1", "202_CT_2",
            "301_CT_1", "301_CT_2", "302_CT_1", "302_CT_2",
            "207_CT_1", "307_CT_1", "101_STEAM_4",
            "123_STEAM_3", "223_STEAM_1", "223_STEAM_3"]
    pruned_unit[PSY.RenewableNonDispatch] =
        ["308_RTPV_1", "313_RTPV_1", "313_RTPV_2", "313_RTPV_3", "313_RTPV_4", "313_RTPV_5",
            "313_RTPV_6", "313_RTPV_7",
            "313_RTPV_8", "313_RTPV_9", "313_RTPV_10", "313_RTPV_11", "313_RTPV_12",
            "320_RTPV_1", "320_RTPV_2", "320_RTPV_3",
            "313_RTPV_13", "320_RTPV_4", "320_RTPV_5", "118_RTPV_1", "118_RTPV_2",
            "118_RTPV_3", "118_RTPV_4", "118_RTPV_5",
            "118_RTPV_6", "320_RTPV_6", "118_RTPV_7", "118_RTPV_8", "118_RTPV_9",
            "118_RTPV_10", "213_RTPV_1"]
    pruned_unit[PSY.RenewableDispatch] = ["309_WIND_1", "212_CSP_1"]
    pruned_unit[PSY.Generator] = ["114_SYNC_COND_1", "314_SYNC_COND_1", "214_SYNC_COND_1"]
    pruned_unit[PSY.EnergyReservoirStorage] = ["313_STORAGE_1"]

    return pruned_unit
end

function create_rts_sys(rts_dir::String,
    base_power::Float64,
    simulation_dir::String,
    scratch_dir::String,
    timeseries_data_dir::String,
    scenarios::Vector{String},
    pcm_scenario::String,
    simulation_years::Int64,
    da_resolution::Int64,
    rt_resolution::Int64,
    MD_horizon::Int64,
    MD_interval::Int64,
    UC_horizon::Int64,
    UC_interval::Int64,
    ED_horizon::Int64,
    ED_interval::Int64,
    outage_dir::Union{Nothing, String} = nothing,
)
    system_config = load_system_config(simulation_dir)
    default_rts_load = system_config.default_rts_load
    ntp_ts_data_dir = joinpath(timeseries_data_dir, "input_processing")
    # ntp_ts_data_dir = POINTER_FILE[:NTPS_TS_DATA_DIR]
    runchecks = false

    ##TODO: revert back to original system after checking storage capacities
    initial_Sienna_system_name = "DA_sys_zonal_with_storage_capacities.json"
    # initial_Sienna_system_name = "DA_sys_zonal.json"

    sys_MDs = Vector{PSY.System}()
    sys_UCs = Vector{PSY.System}()
    sys_EDs = Vector{PSY.System}()
    sys_EDs_dict = Dict(scenario => Vector{PSY.System}() for scenario in scenarios)

    supercc_scenario = get_scenario_pcm_label(system_config, pcm_scenario)

    for sim_year in 1:simulation_years
        MD_sys_filename = joinpath(
            rts_dir,
            "constructed_systems",
            pcm_scenario,
            "sim_year_$(sim_year)",
            "MD_sys_EMIS_$(MD_horizon)hor_$(MD_interval)int.json",
        )
        MD_num_forecast_filename = joinpath(
            rts_dir,
            "constructed_systems",
            pcm_scenario,
            "sim_year_$(sim_year)",
            "MD_num_forecast_$(MD_horizon)hor_$(MD_interval)int.txt",
        )

        loadyear = DEFAULT_LOAD_YEAR + sim_year
        weatheryear = loadyear

        if !(isfile(MD_sys_filename) && isfile(MD_num_forecast_filename))
            @info "Simulation year $(sim_year): MD json file doesn't exist, creating PSY system."
            dir_exists(dirname(MD_sys_filename))
            sys_MD_initial = PSY.System(
                joinpath(rts_dir, initial_Sienna_system_name);
                time_series_directory = scratch_dir,
            );
            # create MD system
            create_sys_w_updated_ts(
                ntp_ts_data_dir,
                sys_MD_initial,
                weatheryear,
                loadyear,
                "dayahead",
                supercc_scenario,
                default_rts_load, # GW
                MD_horizon, # hours
                MD_interval, # hours
                MD_sys_filename,
                true,
                nothing,
                nothing,
                MD_num_forecast_filename,
                outage_dir,
            )
        end
        sys_MD = PSY.System(
            MD_sys_filename;
            time_series_directory = scratch_dir,
            runchecks = runchecks,
        )
        fix_multistart_cost_curves!(sys_MD)
        push!(sys_MDs, sys_MD);

        UC_filename = joinpath(
            rts_dir,
            "constructed_systems",
            pcm_scenario,
            "sim_year_$(sim_year)",
            "DA_sys_EMIS_$(UC_horizon)hor_$(UC_interval)int_$(MD_horizon)mdhor_$(MD_interval)mdint.json",
        )
        if !(isfile(UC_filename))
            @info "Simulation year $(sim_year): UC json file doesn't exist, creating PSY system."
            dir_exists(dirname(UC_filename))
            sys_UC_initial = PSY.System(
                joinpath(rts_dir, initial_Sienna_system_name);
                time_series_directory = scratch_dir,
            );
            create_sys_w_updated_ts(
                ntp_ts_data_dir,
                sys_UC_initial,
                weatheryear,
                loadyear,
                "dayahead",
                supercc_scenario,
                default_rts_load, # GW
                UC_horizon, # hours
                UC_interval, # hours
                UC_filename,
                false,
                MD_horizon,
                MD_interval,
                MD_num_forecast_filename,
                outage_dir,
            )
        end

        sys_UC = PSY.System(
            UC_filename;
            time_series_directory = scratch_dir,
            runchecks = runchecks,
        )
        fix_multistart_cost_curves!(sys_UC)
        push!(sys_UCs, sys_UC);

        for scenario in scenarios
            supercc_scenario_ed = get_scenario_pcm_label(system_config, scenario)

            ED_filename = joinpath(
                rts_dir,
                "constructed_systems",
                scenario,
                "sim_year_$(sim_year)",
                "RT_sys_EMIS_$(ED_horizon)hor_$(ED_interval)int_$(MD_horizon)mdhor_$(MD_interval)mdint.json",
            )
            if !isfile(ED_filename)
                @info "Simulation year $(sim_year): ED json file doesn't exist, creating PSY system."
                dir_exists(dirname(ED_filename))
                sys_ED_initial = PSY.System(
                    joinpath(rts_dir, initial_Sienna_system_name);
                    time_series_directory = scratch_dir,
                );
                create_sys_w_updated_ts(
                    ntp_ts_data_dir,
                    sys_ED_initial,
                    weatheryear,
                    loadyear,
                    "realtime",
                    supercc_scenario_ed,
                    default_rts_load, # GW
                    ED_horizon, # hours
                    ED_interval, # hours
                    ED_filename,
                    false,
                    MD_horizon,
                    MD_interval,
                    MD_num_forecast_filename,
                    outage_dir,
                )
            end

            sys_ED = PSY.System(
                ED_filename;
                time_series_directory = scratch_dir,
                runchecks = runchecks,
            )
            fix_multistart_cost_curves!(sys_ED)
            push!(sys_EDs_dict[scenario], sys_ED);
        end
        sys_EDs = sys_EDs_dict[pcm_scenario]

        # ERCOT specific items
        removegen_name = ["AUSTIN_1", "AUSTIN_2"]
        for sys in [sys_MDs[sim_year], sys_UCs[sim_year], sys_EDs[sim_year]]
            for d in PSY.get_components(PSY.Generator, sys)
                if d.name in removegen_name
                    PSY.remove_component!(sys, d)
                end
            end
        end

        for sys in [sys_MDs[sim_year], sys_UCs[sim_year], sys_EDs[sim_year]]
            d = PSY.get_component(PSY.VariableReserve, sys, "SPIN")
            PSY.remove_component!(sys, d)
        end

        for sys in [sys_MDs[sim_year], sys_UCs[sim_year], sys_EDs[sim_year]]
            d = PSY.get_component(PSY.VariableReserveNonSpinning, sys, "NONSPIN")
            PSY.remove_component!(sys, d)
        end

        for sys in [sys_MDs[sim_year], sys_UCs[sim_year], sys_EDs[sim_year]]
            d = PSY.get_component(PSY.VariableReserve, sys, "REG_DN")
            PSY.set_name!(sys, d, "Reg_Down")
        end

        for sys in [sys_MDs[sim_year], sys_UCs[sim_year], sys_EDs[sim_year]]
            d = PSY.get_component(PSY.VariableReserve, sys, "REG_UP")
            PSY.set_name!(sys, d, "Reg_Up")
        end

        PSY.set_units_base_system!(sys_MDs[sim_year], PSY.IS.UnitSystem.DEVICE_BASE)
        PSY.set_units_base_system!(sys_UCs[sim_year], PSY.IS.UnitSystem.DEVICE_BASE)
        PSY.set_units_base_system!(sys_EDs[sim_year], PSY.IS.UnitSystem.DEVICE_BASE)
    end

    sys_PRAS = Dict{String, PSY.System}()
    for scenario in scenarios
        PRAS_filename = joinpath(
            rts_dir,
            "constructed_systems",
            scenario,
            "sim_year_1",
            "PRAS_sys_EMIS_$(ED_horizon)hor_$(ED_interval)int_$(MD_horizon)mdhor_$(MD_interval)mdint.json",
        )
        if !isfile(PRAS_filename)
            @info "Scenario $(scenario): PRAS json file doesn't exist, creating PSY system."
            create_PRAS_sys_json(sys_EDs_dict[scenario], PRAS_filename, outage_dir)
        end
        sys_PRAS[scenario] = PSY.System(
            PRAS_filename;
            time_series_directory = scratch_dir,
            runchecks = runchecks,
        )
    end

    return sys_MDs,
    sys_UCs,
    sys_EDs,
    sys_PRAS,
    MD_horizon,
    MD_interval,
    UC_horizon,
    UC_interval,
    ED_horizon,
    ED_interval
end

function remove_vre_gens!(sys::PSY.System)
    for gen in get_all_techs(sys)
        if typeof(gen) == PSY.RenewableDispatch
            PSY.remove_component!(sys, gen)
        end
    end
end

function add_outages_to_system!(
    sys::PSY.System,
    outages_dir::Union{Nothing, String},
    dates,
)
    if outages_dir === nothing
        @info "No outages CSV provided — using nominal outage values for all generators."
    end

    for gen in PSY.get_components(PSY.Generator, sys)
        add_nominal_outage_to_component!(sys, gen)
    end
    for storage in PSY.get_components(PSY.Storage, sys)
        add_nominal_outage_to_component!(sys, storage)
    end

    attach_outage_data_from_csv!(
        sys,
        outages_dir;
        mttr_hours = DEFAULT_THERMAL_MTTR_HOURS,
        timestamps = dates,
    )

    return
end

function create_sys_w_updated_ts(
    data_dir::String,
    initial_sys::PSY.System,
    weatheryear::Int64,
    loadyear::Int64,
    market_stage::String, # dayahead, realtime
    scenario::String,
    loadscaler_base::Float64, # GW
    horizon::Int64, # hours
    interval::Int64, # hours
    output_file::String,
    first_stage::Bool = false,
    first_stage_horizon::Union{Nothing, Integer} = nothing, # hour (only need input if first_stage is false)
    first_stage_interval::Union{Nothing, Integer} = nothing, # hour (only need input if first_stage is false)
    first_stage_number_of_forecast_filename::Union{Nothing, String} = nothing, # (only need input if first_stage is false)
    outages_dir::Union{Nothing, String} = nothing,
)

    #--------------------------------------------
    # Calculate load scaling factor: scale 2021 load to 75 GW
    #--------------------------------------------
    loadscaler_profile_rt = DataFrame(
        CSV.File(
            joinpath(
                data_dir,
                "load_actuals_processed",
                "sup3rcc_ecearth3_load_gwh_ercot_$(scenario)_e2021_w2021_cst.csv",
            ),
        ),
    ) #in GW
    loadscaler_profile_da = DataFrame(
        CSV.File(
            joinpath(
                data_dir,
                "load_forecasts_processed",
                "preds_20210101_$(scenario)_365days.csv",
            ),
        ),
    )
    loadscaler_peak_da =
        maximum(sum(eachcol(select(loadscaler_profile_da, Not([:Column1])))))
    loadscaler_peak_rt =
        maximum(sum(eachcol(select(loadscaler_profile_rt, Not([:year, :timestamp])))))
    loadscaler_da = loadscaler_peak_da/loadscaler_base
    loadscaler_rt = loadscaler_peak_rt/loadscaler_base

    sys_MD = initial_sys
    PSY.set_units_base_system!(sys_MD, PSY.IS.UnitSystem.NATURAL_UNITS)

    # Remove the generator that availability is false
    # for sys in [sys_MD]
    # 	d=PSY.get_component(Generator, sys, "Glove Solar")
    # 	PSY.remove_component!(sys, d)
    # end
    # 
    # PSY.get_available(PSY.get_component(Generator, sys_DA, "Glove Solar"))
    # PSY.get_components(PSY.Generator, sys_DA, x -> x.available == false)

    #-----------------------------------------------------------------------------------
    # Construct hydro timeseries with new horizon from existing Deterministic timeseries
    #-----------------------------------------------------------------------------------
    # first_ts_temp_MD = first(PSY.get_time_series_multiple(sys_MD))
    # start_datetime_MD = PSY.IS.get_initial_timestamp(first_ts_temp_MD)
    start_datetime_MD = PSY.get_forecast_initial_timestamp(sys_MD)
    sys_MD_res = PSY.get_time_series_resolutions(sys_MD)[1]
    sys_MD_initial_interval =
        Int(PSY.get_forecast_interval(sys_MD).value / sys_MD_res.value)

    # if not the first stage, I think finish_datetime_MD needs to be first stage's number_of_forecast * interval + (horizon - interval) -- all from first stage (need to pass down these parameters)
    if first_stage == true
        finish_datetime_MD =
            start_datetime_MD + Dates.Hour((DEFAULT_HOURS_PER_YEAR - 1) * sys_MD_res)
    else
        first_stage_number_of_forecast = 0
        open(first_stage_number_of_forecast_filename, "r") do file
            first_stage_number_of_forecast = parse(Int, readline(file))
        end
        finish_datetime_MD =
            start_datetime_MD + Dates.Hour(
                (
                    first_stage_number_of_forecast * first_stage_interval +
                    (first_stage_horizon - first_stage_interval) - 1
                ) * sys_MD_res,
            )
    end

    # timestep here indicate how many MD periods are being constructed
    timestep = StepRange(start_datetime_MD, sys_MD_res * interval, finish_datetime_MD);
    first_stage_total = Int((finish_datetime_MD - start_datetime_MD) / sys_MD_res + 1)
    additional_timestep = Int(
        horizon - (first_stage_total-(interval*(length(timestep)-1))) +
        (first_stage_total - DEFAULT_HOURS_PER_YEAR),
    )
    dates = range(
        SIM_START_DATE;
        step = sys_MD_res,
        length = DEFAULT_HOURS_PER_YEAR + additional_timestep,
    )
    # @info "Line 294 - Constructing time series with timestep length: $(length(timestep)) and additional timestep: $(additional_timestep)."
    for component in
        get_components(x -> has_time_series(x, Deterministic), HydroDispatch, sys_MD)
        forecast = get_time_series(Deterministic, component, "max_active_power")

        reconstruct_single_ts = Float64[]
        for (key, value) in forecast.data
            append!(reconstruct_single_ts, value[1:sys_MD_initial_interval])
        end
        append!(reconstruct_single_ts, reconstruct_single_ts[1:additional_timestep])
        data = TS.TimeArray(dates, reconstruct_single_ts)
        single_time_series = SingleTimeSeries("max_active_power", data)

        add_time_series!(sys_MD, component, single_time_series)
    end

    remove_time_series!(sys_MD, Deterministic)
    remove_time_series!(sys_MD, DeterministicSingleTimeSeries)

    for component in
        get_components(x -> has_time_series(x, SingleTimeSeries), HydroDispatch, sys_MD)
        revisedts = DataStructures.SortedDict{DateTime, Vector{Float64}}()
        newtsdata =
            values(get_time_series(SingleTimeSeries, component, "max_active_power").data)

        for t in 1:length(timestep)
            rtseries=Vector{Float64}()
            datetimeindex = timestep[t]
            rtseries = newtsdata[(interval * (t - 1) + 1):(interval * (t - 1) + horizon)]
            push!(revisedts, datetimeindex => rtseries)
        end

        # convert to deterministic time series
        revisedts_deterministic = PSY.Deterministic(;
            name = "max_active_power",
            data = revisedts,
            resolution = Dates.Hour(1),
            scaling_factor_multiplier = get_max_active_power,
        )

        add_time_series!(sys_MD, component, revisedts_deterministic)
    end

    remove_time_series!(sys_MD, SingleTimeSeries)

    #-----------------------------------------------------------
    # Replacing wind & solar time series  !! In NATURAL_UNITS !!
    #-----------------------------------------------------------
    namemapping = DataFrame(CSV.File(joinpath(data_dir, "GeneratorNameMapping.csv")))
    # To get raw DA data time stamps
    # first_ts_temp_MD = first(PSY.get_time_series_multiple(sys_MD))
    # start_datetime_MD = PSY.IS.get_initial_timestamp(first_ts_temp_MD);
    # sys_MD_res = PSY.get_time_series_resolution(sys_MD)
    # finish_datetime_MD = start_datetime_MD + Dates.Hour(8759*sys_MD_res);
    # timestep here indicate how many MD periods are being constructed
    # timestep = StepRange(start_datetime_MD, sys_MD_res*interval, finish_datetime_MD);
    # hourlytimestep  = StepRange(start_datetime_DA, sys_DA_res, finish_datetime_DA);

    # remove_time_series!(sys_MD, Deterministic)

    for d in get_components(
        x -> get_prime_mover_type(x) in [PrimeMovers.PVe, PrimeMovers.WT],
        Generator,
        sys_MD,
    )
        # println("Processing generator: $(get_name(d))")
        #create dictionary
        # revisedts = Dict{DateTime, Array{Float64}}()
        revisedts = DataStructures.SortedDict{DateTime, Vector{Float64}}()

        # tstype = "1dayahead" ## actuals/, 1dayahead/, intraday/
        if market_stage == "dayahead"
            tstype = "1dayahead"
        elseif market_stage == "realtime"
            # TODO: the real-time timeseries actually uses both "actuals" and "intraday", but it limits the RT horizon to 2-hour.
            tstype = "actuals"
        end

        if get_prime_mover_type(d) == PrimeMovers.WT
            technology = "lbw"  # land-based wind 
        else
            technology = "upv"  # utility pv 
        end

        profile = DataFrame(
            CSV.File(
                joinpath(
                    data_dir,
                    "plant_profiles_processed",
                    "$(tstype)/ercot_$(technology)_build-$(weatheryear)_cst.csv",
                ),
            ),
        )
        newtsdata = profile[
            !,
            namemapping[in([PSY.get_name(d)]).(namemapping.jsonname), :csvname][1],
        ]
        basepower = get_rating(d)
        newtsdata = newtsdata ./ basepower

        for t in 1:length(timestep)
            rtseries = Vector{Float64}()
            datetimeindex = timestep[t]
            if t < DEFAULT_HOURS_PER_YEAR/interval
                rtseries =
                    newtsdata[(interval * (t - 1) + 1):(interval * (t - 1) + horizon)]
            elseif t == ceil(Int, DEFAULT_HOURS_PER_YEAR/interval)
                rtseries = [
                    newtsdata[(interval * (t - 1) + 1):DEFAULT_HOURS_PER_YEAR];
                    newtsdata[1:(horizon - (DEFAULT_HOURS_PER_YEAR - (interval * (t - 1))))]
                ]
            else
                # simply use the end of the previous timeseries as the new start
                new_start =
                    horizon-(
                        DEFAULT_HOURS_PER_YEAR-(
                            interval*(ceil(Int, DEFAULT_HOURS_PER_YEAR/interval)-1)
                        )
                    )
                rtseries = newtsdata[(new_start + (t - ceil(
                    Int,
                    DEFAULT_HOURS_PER_YEAR / interval,
                ) - 1) * interval + 1):(new_start + (t - ceil(
                    Int,
                    DEFAULT_HOURS_PER_YEAR / interval,
                ) - 1) * interval + horizon)]
            end
            push!(revisedts, datetimeindex => rtseries)
        end

        # convert to deterministic time series
        revisedts_deterministic = PSY.Deterministic(;
            name = "max_active_power",
            data = revisedts,
            resolution = Dates.Hour(1),
            scaling_factor_multiplier = get_max_active_power,
        )

        # remove old time series
        # remove_time_series!(sys_MD, Deterministic, d, "max_active_power")
        # add new time series to dataset
        add_time_series!(sys_MD, d, revisedts_deterministic)
    end

    #--------------------------------------------
    # Load Forecasts !!!!!!!!!!!!! CHECK SYSTEM BASE !!!!!!!!!
    #--------------------------------------------
    if market_stage == "dayahead"
        profile=DataFrame(
            CSV.File(
                joinpath(
                    data_dir,
                    "load_forecasts_processed",
                    "preds_$(loadyear)0101_$(scenario)_365days.csv",
                ),
            ),
        ) #in GW; original command: profile=DataFrame(CSV.File(joinpath(data_dir, "Load_Forecast_from_Local", "$(scenario)", "$(loadyear)", "preds_$(weatheryear)0101_365days.csv"))) #in GW
        loadscaler = loadscaler_da
    elseif market_stage == "realtime"
        profile=DataFrame(
            CSV.File(
                joinpath(
                    data_dir,
                    "load_actuals_processed",
                    "sup3rcc_ecearth3_load_gwh_ercot_$(scenario)_e$(loadyear)_w$(weatheryear)_cst.csv",
                ),
            ),
        ) #in GW
        loadscaler = loadscaler_rt
    end

    for d in get_components(PSY_LOADS, sys_MD)
        # println("Processing region: $(get_name(get_bus(d)))")
        #create dictionary
        # revisedts = Dict{DateTime, Array{Float64}}()
        bus = PSY.get_bus(d)
        zone = PSY.get_area(bus)
        revisedts = DataStructures.SortedDict{DateTime, Vector{Float64}}()

        newtsdata = profile[!, lowercase(PSY.get_name(zone))]
        baseload = get_max_active_power(d)

        # Get the loads in the zone and caculate the proportion that this load represents
        loads_in_zone =
            PSY.get_components(x -> PSY.get_area(PSY.get_bus(x)) == zone, PSY_LOADS, sys_MD)
        zonal_load_sum = sum(get_max_active_power.(loads_in_zone))
        proportion = baseload / zonal_load_sum

        newtsdata = newtsdata .* proportion
        newtsdata = newtsdata .* 1000
        newtsdata = newtsdata ./ loadscaler # scale 2021 to 75GW peak

        for t in 1:length(timestep)
            rtseries = Vector{Float64}()
            datetimeindex = timestep[t]
            if t < DEFAULT_HOURS_PER_YEAR/interval
                rtseries =
                    newtsdata[(interval * (t - 1) + 1):(interval * (t - 1) + horizon)]
            elseif t == ceil(Int, DEFAULT_HOURS_PER_YEAR/interval)
                rtseries = [
                    newtsdata[(interval * (t - 1) + 1):DEFAULT_HOURS_PER_YEAR];
                    newtsdata[1:(horizon - (DEFAULT_HOURS_PER_YEAR - (interval * (t - 1))))]
                ]
            else
                # simply use the end of the previous timeseries as the new start
                new_start =
                    horizon-(
                        DEFAULT_HOURS_PER_YEAR-(
                            interval*(ceil(Int, DEFAULT_HOURS_PER_YEAR/interval)-1)
                        )
                    )
                rtseries = newtsdata[(new_start + (t - ceil(
                    Int,
                    DEFAULT_HOURS_PER_YEAR / interval,
                ) - 1) * interval + 1):(new_start + (t - ceil(
                    Int,
                    DEFAULT_HOURS_PER_YEAR / interval,
                ) - 1) * interval + horizon)]
            end
            push!(revisedts, datetimeindex => rtseries)
        end

        # convert to deterministic time series
        revisedts_deterministic = PSY.Deterministic(;
            name = "max_active_power",
            data = revisedts,
            resolution = Dates.Hour(1),
            scaling_factor_multiplier = get_max_active_power,
        )

        # remove old time series
        # remove_time_series!(sys_DA, Deterministic, d, "max_active_power")
        # add new time series to dataset
        add_time_series!(sys_MD, d, revisedts_deterministic)
    end

    #-----------------------------------------------------------------
    # Regulation time series update
    #-----------------------------------------------------------------
    if market_stage == "dayahead"
        regdown_ts = DataFrame(
            CSV.File(
                joinpath(
                    data_dir,
                    "regulation_reserves",
                    scenario,
                    "DA_regDown_$(scenario)_gmlc$(loadyear).csv",
                ),
            ),
        ) #in MW; previous command: regdown_ts = DataFrame(CSV.File(joinpath(data_dir, "TS_for_Regulation_Req_Calc", "RegulationTS_Bethany", "DA_regDown_baseline_gmlc$(weatheryear).csv")))
        regup_ts = DataFrame(
            CSV.File(
                joinpath(
                    data_dir,
                    "regulation_reserves",
                    scenario,
                    "DA_regUp_$(scenario)_gmlc$(loadyear).csv",
                ),
            ),
        ) #in MW; previous command: regup_ts = DataFrame(CSV.File(joinpath(data_dir, "TS_for_Regulation_Req_Calc", "RegulationTS_Bethany", "DA_regUp_baseline_gmlc$(weatheryear).csv")))
    else
        regdown_ts = DataFrame(
            CSV.File(
                joinpath(
                    data_dir,
                    "regulation_reserves",
                    scenario,
                    "RT_regDown_$(scenario)_gmlc$(loadyear).csv",
                ),
            ),
        ) #in MW; previous command: regdown_ts = DataFrame(CSV.File(joinpath(data_dir, "TS_for_Regulation_Req_Calc", "RegulationTS_Bethany", "DA_regDown_baseline_gmlc$(weatheryear).csv")))
        regup_ts = DataFrame(
            CSV.File(
                joinpath(
                    data_dir,
                    "regulation_reserves",
                    scenario,
                    "RT_regUp_$(scenario)_gmlc$(loadyear).csv",
                ),
            ),
        ) #in MW; previous command: regup_ts = DataFrame(CSV.File(joinpath(data_dir, "TS_for_Regulation_Req_Calc", "RegulationTS_Bethany", "DA_regUp_baseline_gmlc$(weatheryear).csv")))
    end

    regdown_ts = select!(regdown_ts, Not(:DATETIME))
    regup_ts = select!(regup_ts, Not(:DATETIME))
    regdown_ts[!, "RegDown"] = sum(eachcol(regdown_ts))
    regup_ts[!, "RegUp"] = sum(eachcol(regup_ts))
    reg_profile =
        DataFrame(; REG_DN = regdown_ts[!, "RegDown"], REG_UP = regup_ts[!, "RegUp"])

    for d in get_components(x -> PSY.get_name(x) in ["REG_DN", "REG_UP"], Service, sys_MD)
        # println("Processing serve: $(get_name(d))")
        #create dictionary
        # revisedts = Dict{DateTime, Array{Float64}}()
        revisedts = DataStructures.SortedDict{DateTime, Vector{Float64}}()

        newtsdata = reg_profile[!, PSY.get_name(d)]
        basereq = get_requirement(d)
        newtsdata = newtsdata ./ basereq/100

        for t in 1:length(timestep)
            rtseries = Vector{Float64}()
            datetimeindex = timestep[t]
            if t < DEFAULT_HOURS_PER_YEAR/interval
                rtseries =
                    newtsdata[(interval * (t - 1) + 1):(interval * (t - 1) + horizon)]
            elseif t == ceil(Int, DEFAULT_HOURS_PER_YEAR/interval)
                rtseries = [
                    newtsdata[(interval * (t - 1) + 1):DEFAULT_HOURS_PER_YEAR];
                    newtsdata[1:(horizon - (DEFAULT_HOURS_PER_YEAR - (interval * (t - 1))))]
                ]
            else
                # simply use the end of the previous timeseries as the new start
                new_start =
                    horizon-(
                        DEFAULT_HOURS_PER_YEAR-(
                            interval*(ceil(Int, DEFAULT_HOURS_PER_YEAR/interval)-1)
                        )
                    )
                rtseries = newtsdata[(new_start + (t - ceil(
                    Int,
                    DEFAULT_HOURS_PER_YEAR / interval,
                ) - 1) * interval + 1):(new_start + (t - ceil(
                    Int,
                    DEFAULT_HOURS_PER_YEAR / interval,
                ) - 1) * interval + horizon)]
            end
            push!(revisedts, datetimeindex => rtseries)
        end

        # convert to deterministic time series
        revisedts_deterministic = PSY.Deterministic(;
            name = "requirement",
            data = revisedts,
            resolution = Dates.Hour(1),
            scaling_factor_multiplier = get_requirement,
        )

        # remove old time series
        # remove_time_series!(sys_MD, Deterministic, d, "requirement")
        # add new time series to dataset
        add_time_series!(sys_MD, d, revisedts_deterministic)
    end

    # Add outage supplemental attributes and time series
    add_outages_to_system!(sys_MD, outages_dir, dates)

    PSY.set_units_base_system!(sys_MD, PSY.IS.UnitSystem.SYSTEM_BASE)
    to_json(sys_MD, output_file; force = true)

    number_of_forecast = get_forecast_window_count(sys_MD)
    if first_stage
        open(first_stage_number_of_forecast_filename, "w") do file
            write(first_stage_number_of_forecast_filename, string(number_of_forecast))
        end
    end

    return
end

function create_sys_w_updated_ts_ny(
    data_dir::String,
    initial_sys::PSY.System,
    weatheryear::Int64,
    loadyear::Int64,
    market_stage::String, # dayahead, realtime
    scenario::String,
    loadscaler_base::Float64, #GW
    horizon::Int64, # hours
    interval::Int64, # hours
    output_file::String,
    first_stage::Bool = false,
    first_stage_horizon::Union{Nothing, Integer} = nothing, # hour (only need input if first_stage is false)
    first_stage_interval::Union{Nothing, Integer} = nothing, # hour (only need input if first_stage is false)
    first_stage_number_of_forecast_filename::Union{Nothing, String} = nothing, # (only need input if first_stage is false)
)
    sys_MD = initial_sys
    PSY.set_units_base_system!(sys_MD, PSY.IS.UnitSystem.SYSTEM_BASE)
    to_json(sys_MD, output_file; force = true)

    return
end

function create_PRAS_sys_json(
    sys_EDs::Vector{PSY.System},
    output_file::String,
    outages_dir::Union{Nothing, String} = nothing,
)
    sys_PRAS = deepcopy(first(sys_EDs))

    start_datetime_PRAS = PSY.get_forecast_initial_timestamp(sys_PRAS)
    sys_PRAS_res = PSY.get_time_series_resolutions(sys_PRAS)[1]
    sys_PRAS_interval = Int(PSY.get_forecast_interval(sys_PRAS).value / sys_PRAS_res.value)
    sys_PRAS_horizon = Int(PSY.get_forecast_horizon(sys_PRAS).value / sys_PRAS_res.value)
    finish_datetime_PRAS =
        start_datetime_PRAS +
        Dates.Hour(DEFAULT_HOURS_PER_YEAR * sys_PRAS_res * length(sys_EDs) - sys_PRAS_res)

    timestep = StepRange(
        start_datetime_PRAS,
        sys_PRAS_res * sys_PRAS_interval,
        finish_datetime_PRAS,
    );

    ts_objects = Dict{String, Any}()
    ts_objects["gens"] = collect(
        PSY.get_components(
            x ->
                has_time_series(x, Deterministic) && PSY.get_prime_mover_type(x) in
                [PrimeMovers.PVe, PrimeMovers.WT, PrimeMovers.HY],
            Generator,
            sys_PRAS,
        ),
    )
    ts_objects["loads"] = collect(
        PSY.get_components(x -> has_time_series(x, Deterministic), PSY_LOADS, sys_PRAS),
    )

    for (key, devices) in ts_objects
        for component_PRAS in devices
            name = PSY.get_name(component_PRAS)
            reconstruct_single_ts = Float64[]
            for sys in sys_EDs
                component = first(PSY.get_components_by_name(PSY.Device, sys, name))
                forecast = PSY.get_time_series(Deterministic, component, "max_active_power")
                weather_year = Dates.year(first(keys(forecast.data)))
                for (key, value) in forecast.data
                    if Dates.year(key) == weather_year
                        append!(reconstruct_single_ts, value[1:sys_PRAS_interval])
                    end
                end
            end
            # append!(reconstruct_single_ts, reconstruct_single_ts[1:sys_PRAS_horizon - 1])
            # dates = range(DateTime("2018-01-01T00:00:00"), step = sys_PRAS_res, length = length(timestep) + sys_PRAS_horizon - 1)
            dates = range(SIM_START_DATE; step = sys_PRAS_res, length = length(timestep))
            # @info "Line 633 - Reconstructing time series for component $(name) with length: $(length(reconstruct_single_ts))."
            data = TS.TimeArray(dates, reconstruct_single_ts)
            single_time_series = SingleTimeSeries("max_active_power", data)
            add_time_series!(sys_PRAS, component_PRAS, single_time_series)
        end
    end

    services = collect(
        PSY.get_components(
            x -> PSY.get_name(x) in ["Reg_Down", "Reg_Up"],
            Service,
            sys_PRAS,
        ),
    )
    for service_PRAS in services
        name = PSY.get_name(service_PRAS)
        reconstruct_single_ts = Float64[]
        for sys in sys_EDs
            service = first(PSY.get_components_by_name(PSY.Service, sys, name))
            forecast = PSY.get_time_series(Deterministic, service, "requirement")
            weather_year = Dates.year(first(keys(forecast.data)))
            for (key, value) in forecast.data
                if Dates.year(key) == weather_year
                    append!(reconstruct_single_ts, value[1:sys_PRAS_interval])
                end
            end
        end
        # append!(reconstruct_single_ts, reconstruct_single_ts[1:sys_PRAS_horizon - 1])
        # dates = range(DateTime("2018-01-01T00:00:00"), step = sys_PRAS_res, length = length(timestep) + sys_PRAS_horizon - 1)
        dates = range(SIM_START_DATE; step = sys_PRAS_res, length = length(timestep))
        # @info "Line 644 - Creating time series for service $(name) with length: $(length(reconstruct_single_ts))."
        data = TS.TimeArray(dates, reconstruct_single_ts)
        single_time_series = SingleTimeSeries("requirement", data)
        add_time_series!(sys_PRAS, service_PRAS, single_time_series)
    end
    remove_time_series!(sys_PRAS, Deterministic)
    remove_time_series!(sys_PRAS, DeterministicSingleTimeSeries)

    for (key, devices) in ts_objects
        for component in devices
            revisedts = DataStructures.SortedDict{DateTime, Vector{Float64}}()
            newtsdata = values(
                get_time_series(SingleTimeSeries, component, "max_active_power").data,
            )
            # @info "Line 668 - Processing component $(PSY.get_name(component)) with original time series length: $(length(newtsdata))."

            for t in 1:(length(timestep) - 1)
                rtseries=Vector{Float64}()
                datetimeindex = timestep[t]
                rtseries =
                    newtsdata[(sys_PRAS_interval * (t - 1) + 1):(sys_PRAS_interval * (t - 1) + sys_PRAS_horizon)]
                push!(revisedts, datetimeindex => rtseries)
            end

            # @info "Line 677 - Processing component $(PSY.get_name(component)) with revised time series length: $(length(revisedts))."
            # convert to deterministic time series
            revisedts_deterministic = PSY.Deterministic(;
                name = "max_active_power",
                data = revisedts,
                resolution = Dates.Hour(1),
                scaling_factor_multiplier = get_max_active_power,
            )

            add_time_series!(sys_PRAS, component, revisedts_deterministic)
        end
    end

    for service in services
        revisedts = DataStructures.SortedDict{DateTime, Vector{Float64}}()
        newtsdata = values(get_time_series(SingleTimeSeries, service, "requirement").data)

        for t in 1:(length(timestep) - 1)
            rtseries=Vector{Float64}()
            datetimeindex = timestep[t]
            rtseries =
                newtsdata[(sys_PRAS_interval * (t - 1) + 1):(sys_PRAS_interval * (t - 1) + sys_PRAS_horizon)]
            push!(revisedts, datetimeindex => rtseries)
        end

        # @info "Line 701 - Processing service $(PSY.get_name(service)) with revised time series length: $(length(revisedts))."
        # convert to deterministic time series
        revisedts_deterministic = PSY.Deterministic(;
            name = "requirement",
            data = revisedts,
            resolution = Dates.Hour(1),
            scaling_factor_multiplier = get_requirement,
        )

        add_time_series!(sys_PRAS, service, revisedts_deterministic)
    end

    remove_time_series!(sys_PRAS, Deterministic)

    # sys_PRAS was deepcopied from an ED system that already has GeometricDistributionForcedOutage
    # supplemental attributes (with single-year time series) embedded. Strip them here at the
    # call site before re-adding with the full multi-year PRAS date range.
    # Note: remove_supplemental_attribute! calls prepare_for_removal! which clears time series too.
    for comp in PSY.get_components(PSY.Generator, sys_PRAS)
        for attr in collect(
            PSY.get_supplemental_attributes(PSY.GeometricDistributionForcedOutage, comp),
        )
            PSY.remove_supplemental_attribute!(sys_PRAS, comp, attr)
        end
    end
    for comp in PSY.get_components(PSY.Storage, sys_PRAS)
        for attr in collect(
            PSY.get_supplemental_attributes(PSY.GeometricDistributionForcedOutage, comp),
        )
            PSY.remove_supplemental_attribute!(sys_PRAS, comp, attr)
        end
    end

    dates = range(SIM_START_DATE; step = sys_PRAS_res, length = length(timestep))
    add_outages_to_system!(sys_PRAS, outages_dir, dates)

    to_json(sys_PRAS, output_file; force = true)

    return
end

function create_PRAS_sys_NY_json(
    sys_EDs::Vector{PSY.System},
    output_file::String,
)
    sys_PRAS = deepcopy(first(sys_EDs))

    first_ts_temp_PRAS = first(PSY.get_time_series_multiple(sys_PRAS))
    start_datetime_PRAS = PSY.IS.get_initial_timestamp(first_ts_temp_PRAS)
    sys_PRAS_res = PSY.get_time_series_resolutions(sys_PRAS)[1]
    finish_datetime_PRAS =
        start_datetime_PRAS +
        Dates.Hour(DEFAULT_HOURS_PER_YEAR * sys_PRAS_res * length(sys_EDs) - sys_PRAS_res)
    timestep = StepRange(start_datetime_PRAS, sys_PRAS_res, finish_datetime_PRAS);

    ts_objects = Dict{String, Any}()

    # ts_objects["gens"] = collect(PSY.get_components(x -> PSY.get_prime_mover_type(x) in [PrimeMovers.PVe, PrimeMovers.WT, PrimeMovers.HY], Generator, sys_PRAS))
    ts_objects["gens"] = collect(
        PSY.get_components(x -> PSY.has_time_series(x) == true, Generator, sys_PRAS),
    )
    ts_objects["loads"] = collect(PSY.get_components(PSY_LOADS, sys_PRAS))

    for (key, devices) in ts_objects
        if key == "gens"
            ts_name = "fuel_cost"
        elseif key == "loads"
            ts_name = "max_active_power"
        end

        for component_PRAS in devices
            name = PSY.get_name(component_PRAS)
            reconstruct_single_ts = Float64[]
            for sys in sys_EDs
                component = first(PSY.get_components_by_name(PSY.Device, sys, name))
                ts_year = PSY.get_time_series(SingleTimeSeries, component, ts_name)
                append!(reconstruct_single_ts, values(ts_year.data))
            end
            dates = timestep
            data = TS.TimeArray(dates, reconstruct_single_ts)
            single_time_series = SingleTimeSeries(ts_name, data)
            remove_time_series!(sys_PRAS, SingleTimeSeries, component_PRAS, ts_name)
            add_time_series!(sys_PRAS, component_PRAS, single_time_series)
        end
    end

    to_json(sys_PRAS, output_file; force = true)

    return
end

function fix_multistart_cost_curves!(sys::PSY.System)
    base_power = PSY.get_base_power(sys)
    for unit in PSY.get_components(PSY.ThermalMultiStart, sys)
        # @info "Checking cost curve for multi-start unit $(PSY.get_name(unit))"
        cost = PSY.get_operation_cost(unit)
        vc = PSY.get_variable(cost)
        curve = PSY.get_value_curve(vc)
        pmax_mw = PSY.get_active_power_limits(unit).max * base_power

        if curve isa PSY.PiecewisePointCurve
            points = PSY.get_points(curve)
            if points[end].x < pmax_mw
                # @warn "$(PSY.get_name(unit)): PiecewisePointCurve ends at $(points[end].x) MW < Pmax $(pmax_mw) MW — extending"
                # @info "Original points: $(points)"
                slope = points[end].x > 0 ? points[end].y / points[end].x : 1.0
                new_points = vcat(points, [(x = pmax_mw, y = pmax_mw * slope)])
                new_curve = PSY.PiecewisePointCurve(new_points)
                new_vc =
                    PSY.CostCurve(new_curve, PSY.get_power_units(vc), PSY.get_vom_cost(vc))
                PSY.set_operation_cost!(
                    unit,
                    PSY.ThermalGenerationCost(
                        new_vc, PSY.get_fixed(cost), PSY.get_start_up(cost),
                        PSY.get_shut_down(cost),
                    ),
                )
            end

        elseif curve isa PSY.PiecewiseIncrementalCurve
            x_coords = PSY.get_x_coords(curve)
            if x_coords[end] < pmax_mw - 1e-3 # allow small numerical tolerance
                # @warn "$(PSY.get_name(unit)): PiecewiseIncrementalCurve ends at $(x_coords[end]) MW < Pmax $(pmax_mw) MW — extending"
                # @info "Original x coords: $(x_coords)"
                slopes = PSY.get_slopes(curve)
                # @info "Original slopes: $(slopes)"
                new_x = vcat(x_coords, pmax_mw)
                new_slopes = vcat(slopes, slopes[end])  # extend with last slope
                new_curve = PSY.PiecewiseIncrementalCurve(
                    PSY.get_initial_input(curve), new_x, new_slopes,
                )
                new_vc =
                    PSY.CostCurve(new_curve, PSY.get_power_units(vc), PSY.get_vom_cost(vc))
                PSY.set_operation_cost!(
                    unit,
                    PSY.ThermalGenerationCost(
                        new_vc, PSY.get_fixed(cost), PSY.get_start_up(cost),
                        PSY.get_shut_down(cost),
                    ),
                )
            end
        else
            # @debug "$(PSY.get_name(unit)): $(typeof(curve)) — no extension needed"
        end
    end
end
