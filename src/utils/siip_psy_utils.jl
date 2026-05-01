function get_name(device::D) where {D <: PSY.Device}
    return device.name
end

"""
This function returns all generation and storage technologies included in the PSY System.
"""
function get_all_techs(sys::PSY.System)
    # sys_gens = PSY.get_components(PSY.Generator, sys)
    sys_gens = PSY.get_components(x -> PSY.get_available(x) == true, Generator, sys)
    sys_storage = PSY.get_components(PSY.Storage, sys)
    sys_techs = union(sys_gens, sys_storage)
    return sys_techs
end

"""
This function does nothing if the PSY System is not defined.
"""
function get_system_services(sys::Nothing)
    services = PSY.Service[]
    return services
end

"""
This function returns all Services modeled in the PSY System.
"""
function get_system_services(sys::PSY.System)
    sys_techs = get_all_techs(sys)

    services = PSY.Service[]
    for tech in sys_techs
        for service in PSY.get_services(tech)
            if !in(service, services)
                push!(services, service)
            end
        end
    end

    return services
end

"""
This function does nothing if Service is not of ReserveUp type.
"""
function add_device_services!(sys::PSY.System,
    device::PSY.Device,
    product::P) where {P <: Product}
    return
end

"""
This function adds ReserveUp service for PSY Devices.
"""
function add_device_services!(sys::PSY.System,
    device::PSY.Device,
    product::Union{
        OperatingReserve{ReserveUpEMIS},
        OperatingReserve{ReserveDownEMIS},
        Inertia,
    })
    for service in get_system_services(sys)
        if PSY.get_name(service) == String(get_name(product))
            PSY.add_service!(device, service, sys)
        end
    end

    return
end

function add_clean_energy_contribution!(sys_UC::PSY.System,
    PSY_project_UC::PSY.Device)
    return
end

function add_clean_energy_contribution!(sys::PSY.System,
    device::T,
) where {T <: Union{
    PSY.ThermalStandard,
    PSY.RenewableDispatch,
    PSY.HydroTurbine,
    PSY.HydroDispatch,
}}
    services = get_system_services(sys)
    for service in services
        if PSY.get_name(service) == "Clean_Energy"
            PSY.add_service!(device, service, sys)
        end
    end
    return
end

"""
This function does nothing if Device is not of RenewableGen type.
"""
function add_device_forecast!(sys_MD::PSY.System,
    sys_UC::PSY.System,
    sys_ED::PSY.System,
    device_MD::D,
    device_UC::D,
    device_ED::D,
    availability_df::Vector{Float64},
    availability_raw_rt::Vector{Float64},
    da_resolution::Int64,
    rt_resolution::Int64,
) where {D <: Union{PSY.ThermalGen, PSY.HydroGen, PSY.Storage}}
    return
end

"""
This function adds forecast timeseries if Device is of RenewableGen type.
"""
function add_device_forecast!(sys_MD::PSY.System,
    sys_UC::PSY.System,
    sys_ED::PSY.System,
    device_MD::D,
    device_UC::D,
    device_ED::D,
    availability_raw::Vector{Float64},
    availability_raw_rt::Vector{Float64},
    da_resolution::Int64,
    rt_resolution::Int64) where {D <: PSY.RenewableGen}

    for (sys, device, avail_data, resolution_minutes) in [
        (sys_MD, device_MD, availability_raw,    da_resolution),
        (sys_UC, device_UC, availability_raw,    da_resolution),
        (sys_ED, device_ED, availability_raw_rt, rt_resolution),
    ]
        # Read forecast parameters from an existing component in this system
        any_load = first(PSY.get_components(
            x -> PSY.has_time_series(x, PSY.Deterministic, "max_active_power"),
            PSY.PowerLoad, sys))
        existing_ts = PSY.get_time_series(PSY.Deterministic, any_load, "max_active_power")
        forecast_count = PSY.get_count(existing_ts)
        sys_horizon = PSY.get_horizon(existing_ts)
        sys_horizon = Dates.Hour(sys_horizon).value
        sys_interval = PSY.get_forecast_interval(sys)
        sys_resolution = PSY.get_time_series_resolutions(sys)[1]
        start_datetime = SIM_START_DATE

        # Build timestamps to match existing system length exactly
        finish_datetime =
        start_datetime + Dates.Hour((
            forecast_count * sys_interval/sys_resolution +
            (sys_horizon - sys_interval/sys_resolution) - 1
        ))
        time_stamps = StepRange(start_datetime, Dates.Hour(1), finish_datetime)

        # Pad availability data by wrapping from the start of the year
        additional_timestep = length(time_stamps) - DEFAULT_HOURS_PER_YEAR
        avail_padded = additional_timestep > 0 ?
            [avail_data; avail_data[1:additional_timestep]] :
            avail_data[1:length(time_stamps)]

        sys_interval_steps = Int(sys_interval / sys_resolution)
        data = Dict(
            time_stamps[i] => avail_padded[i:(i + sys_horizon - 1)]
            for i in 1:sys_interval_steps:(length(time_stamps) - sys_horizon + 1)
        )
        forecast = PSY.Deterministic(
            "max_active_power", data, Dates.Minute(resolution_minutes))
        PSY.add_time_series!(sys, device, forecast)
    end

    #  ###### MD ######
    # start_datetime = Dates.DateTime("2018-01-01T00:00:00")                          
    # timestamp = range(start_datetime, step = Hour(1), length = 8760)
    # value_ts = availability_raw
    # data = TS.TimeArray(timestamp, value_ts)
    # single_time_series = PSY.SingleTimeSeries("max_active_power", data, scaling_factor_multiplier=get_max_active_power)
    # PSY.add_time_series!(sys_MD, device_MD, single_time_series)

    # ###### UC ######
    # value_ts = availability_raw
    # data = TS.TimeArray(timestamp, value_ts)
    # single_time_series = PSY.SingleTimeSeries("max_active_power", data, scaling_factor_multiplier=get_max_active_power)
    # PSY.add_time_series!(sys_UC, device_UC, single_time_series)

    # ###### ED ######
    # value_ts = availability_raw_rt
    # data = TS.TimeArray(timestamp, value_ts)
    # single_time_series = PSY.SingleTimeSeries("max_active_power", data, scaling_factor_multiplier=get_max_active_power)
    # PSY.add_time_series!(sys_ED, device_ED, single_time_series)


    return
end

"""
This function does nothing if Device is not of RenewableGen type.
"""
function add_device_forecast_PRAS!(sys::PSY.System,
    device::D,
    availability_raw_rt::Vector{Float64},
    rt_resolution::Int64,
    simulation_years::Int64) where {D <: Union{PSY.ThermalGen, PSY.HydroGen, PSY.Storage}}
    return
end

"""
This function adds forecast timeseries to PRAS System if Device is of RenewableGen type.
"""
function add_device_forecast_PRAS!(sys::PSY.System,
    device::D,
    availability_raw_rt::Vector{Float64},
    rt_resolution::Int64,
    simulation_years::Int64
) where {D <: Union{PSY.RenewableGen, PSY.RenewableDispatch}}
    start_datetime = SIM_START_DATE
    timestamp = range(start_datetime; step = Hour(1), length = DEFAULT_HOURS_PER_YEAR * simulation_years)
    value_ts = availability_raw_rt
    data = TS.TimeArray(timestamp, value_ts)
    single_time_series = PSY.SingleTimeSeries(
        "max_active_power",
        data;
        scaling_factor_multiplier = get_max_active_power,
    )
    PSY.add_time_series!(sys, device, single_time_series)

    ######### Adding to PRAS##########
    # sys_interval = sys.data.time_series_params.forecast_params.interval
    # sys_horizon = sys.data.time_series_params.forecast_params.horizon
    # forecast_count = sys.data.time_series_params.forecast_params.count
    # sys_resolution = sys.data.time_series_params.resolution
    # start_datetime = sys.data.time_series_params.forecast_params.initial_timestamp
    # finish_datetime = start_datetime + Dates.Hour((forecast_count * sys_interval/sys_resolution + (sys_horizon - sys_interval/sys_resolution) - 1))
    # time_stamps = StepRange(start_datetime, Dates.Hour(1), finish_datetime);

    # append!(availability_raw_rt, availability_raw_rt[1:sys_horizon - 1])

    # data = Dict(time_stamps[i] => availability_raw_rt[i:(i + sys_horizon - 1)] for i in 1:Int(sys_interval/sys_resolution):(length(time_stamps)-sys_horizon + 1))
    # forecast = PSY.Deterministic("max_active_power", data, Dates.Minute(rt_resolution))
    # PSY.add_time_series!(sys, device, forecast)
    return
end

"""
This function does nothing if Device is not of RenewableGen type.
"""
function write_vg_data(sys::PSY.System,
    sys_ED::PSY.System,
    device_UC::D,
    device_ED::D,
    simulation_dir::String,
    availability_raw::Vector{Float64},
    availability_raw_rt::Vector{Float64},
    scenario_names::Vector{String},
    year::Int64,
) where {D <: Union{PSY.ThermalGen, PSY.HydroGen, PSY.Storage}}
end

"""
This function adds forecast timeseries to load_n_vg data files if Device is of RenewableGen type.
"""
function write_vg_data(sys_UC::PSY.System,
    sys_ED::PSY.System,
    device_UC::D,
    device_ED::D,
    simulation_dir::String,
    availability_raw::Vector{Float64},
    availability_raw_rt::Vector{Float64},
    scenario_names::Vector{String},
    year::Int64,
) where {D <: PSY.RenewableGen}

    ########## Adding to Net Load Data File ##############
    for scenario in scenario_names
        load_n_vg_df = read_data(
            joinpath(
                simulation_dir,
                "timeseries_data_files",
                scenario,
                "sim_year_$(year)",
                "Net Load Data",
                "load_n_vg_data.csv",
            ),
        )
        load_n_vg_df[:, get_name(device_UC)] =
            availability_raw[1:DataFrames.nrow(load_n_vg_df)] * get_device_size(device_UC) *
            PSY.get_base_power(sys_UC)
        load_n_vg_df_rt = read_data(
            joinpath(
                simulation_dir,
                "timeseries_data_files",
                scenario,
                "sim_year_$(year)",
                "Net Load Data",
                "load_n_vg_data_rt.csv",
            ),
        )
        load_n_vg_df_rt[:, get_name(device_ED)] =
            availability_raw_rt[1:DataFrames.nrow(load_n_vg_df_rt)] *
            get_device_size(device_ED) * PSY.get_base_power(sys_ED)

        write_data(
            joinpath(
                simulation_dir,
                "timeseries_data_files",
                scenario,
                "sim_year_$(year)",
                "Net Load Data",
            ),
            "load_n_vg_data.csv",
            load_n_vg_df,
        )
        write_data(
            joinpath(
                simulation_dir,
                "timeseries_data_files",
                scenario,
                "sim_year_$(year)",
                "Net Load Data",
            ),
            "load_n_vg_data_rt.csv",
            load_n_vg_df_rt,
        )
    end

    return
end

"""
This function does nothing if PSY System is not defined.
"""
function update_PSY_timeseries!(sys_UC::Nothing,
    rec_requirement::Float64,
    simulation_dir::String,
    type::String,
    pcm_scenario::String,
    iteration_year::Int64,
    da_resolution::Int64,
    rt_resolution::Int64)
    return
end

"""
This function updates the PSY load and reserve requirment timeseries each year.
"""
function update_PSY_timeseries!(
    simulation::AgentSimulation,
    sys::PSY.System,
    rec_requirement::Float64,
    simulation_dir::String,
    type::String,
    pcm_scenario::String,
    iteration_year::Int64,
    da_resolution::Int64,
    rt_resolution::Int64)
    total_active_power = 0.0

    first_ts_temp = first(PSY.get_time_series_multiple(sys))
    start_datetime = PSY.IS.get_initial_timestamp(first_ts_temp)
    sys_res = PSY.get_time_series_resolutions(sys)[1]
    finish_datetime = start_datetime + Dates.Hour(8760 * sys_res - sys_res)
    time_stamps = StepRange(start_datetime, sys_res, finish_datetime)

    # Only update ORDC product time series.
    services = get_system_services(sys)
    ordc_products = split(
        read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[
            1,
            "ordc_products",
        ],
        "; ",
    )
    for service in services
        service_name = PSY.get_name(service)
        println("Current service is $(service_name)")
        if service_name in ordc_products
            ### TODO: CHange when ORDC is available
            ### NY_change: ORDC timeseries is not available now
            # @warn "UPDATE ORDC TIMESERIES REMOVAL BASED ON NEW SIENNA SYSTEM"

            # PSY.remove_time_series!(sys, PSY.Deterministic, service, "variable_cost")
            # PSY.remove_time_series!(sys, PSY.SingleTimeSeries, service, "variable_cost")

            # time_stamps = StepRange(start_datetime, Dates.Hour(1), finish_datetime);
            # if type == "ED"
            #     product_ts_raw = read_data(joinpath(simulation_dir, "timeseries_data_files", pcm_scenario, "sim_year_$(iteration_year)", "Reserves", "$(service_name)_REAL_TIME.csv"))[:, service_name]
            # else
            #     product_ts_raw = read_data(joinpath(simulation_dir, "timeseries_data_files", pcm_scenario, "sim_year_$(iteration_year)", "Reserves", "$(service_name).csv"))[:, service_name]
            # end
            # product_data_ts = process_ordc_data_for_siip(product_ts_raw)
            # dates = time_stamps
            # data = TS.TimeArray(dates, product_data_ts)
            # single_time_series = SingleTimeSeries("variable_cost", data)
            # add_time_series!(sys, service, single_time_series)
        elseif service_name == "Clean_Energy"
            @warn "Service is clean energy"
        else
            @warn "Update service timeseries (other than ORDC)"
            # reserve_scaling_factor = calculate_reserve_scaling_factor(simulation)
            load_initial = read_data(
                joinpath(
                    simulation_dir,
                    "timeseries_data_files",
                    pcm_scenario,
                    "sim_year_1",
                    "Load",
                    "load.csv",
                ),
            )
            load_current = read_data(
                joinpath(
                    simulation_dir,
                    "timeseries_data_files",
                    pcm_scenario,
                    "sim_year_$(iteration_year)",
                    "Load",
                    "load.csv",
                ),
            )
            load_initial_total =
                sum(sum(eachcol(load_initial[:, Not(:Year, :Month, :Day, :Period)])))
            load_current_total =
                sum(sum(eachcol(load_current[:, Not(:Year, :Month, :Day, :Period)])))

            reserve_scaling_factor =
                (load_current_total - load_initial_total) / load_initial_total

            @info "Reserve scaling factor is $(reserve_scaling_factor)."
            scaled_requirement =
                deepcopy(PSY.get_requirement(service)) * (1 + reserve_scaling_factor)
            PSY.set_requirement!(service, scaled_requirement)
        end
    end

    return
end

"""
This function updates thermal generator outagestimeseries each year.
"""
function update_PSY_outage_timeseries!(sys_UC::PSY.System,
    sys_ED::PSY.System,
    resultsfolder::String,
    base_dir::String,
    iteration_year::Int64)
    if iteration_year == 1
        outagefolder = joinpath(resultsfolder, "GeneratorOutage")
        dir_exists(outagefolder)
        cp(joinpath(base_dir, "GeneratorOutage"), outagefolder; force = true)
    end
    # to update outage profile at the beginning of each iteration year
    outage_csv_location=joinpath(resultsfolder, "GeneratorOutage") #get_base_dir(case)

    # Remove existing outage time series for thermal generators
    for sys in [sys_UC, sys_ED]
        for d in PSY.get_components(ThermalGen, sys)  #including ThermalStandard & ThermalFastStartSIIP
            if ("outage" in get_time_series_names(Deterministic, d))
                PSY.remove_time_series!(sys, Deterministic, d, "outage")
            end
        end
    end

    # Updating outage time series for thermal generators from last year's PRAS samples
    outagescenario = 1
    sys_UC, sys_ED = add_csv_time_series!(
        sys_UC,
        sys_ED,
        outage_csv_location;
        add_scenario = outagescenario,
        iteration_year = iteration_year,
    ); # align outage structure to be the same as PRAS-ED (same value persist 36 times)

    # Check to see if all ThermalStandard has time series data
    # For generators that don't have time series, add artifical outage time series
    file_location = joinpath(outage_csv_location, "1/Generator_year$(iteration_year).csv")
    df_outage_profile = DataFrames.DataFrame(CSV.File(file_location));
    for d in PSY.get_components(ThermalGen, sys_UC)
        if ~("outage" in get_time_series_names(Deterministic, d))
            @info "Adding outage time series to new generators..."
            first_ts_temp_DA = first(PSY.get_time_series_multiple(sys_UC));
            len_ts_data_DA = length(first_ts_temp_DA.data);
            days_of_interest = range(1; length = len_ts_data_DA)
            num_days = length(days_of_interest);
            period_of_interest = range(
                days_of_interest.start==1 ? 1 : (24*days_of_interest.start)+1;
                length = num_days*24,
            )
            N = length(period_of_interest);
            start_datetime_DA = PSY.IS.get_initial_timestamp(first_ts_temp_DA);
            sys_DA_res_in_hour = PSY.get_time_series_resolution(sys_UC)
            start_datetime_DA =
                start_datetime_DA +
                Dates.Hour((period_of_interest.start-1)*sys_DA_res_in_hour);
            finish_datetime_DA = start_datetime_DA + Dates.Hour((N-1)*sys_DA_res_in_hour);
            all_timestamps =
                StepRange(start_datetime_DA, sys_DA_res_in_hour, finish_datetime_DA);
            allnames = PSY.get_name.(
                PSY.get_components(
                    ThermalFastStartSIIP,
                    sys_UC,
                    x ->
                        get_prime_mover_type(x) in [PrimeMovers.CT, PrimeMovers.GT] &&
                        PSY.get_fuel(x) in [ThermalFuels.NATURAL_GAS],
                ),
            )
            allnames = allnames[findall(x->x in names(df_outage_profile), allnames)]
            asset_name = allnames[rand(1:length(allnames))]
            @info "Use $(asset_name) time series for $(PSY.get_name(d))..."

            DA_data = Dict()
            RT_data = Dict()

            for (idx, timestamp) in enumerate(all_timestamps) #
                if (rem(idx, 24) == 1)
                    if (df_outage_profile[idx, asset_name] == 0)
                        push!(DA_data, timestamp => zeros(36))
                    else
                        push!(DA_data, timestamp => ones(36))
                    end
                end

                if (df_outage_profile[idx, asset_name] == 0)
                    # If PRAS actual status is actually 0, just use 0 for RT, ReliablePSY can handle this transition
                    push!(RT_data, timestamp => zeros(2))
                else
                    # If PRAS actual status is actually 1 and DA says 1, use 1 for RT, if DA says 0, use O because ReliablePSY currently can't handle this
                    if (rem(idx, 24) == 0)
                        offset = 24
                    else
                        offset = rem(idx, 24)
                    end

                    if (df_outage_profile[(idx - offset + 1), asset_name] == 1) # If PRAS Actual == 1 and DA ==1
                        push!(RT_data, timestamp => ones(2))
                    else
                        push!(RT_data, timestamp => ones(2)) # # If PRAS Actual == 1 and DA ==0 - Special case
                    end
                end
            end

            DA_availability_forecast =
                PSY.Deterministic("outage", DA_data, sys_DA_res_in_hour)
            RT_availability_forecast =
                PSY.Deterministic("outage", RT_data, sys_DA_res_in_hour)

            # Adding TimeSeries Data to PSY System
            PSY.add_time_series!(
                sys_UC,
                PSY.get_component(PSY.Generator, sys_UC, PSY.get_name(d)),
                DA_availability_forecast,
            )
            PSY.add_time_series!(
                sys_ED,
                PSY.get_component(PSY.Generator, sys_ED, PSY.get_name(d)),
                RT_availability_forecast,
            )
        end
    end
end

"""
This function finds the buses located in each zone.
"""
function find_zonal_bus(zone::String, sys::PSY.System)
    # buses = PSY.get_components(PSY.Bus, sys)

    # zonal_bus = nothing

    # for b in buses
    #     if "zone_$(PSY.get_name(PSY.get_area(b)))" == zone
    #         zonal_bus = b
    #     end
    # end

    zonal_bus = PSY.get_component(PSY.ACBus, sys, zone)

    return zonal_bus
end

function find_zonal_area(zone::String, sys::PSY.System)
    zone_number_map = Dict(
        "FarWest" => 1,
        "North" => 2,
        "West" => 3,
        "Southern" => 4,
        "NorthCentral" => 5,
        "SouthCentral" => 6,
        "Coast" => 7,
        "East" => 8,
    )
    for area in PSY.get_components(PSY.Area, sys)
        name = PSY.get_name(area)
        if zone == "zone_$(zone_number_map[name])"
            return area
        end
    end
    return nothing
end

"""
This function transforms the timeseries of PSY Systems.
"""

function transform_psy_timeseries!(sys_MD::Nothing,
    sys_UC::Nothing,
    sys_ED::Nothing,
    da_resolution::Int64,
    rt_resolution::Int64,
    da_horizon::Int64,
    rt_horizon::Int64,
    md_interval::Int64,
    da_interval::Int64,
    rt_interval::Int64)
    return
end

function transform_psy_timeseries!(sys_MD::PSY.System,
    sys_UC::PSY.System,
    sys_ED::PSY.System,
    da_resolution::Int64,
    rt_resolution::Int64,
    md_horizon::Int64,
    da_horizon::Int64,
    rt_horizon::Int64,
    md_interval::Int64,
    da_interval::Int64,
    rt_interval::Int64)

    @info "Updating timeseries in transform_psy_timeseries!"
    
    # TODO: may want to add md_resolution
    @info "Updating MD timeseries with horizon $(md_horizon) and interval $(md_interval) and da_resolution $(da_resolution)..."
    PSY.transform_single_time_series!(
        sys_MD,
        Dates.Hour(Int(md_horizon * 60 / da_resolution)),
        Dates.Hour(md_interval),
    )
    @info "Updating UC timeseries with horizon $(da_horizon) and interval $(da_interval) and da_resolution $(da_resolution)..."
    PSY.transform_single_time_series!(
        sys_UC,
        Dates.Hour(Int(da_horizon * 60 / da_resolution)),
        Dates.Hour(da_interval),
    )
    @info "Updating ED timeseries with horizon $(rt_horizon) and interval $(rt_interval) and rt_resolution $(rt_resolution)..."
    PSY.transform_single_time_series!(
        sys_ED,
        Dates.Hour(Int(rt_horizon * 60 / rt_resolution)),
        Dates.Hour(rt_interval),
    )
    return
end

function add_psy_inertia!(simulation_dir::String,
    sys::Nothing,
    reserve_penalty::String,
    system_peak_load::Float64)
    return
end

function add_psy_inertia!(simulation_dir::String,
    sys::PSY.System,
    type::String,
    reserve_penalty::String,
    system_peak_load::Float64)
    inertia_data = read_data(
        joinpath(
            simulation_dir,
            "markets_data",
            "$(reserve_penalty)_reserve_penalty",
            "Inertia.csv",
        ),
    )

    inertia_requirement =
        system_peak_load * inertia_data[1, "requirement_multiplier"] /
        PSY.get_base_power(sys)

    ####### Adding Inertia reserve

    inertia_reserve = PSY.VariableReserve{PSY.ReserveUp}(
        "Inertia",
        true,
        3600,
        inertia_requirement,
    )
    contri_devices =
        vcat(collect(PSY.get_components(PSY.ThermalStandard, sys)),
            collect(PSY.get_components(ThermalFastStartSIIP, sys)),
            collect(PSY.get_components(PSY.RenewableDispatch, sys)),
            collect(PSY.get_components(PSY.HydroDispatch, sys)),
            collect(PSY.get_components(PSY.HydroTurbine, sys)),
            collect(PSY.get_components(PSY.EnergyReservoirStorage, sys)),
        );

    PSY.add_service!(sys, inertia_reserve, contri_devices)

    # time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
    #                 PSY.SingleTimeSeries,
    #                 first(PSY.get_components(PSY.ElectricLoad, sys)),
    #                 "max_active_power"
    #                 )))

    if type in ["UC", "ED"]
        if type == "UC"
            time_stamps = StepRange(
                SIM_START_DATE,
                Dates.Hour(1),
                SIM_END_DATE,
            );
        elseif type == "ED"
            time_stamps = StepRange(
                SIM_START_DATE,
                Dates.Hour(1),
                SIM_END_DATE,
            );
        else
            error("Type should be UC or ED")
        end

        ts_data = ones(length(time_stamps))
        ts = TimeSeries.TimeArray(time_stamps, ts_data);
        forecast = PSY.SingleTimeSeries("requirement", ts)
        PSY.add_time_series!(sys, inertia_reserve, forecast)
    end
end

function add_psy_clean_energy_constraint!(simulation_dir::String,
    sys::Nothing,
    requirement::Float64)
    return
end

function add_psy_clean_energy_constraint!(sys::PSY.System,
    requirement::Float64)

    ####### Adding Clean Energy Constraint
    total_active_power = 0.0
    # time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
    #                 PSY.SingleTimeSeries,
    #                 first(PSY.get_components(PSY.ElectricLoad, sys)),
    #                 "max_active_power"
    #                 )))
    sys_interval = PSY.get_forecast_interval(sys)
    sys_horizon = PSY.get_forecast_horizon(sys)
    sys_horizon = Dates.Hour(sys_horizon).value
    forecast_count = PSY.get_forecast_window_count(sys)
    sys_resolution = PSY.get_time_series_resolutions(sys)[1]
    start_datetime = PSY.get_forecast_initial_timestamp(sys)
    finish_datetime =
        start_datetime + Dates.Hour((
            forecast_count * sys_interval/sys_resolution +
            (sys_horizon - sys_interval/sys_resolution) - 1
        ))
    time_stamps = StepRange(start_datetime, Dates.Hour(1), finish_datetime);
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    clean_energy_ts_data = zeros(length(time_stamps))

    nodal_loads = PSY.get_components(PSY_LOADS, sys)
    for load in nodal_loads
        load_active_power = PSY.get_max_active_power(load)
        total_active_power += load_active_power
        # ts_data = deepcopy(PSY.get_data(PSY.get_time_series(PSY.SingleTimeSeries,
        #                                        load,
        #                                        "max_active_power"
        #                                         )))
        loadts_raw=PSY.get_data(
            PSY.get_time_series(PSY.Deterministic,
                load,
                "max_active_power",
            ),
        )
        ts_data = []
        for timestep in keys(loadts_raw)
            ts_data=[ts_data; loadts_raw[timestep][1:Int(sys_interval / sys_resolution)]]
        end
        ts_data=[ts_data; loadts_raw[end][(Int(sys_interval / sys_resolution) + 1):end]]
        clean_energy_ts_data .+= (TS.values(ts_data) .* load_active_power)
    end

    clean_energy_ts_data = clean_energy_ts_data / total_active_power

    clean_energy_reserve = PSY.VariableReserve{PSY.ReserveUp}(
        "Clean_Energy",
        true,
        3600,
        total_active_power * requirement,
    )
    contri_devices =
        vcat(
            collect(
                PSY.get_components(
                    x -> (
                        occursin("NUC", PSY.get_name(x)) ||
                        occursin("RECT", PSY.get_name(x))
                    ),
                    PSY.ThermalStandard,
                    sys,
                ),
            ),
            collect(PSY.get_components(PSY.RenewableDispatch, sys)),
            collect(PSY.get_components(PSY.HydroDispatch, sys)),
            collect(PSY.get_components(PSY.HydroTurbine, sys)),
        );

    PSY.add_service!(sys, clean_energy_reserve, contri_devices)

    ts_data = TimeSeries.TimeArray(time_stamps, clean_energy_ts_data)
    forecast = PSY.SingleTimeSeries("requirement", ts_data)
    PSY.add_time_series!(sys, clean_energy_reserve, forecast)
end

function calculate_total_load(
    sys::PSY.System,
    time_resolution::Int64,
    period_of_interest::UnitRange{Int64},
)
    @info "Calculating total load for period of interest: $(period_of_interest)"
    total_load = 0.0

    nodal_loads = PSY.get_components(PSY_LOADS, sys)
    for load in nodal_loads
        ts_data = PSY.get_data(
            PSY.get_time_series(
                PSY.SingleTimeSeries,
                load,
                "max_active_power"),
        )

        # @info "Calculating load for $(PSY.get_name(load)), load length is $(length(ts_data))"

        total_load +=
            sum(TS.values(ts_data[period_of_interest]) * PSY.get_max_active_power(load)) *
            time_resolution / DEFAULT_TIME_RESOLUTION
    end
    return total_load
end

function convert_to_thermal_fast_start!(d::PSY.ThermalStandard, system::PSY.System)
    new = ThermalFastStartSIIP(;
        name = d.name,
        available = d.available,
        status = d.status,
        bus = d.bus,
        active_power = d.active_power,
        reactive_power = d.reactive_power,
        rating = d.rating,
        active_power_limits = d.active_power_limits,
        reactive_power_limits = d.reactive_power_limits,
        ramp_limits = d.ramp_limits,
        operation_cost = d.operation_cost,
        base_power = d.base_power,
        time_limits = d.time_limits,
        prime_mover_type = d.prime_mover_type,
        fuel = d.fuel,
        ext = d.ext,
    )

    PSY.add_component!(system, new)
    for service in PSY.get_services(d)
        PSY.add_service!(new, service, system)
    end

    ts_names = []
    all_ts = PSY.get_time_series_multiple(d)
    for i in all_ts
        push!(ts_names, i.name)
    end

    if ("outage" in ts_names)
        newts = PSY.get_time_series(PSY.Deterministic, d, "outage")
        PSY.add_time_series!(system, new, newts)
    end

    PSY.remove_component!(system, d)
    return
end

function convert_thermal_fast_start!(system::PSY.System)
    for gen in PSY.get_components(PSY.ThermalStandard, system)
        prime_mover_type = PSY.get_prime_mover_type(gen)
        fuel = PSY.get_fuel(gen)

        target_prime_mover = [PSY.PrimeMovers.CT, PSY.PrimeMovers.GT]
        target_fuel = PSY.ThermalFuels.NATURAL_GAS

        if prime_mover_type in target_prime_mover && fuel == target_fuel
            convert_to_thermal_fast_start!(gen, system)
        end
    end
end

function apply_PSY_past_load_growth!(sys::PSY.System,
    load_growth::AxisArrays.AxisArray{Float64, 1},
    simulation_dir::String)

    # update load timeseries.
    nodal_loads = PSY.get_components(PSY_LOADS, sys)

    for load in nodal_loads
        zone = "load_zone_$(PSY.get_name(PSY.get_area(PSY.get_bus(load))))"
        scaled_active_power =
            deepcopy(PSY.get_max_active_power(load)) * (1 + load_growth[zone])
        PSY.set_max_active_power!(load, scaled_active_power)
    end

    average_load_growth = Statistics.mean(load_growth)

    # update service requirement timeseries.
    services = get_system_services(sys)
    ordc_products = split(
        read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[
            1,
            "ordc_products",
        ],
        "; ",
    )
    for service in services
        service_name = PSY.get_name(service)
        if !(service_name in ordc_products || service_name == "Clean_Energy")
            scaled_requirement =
                deepcopy(PSY.get_requirement(service)) * (1 + average_load_growth)
            PSY.set_requirement!(service, scaled_requirement)
        end
    end

    return
end
