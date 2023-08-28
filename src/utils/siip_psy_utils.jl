function get_name(device::D) where D <: PSY.Device
    return device.name
end

"""
This function returns all generation and storage technologies included in the PSY System.
"""
function get_all_techs(sys::PSY.System)
    # sys_gens = PSY.get_components(PSY.Generator, sys)
    sys_gens = PSY.get_components(PSY.Generator, sys, x -> PSY.get_available(x) == true)
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
                              product::P) where P <: Product
    return
end

"""
This function adds ReserveUp service for PSY Devices.
"""
function add_device_services!(sys::PSY.System,
                              device::PSY.Device,
                              product::Union{OperatingReserve{ReserveUpEMIS}, OperatingReserve{ReserveDownEMIS}, Inertia})


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
                                        device::T) where T <:Union{PSY.ThermalStandard, PSY.RenewableDispatch, PSY.HydroEnergyReservoir, PSY.HydroDispatch}

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
function add_device_forecast!(simulation_dir::String,
                              sys_UC::PSY.System,
                              sys_ED::PSY.System,
                              device_UC::D,
                              device_ED::D,
                              availability_df::Vector{Float64},
                              availability_raw_rt::Vector{Float64},
                              da_resolution::Int64,
                              rt_resolution::Int64
                              ) where D <: Union{PSY.ThermalGen, PSY.HydroGen, PSY.Storage}
    return
end

"""
This function adds forecast timeseries if Device is of RenewableGen type.
"""
function add_device_forecast!(simulation_dir::String,
                              sys_UC::PSY.System,
                              sys_ED::PSY.System,
                              device_UC::D,
                              device_ED::D,
                              availability_raw::Vector{Float64},
                              availability_raw_rt::Vector{Float64},
                              da_resolution::Int64,
                              rt_resolution::Int64) where D <: PSY.RenewableGen



    ######### Adding to UC##########
    # time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
    #                                                 PSY.SingleTimeSeries,
    #                                                 first(PSY.get_components(PSY.ElectricLoad, sys_UC)),
    #                                                 "max_active_power"
    #                                                 )))
    # time_stamps = StepRange(Dates.DateTime("2018-01-01T00:00:00"), Dates.Hour(1), Dates.DateTime("2018-12-31T23:00:00"));
    # intervals = Int(24 * 60 / da_resolution)
    # append!(availability_raw, availability_raw[(length(availability_raw) - intervals + 1):end])
    # data = Dict(time_stamps[i] => availability_raw[i:(i + intervals - 1)] for i in 1:intervals:length(time_stamps))
    # forecast = PSY.Deterministic("max_active_power", data, Dates.Minute(da_resolution))
    # PSY.add_time_series!(sys_UC, device_UC, forecast)

    time_stamps = StepRange(Dates.DateTime("2018-01-01T00:00:00"), Dates.Hour(1), Dates.DateTime("2019-01-01T11:00:00"));
    intervals = Int(36 * 60 / da_resolution)
    append!(availability_raw, availability_raw[(length(availability_raw) - intervals + 25):end])
    data = Dict(time_stamps[i] => availability_raw[i:(i + intervals - 1)] for i in 1:Int(24 * 60 / da_resolution):8760)
    forecast = PSY.Deterministic("max_active_power", data, Dates.Minute(da_resolution))
    PSY.add_time_series!(sys_UC, device_UC, forecast)

    ######### Adding to ED##########
    # time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
    #                                                 PSY.SingleTimeSeries,
    #                                                 first(PSY.get_components(PSY.ElectricLoad, sys_ED)),
    #                                                 "max_active_power"
    #                                                 )))
    time_stamps = StepRange(Dates.DateTime("2018-01-01T00:00:00"), Dates.Hour(1), Dates.DateTime("2019-01-01T01:00:00"));
    intervals =  Int(2*60 / rt_resolution)
    append!(availability_raw_rt, availability_raw_rt[(length(availability_raw_rt) - intervals + 1):end])
    data = Dict(time_stamps[i] => availability_raw_rt[i:(i + intervals - 1)] for i in 1:Int(60 / rt_resolution):8760)
    forecast = PSY.Deterministic("max_active_power", data, Dates.Minute(rt_resolution))
    PSY.add_time_series!(sys_ED, device_ED, forecast)

    ########## Adding to Net Load Data ##############
    load_n_vg_df =  read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))
    load_n_vg_df[:, get_name(device_UC)] = availability_raw[1:DataFrames.nrow(load_n_vg_df)] * get_device_size(device_UC) * PSY.get_base_power(sys_UC)

    load_n_vg_df_rt =  read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data_rt.csv"))
    load_n_vg_df_rt[:, get_name(device_UC)] = availability_raw_rt[1:DataFrames.nrow(load_n_vg_df_rt)] * get_device_size(device_UC) * PSY.get_base_power(sys_UC)

    write_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data"), "load_n_vg_data.csv", load_n_vg_df)
    write_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data"), "load_n_vg_data_rt.csv", load_n_vg_df_rt)

    return
end


"""
This function does nothing if PSY System is not defined.
"""
function update_PSY_timeseries!(sys_UC::Nothing,
                               load_growth::AxisArrays.AxisArray{Float64, 1},
                               rec_requirement::Float64,
                               simulation_dir::String,
                               type::String,
                               iteration_year::Int64,
                               da_resolution::Int64,
                               rt_resolution::Int64)
    return
end

"""
This function updates the PSY load and reserve requirment timeseries each year.
"""
function update_PSY_timeseries!(sys::PSY.System,
                               load_growth::AxisArrays.AxisArray{Float64, 1},
                               rec_requirement::Float64,
                               simulation_dir::String,
                               type::String,
                               iteration_year::Int64,
                               da_resolution::Int64,
                               rt_resolution::Int64)

    total_active_power = 0.0

    # update load timeseries.
    nodal_loads = PSY.get_components(PSY.PowerLoad, sys)

    for load in nodal_loads
        zone = "zone_$(PSY.get_name(PSY.get_area(PSY.get_bus(load))))"
        scaled_active_power = deepcopy(PSY.get_max_active_power(load)) * (1 + load_growth[zone])
        println(PSY.get_name(load))
        println(scaled_active_power)

        PSY.set_max_active_power!(load, scaled_active_power)

        total_active_power += scaled_active_power
    end

    println(total_active_power)

    average_load_growth = Statistics.mean(load_growth)

    # update service requirement timeseries.
    services = get_system_services(sys)
    ordc_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"ordc_products"], "; ")
    for service in services
        service_name = PSY.get_name(service)
        println("Current service is $(service_name)")
        if service_name in ordc_products
            # time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
            #                                         PSY.SingleTimeSeries,
            #                                         first(PSY.get_components(PSY.ElectricLoad, sys)),
            #                                         "max_active_power"
            #                                         )))
            # time_stamps = StepRange(Dates.DateTime("2018-01-01T00:00:00"), Dates.Hour(1), Dates.DateTime("2018-12-31T23:00:00"));

            PSY.remove_time_series!(sys, PSY.Deterministic, service, "variable_cost")
            if type == "UC"
                # product_ts_raw = read_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "$(service_name)_$(iteration_year - 1).csv"))[:, service_name]
                # product_data_ts = process_ordc_data_for_siip(product_ts_raw)
                # intervals = Int(24 * 60 / da_resolution)
                # append!(product_data_ts, product_data_ts[(length(product_data_ts) - intervals + 1):end])
                # data = Dict(time_stamps[i] => product_data_ts[i:(i + intervals - 1)] for i in 1:intervals:length(time_stamps))
                # forecast = PSY.Deterministic("variable_cost", data, Dates.Minute(da_resolution))

                time_stamps = StepRange(Dates.DateTime("2018-01-01T00:00:00"), Dates.Hour(1), Dates.DateTime("2019-01-01T11:00:00"));
                product_ts_raw = read_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "$(service_name)_$(iteration_year - 1).csv"))[:, service_name]
                product_data_ts = process_ordc_data_for_siip(product_ts_raw)
                intervals = Int(36 * 60 / da_resolution)
                append!(product_data_ts, product_data_ts[(length(product_data_ts) - intervals + 25):end])
                data = Dict(time_stamps[i] => product_data_ts[i:(i + intervals - 1)] for i in 1:Int(24 * 60 / da_resolution):8760)
                forecast = PSY.Deterministic("variable_cost", data, Dates.Minute(da_resolution))

            elseif type == "ED"
                # product_ts_raw = read_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "$(service_name)_REAL_TIME_$(iteration_year - 1).csv"))[:, service_name]
                # product_data_ts = process_ordc_data_for_siip(product_ts_raw)
                # intervals =  Int(60 / rt_resolution)
                # append!(product_data_ts, product_data_ts[(length(product_data_ts) - intervals + 1):end])
                # data = Dict(time_stamps[i] => product_data_ts[i:(i + intervals  - 1)] for i in 1:intervals:length(time_stamps))
                # forecast = PSY.Deterministic("variable_cost", data, Dates.Minute(rt_resolution))

                time_stamps = StepRange(Dates.DateTime("2018-01-01T00:00:00"), Dates.Hour(1), Dates.DateTime("2019-01-01T01:00:00"));
                product_ts_raw = read_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "$(service_name)_REAL_TIME_$(iteration_year - 1).csv"))[:, service_name]
                product_data_ts = process_ordc_data_for_siip(product_ts_raw)
                intervals =  Int(2*60 / rt_resolution)
                append!(product_data_ts, product_data_ts[(length(product_data_ts) - intervals + 1):end])
                data = Dict(time_stamps[i] => product_data_ts[i:(i + intervals  - 1)] for i in 1:Int(60 / rt_resolution):8760)
                forecast = PSY.Deterministic("variable_cost", data, Dates.Minute(rt_resolution))
            else
                error("Type should be UC or ED")
            end

            PSY.add_time_series!(sys, service, forecast)
        elseif service_name == "Clean_Energy"
            clean_energy_requirement = total_active_power * rec_requirement * 1.0
            println(rec_requirement)
            println(clean_energy_requirement)
            PSY.set_requirement!(service, clean_energy_requirement)
        else
            scaled_requirement = deepcopy(PSY.get_requirement(service)) * (1 + average_load_growth)
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

    if iteration_year ==1
        outagefolder = joinpath(resultsfolder,"GeneratorOutage")
        dir_exists(outagefolder)
        cp(joinpath(base_dir,"GeneratorOutage"), outagefolder, force = true)
    end
    # to update outage profile at the beginning of each iteration year
    outage_csv_location=joinpath(resultsfolder,"GeneratorOutage") #get_base_dir(case)

    # Remove existing outage time series for thermal generators
    for sys in [sys_UC,sys_ED]
        for d in PSY.get_components(ThermalGen,sys)  #including ThermalStandard & ThermalFastStartSIIP
            if ("outage" in get_time_series_names(Deterministic, d))
                PSY.remove_time_series!(sys, Deterministic, d, "outage")
            end
        end
    end

    # Updating outage time series for thermal generators from last year's PRAS samples
    outagescenario =1
    sys_UC,sys_ED = add_csv_time_series!(sys_UC,sys_ED,outage_csv_location,add_scenario = outagescenario, iteration_year=iteration_year); # align outage structure to be the same as PRAS-ED (same value persist 36 times)

    # Check to see if all ThermalStandard has time series data
    # For generators that don't have time series, add artifical outage time series
    file_location = joinpath(outage_csv_location,"1/Generator_year$(iteration_year).csv")
    df_outage_profile = DataFrames.DataFrame(CSV.File(file_location));
    for d in PSY.get_components(ThermalGen,sys_UC)
        if ~("outage" in get_time_series_names(Deterministic, d))
            @info "Adding outage time series to new generators..."
            first_ts_temp_DA = first(PSY.get_time_series_multiple(sys_UC));
            len_ts_data_DA = length(first_ts_temp_DA.data);
            days_of_interest = range(1,length = len_ts_data_DA)
            num_days = length(days_of_interest);
            period_of_interest = range(days_of_interest.start==1 ? 1 : (24*days_of_interest.start)+1,length=num_days*24)
            N = length(period_of_interest);
            start_datetime_DA = PSY.IS.get_initial_timestamp(first_ts_temp_DA);
            sys_DA_res_in_hour = PSY.get_time_series_resolution(sys_UC)
            start_datetime_DA = start_datetime_DA + Dates.Hour((period_of_interest.start-1)*sys_DA_res_in_hour);
            finish_datetime_DA = start_datetime_DA +  Dates.Hour((N-1)*sys_DA_res_in_hour);
            all_timestamps = StepRange(start_datetime_DA, sys_DA_res_in_hour, finish_datetime_DA);
            allnames = PSY.get_name.(PSY.get_components(ThermalFastStartSIIP,sys_UC,x -> get_prime_mover(x) in [PrimeMovers.CT, PrimeMovers.GT] && PSY.get_fuel(x) in [ThermalFuels.NATURAL_GAS]))
            allnames = allnames[findall(x->x in names(df_outage_profile), allnames)]
            asset_name = allnames[rand(1:length(allnames))]
            @info "Use $(asset_name) time series for $(PSY.get_name(d))..."

            DA_data = Dict()
            RT_data = Dict()

            for (idx,timestamp) in enumerate(all_timestamps) #
                if (rem(idx,24) == 1)
                    if (df_outage_profile[idx,asset_name] ==0)
                        push!(DA_data,timestamp => zeros(36))
                    else
                        push!(DA_data,timestamp => ones(36))
                    end
                end

                if (df_outage_profile[idx,asset_name] ==0)
                    # If PRAS actual status is actually 0, just use 0 for RT, ReliablePSY can handle this transition
                    push!(RT_data,timestamp => zeros(2))
                else
                    # If PRAS actual status is actually 1 and DA says 1, use 1 for RT, if DA says 0, use O because ReliablePSY currently can't handle this
                    if (rem(idx,24) == 0)
                        offset = 24
                    else
                        offset = rem(idx,24)
                    end

                    if (df_outage_profile[(idx-offset+1),asset_name] ==1) # If PRAS Actual == 1 and DA ==1
                        push!(RT_data,timestamp => ones(2))
                    else
                        push!(RT_data,timestamp => ones(2)) # # If PRAS Actual == 1 and DA ==0 - Special case
                    end
                end
            end

            DA_availability_forecast = PSY.Deterministic("outage", DA_data,sys_DA_res_in_hour)
            RT_availability_forecast = PSY.Deterministic("outage", RT_data,sys_DA_res_in_hour)

            # Adding TimeSeries Data to PSY System
            PSY.add_time_series!(sys_UC,PSY.get_component(PSY.Generator,sys_UC, PSY.get_name(d)), DA_availability_forecast)
            PSY.add_time_series!(sys_ED,PSY.get_component(PSY.Generator,sys_ED, PSY.get_name(d)), RT_availability_forecast)
        end
    end
end

"""
This function finds the buses located in each zone.
"""
function find_zonal_bus(zone::String, sys::PSY.System)
    buses = PSY.get_components(PSY.Bus, sys)

    zonal_bus = nothing

    for b in buses
        if "zone_$(PSY.get_name(PSY.get_area(b)))" == zone
            zonal_bus = b
        end
    end

    return zonal_bus

end

"""
This function transforms the timeseries of PSY Systems.
"""

function transform_psy_timeseries!(sys_UC::Nothing,
                                   sys_ED::Nothing,
                                   da_resolution::Int64,
                                   rt_resolution::Int64,
                                   da_horizon::Int64,
                                   rt_horizon::Int64)
    return
end

function transform_psy_timeseries!(sys_UC::PSY.System,
                                   sys_ED::PSY.System,
                                   da_resolution::Int64,
                                   rt_resolution::Int64,
                                   da_horizon::Int64,
                                   rt_horizon::Int64)

    PSY.transform_single_time_series!(sys_UC, Int(da_horizon * 60 / da_resolution), Dates.Hour(24))
    PSY.transform_single_time_series!(sys_ED, Int(rt_horizon * 60/ rt_resolution), Dates.Hour(1))
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

    inertia_data = read_data(joinpath(simulation_dir, "markets_data", "$(reserve_penalty)_reserve_penalty", "Inertia.csv"))

    inertia_requirement = system_peak_load * inertia_data[1, "requirement_multiplier"] / PSY.get_base_power(sys)

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
        collect(PSY.get_components(PSY.HydroEnergyReservoir, sys)),
        collect(PSY.get_components(PSY.GenericBattery, sys))
        );

    PSY.add_service!(sys, inertia_reserve, contri_devices)

    # time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
    #                 PSY.SingleTimeSeries,
    #                 first(PSY.get_components(PSY.ElectricLoad, sys)),
    #                 "max_active_power"
    #                 )))
    if type == "UC"
        time_stamps = StepRange(Dates.DateTime("2018-01-01T00:00:00"), Dates.Hour(1), Dates.DateTime("2019-01-01T11:00:00"));
    elseif type == "ED"
        time_stamps = StepRange(Dates.DateTime("2018-01-01T00:00:00"), Dates.Hour(1), Dates.DateTime("2019-01-01T00:00:00"));
    else
        error("Type should be UC or ED")
    end

    ts_data = ones(length(time_stamps))
    ts = TimeSeries.TimeArray(time_stamps, ts_data);
    forecast = PSY.SingleTimeSeries("requirement", ts)
    PSY.add_time_series!(sys, inertia_reserve, forecast)

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
    time_stamps = StepRange(Dates.DateTime("2018-01-01T00:00:00"), Dates.Hour(1), Dates.DateTime("2019-01-01T11:00:00"));
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    clean_energy_ts_data = zeros(length(time_stamps))

    nodal_loads = PSY.get_components(PSY.PowerLoad, sys)
    for load in nodal_loads
        load_active_power = PSY.get_max_active_power(load)
        total_active_power += load_active_power
        # ts_data = deepcopy(PSY.get_data(PSY.get_time_series(PSY.SingleTimeSeries,
        #                                        load,
        #                                        "max_active_power"
        #                                         )))
        loadts_raw=PSY.get_data(PSY.get_time_series(PSY.Deterministic,
                                               load,
                                               "max_active_power"
                                                ))
        ts_data =[]
        for timestep in keys(loadts_raw)
            ts_data=[ts_data;loadts_raw[timestep][1:24]]
        end
        ts_data=[ts_data;loadts_raw[Dates.DateTime("2018-12-31T00:00:00")][25:36]]
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
        vcat(collect(PSY.get_components(PSY.ThermalStandard, sys, x -> (occursin("NUC", PSY.get_name(x)) || occursin("RECT", PSY.get_name(x))))),
        collect(PSY.get_components(PSY.RenewableDispatch, sys)),
        collect(PSY.get_components(PSY.HydroDispatch, sys)),
        collect(PSY.get_components(PSY.HydroEnergyReservoir, sys)),
    );

    PSY.add_service!(sys, clean_energy_reserve, contri_devices)

    ts_data = TimeSeries.TimeArray(time_stamps, clean_energy_ts_data)
    forecast = PSY.SingleTimeSeries("requirement", ts_data)
    PSY.add_time_series!(sys, clean_energy_reserve, forecast)

end

function calculate_total_load(sys::PSY.System, time_resolution::Int64)
    total_load = 0.0

    nodal_loads = PSY.get_components(PSY.PowerLoad, sys)
    for load in nodal_loads
        zone = "zone_$(PSY.get_name(PSY.get_area(PSY.get_bus(load))))"
        # ts_data = PSY.get_data(PSY.get_time_series(
        #                                           PSY.SingleTimeSeries,
        #                                           load,
        #                                           "max_active_power"
        #                                           )
        #                       )
        loadts_raw=PSY.get_data(PSY.get_time_series(PSY.Deterministic,
                                               load,
                                               "max_active_power"
                                                ))
        ts_data =[]
        for timestep in keys(loadts_raw)
            ts_data=[ts_data;loadts_raw[timestep][1:Int(8760/length(keys(loadts_raw)))]]
        end
        total_load += sum(TS.values(ts_data) * PSY.get_max_active_power(load)) * time_resolution / 60

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
        prime_mover = d.prime_mover,
        fuel = d.fuel,
        ext = d.ext
    )

    PSY.add_component!(system, new)
    for service in PSY.get_services(d)
        PSY.add_service!(new, service, system)
    end

    if ("outage" in get_time_series_names(Deterministic, d))
        newts =PSY.get_time_series(PSY.Deterministic, d, "outage")
        PSY.add_time_series!(system,new, newts)
    end

    PSY.remove_component!(system, d)
    return
end

function convert_thermal_fast_start!(system::PSY.System)
    for gen in PSY.get_components(PSY.ThermalStandard, system)
        prime_mover = PSY.get_prime_mover(gen)
        fuel = PSY.get_fuel(gen)

        target_prime_mover = [PSY.PrimeMovers.CT,PSY.PrimeMovers.GT]
        target_fuel = PSY.ThermalFuels.NATURAL_GAS

        if prime_mover in target_prime_mover && fuel == target_fuel
            convert_to_thermal_fast_start!(gen, system)
        end

    end
end

function apply_PSY_past_load_growth!(sys::PSY.System,
                               load_growth::AxisArrays.AxisArray{Float64, 1},
                               simulation_dir::String)

    # update load timeseries.
    nodal_loads = PSY.get_components(PSY.PowerLoad, sys)

    for load in nodal_loads
        zone = "load_zone_$(PSY.get_name(PSY.get_area(PSY.get_bus(load))))"
        scaled_active_power = deepcopy(PSY.get_max_active_power(load)) * (1 + load_growth[zone])
        PSY.set_max_active_power!(load, scaled_active_power)
    end

    average_load_growth = Statistics.mean(load_growth)

    # update service requirement timeseries.
    services = get_system_services(sys)
    ordc_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"ordc_products"], "; ")
    for service in services
        service_name = PSY.get_name(service)
        if !(service_name in ordc_products || service_name == "Clean_Energy")
            scaled_requirement = deepcopy(PSY.get_requirement(service)) * (1 + average_load_growth)
            PSY.set_requirement!(service, scaled_requirement)
        end

    end

    return
end