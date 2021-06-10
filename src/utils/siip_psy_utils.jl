function get_name(device::D) where D <: PSY.Device
    return device.name
end

"""
This function returns all generation and storage technologies included in the PSY System.
"""
function get_all_techs(sys::PSY.System)
    sys_gens = PSY.get_components(PSY.Generator, sys)
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
    time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
                                                    PSY.SingleTimeSeries,
                                                    first(PSY.get_components(PSY.ElectricLoad, sys_UC)),
                                                    "max_active_power"
                                                    )))

    intervals = Int(24 * 60 / da_resolution)
    append!(availability_raw, availability_raw[(length(availability_raw) - intervals + 1):end])
    data = Dict(time_stamps[i] => availability_raw[i:(i + intervals - 1)] for i in 1:intervals:length(time_stamps))
    forecast = PSY.Deterministic("max_active_power", data, Dates.Minute(da_resolution))
    PSY.add_time_series!(sys_UC, device_UC, forecast)

    ######### Adding to ED##########
    time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
                                                    PSY.SingleTimeSeries,
                                                    first(PSY.get_components(PSY.ElectricLoad, sys_ED)),
                                                    "max_active_power"
                                                    )))

    intervals = Int(60 / rt_resolution)
    append!(availability_raw_rt, availability_raw_rt[(length(availability_raw) - intervals + 1):end])
    data = Dict(time_stamps[i] => availability_raw_rt[i:(i + intervals - 1)] for i in 1:intervals:length(time_stamps))
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
                               simulation_dir::String,
                               type::String,
                               iteration_year::Int64,
                               da_resolution::Int64,
                               rt_resolution::Int64)

    # update load timeseries.
    nodal_loads = PSY.get_components(PSY.ElectricLoad, sys)

    for load in nodal_loads
        zone = "zone_$(PSY.get_name(PSY.get_area(PSY.get_bus(load))))"


        ts_data = PSY.get_data(PSY.get_time_series(
                                                  PSY.SingleTimeSeries,
                                                  load,
                                                  "max_active_power"
                                                  )
                             )
        #=
        ts_timestamps = TS.timestamp(ts_data)
        ts_values = TS.values(ts_data)
        println(ts_values)

        PSY.remove_time_series!(sys,
                                PSY.SingleTimeSeries,
                                load,
                                "max_active_power"
                                )

        new_ts = TS.TimeArray(ts_timestamps, ts_values * (1 + load_growth[zone]))
        PSY.add_time_series!(sys, load, PSY.SingleTimeSeries("max_active_power", new_ts))
        =#
        scaled_active_power = deepcopy(PSY.get_max_active_power(load)) * (1 + load_growth[zone])

        PSY.set_max_active_power!(load, scaled_active_power)

    end

    average_load_growth = Statistics.mean(load_growth)

    # update service requirement timeseries.
    services = get_system_services(sys)
    ordc_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"ordc_products"], "; ")
    for service in services
        service_name = PSY.get_name(service)
        if service_name in ordc_products
            time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
                                                    PSY.SingleTimeSeries,
                                                    first(PSY.get_components(PSY.ElectricLoad, sys)),
                                                    "max_active_power"
                                                    )))

            PSY.remove_time_series!(sys, PSY.Deterministic, service, "variable_cost")
            if type == "UC"
                product_ts_raw = read_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "$(service_name)_$(iteration_year - 1).csv"))[:, service_name]
                product_data_ts = process_ordc_data_for_siip(product_ts_raw)
                intervals = Int(24 * 60 / da_resolution)
                append!(product_data_ts, product_data_ts[(length(product_data_ts) - intervals + 1):end])
                data = Dict(time_stamps[i] => product_data_ts[i:(i + intervals - 1)] for i in 1:intervals:length(time_stamps))
                forecast = PSY.Deterministic("variable_cost", data, Dates.Minute(da_resolution))
            elseif type == "ED"
                product_ts_raw = read_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "$(service_name)_REAL_TIME_$(iteration_year - 1                                       ).csv"))[:, service_name]
                product_data_ts = process_ordc_data_for_siip(product_ts_raw)
                intervals =  Int(60 / rt_resolution)
                append!(product_data_ts, product_data_ts[(length(product_data_ts) - intervals + 1):end])
                data = Dict(time_stamps[i] => product_data_ts[i:(i + intervals  - 1)] for i in 1:intervals:length(time_stamps))
                forecast = PSY.Deterministic("variable_cost", data, Dates.Minute(rt_resolution))
            else
                error("Type should be UC or ED")
            end

            PSY.add_time_series!(sys, service, forecast)
        else
            scaled_requirement = deepcopy(PSY.get_requirement(service)) * (1 + average_load_growth)
            PSY.set_requirement!(service, scaled_requirement)
        end

    end

    return
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
                                   rt_resolution::Int64)
    return
end

function transform_psy_timeseries!(sys_UC::PSY.System,
                                   sys_ED::PSY.System,
                                   da_resolution::Int64,
                                   rt_resolution::Int64)

    PSY.transform_single_time_series!(sys_UC, Int(24 * 60 / da_resolution), Dates.Hour(24))
    PSY.transform_single_time_series!(sys_ED, Int(60 / rt_resolution), Dates.Hour(1))
    return
end

function add_psy_inertia!(simulation_dir::String,
                          sys::Nothing,
                          system_peak_load::Float64)

    return

end

function add_psy_inertia!(simulation_dir::String,
                          sys::PSY.System,
                          system_peak_load::Float64)

    inertia_data = read_data(joinpath(simulation_dir, "markets_data", "Inertia.csv"))

    inertia_requirement = system_peak_load * inertia_data[1, "requirement_multiplier"]

    ####### Adding Inertia reserve
    inertia_reserve = PSY.VariableReserve{PSY.ReserveSymmetric}(
        "Inertia",
        true,
        60,
        inertia_requirement,
    )
    contri_devices =
        vcat(collect(PSY.get_components(PSY.ThermalStandard, sys)),
        collect(PSY.get_components(PSY.RenewableDispatch, sys)),
        collect(PSY.get_components(PSY.HydroDispatch, sys)),
        collect(PSY.get_components(PSY.HydroEnergyReservoir, sys)),
        collect(PSY.get_components(PSY.GenericBattery, sys))
        );

    PSY.add_service!(sys, inertia_reserve, contri_devices)

    time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
                    PSY.SingleTimeSeries,
                    first(PSY.get_components(PSY.ElectricLoad, sys)),
                    "max_active_power"
                    )))

    ts_data = ones(length(time_stamps))
    ts = TimeSeries.TimeArray(time_stamps, ts_data);
    forecast = PSY.SingleTimeSeries("requirement", ts)
    PSY.add_time_series!(sys, inertia_reserve, forecast)

end
