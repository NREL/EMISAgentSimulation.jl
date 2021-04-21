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
function add_device_services!(services::Vector{PSY.Service},
                              device::PSY.Device,
                              product::P) where P <: Product
    return
end

"""
This function adds ReserveUp service for PSY Devices.
"""
function add_device_services!(services::Vector{PSY.Service},
                              device::PSY.Device,
                              product::OperatingReserve{ReserveUpEMIS})

    device_zone = PSY.get_name(PSY.get_load_zone(PSY.get_bus(device)))

    for service in services
        if typeof(service) == PSY.VariableReserve{PSY.ReserveUp}
            service_zone = "zone_$(split(PSY.get_name(service), "_")[end])"
            if device_zone == service_zone
                PSY.add_service!(device, service)
            end
        end
    end

    return
end

"""
This function adds ReserveDown service for PSY Devices.
"""
function add_device_services!(services::Vector{PSY.Service},
                              device::PSY.Device,
                              product::OperatingReserve{ReserveDownEMIS})

    device_zone = PSY.get_name(PSY.get_load_zone(PSY.get_bus(device)))

    for service in services
        if typeof(service) == PSY.VariableReserve{PSY.ReserveDown}
            service_zone = "zone_$(split(PSY.get_name(service), "_")[end])"
            if device_zone == service_zone
                PSY.add_service!(device, service)
            end
        end
    end

    return
end

"""
This function does nothing if Device is not of RenewableGen type.
"""
function add_device_forecast!(simulation_dir::String,
                              sys::PSY.System,
                              device::D,
                              availability_df::Vector{Float64},
                              availability_raw_rt::Vector{Float64},
                              start_year::Int64
                              ) where D <: Union{PSY.ThermalGen, PSY.HydroGen, PSY.Storage}
    return
end

"""
This function adds forecast timeseries if Device is of RenewableGen type.
"""
function add_device_forecast!(simulation_dir::String,
                              sys::PSY.System,
                              device::D,
                              availability_raw::Vector{Float64},
                              availability_raw_rt::Vector{Float64},
                              start_year::Int64) where D <: PSY.RenewableGen


    YearAhead = collect(
    Dates.DateTime("1/1/$(start_year)  0:00:00", "d/m/y  H:M:S"):Dates.Hour(1):Dates.DateTime(
        "31/12/$(start_year)  23:00:00", "d/m/y  H:M:S",),)

        timeseries = TS.TimeArray(YearAhead, availability_raw)

        PSY.add_forecast!(sys, device, IS.Deterministic("get_rating", timeseries))

        load_n_vg_df =  read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))
        load_n_vg_df[:, get_name(device)] = availability_raw * get_device_size(device) * PSY.get_base_power(sys)

        load_n_vg_df_rt =  read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data_rt.csv"))
        load_n_vg_df_rt[:, get_name(device)] = availability_raw_rt * get_device_size(device) * PSY.get_base_power(sys)

        write_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data"), "load_n_vg_data.csv", load_n_vg_df)
        write_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data"), "load_n_vg_data_rt.csv", load_n_vg_df_rt)

    return
end

"""
This function does nothing if PSY System is not defined.
"""
function update_PSY_timeseries!(sys_UC::Nothing,
                               load_growth::AxisArrays.AxisArray{Float64, 1})
    return
end

"""
This function updates the PSY load and reserve requirment timeseries each year.
"""
function update_PSY_timeseries!(sys::PSY.System,
                               load_growth::AxisArrays.AxisArray{Float64, 1},
                               simulation_dir::String)

    # update load timeseries.
    nodal_loads = PSY.get_components(PSY.ElectricLoad, sys)
    println("FDSFSdfsd")
    quit()

    for load in nodal_loads
        zone = "zone_$(PSY.get_name(PSY.get_area(PSY.get_bus(load))))"

        #=
        ts_data = PSY.get_data(PSY.get_time_series(
                                                  PSY.SingleTimeSeries,
                                                  load,
                                                  "max_active_power"
                                                  )
                             )

                             println(ts_data)

        ts_timestamps = TS.timestamp(ts_data)
        ts_values = TS.values(ts_data)

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
        #=
        try
            ts_data= PSY.get_data(PSY.get_time_series(
                                                    PSY.SingleTimeSeries,
                                                    service,
                                                    "requirement"
                                                    )
                            )
            ts_timestamps = TS.timestamp(ts_data)
            ts_values = TS.values(ts_data)

            PSY.remove_time_series!(
                                    sys,
                                    PSY.SingleTimeSeries,
                                    service,
                                    "requirement",
                                    )

            new_ts = TS.TimeArray(ts_timestamps, ts_values * (1 + average_load_growth))

            PSY.add_time_series!(sys, service, PSY.SingleTimeSeries("requirement", new_ts))
        catch

        end
        =#

        if !(PSY.get_name(service) in ordc_products)
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
