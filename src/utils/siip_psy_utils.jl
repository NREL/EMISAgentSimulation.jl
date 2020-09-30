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
                              start_year::Int64) where D <: PSY.RenewableGen


    YearAhead = collect(
    Dates.DateTime("1/1/$(start_year)  0:00:00", "d/m/y  H:M:S"):Dates.Hour(1):Dates.DateTime(
        "31/12/$(start_year)  23:00:00", "d/m/y  H:M:S",),)

        timeseries = TS.TimeArray(YearAhead, availability_raw)

        PSY.add_forecast!(sys, device, IS.Deterministic("get_rating", timeseries))

        load_n_vg_df =  read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))
        load_n_vg_df[:, get_name(device)] = availability_raw * get_device_size(device) * PSY.get_base_power(sys)

        write_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data"), "load_n_vg_data.csv", load_n_vg_df)

    return
end

"""
This function does nothing if PSY System is not defined.
"""
function update_PSY_timeseries!(sys_UC::Nothing,
                               load_growth::AxisArrays.AxisArray{Float64, 1},
                               days::Int64)
    return
end

"""
This function updates the PSY load and reserve requirment timeseries each year.
"""
function update_PSY_timeseries!(sys::PSY.System,
                               load_growth::AxisArrays.AxisArray{Float64, 1},
                               num_days::Int64)

    num_hours = 24
    initial_time = PSY.get_forecasts_initial_time(sys)

    # update load forecasts.
    nodal_loads = PSY.get_components(PSY.ElectricLoad, sys)
    for load in nodal_loads
        zone = "zone_$(PSY.get_name(PSY.get_area(PSY.get_bus(load))))"
        ts_data= PSY.get_data(PSY.get_forecast(
            PSY.Deterministic,
            load,
            initial_time,
            "get_max_active_power",
            num_days * num_hours,
        ))
        ts_timestamps = TS.timestamp(ts_data)
        ts_values = TS.values(ts_data)

        PSY.remove_forecast!(PSY.Deterministic,
            sys,
            load,
            initial_time,
            "get_max_active_power",
            )

        new_ts = TS.TimeArray(ts_timestamps, ts_values * (1 + load_growth[zone]))
        PSY.add_forecast!(sys, load, IS.Deterministic("get_max_active_power", new_ts))

    end

    # update service requirement forecasts.
    services = get_system_services(sys)
    for service in services
        zone = "zone_$(last(split(PSY.get_name(service), "_")[end], 1))"

        ts_data= PSY.get_data(PSY.get_forecast(
            PSY.Deterministic,
            service,
            initial_time,
            "get_requirement",
            num_days * num_hours,
        ))
        ts_timestamps = TS.timestamp(ts_data)
        ts_values = TS.values(ts_data)

        PSY.remove_forecast!(PSY.Deterministic,
            sys,
            service,
            initial_time,
            "get_requirement",
            )

        new_ts = TS.TimeArray(ts_timestamps, ts_values * (1 + load_growth[zone]))

        PSY.add_forecast!(sys, service, IS.Deterministic("get_requirement", new_ts))
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
