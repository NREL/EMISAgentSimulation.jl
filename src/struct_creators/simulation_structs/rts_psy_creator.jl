function prune_system_devices!(sys::PSY.System, prune_dict::Dict{Type{<:PSY.Component}, Array{AbstractString}})
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
    pruned_unit[PSY.ThermalStandard] = ["115_STEAM_1", "115_STEAM_2", "315_STEAM_1", "315_STEAM_2", "315_STEAM_3", "315_STEAM_4", "315_STEAM_5",
                                         "101_CT_1", "101_CT_2", "102_CT_1", "102_CT_2",
                                         "201_CT_1", "201_CT_2", "202_CT_1", "202_CT_2",
                                         "301_CT_1", "301_CT_2", "302_CT_1", "302_CT_2",
                                         "207_CT_1", "307_CT_1", "101_STEAM_4",
                                         "123_STEAM_3", "223_STEAM_1", "223_STEAM_3"]
    pruned_unit[PSY.RenewableFix] = ["308_RTPV_1", "313_RTPV_1", "313_RTPV_2", "313_RTPV_3", "313_RTPV_4", "313_RTPV_5", "313_RTPV_6", "313_RTPV_7",
                                        "313_RTPV_8", "313_RTPV_9", "313_RTPV_10", "313_RTPV_11", "313_RTPV_12", "320_RTPV_1", "320_RTPV_2", "320_RTPV_3",
                                        "313_RTPV_13", "320_RTPV_4", "320_RTPV_5", "118_RTPV_1", "118_RTPV_2", "118_RTPV_3", "118_RTPV_4", "118_RTPV_5",
                                        "118_RTPV_6", "320_RTPV_6", "118_RTPV_7", "118_RTPV_8", "118_RTPV_9", "118_RTPV_10", "213_RTPV_1"]
    pruned_unit[PSY.RenewableDispatch] = ["309_WIND_1", "212_CSP_1"]
    pruned_unit[PSY.Generator] = ["114_SYNC_COND_1", "314_SYNC_COND_1", "214_SYNC_COND_1"]
    pruned_unit[PSY.GenericBattery] = ["313_STORAGE_1"]

    return pruned_unit
    end


function create_rts_sys(rts_dir::String,
                        base_power::Float64,
                        simulation_dir::String,
                        da_resolution::Int64,
                        rt_resolution::Int64)

    da_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"da_products"], "; ")

    rts_src_dir = joinpath(rts_dir, "RTS_Data", "SourceData")
    rts_siip_dir = joinpath(rts_dir, "RTS_Data", "FormattedData", "SIIP");

    rawsys = PSY.PowerSystemTableData(
            rts_src_dir,
            base_power,
            joinpath(rts_siip_dir, "user_descriptors.yaml"),
            timeseries_metadata_file = joinpath(rts_siip_dir, "timeseries_pointers.json"),
            );

    sys_UC = PSY.System(rawsys; time_series_resolution = Dates.Minute(da_resolution));

    services_UC = get_system_services(sys_UC)

    for service in services_UC
        if !(PSY.get_name(service) in da_products)
            PSY.remove_component!(sys_UC, service)
        end
    end

    rt_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"rt_products"], "; ")

    sys_ED = PSY.System(rawsys; time_series_resolution = Dates.Minute(rt_resolution));

    services_ED = get_system_services(sys_ED)

    for service in services_ED
        if !(PSY.get_name(service) in rt_products)
            PSY.remove_component!(sys_ED, service)
        end
    end

    pruned_unit = specify_pruned_units()
    prune_system_devices!(sys_UC, pruned_unit)
    prune_system_devices!(sys_ED, pruned_unit)

    return sys_UC, sys_ED
end

function remove_vre_gens!(sys::PSY.System)
    for gen in get_all_techs(sys)
        if typeof(gen) == PSY.RenewableDispatch
            println(PSY.get_name(gen))
            println(PSY.get_ext(gen))
            PSY.remove_component!(sys, gen)
        end
    end
end
