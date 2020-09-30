function prune_system_devices!(sys::PSY.System, prune_dict::Dict{Type{<:PSY.Component}, Array{AbstractString}})
    for (type, device_names) in prune_dict
        for name in device_names
            device = PSY.get_component(type, sys, name)
            PSY.remove_component!(sys, device)
        end
    end
    return"""
    THis function removes the pruned devices from PSY System.
    """
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
                                         "207_CT_1", "307_CT_1", "101_STEAM_4"]
        pruned_unit[PSY.RenewableFix] = ["308_RTPV_1", "313_RTPV_1", "313_RTPV_2", "313_RTPV_3", "313_RTPV_4", "313_RTPV_5", "313_RTPV_6", "313_RTPV_7",
                                        "313_RTPV_8", "313_RTPV_9", "313_RTPV_10", "313_RTPV_11", "313_RTPV_12", "320_RTPV_1", "320_RTPV_2", "320_RTPV_3",
                                        "313_RTPV_13", "320_RTPV_4", "320_RTPV_5", "118_RTPV_1", "118_RTPV_2", "118_RTPV_3", "118_RTPV_4", "118_RTPV_5",
                                        "118_RTPV_6", "320_RTPV_6", "118_RTPV_7", "118_RTPV_8", "118_RTPV_9", "118_RTPV_10", "213_RTPV_1"]
        pruned_unit[PSY.RenewableDispatch] = ["309_WIND_1", "212_CSP_1"]
        pruned_unit[PSY.GenericBattery] = ["313_HEAD_STORAGE"]

        return pruned_unit
    end

    """
    This function creates a PSY System for Unit Commitment.
    """
    function create_rts_sysUC(rts_dir::String,
                            base_power::Float64)
        rawsys1 = PSY.PowerSystemTableData(rts_dir, base_power, joinpath(rts_dir,"../FormattedData/SIIP/user_descriptors.yaml"),
        generator_mapping_file = joinpath(rts_dir,"../FormattedData/SIIP/generator_mapping.yaml"),
        timeseries_metadata_file = joinpath(rts_dir,"../FormattedData/SIIP/timeseries_pointers.json"))
        sys_UC = PSY.System(rawsys1; forecast_resolution = Dates.Hour(1))

        pruned_unit = specify_pruned_units()
        device_list = Dict{Type{<:PSY.Component}, Array{AbstractString}}()
        prune_system_devices!(sys_UC, pruned_unit)

        return sys_UC
    end

    """
    This function creates a PSY System for Economic Dispatch.
    """
    function create_rts_sysED(rts_dir::String,
                            base_power::Float64)
        rawsys1 = PSY.PowerSystemTableData(rts_dir, base_power, joinpath(rts_dir,"../FormattedData/SIIP/user_descriptors.yaml"),
        generator_mapping_file = joinpath(rts_dir,"../FormattedData/SIIP/generator_mapping.yaml"),
        timeseries_metadata_file = joinpath(rts_dir,"../FormattedData/SIIP/timeseries_pointers.json"))
        sys_ED = PSY.System(rawsys1; forecast_resolution = Dates.Hour(1))

        pruned_unit = specify_pruned_units()

        device_list = Dict{Type{<:PSY.Component}, Array{AbstractString}}()

        prune_system_devices!(sys_ED, pruned_unit)

        return sys_ED
    end

end

function specify_pruned_units()
    pruned_unit = Dict{Type{<:PSY.Component}, Array{AbstractString}}()
    pruned_unit[PSY.ThermalStandard] = ["115_STEAM_1", "115_STEAM_2", "315_STEAM_1", "315_STEAM_2", "315_STEAM_3", "315_STEAM_4", "315_STEAM_5",
                                     "101_CT_1", "101_CT_2", "102_CT_1", "102_CT_2",
                                     "201_CT_1", "201_CT_2", "202_CT_1", "202_CT_2",
                                     "301_CT_1", "301_CT_2", "302_CT_1", "302_CT_2",
                                     "207_CT_1", "307_CT_1", "101_STEAM_4"]
    pruned_unit[PSY.RenewableFix] = ["308_RTPV_1", "313_RTPV_1", "313_RTPV_2", "313_RTPV_3", "313_RTPV_4", "313_RTPV_5", "313_RTPV_6", "313_RTPV_7",
                                    "313_RTPV_8", "313_RTPV_9", "313_RTPV_10", "313_RTPV_11", "313_RTPV_12", "320_RTPV_1", "320_RTPV_2", "320_RTPV_3",
                                    "313_RTPV_13", "320_RTPV_4", "320_RTPV_5", "118_RTPV_1", "118_RTPV_2", "118_RTPV_3", "118_RTPV_4", "118_RTPV_5",
                                    "118_RTPV_6", "320_RTPV_6", "118_RTPV_7", "118_RTPV_8", "118_RTPV_9", "118_RTPV_10", "213_RTPV_1"]
    pruned_unit[PSY.RenewableDispatch] = ["309_WIND_1", "212_CSP_1"]
    pruned_unit[PSY.GenericBattery] = ["313_HEAD_STORAGE"]

    return pruned_unit
end

function create_rts_sysUC(rts_dir::String,
                        base_power::Float64)
    rawsys1 = PSY.PowerSystemTableData(rts_dir, base_power, joinpath(rts_dir,"../FormattedData/SIIP/user_descriptors.yaml"),
    generator_mapping_file = joinpath(rts_dir,"../FormattedData/SIIP/generator_mapping.yaml"),
    timeseries_metadata_file = joinpath(rts_dir,"../FormattedData/SIIP/timeseries_pointers.json"))
    sys_UC = PSY.System(rawsys1; forecast_resolution = Dates.Hour(1))

    pruned_unit = specify_pruned_units()
    device_list = Dict{Type{<:PSY.Component}, Array{AbstractString}}()
    prune_system_devices!(sys_UC, pruned_unit)

    return sys_UC
end

function create_rts_sysED(rts_dir::String,
                        base_power::Float64)
    rawsys1 = PSY.PowerSystemTableData(rts_dir, base_power, joinpath(rts_dir,"../FormattedData/SIIP/user_descriptors.yaml"),
    generator_mapping_file = joinpath(rts_dir,"../FormattedData/SIIP/generator_mapping.yaml"),
    timeseries_metadata_file = joinpath(rts_dir,"../FormattedData/SIIP/timeseries_pointers.json"))
    sys_ED = PSY.System(rawsys1; forecast_resolution = Dates.Hour(1))

    pruned_unit = specify_pruned_units()

    device_list = Dict{Type{<:PSY.Component}, Array{AbstractString}}()

    prune_system_devices!(sys_ED, pruned_unit)

    return sys_ED
end
