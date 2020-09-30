"""
This function reads the RTS-GMLC timeseries data and parses it into the format
compatible with AgentSimulation
"""

function read_rts(data_dir::String,
                  test_system_dir::String,
                  base_dir::String,
                  annual_growth::AxisArrays.AxisArray{Float64, 2},
                  start_year::Int64,
                  n_rep_days::Int64)

    test_system_load = DataFrames.DataFrame(CSV.File(joinpath(test_system_dir, "RTS_Data", "timeseries_data_files", "Load", "DAY_AHEAD_regional_Load.csv")))

    base_year = test_system_load[1, "Year"]

    test_sys_num_hours = DataFrames.nrow(test_system_load)
    if test_sys_num_hours >= 8760
        test_sys_hour_weight = ones(test_sys_num_hours)
    else
        test_sys_hour_weight = ones(test_sys_num_hours) * 8760 / test_sys_num_hours
    end

    zone_numbers = names(test_system_load)[5:end]
    zones = ["zone_$(i)" for i in zone_numbers]

    test_system_reserve_up = test_system_load[:, 1:4]
    test_system_reserve_down = test_system_load[:, 1:4]

    for zone in zone_numbers
        reserve_up_df = DataFrames.DataFrame(CSV.File(joinpath(test_system_dir, "RTS_Data", "timeseries_data_files", "Reserves", "DAY_AHEAD_regional_Reg_Up_R$(zone).csv")))

        test_system_reserve_up[:, zone] = reserve_up_df[:, "Reg_Up_R$(zone)"]

        reserve_down_df = DataFrames.DataFrame(CSV.File(joinpath(test_system_dir, "RTS_Data", "timeseries_data_files", "Reserves", "DAY_AHEAD_regional_Reg_Down_R$(zone).csv")))
        test_system_reserve_down[:, zone] = zeros(DataFrames.nrow(test_system_reserve_down))
        for d in 1:Int(DataFrames.nrow(test_system_reserve_down)/24)
            for h in 1:24
                test_system_reserve_down[(d - 1) * 24 + h, zone] = reserve_down_df[d, Symbol(h)]
            end
        end
    end

    scaled_test_sys_load = deepcopy(test_system_load)
    scaled_test_sys_reserve_up = deepcopy(test_system_reserve_up)
    scaled_test_sys_reserve_down = deepcopy(test_system_reserve_down)

    scaled_test_sys_load[:, "Year"] = fill(start_year, DataFrames.nrow(scaled_test_sys_load))
    remove_leap_day!(scaled_test_sys_load, start_year)

    scaled_test_sys_reserve_up[:, "Year"] = fill(start_year, DataFrames.nrow(scaled_test_sys_reserve_up))
    remove_leap_day!(scaled_test_sys_reserve_up, start_year)

    scaled_test_sys_reserve_down[:, "Year"] = fill(start_year, DataFrames.nrow(scaled_test_sys_reserve_down))
    remove_leap_day!(scaled_test_sys_reserve_down, start_year)

    for y in 1:(start_year - base_year)
        for zone in zone_numbers
            scaled_test_sys_load[:, "$(zone)"] =  scaled_test_sys_load[:, "$(zone)"] * (1 + annual_growth[y, Symbol("load_zone_$(zone)")])
            scaled_test_sys_reserve_up[:, "$(zone)"] =  scaled_test_sys_reserve_up[:, "$(zone)"] * (1 + annual_growth[y, Symbol("load_zone_$(zone)")])
            scaled_test_sys_reserve_down[:, "$(zone)"] =  scaled_test_sys_reserve_down[:, "$(zone)"] * (1 + annual_growth[y, Symbol("load_zone_$(zone)")])
        end
    end

    net_load_df = scaled_test_sys_load[:, 1:4]

    for zone in zone_numbers
        net_load_df[:, "load_zone_$(zone)"] = scaled_test_sys_load[:, zone]
    end

    existing_generator_data = DataFrames.DataFrame(CSV.File(joinpath(test_system_dir, "RTS_Data", "SourceData", "gen.csv")))

    wind_timeseries_file = joinpath(test_system_dir, "RTS_Data", "timeseries_data_files", "WIND", "DAY_AHEAD_wind.csv")
    if isfile(wind_timeseries_file)
        wind_timeseries_data = DataFrames.DataFrame(CSV.File(wind_timeseries_file))
        remove_leap_day!(wind_timeseries_data, start_year)
    end

    pv_timeseries_file = joinpath(test_system_dir, "RTS_Data", "timeseries_data_files", "PV", "DAY_AHEAD_pv.csv")
    if isfile(pv_timeseries_file)
        pv_timeseries_data = DataFrames.DataFrame(CSV.File(pv_timeseries_file))
        remove_leap_day!(pv_timeseries_data, start_year)
    end

    rtpv_timeseries_file = joinpath(test_system_dir, "RTS_Data", "timeseries_data_files", "RTPV", "DAY_AHEAD_rtpv.csv")
    if isfile(rtpv_timeseries_file)
        rtpv_timeseries_data = DataFrames.DataFrame(CSV.File(rtpv_timeseries_file))
        remove_leap_day!(rtpv_timeseries_data, start_year)
    end

    hydro_timeseries_file = joinpath(test_system_dir, "RTS_Data", "timeseries_data_files", "Hydro", "DAY_AHEAD_hydro.csv")
    if isfile(hydro_timeseries_file)
        hydro_timeseries_data = DataFrames.DataFrame(CSV.File(hydro_timeseries_file))
        remove_leap_day!(hydro_timeseries_data, start_year)
    end

    gen_availability_df = scaled_test_sys_load[:, 1:4]

    for i in 1:DataFrames.nrow(existing_generator_data)
        if existing_generator_data[i, "Unit Type"] == "WIND"
            gen_availability_df[:, existing_generator_data[i, "GEN UID"]] = wind_timeseries_data[:, existing_generator_data[i, "GEN UID"]] / existing_generator_data[i, "PMax MW"]
            net_load_df[:, existing_generator_data[i, "GEN UID"]] = wind_timeseries_data[:, existing_generator_data[i, "GEN UID"]]
        elseif existing_generator_data[i, "Unit Type"] == "PV"
            gen_availability_df[:, existing_generator_data[i, "GEN UID"]] = pv_timeseries_data[:, existing_generator_data[i, "GEN UID"]] / existing_generator_data[i, "PMax MW"]
            net_load_df[:, existing_generator_data[i, "GEN UID"]] = pv_timeseries_data[:, existing_generator_data[i, "GEN UID"]]
        elseif existing_generator_data[i, "Unit Type"] == "RTPV"
            gen_availability_df[:, existing_generator_data[i, "GEN UID"]] = rtpv_timeseries_data[:, existing_generator_data[i, "GEN UID"]] / existing_generator_data[i, "PMax MW"]
        elseif existing_generator_data[i, "Unit Type"] == "HYDRO"
            gen_availability_df[:, existing_generator_data[i, "GEN UID"]] = hydro_timeseries_data[:, existing_generator_data[i, "GEN UID"]] / existing_generator_data[i, "PMax MW"]
        else
            gen_availability_df[:, existing_generator_data[i, "GEN UID"]] = ones(DataFrames.nrow(gen_availability_df))
        end
    end

    write_data(joinpath(data_dir, "timeseries_data_files", "Availability"), "DAY_AHEAD_availability.csv", gen_availability_df)
    write_data(joinpath(data_dir, "timeseries_data_files", "Net Load Data"), "load_n_vg_data.csv", net_load_df)

    representative_days = find_representative_days(data_dir, base_dir, n_rep_days)
    rep_load_data = filter(row -> in(Dates.Date(row[:Year], row[:Month], row[:Day]), keys(representative_days)), scaled_test_sys_load)

    rep_hour_weight = zeros(DataFrames.nrow(rep_load_data))

    for i in 1:length(rep_hour_weight)
        rep_hour_weight[i] = representative_days[Dates.Date(rep_load_data[i, :Year], rep_load_data[i, :Month], rep_load_data[i, :Day])]
    end

    write_data(joinpath(data_dir, "timeseries_data_files", "Load"), "load_0.csv", scaled_test_sys_load)
    write_data(joinpath(data_dir, "timeseries_data_files", "Load"), "rep_load_0.csv", rep_load_data)

    system_peak_load = maximum(sum(scaled_test_sys_load[:, zone] for zone in zone_numbers))

    rep_reserve_up_data = filter(row -> in(Dates.Date(row[:Year], row[:Month], row[:Day]), keys(representative_days)), scaled_test_sys_reserve_up)
    write_data(joinpath(data_dir, "timeseries_data_files", "Reserves"), "reserve_up_0.csv", scaled_test_sys_reserve_up)
    write_data(joinpath(data_dir, "timeseries_data_files", "Reserves"), "rep_reserve_up_0.csv", rep_reserve_up_data)

    rep_reserve_down_data = filter(row -> in(Dates.Date(row[:Year], row[:Month], row[:Day]), keys(representative_days)), scaled_test_sys_reserve_down)
    write_data(joinpath(data_dir, "timeseries_data_files", "Reserves"), "reserve_down_0.csv", scaled_test_sys_reserve_down)
    write_data(joinpath(data_dir, "timeseries_data_files", "Reserves"), "rep_reserve_down_0.csv", rep_reserve_down_data)

    reserve_params_df = read_data(joinpath(test_system_dir, "RTS_Data", "SourceData", "reserves.csv"))

    reserve_up_eligible_categories = String[]
    reserve_down_eligible_categories = String[]

    for zone in zone_numbers
        push!(reserve_up_eligible_categories, reserve_params_df[findfirst(x-> x == "Reg_Up_R$(zone)", reserve_params_df[:, "Reserve Product"]), "Eligible Device SubCategories"])
        push!(reserve_down_eligible_categories, reserve_params_df[findfirst(x-> x == "Reg_Down_R$(zone)", reserve_params_df[:, "Reserve Product"]), "Eligible Device SubCategories"])
    end

    reserve_up_mkt_df = read_data(joinpath(data_dir, "markets_data", "reserve_up_mkt_param.csv"))
    reserve_down_mkt_df = read_data(joinpath(data_dir, "markets_data", "reserve_down_mkt_param.csv"))

    reserve_up_mkt_df[!, "eligible categories"] = reserve_up_eligible_categories
    reserve_down_mkt_df[!, "eligible categories"] = reserve_down_eligible_categories

    write_data(joinpath(data_dir, "markets_data"), "reserve_up_mkt_param.csv", reserve_up_mkt_df)
    write_data(joinpath(data_dir, "markets_data"), "reserve_down_mkt_param.csv", reserve_down_mkt_df)


    # Create zonal lines

    zonal_lines = ZonalLine[]

    branches = read_data(joinpath(test_system_dir, "RTS_Data", "SourceData", "branch.csv"))
    dc_branches = read_data(joinpath(test_system_dir, "RTS_Data", "SourceData", "dc_branch.csv"))

    for b in 1:DataFrames.nrow(branches)
        from_bus = "$(branches[b, "From Bus"])"
        from_zone = "zone_$(first(from_bus, 1))"

        to_bus = "$(branches[b, "To Bus"])"
        to_zone = "zone_$(first(to_bus, 1))"

        similar_line = filter(l -> (in(from_zone, [get_from_zone(l), get_to_zone(l)]) && in(to_zone, [get_from_zone(l), get_to_zone(l)])), zonal_lines)

        if length(similar_line) < 1
            push!(zonal_lines, ZonalLine(branches[b, "UID"], from_zone, to_zone, branches[b, "Cont Rating"]))
        else
            set_active_power_limit!(similar_line[1], get_active_power_limit(similar_line[1]) + branches[b, "Cont Rating"])
        end
    end

    for b in 1:DataFrames.nrow(dc_branches)
        from_bus = "$(dc_branches[b, "From Bus"])"
        from_zone = "zone_$(first(from_bus, 1))"

        to_bus = "$(dc_branches[b, "To Bus"])"
        to_zone = "zone_$(first(to_bus, 1))"

        similar_line = filter(l -> (in(from_zone, [get_from_zone(l), get_to_zone(l)]) && in(to_zone, [get_from_zone(l), get_to_zone(l)])), zonal_lines)

        if length(similar_line) < 1
            push!(zonal_lines, ZonalLine(dc_branches[b, "UID"], from_zone, to_zone, dc_branches[b, "MW Load"]))
        else
            set_active_power_limit!(similar_line[1], get_active_power_limit(similar_line[1]) + dc_branches[b, "MW Load"])
        end
    end

    return zones, representative_days, rep_hour_weight, system_peak_load, test_sys_hour_weight, zonal_lines
end
