"""
This function reads the test system time series and returns representative hour weight.
"""
function read_test_system(data_dir::String,
                          test_system_dir::String,
                          base_dir::String,
                          annual_growth::AxisArrays.AxisArray{Float64, 2},
                          start_year::Int64,
                          n_rep_days::Int64)

    test_sys_hour_weight = nothing
    zones = nothing
    representative_days = nothing
    rep_hour_weight = nothing
    system_peak_load = nothing
    test_sys_hour_weight = nothing
    zonal_lines = nothing

    if occursin("RTS", test_system_dir)
        zones,
        representative_days,
        rep_hour_weight,
        system_peak_load,
        test_sys_hour_weight,
        zonal_lines = read_rts(data_dir,
                               test_system_dir,
                               base_dir,
                               annual_growth,
                               start_year,
                               n_rep_days)
    end

    return zones, representative_days, rep_hour_weight, system_peak_load, test_sys_hour_weight, zonal_lines
end
