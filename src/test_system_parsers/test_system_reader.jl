"""
This function reads the test system time series and returns representative hour weight.
"""
function read_test_system(data_dir::String,
                          test_system_dir::String,
                          base_dir::String,
                          test_system_load_da::DataFrames.DataFrame,
                          test_system_load_rt::DataFrames.DataFrame,
                          base_year::Int64,
                          annual_growth_past::AxisArrays.AxisArray{Float64, 2},
                          start_year::Int64,
                          rep_period_interval::Int64,
                          n_rep_periods::Int64,
                          rep_checkpoint::Int64)

    test_sys_hour_weight = nothing
    zones = nothing
    representative_days = nothing
    chron_weights = nothing
    rep_hour_weight = nothing
    system_peak_load = nothing
    test_sys_hour_weight = nothing
    zonal_lines = nothing

    if occursin("RTS", test_system_dir)
        zones,
        representative_periods,
        rep_hour_weight,
        chron_weights,
        system_peak_load,
        test_sys_hour_weight,
        zonal_lines = read_rts(data_dir,
                               test_system_dir,
                               base_dir,
                               test_system_load_da,
                               test_system_load_rt,
                               base_year,
                               annual_growth_past,
                               start_year,
                               rep_period_interval,
                               n_rep_periods,
                               rep_checkpoint)
    end

    return zones, representative_periods, rep_hour_weight, chron_weights, system_peak_load, test_sys_hour_weight, zonal_lines
end
