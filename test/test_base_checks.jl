@testset "Test base years" begin
    start_year = get_start_year(case)
    data_dir = get_data_dir(case)
    scenarios = string.(get_all_scenario_names(data_dir))

    for scenario in scenarios
        for sim_year in collect(1:simulation_years)
            test_system_load_da = DataFrames.DataFrame(
                CSV.File(
                    joinpath(
                        test_system_dir,
                        "RTS_Data",
                        "timeseries_data_files",
                        scenario,
                        "sim_year_$(sim_year)",
                        "Load",
                        "DAY_AHEAD_regional_Load.csv",
                    ),
                ),
            )

            base_year = test_system_load_da[1, "Year"]
            @assert base_year <= start_year
        end
    end
end
