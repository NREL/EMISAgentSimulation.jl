using Test
using DataFrames
using CSV
using EMISAgentSimulation

@testset "SystemConfig" begin
    tmpdir = mktempdir()
    config_dir = joinpath(tmpdir, "system_config")
    mkpath(config_dir)

    DataFrames.DataFrame(
        zone_id = [1, 2],
        zone_name = ["zone_a", "zone_b"],
        load_column = ["load_zone_1", "load_zone_2"],
    ) |> df -> CSV.write(joinpath(config_dir, "zones.csv"), df)

    DataFrames.DataFrame(
        scenario_name = ["scenario_1", "scenario_2"],
        probability = [0.8, 0.2],
        pcm_label = ["baseline", "central"],
        weather_year = [2020, 2021],
    ) |> df -> CSV.write(joinpath(config_dir, "scenarios.csv"), df)

    DataFrames.DataFrame(
        unit_type = ["WT", "CC"],
        category = ["Wind", "Gas CC"],
        class = ["renewable", "thermal"],
        capacity_eligible = [true, true],
        rec_eligible = [true, false],
        mopr_exempt = [true, false],
        macrs_years = [20, 25],
        duration_hr = [missing, missing],
    ) |> df -> CSV.write(joinpath(config_dir, "technologies.csv"), df)

    DataFrames.DataFrame(
        device_type = ["ThermalStandard"],
        device_name = ["AUSTIN_1"],
    ) |> df -> CSV.write(joinpath(config_dir, "devices_to_remove.csv"), df)

    cfg = load_system_config(tmpdir)
    @test cfg.zone_names == ["zone_a", "zone_b"]
    @test get_zone_name(cfg, 1) == "zone_a"
    @test get_zone_name(cfg, 2) == "zone_b"
    @test cfg.scenario_names == ["scenario_1", "scenario_2"]
    @test get_default_scenario(cfg) == "scenario_1"
    @test cfg.default_rts_load == EMISAgentSimulation.DEFAULT_RTS_LOAD

    fallback = load_system_config(joinpath(tmpdir, "missing"))
    @test fallback.zone_names == ["zone_$(index)" for index in 1:8]
    @test get_default_scenario(fallback) == "scenario_1"
    @test EMISAgentSimulation.get_scenario_pcm_label(fallback, "scenario_1") == "baseline"
    @test EMISAgentSimulation.get_scenario_pcm_label(fallback, "scenario_2") == "central"
    @test EMISAgentSimulation.get_scenario_pcm_label(fallback, "scenario_3") == "ira"
    @test fallback.default_rts_load == EMISAgentSimulation.DEFAULT_RTS_LOAD
end
