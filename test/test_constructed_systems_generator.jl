using Test
using EMISAgentSimulation

@testset "constructed system generator contract" begin
    base_dir = mktempdir()
    scenario = "scenario_1"
    sim_year = 1
    market_stage = "dayahead"
    horizon = 24
    interval = 1

    md_path = canonical_constructed_system_path(base_dir, scenario, sim_year, market_stage, horizon, interval; kind="md")
    sys_path = canonical_constructed_system_path(base_dir, scenario, sim_year, market_stage, horizon, interval; kind="sys")
    forecast_path = canonical_constructed_system_path(base_dir, scenario, sim_year, market_stage, horizon, interval; kind="forecast")

    create_sys_with_timeseries(
        nothing,
        base_dir,
        scenario,
        sim_year,
        market_stage,
        horizon,
        interval,
        sys_path;
        forecast_count = 3,
    )

    bundle = canonical_constructed_system_paths(base_dir, scenario, sim_year, market_stage, horizon, interval)
    @test endswith(bundle.sys, "DA_sys_EMIS_24hor_1int.json")
    @test endswith(bundle.md, "MD_sys_EMIS_24hor_1int.json")
    @test endswith(bundle.forecast, "MD_num_forecast_24hor_1int.txt")

    @test isfile(sys_path)
    @test isfile(forecast_path)
    @test occursin("constructed", read(sys_path, String))
    @test read(forecast_path, String) == "3\n"
end

@testset "canonical profile input contract" begin
    profile_file = joinpath(mktempdir(), "load.csv")
    open(profile_file, "w") do io
        println(io, "load_a,load_b")
        println(io, "1.0,2.0")
        println(io, "3.0,4.0")
    end

    @test EMISAgentSimulation._profile_values(profile_file, "load_a") == [1.0, 3.0]
    @test EMISAgentSimulation._profile_values(profile_file, :load_b) == [2.0, 4.0]
    @test EMISAgentSimulation._cyclic_profile([1.0, 2.0], 2, 5) == [2.0, 1.0, 2.0, 1.0, 2.0]
    @test_throws ErrorException EMISAgentSimulation._profile_values(profile_file, nothing)
    @test_throws ErrorException EMISAgentSimulation._profile_values(profile_file, "missing")

    @test endswith(canonical_timeseries_path("/tmp/input", "scenario_1", 2, "da", "load"),
        "scenario_1/sim_year_2/Load/DAY_AHEAD_regional_Load.csv")
    @test endswith(canonical_timeseries_path("/tmp/input", "scenario_1", 2, "rt", "wind"),
        "scenario_1/sim_year_2/WIND/REAL_TIME_wind.csv")
    @test endswith(canonical_timeseries_path("/tmp/input", "scenario_1", 2, "da", "pv"),
        "scenario_1/sim_year_2/PV/DAY_AHEAD_pv.csv")
    @test endswith(canonical_timeseries_path("/tmp/input", "scenario_1", 2, "da", "reserve"; product="SPIN"),
        "scenario_1/sim_year_2/Reserves/DAY_AHEAD_regional_SPIN.csv")

    defaults_file = joinpath(mktempdir(), "timeseries_defaults.csv")
    open(defaults_file, "w") do io
        println(io, "kind,market_stage,directory,filename")
        println(io, "load,dayahead,Demand,DA_load.csv")
    end
    @test endswith(canonical_timeseries_path(
        "/tmp/input", "scenario_1", 2, "da", "load"; defaults_file=defaults_file
    ), "scenario_1/sim_year_2/Demand/DA_load.csv")
end
