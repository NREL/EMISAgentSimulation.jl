using Test
using EMISAgentSimulation
using CSV
using DataFrames

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))

@testset "Project init validation" begin
    spec_dir = joinpath(PROJECT_ROOT, "config", "project_templates", "project_spec")
    init_dir = joinpath(mktempdir(), "project_init_validation")
    initialize_emis_project(spec_dir; output_dir=init_dir, reference_case_dir=nothing)

    @test validate_emis_project(init_dir) == init_dir
end

@testset "Generated reserves only include products with matching time series" begin
    tmp = mktempdir()
    ts_root = joinpath(tmp, "timeseries")
    scenario_dir = joinpath(ts_root, "scenario_1", "sim_year_1", "Reserves")
    mkpath(scenario_dir)
    CSV.write(joinpath(scenario_dir, "DAY_AHEAD_regional_SPIN.csv"), DataFrame(Year=[2020], Month=[1], Day=[1], Period=[1], SPIN=[1.0]))
    CSV.write(joinpath(scenario_dir, "REAL_TIME_regional_SPIN.csv"), DataFrame(Year=[2020], Month=[1], Day=[1], Period=[1], SPIN=[1.0]))

    reserves_path = joinpath(tmp, "reserves.csv")
    CSV.write(reserves_path, DataFrame(
        "Reserve Product" => ["SPIN", "NONSPIN"],
        "Timeframe (sec)" => [Int(300), Int(300)],
        "Requirement (MW)" => [10.0, 5.0],
        "Eligible Regions" => ["", ""],
        "Eligible Device Categories" => ["", ""],
        "Eligible Device SubCategories" => ["", ""],
        "Direction" => ["Up", "Up"],
    ))

    filtered = EMISAgentSimulation._filter_reserve_products_by_timeseries(reserves_path, ts_root)
    @test String.(filtered[:, "Reserve Product"]) == ["SPIN"]

    empty_ts_root = joinpath(tmp, "empty_timeseries")
    empty_reserves_path = joinpath(tmp, "empty_reserves.csv")
    CSV.write(empty_reserves_path, DataFrame(
        "Reserve Product" => ["SPIN", "NONSPIN"],
        "Timeframe (sec)" => [Int(300), Int(300)],
        "Requirement (MW)" => [10.0, 5.0],
        "Eligible Regions" => ["", ""],
        "Eligible Device Categories" => ["", ""],
        "Eligible Device SubCategories" => ["", ""],
        "Direction" => ["Up", "Up"],
    ))

    empty_filtered = EMISAgentSimulation._filter_reserve_products_by_timeseries(empty_reserves_path, empty_ts_root)
    @test nrow(empty_filtered) == 0
    @test names(empty_filtered) == names(DataFrame(CSV.File(empty_reserves_path; stringtype=String)))
end

@testset "Reserve filtering matches case-insensitively and canonicalizes casing" begin
    tmp = mktempdir()
    ts_root = joinpath(tmp, "timeseries")
    scenario_dir = joinpath(ts_root, "scenario_1", "sim_year_1", "Reserves")
    mkpath(scenario_dir)
    # Real-world data uses Title-cased file names (e.g. "Spin", "Reg_Up") while the
    # PSY reserve service names are upper-cased (e.g. "SPIN", "REG_UP").
    CSV.write(joinpath(scenario_dir, "DAY_AHEAD_regional_Spin.csv"), DataFrame(Year=[2020], Month=[1], Day=[1], Period=[1], Spin=[1.0]))
    CSV.write(joinpath(scenario_dir, "DAY_AHEAD_regional_Reg_Up.csv"), DataFrame(Year=[2020], Month=[1], Day=[1], Period=[1], Reg_Up=[1.0]))

    reserves_path = joinpath(tmp, "case_reserves.csv")
    CSV.write(reserves_path, DataFrame(
        "Reserve Product" => ["SPIN", "REG_UP", "REG_DN", "NONSPIN"],
        "Timeframe (sec)" => [Int(300), Int(300), Int(300), Int(300)],
        "Requirement (MW)" => [10.0, 5.0, 5.0, 5.0],
        "Eligible Regions" => ["", "", "", ""],
        "Eligible Device Categories" => ["", "", "", ""],
        "Eligible Device SubCategories" => ["", "", "", ""],
        "Direction" => ["Up", "Up", "Down", "Up"],
    ))

    filtered = EMISAgentSimulation._filter_reserve_products_by_timeseries(reserves_path, ts_root)
    # SPIN/REG_UP match case-insensitively and are canonicalized to the on-disk casing;
    # REG_DN has no matching file (only "Reg_Down" would, which is a different name) and
    # NONSPIN has no file at all, so both are dropped.
    @test String.(filtered[:, "Reserve Product"]) == ["Spin", "Reg_Up"]
end
