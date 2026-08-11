# Seasonal Capacity Markets — feature validation tests
#
# Validates every new struct field, accessor, and utility function introduced by
# Capacity_markets_implementation.md, organized under the same phase numbering as
# that document (see /projects/gmlcmarkets/GS_EMIS_packages_AAYAD/Capacity_markets_implementation.md).
#
# Written before implementation (TDD): as of this writing, none of the phases have
# landed. Testsets below will report "Error" (UndefVarError / MethodError from a
# function or field that doesn't exist yet) or "Fail" (wrong type/value) until the
# corresponding phase is implemented — that is expected and intentional, not a sign
# this file is broken. Run this file after each phase lands and record pass/fail
# per testset in implementation_cap_markets.md.
#
# Self-contained — no RTS/test_system_dir fixture required. Follows the fixture-
# building convention in test_hdf5_serialization.jl (small hand-built structs) and
# the EMISAgentSimulation.<fn> convention in test_PSI.jl for calling unexported
# internal functions.
#
# Deliberately NOT covered here (see comments at each phase boundary for why, and
# Phase 5 in the implementation guide for where they belong instead):
#   - calculate_derating_data / calculate_derating_factors (Phase 2a): full TopNetLoad
#     / PRAS methodology needs real net-load timeseries / PRAS systems; numerical
#     correctness of the seasonal split belongs in an integration test.
#   - cem() (Phase 3b) and actual_market_simulation's capacity section (Phase 4b):
#     full multi-market JuMP models; the season-indexing *pattern* they rely on
#     (make_capacity_demand_vectors, the Dict-per-season structs) is covered here,
#     but solving the real model needs a full MarketClearingProblem/System fixture.
#   - calculate_realized_profit (Phase 4d): ~15-argument function entangled with the
#     whole realized-profit pipeline; the season-summation logic it depends on
#     (get_prices(...)["realized"] being season-keyed) is covered here via Phase 1e.
#
# Run from REPL after loading the package:
#   using EMISAgentSimulation
#   include("test/test_seasonal_capacity_markets.jl")

using Test
using EMISAgentSimulation
using AxisArrays
using DataFrames
using CSV
using HDF5
import JuMP
import GLPK
import HiGHS

# ─────────────────────────────────────────────────────────────────────────────
# Shared fixtures
# ─────────────────────────────────────────────────────────────────────────────

"""
Builds a minimal, disk-safe CaseDefinition via the keyword-default outer
constructor (src/structs/CaseDefinition.jl ~line 496). Pre-creates the case
data directory so the constructor's `make_case_data_dir` call sees it already
exists and skips its file-copy side effect.
"""
function make_case(; seasonal_capacity_market::Bool = false)
    base_dir = mktempdir()
    name = "TestCase"
    mkpath(joinpath(base_dir, name))
    solver = JuMP.optimizer_with_attributes(GLPK.Optimizer)
    return CaseDefinition(
        name, base_dir, base_dir, base_dir, base_dir, base_dir, solver;
        seasonal_capacity_market = seasonal_capacity_market,
    )
end

function make_capacity_product(;
    derating = Dict{String, Dict{String, Float64}}(),
    accepted_perc = Dict{String, Dict{String, Vector{Float64}}}(),
    capacity_bid = 100.0,
)
    return Capacity(:Capacity, derating, accepted_perc, capacity_bid)
end

function make_energy_product()
    return Energy(:Energy, Dict{String, Array{Float64, 2}}(), 3.5, 0.0)
end

function make_thermal_tech(type::String = "CC")
    return ThermalTech(
        type, "Gas",
        (min = 50.0, max = 400.0),
        (up = 100.0, down = 100.0),
        (up = 2.0, down = 4.0),
        nothing,
        4.5,
        [(0.0, 8.0), (400.0, 9.5)],
        "BusA", "ZoneA", 0.05, 24,
    )
end

function make_hydro_tech(type::String = "HY")
    return HydroTech(
        type,
        (min = 50.0, max = 400.0),
        nothing,
        nothing,
        nothing,
        "BusA",
        "ZoneA",
        0.05,
        24,
    )
end

function make_renewable_tech(type::String = "WT")
    return RenewableTech(
        type,
        (min = 50.0, max = 400.0),
        nothing,
        nothing,
        "BusA",
        "Z1",
        0.0,
        0,
    )
end

function make_battery_tech()
    return BatteryTech(
        "BA",
        (min = 0.0, max = 100.0),
        (min = 0.0, max = 100.0),
        nothing,
        (min = 0.0, max = 200.0),
        (min = 0.0, max = 1.0),
        0.5,
        100.0,
        0.5,
        (in = 0.9, out = 0.9),
        "BusA",
        "Z1",
        0.05,
        24,
        100.0,
    )
end

function make_finance_stub()
    return Finance(
        [0.0],                                     # investment_cost
        0.0,                                        # effective_investment_cost
        [1.0],                                       # preference_multiplier
        0, 30, 20,                                    # lag_time, life_time, capex_years
        0.0,                                         # fixed_OM_cost
        [0.0],                                       # queue_cost
        Dict{String, Array{Float64, 2}}(),            # scenario_total_utilization
        Dict{String, Vector{AxisArrays.AxisArray{Float64, 2}}}(),  # scenario_profit
        AxisArrays.AxisArray(zeros(1, 1), AxisArrays.Axis{:year}(1:1), AxisArrays.Axis{:product}([:Capacity])),  # realized_profit
        0.08,                                        # discount_rate
        Dict{String, Vector{Float64}}(),              # scenario_npv
        Float64[],                                    # expected_npv
        Dict{String, Vector{Float64}}(),              # scenario_utility
        Float64[],                                    # expected_utility
        Float64[],                                    # annual_cashflow
        "TestInvestor",                               # ownedby
    )
end

function make_thermal_project(products::Vector{<:Any})
    return ThermalGenEMIS{Existing}(
        "Gen1", make_thermal_tech(), 1, 1, 30, 30, Product[products...], make_finance_stub(),
    )
end

function make_hydro_project(products::Vector{<:Any})
    return HydroGenEMIS{Existing}(
        "Hyd1",
        make_hydro_tech(),
        1,
        1,
        30,
        30,
        Product[products...],
        make_finance_stub(),
    )
end

function make_renewable_existing_project(products::Vector{<:Any}; type::String = "WT")
    return RenewableGenEMIS{Existing}(
        "RenEx1",
        make_renewable_tech(type),
        1,
        1,
        30,
        30,
        Product[products...],
        make_finance_stub(),
    )
end

function make_renewable_option_project(products::Vector{<:Any}; type::String = "WT")
    return RenewableGenEMIS{Option}(
        "RenOpt1",
        make_renewable_tech(type),
        1,
        1,
        30,
        30,
        Product[products...],
        make_finance_stub(),
    )
end

function make_battery_existing_project(products::Vector{<:Any})
    return BatteryEMIS{Existing}(
        "BattEx1",
        make_battery_tech(),
        1,
        1,
        30,
        30,
        Product[products...],
        make_finance_stub(),
    )
end

function make_battery_option_project(products::Vector{<:Any})
    return BatteryEMIS{Option}(
        "BattOpt1",
        make_battery_tech(),
        1,
        1,
        30,
        30,
        Product[products...],
        make_finance_stub(),
    )
end

"""
Mirrors the ~40-field MarketProject inner constructor
(src/structs/market_structs/MarketProject.jl). Only `derating_factor` varies
between tests — that's the field under test in Phase 1c.
"""
function make_market_project(; derating_factor)
    return MarketProject(
        "TestGen", "Generator", "CC",
        1000.0,                          # fixed_cost
        [0.0, 0.0],                      # queue_cost
        20.0,                            # marginal_energy_cost
        Dict("reg_up" => 5.0),           # marginal_reserve_cost
        0.4,                             # emission_intensity
        [500_000.0],                     # expansion_cost
        0.08,                            # discount_rate
        50.0, 400.0,                     # min_gen, max_gen
        0.0, 0.0,                        # min_input, max_input
        1.0, 1.0,                        # efficiency_in, efficiency_out
        0.0, 0.0, 0.0,                   # min_storage, max_storage, init_storage
        ones(1, 24),                     # availability
        derating_factor,                 # derating_factor  <-- under test
        (up = 100.0, down = 100.0),      # ramp_limits
        Dict("reg_up" => 50.0),          # max_reserve_limits
        1.0,                             # existing_units
        [0.0],                           # units_in_queue
        3, 0,                            # build_lead_time, remaining_build_time
        1, 0,                            # max_new_options, base_cost_units
        20, 30, 25,                      # capex_years, life_time, remaining_life
        true, false,                     # capacity_eligible, rec_eligible
        0.0,                             # rec_correction
        3.0, true,                       # inertia_constant, synchronous_inertia
        "ZoneA", ["inv1"],               # zone, ownedby
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 1a — CaseDefinition seasonal toggle
# ─────────────────────────────────────────────────────────────────────────────

@testset "Phase 1a — CaseDefinition.seasonal_capacity_market" begin
    @testset "field defaults to false and is readable via getter" begin
        case = make_case()
        @test EMISAgentSimulation.get_seasonal_capacity_market(case) == false
    end

    @testset "toggling true is respected" begin
        case = make_case(; seasonal_capacity_market = true)
        @test EMISAgentSimulation.get_seasonal_capacity_market(case) == true
    end

    @testset "HDF5 round-trip preserves the field (both values)" begin
        for flag in (true, false)
            tmpfile = tempname() * ".h5"
            try
                case = make_case(; seasonal_capacity_market = flag)
                h5open(tmpfile, "w") do f
                    EMISAgentSimulation.save_case_definition!(HDF5.create_group(f, "case"), case)
                end
                loaded = h5open(tmpfile, "r") do f
                    EMISAgentSimulation.load_case_definition(f["case"])
                end
                @test EMISAgentSimulation.get_seasonal_capacity_market(loaded) == flag
            finally
                isfile(tmpfile) && rm(tmpfile)
            end
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 1b — Capacity product: season-nested derating / accepted_perc
# ─────────────────────────────────────────────────────────────────────────────

@testset "Phase 1b — Capacity product season dispatch" begin
    scenario = "scenario_1"

    @testset "set_derating!/get_derating with an explicit season" begin
        prod = make_capacity_product()
        EMISAgentSimulation.set_derating!(prod, scenario, "summer", 0.82)
        EMISAgentSimulation.set_derating!(prod, scenario, "winter", 0.55)

        @test get_derating(prod, scenario, "summer") == 0.82
        @test get_derating(prod, scenario, "winter") == 0.55
    end

    @testset "get_derating(prod, scenario) backward-compat wrapper reads the 'annual' key" begin
        prod = make_capacity_product()
        EMISAgentSimulation.set_derating!(prod, scenario, "annual", 0.9)
        @test get_derating(prod, scenario) == 0.9
    end

    @testset "set_accepted_perc!/get_accepted_perc with an explicit season" begin
        prod = make_capacity_product()
        EMISAgentSimulation.set_accepted_perc!(prod, scenario, "summer", [0.7, 0.8, 0.9])
        @test get_accepted_perc(prod, scenario, "summer") == [0.7, 0.8, 0.9]
    end

    @testset "catch-all methods accept the season arg for non-Capacity products" begin
        # Recurring bug class in this codebase: a catch-all/no-op method's arity or
        # argument types drift out of sync with the specialized method. Season-aware
        # call sites (Phase 2/3/4) loop over ALL products on a project and call these
        # functions unconditionally, so the non-Capacity path must not error.
        energy_prod = make_energy_product()
        @test_nowarn EMISAgentSimulation.set_derating!(energy_prod, scenario, "summer", 0.5)
        @test_nowarn EMISAgentSimulation.set_accepted_perc!(energy_prod, scenario, "summer", [0.1])

        # get_derating(prod::Product) = nothing already exists as a 1-arg catch-all
        # today; confirm an equivalent (season-arg-accepting) catch-all exists once
        # get_derating(::Capacity, scenario, season) is added — this is a second,
        # related catch-all gap not called out in the original plan text.
        @test isnothing(get_derating(energy_prod, scenario, "summer"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 1c — MarketProject.derating_factor as a season Dict
# ─────────────────────────────────────────────────────────────────────────────

@testset "Phase 1c — MarketProject.derating_factor season Dict" begin
    @testset "constructor accepts a season-keyed Dict" begin
        mp = make_market_project(derating_factor = Dict("summer" => 0.85, "winter" => 0.60))
        @test mp.derating_factor isa Dict{String, Float64}
        @test mp.derating_factor["summer"] == 0.85
        @test mp.derating_factor["winter"] == 0.60
    end

    @testset "negative derating in any season still throws AssertionError" begin
        @test_throws AssertionError make_market_project(derating_factor = Dict("summer" => -0.1))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 1d — MarketCollection.capacity as a season Dict
# ─────────────────────────────────────────────────────────────────────────────

@testset "Phase 1d — MarketCollection.capacity season Dict (type-level)" begin
    # Full construction needs EnergyMarket/RECMarket/InertiaMarket fixtures that
    # are out of scope for this feature's tests; the field-type check below is a
    # lightweight structural contract check. Add a full construction test once
    # Phase 1d actually lands, alongside whatever EnergyMarket-family fixtures
    # Phase 3/4's own tests end up building.
    @test fieldtype(MarketCollection, :capacity) == Dict{String, CapacityMarket}
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 1e — MarketPrices: shared expected/realized capacity_price storage
# ─────────────────────────────────────────────────────────────────────────────

@testset "Phase 1e — MarketPrices season-nested capacity_price" begin
    @test fieldtype(MarketPrices, :capacity_price) ==
          Union{Nothing, Dict{String, Dict{String, AxisArrays.AxisArray{Float64, 1}}}}

    @testset "set_capacity_price!/get_capacity_price round-trip a season Dict" begin
        mp = MarketPrices(nothing, nothing, nothing, nothing, nothing)

        price_summer = AxisArrays.AxisArray(rand(3), AxisArrays.Axis{:year}(1:3))
        price_winter = AxisArrays.AxisArray(rand(3), AxisArrays.Axis{:year}(1:3))
        seasonal_prices = Dict("summer" => price_summer, "winter" => price_winter)

        EMISAgentSimulation.set_capacity_price!(mp, "realized", seasonal_prices)

        stored = get_capacity_price(mp)["realized"]
        @test stored["summer"] == price_summer
        @test stored["winter"] == price_winter
    end

    @testset "expected-price call site shape (mirrors investor_iteration.jl:51)" begin
        # investor_iteration.jl feeds CEM/expected prices through the exact same
        # setter under the scenario name as key — not covered by the original
        # plan's file list; this test stands in for that call site.
        mp = MarketPrices(nothing, nothing, nothing, nothing, nothing)
        expected_prices = Dict(
            "summer" => AxisArrays.AxisArray(rand(3), AxisArrays.Axis{:year}(1:3)),
            "winter" => AxisArrays.AxisArray(rand(3), AxisArrays.Axis{:year}(1:3)),
        )
        EMISAgentSimulation.set_capacity_price!(mp, "scenario_1", expected_prices)
        @test Set(keys(get_capacity_price(mp)["scenario_1"])) == Set(["summer", "winter"])
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# capacity_seasons.csv / load_season_months (referenced by Phases 2, 3c, 4b)
# ─────────────────────────────────────────────────────────────────────────────

@testset "load_season_months utility" begin
    mktempdir() do dir
        markets_dir = joinpath(dir, "markets_data")
        mkpath(markets_dir)
        seasons_file = joinpath(markets_dir, "capacity_seasons.csv")

        @testset "wrap-around season parsing (summer/winter)" begin
            CSV.write(seasons_file, DataFrame("name" => ["summer", "winter"], "months" => ["May-Sep", "Oct-Apr"]))
            season_months = EMISAgentSimulation.load_season_months(seasons_file)

            @test Set(keys(season_months)) == Set(["summer", "winter"])
            @test season_months["summer"] == collect(5:9)
            @test season_months["winter"] == vcat(collect(1:4), collect(10:12))
        end

        @testset "seasons partition the year — no gaps, no overlaps" begin
            CSV.write(seasons_file, DataFrame("name" => ["summer", "winter"], "months" => ["May-Sep", "Oct-Apr"]))
            season_months = EMISAgentSimulation.load_season_months(seasons_file)
            all_months = sort(vcat(values(season_months)...))
            @test all_months == collect(1:12)   # full coverage, no duplicate months
            @test length(all_months) == length(unique(all_months))  # no overlap
        end

        @testset "4-season (with shoulders) partition" begin
            CSV.write(seasons_file, DataFrame(
                "name" => ["summer", "spring", "fall", "winter"],
                "months" => ["Jun-Aug", "Mar-May", "Sep-Nov", "Dec-Feb"],
            ))
            season_months = EMISAgentSimulation.load_season_months(seasons_file)
            all_months = sort(vcat(values(season_months)...))
            @test all_months == collect(1:12)
        end

        @testset "annual fallback when the file is absent" begin
            missing_file = joinpath(markets_dir, "does_not_exist.csv")
            @test !isfile(missing_file)
            season_months = EMISAgentSimulation.load_season_months(missing_file)
            @test season_months == Dict("annual" => collect(1:12))
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 2b — update_derating_factor! seasonal dispatch
# ─────────────────────────────────────────────────────────────────────────────

@testset "Phase 2b — update_derating_factor! iterates derating_dict.csv season rows" begin
    scenario = "scenario_1"
    fixture_dir = joinpath(@__DIR__, "test_data", "derating_update_cases")

        """
        Installs fixture CSVs into the per-test temporary simulation directory layout
        expected by update_derating_factor!:
            markets_data/derating_data/{scenario}/derating_dict.csv

        Optionally installs cc_scalar.csv when scalar_file is provided, so tests can
        exercise the scalar-multiplication path without mutating shared fixture files.
        """
        function install_derating_fixture!(
        dir::String;
        derating_file::String,
        scalar_file::Union{Nothing, String} = nothing,
    )
        derating_dir = joinpath(dir, "markets_data", "derating_data", scenario)
        mkpath(derating_dir)
        cp(
            joinpath(fixture_dir, derating_file),
            joinpath(derating_dir, "derating_dict.csv");
            force = true,
        )
        # Only copy cc_scalar.csv for tests that need scalar multiplication.
        # When scalar_file is nothing, the update path should use default scalar 1.0.
        if !isnothing(scalar_file)
            cp(
                joinpath(fixture_dir, scalar_file),
                joinpath(derating_dir, "cc_scalar.csv");
                force = true,
            )
        end
        return
    end

    @testset "seasonal mode — row iteration across update_derating_factor! families" begin
        cases = [
            (
                name = "thermal via type column",
                project = make_thermal_project,
                marginal_cc = false,
                expected = Dict("summer" => 0.81, "winter" => 0.91),
            ),
            (
                name = "hydro via type column",
                project = make_hydro_project,
                marginal_cc = false,
                expected = Dict("summer" => 0.62, "winter" => 0.66),
            ),
            (
                name = "renewable existing",
                project = make_renewable_existing_project,
                marginal_cc = false,
                expected = Dict("summer" => 0.41, "winter" => 0.27),
            ),
            (
                name = "renewable option marginal_cc=true",
                project = make_renewable_option_project,
                marginal_cc = true,
                expected = Dict("summer" => 0.33, "winter" => 0.19),
            ),
            (
                name = "renewable option marginal_cc=false",
                project = make_renewable_option_project,
                marginal_cc = false,
                expected = Dict("summer" => 0.41, "winter" => 0.27),
            ),
            (
                name = "battery existing",
                project = make_battery_existing_project,
                marginal_cc = false,
                expected = Dict("summer" => 0.73, "winter" => 0.65),
            ),
            (
                name = "battery option marginal_cc=true",
                project = make_battery_option_project,
                marginal_cc = true,
                expected = Dict("summer" => 0.69, "winter" => 0.58),
            ),
            (
                name = "battery option marginal_cc=false",
                project = make_battery_option_project,
                marginal_cc = false,
                expected = Dict("summer" => 0.73, "winter" => 0.65),
            ),
        ]

        for tc in cases
            @testset "$(tc.name)" begin
                mktempdir() do dir
                    install_derating_fixture!(dir; derating_file = "derating_dict_seasonal.csv")

                    prod = make_capacity_product()
                    project = tc.project([prod])
                    EMISAgentSimulation.update_derating_factor!(project, dir, scenario, tc.marginal_cc)

                    @test get_derating(prod, scenario, "summer") == tc.expected["summer"]
                    @test get_derating(prod, scenario, "winter") == tc.expected["winter"]
                end
            end
        end
    end

    @testset "seasonal mode — cc_scalar multiplication is applied per season" begin
        cases = [
            (
                name = "thermal scalar",
                project = make_thermal_project,
                marginal_cc = false,
                expected = Dict("summer" => 1.0125, "winter" => 1.1375),
            ),
            (
                name = "renewable existing scalar",
                project = make_renewable_existing_project,
                marginal_cc = false,
                expected = Dict("summer" => 0.615, "winter" => 0.405),
            ),
            (
                name = "battery option scalar",
                project = make_battery_option_project,
                marginal_cc = true,
                expected = Dict("summer" => 0.759, "winter" => 0.638),
            ),
        ]

        for tc in cases
            @testset "$(tc.name)" begin
                mktempdir() do dir
                    install_derating_fixture!(
                        dir;
                        derating_file = "derating_dict_seasonal.csv",
                        scalar_file = "cc_scalar.csv",
                    )

                    prod = make_capacity_product()
                    project = tc.project([prod])
                    EMISAgentSimulation.update_derating_factor!(project, dir, scenario, tc.marginal_cc)

                    @test get_derating(prod, scenario, "summer") ≈ tc.expected["summer"]
                    @test get_derating(prod, scenario, "winter") ≈ tc.expected["winter"]
                end
            end
        end
    end

    @testset "annual mode — single-row derating assigns only the annual key" begin
        cases = [
            (
                name = "thermal annual fallback",
                project = make_thermal_project,
                marginal_cc = false,
                expected = 0.88,
            ),
            (
                name = "renewable annual fallback",
                project = make_renewable_existing_project,
                marginal_cc = false,
                expected = 0.51,
            ),
            (
                name = "battery annual fallback",
                project = make_battery_existing_project,
                marginal_cc = false,
                expected = 0.67,
            ),
        ]

        for tc in cases
            @testset "$(tc.name)" begin
                mktempdir() do dir
                    install_derating_fixture!(dir; derating_file = "derating_dict_annual.csv")

                    prod = make_capacity_product()
                    project = tc.project([prod])
                    EMISAgentSimulation.update_derating_factor!(project, dir, scenario, tc.marginal_cc)

                    # In annual mode, only the annual season key should be assigned.
                    @test get_derating(prod, scenario, "annual") == tc.expected
                    @test_throws KeyError get_derating(prod, scenario, "summer")
                end
            end
        end
    end

    @testset "seasonal mode — missing required derating columns throw" begin
        @testset "renewable option with marginal_cc=true requires new_ column" begin
            mktempdir() do dir
                install_derating_fixture!(dir; derating_file = "derating_dict_missing_new_wind.csv")

                prod = make_capacity_product()
                project = make_renewable_option_project([prod])
                @test_throws Exception EMISAgentSimulation.update_derating_factor!(project, dir, scenario, true)
            end
        end

        @testset "battery option with marginal_cc=true requires new_STOR_N column" begin
            mktempdir() do dir
                install_derating_fixture!(dir; derating_file = "derating_dict_missing_new_stor2.csv")

                prod = make_capacity_product()
                project = make_battery_option_project([prod])
                @test_throws Exception EMISAgentSimulation.update_derating_factor!(project, dir, scenario, true)
            end
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 2d — get_project_derating
# ─────────────────────────────────────────────────────────────────────────────

@testset "Phase 2d — get_project_derating returns a season Dict" begin
    scenario = "scenario_1"

    @testset "project with a Capacity product returns its full season Dict" begin
        prod = make_capacity_product()
        EMISAgentSimulation.set_derating!(prod, scenario, "summer", 0.8)
        EMISAgentSimulation.set_derating!(prod, scenario, "winter", 0.5)
        project = make_thermal_project([prod])

        result = EMISAgentSimulation.get_project_derating(project, scenario)
        @test result isa Dict{String, Float64}
        @test result["summer"] == 0.8
        @test result["winter"] == 0.5
    end

    @testset "project with no Capacity product falls back to an empty/zero Dict, not a scalar" begin
        # Pre-implementation review flagged this: the current fallback is a bare
        # `0.` scalar (see market_project_creator.jl:135-144), which won't type-match
        # a Dict{String,Float64} return. Implementation must pick a convention
        # (empty Dict + "missing key means 0" downstream, most likely) — this test
        # encodes that expectation and should be adjusted if a different convention
        # is chosen instead.
        project = make_thermal_project([make_energy_product()])

        result = EMISAgentSimulation.get_project_derating(project, scenario)
        @test result isa Dict{String, Float64}
        @test isempty(result) || all(values(result) .== 0.0)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 3a — market_utils.jl demand curve vectors
# ─────────────────────────────────────────────────────────────────────────────

@testset "Phase 3a — capacity demand curve vector conversion" begin
    summer_market = CapacityMarket([0.0, 100.0, 200.0], [50.0, 30.0, 0.0])
    winter_market = CapacityMarket([0.0, 150.0, 250.0], [60.0, 25.0, 0.0])

    @testset "make_capacity_demand_vectors accepts Vector{Dict{String,CapacityMarket}}" begin
        capmarkets = [Dict("summer" => summer_market, "winter" => winter_market)]  # 1 investment period

        segsize, seggrad, pricepoints, numsegs = EMISAgentSimulation.make_capacity_demand_vectors(capmarkets)

        @test Set(keys(segsize)) == Set(["summer", "winter"])
        @test numsegs["summer"][1] == 2
        @test segsize["summer"][1] == [100.0, 100.0]
        @test pricepoints["winter"][1] == [60.0, 25.0, 0.0]
        @test seggrad["summer"][1] == [(30.0 - 50.0) / 100.0, (0.0 - 30.0) / 100.0]
    end

    @testset "make_capacity_demand (single-market, actual clearing) is unchanged — regression" begin
        # Plan explicitly states this function is NOT touched; this is a
        # won't-break check, expected to pass both before and after the feature lands.
        segsize, seggrad, pricepoints = EMISAgentSimulation.make_capacity_demand(summer_market)
        @test segsize == [100.0, 100.0]
        @test pricepoints == [50.0, 30.0, 0.0]
        @test seggrad == [(30.0 - 50.0) / 100.0, (0.0 - 30.0) / 100.0]
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 3d — capacity_profit.jl: update_capacity_revenues! sums over seasons
# ─────────────────────────────────────────────────────────────────────────────

@testset "Phase 3d — update_capacity_revenues! sums over seasons" begin
    scenario = "scenario_1"
    price_years = (start_year = 1, end_year = 3)
    max_cap = 100.0

    prod = make_capacity_product()
    EMISAgentSimulation.set_derating!(prod, scenario, "summer", 0.8)
    EMISAgentSimulation.set_derating!(prod, scenario, "winter", 0.6)
    EMISAgentSimulation.set_accepted_perc!(prod, scenario, "summer", [1.0, 1.0, 1.0])
    EMISAgentSimulation.set_accepted_perc!(prod, scenario, "winter", [1.0, 1.0, 1.0])

    capacity_prices = Dict(scenario => Dict(
        "summer" => AxisArrays.AxisArray(fill(10.0, 3), AxisArrays.Axis{:year}(1:3)),
        "winter" => AxisArrays.AxisArray(fill(20.0, 3), AxisArrays.Axis{:year}(1:3)),
    ))

    @testset "revenue is the sum of each season's price x accepted_perc x derated capacity" begin
        revenues = EMISAgentSimulation.update_capacity_revenues!(
            prod, scenario, price_years, zeros(3), capacity_prices, max_cap,
        )
        # summer: 10 * 1.0 * 100 * 0.8 = 800/yr ; winter: 20 * 1.0 * 100 * 0.6 = 1200/yr
        @test revenues ≈ fill(2000.0, 3)
    end

    @testset "catch-all for non-Capacity products matches the new nested capacity_prices type" begin
        # Same catch-all/specialized-method arity/type drift risk flagged for
        # set_derating! and update_capacity_supply_curve! — capacity_prices here is
        # explicitly typed on the catch-all today (Dict{String,AxisArray}); it must
        # be updated to Dict{String,Dict{String,AxisArray}} in the same change.
        energy_prod = make_energy_product()
        @test_nowarn EMISAgentSimulation.update_capacity_revenues!(
            energy_prod, scenario, price_years, zeros(3), capacity_prices, max_cap,
        )
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 4a — actual_capacity_mkt_clearing.jl
# ─────────────────────────────────────────────────────────────────────────────

@testset "Phase 4a — supply curve + demand curve construction" begin
    scenario = "scenario_1"

    @testset "update_capacity_supply_curve! takes a season arg and uses seasonal derating" begin
        prod = make_capacity_product()
        EMISAgentSimulation.set_derating!(prod, scenario, "summer", 0.8)
        project = make_thermal_project([prod])

        curve = Vector{Union{String, Float64}}[]
        curve = EMISAgentSimulation.update_capacity_supply_curve!(curve, prod, project, scenario, "summer")

        @test length(curve) == 1
        @test curve[1][1] == get_name(project)
        @test curve[1][4] == 0.8  # seasonal derating value, not a bare scenario scalar
    end

    @testset "catch-all matches the new 5-arg signature for non-Capacity products" begin
        # This is the sharpest instance of the catch-all/arity-mismatch bug class
        # in this feature: Phase 4b's per-season loop calls this function for EVERY
        # product on every project, not just Capacity products.
        project = make_thermal_project([make_energy_product()])
        curve = Vector{Union{String, Float64}}[]
        @test_nowarn EMISAgentSimulation.update_capacity_supply_curve!(
            curve, make_energy_product(), project, scenario, "summer",
        )
        @test isempty(curve)
    end

    @testset "create_capacity_demand_curve accepts a pre-loaded params row, not a file path" begin
        # Both Phase 3c's and Phase 4b's example code call this with a DataFrameRow
        # already extracted from a multi-row Capacity.csv — the current
        # implementation takes a file path and does read_data(path)[1,:] internally.
        params_df = DataFrame(
            "season" => ["summer"],
            "IRM" => [0.15],
            "EFORd" => [0.08],
            "IRM perc points" => ["0.02;0.04"],
            "Net CONE per day" => [40.0],
            "Net CONE perc points" => ["1.0;0.75"],
            "Gross CONE per day" => [55.0],
            "Max Clear" => [1.2],
        )
        row = first(params_df)

        curve = EMISAgentSimulation.create_capacity_demand_curve(row, 1000.0, 1.0, 0.0, true)
        @test curve isa CapacityMarket
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 4b — actual_market_simulation.jl seasonal capacity clearing
# ─────────────────────────────────────────────────────────────────────────────
#
# The full create_realized_marketdata entry point needs an AgentSimulation + PSY
# system fixture (deferred to an integration test, per this file's header). These
# testsets instead exercise the *units we changed* in the capacity section, with
# simple in-memory mocking:
#   - create_capacity_demand_curve(::DataFrameRow) — seasonal per-row variant
#   - the caller's season-row indexing pattern (cap_params_by_season + annual fallback)
#   - is_annual_only_seasons gating (the toggle that collapses to "annual")
#   - the per-season supply-build -> clear pipeline that produces one price per season

@testset "Phase 4b — seasonal capacity clearing units" begin
    scenario = "scenario_1"
    # capacity_market_clearing builds a quadratic (welfare) objective, so it needs a
    # QP-capable solver. GLPK (used elsewhere in this file) is LP-only; use HiGHS here.
    solver = JuMP.optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)

    # Two-season market params, distinct demand-curve shapes per season.
    cap_params_df = DataFrame(
        "season"               => ["summer", "winter"],
        "introduction_year"    => [2020, 2020],
        "discontinuation_year" => [2050, 2050],
        "IRM"                  => [0.15, 0.17],
        "EFORd"                => [0.08, 0.10],
        "IRM perc points"      => ["0.02;0.04", "0.02;0.04"],
        "Net CONE per day"     => [40.0, 45.0],
        "Net CONE perc points" => ["1.0;0.75", "1.0;0.75"],
        "Gross CONE per day"   => [55.0, 60.0],
        "Max Clear"            => [1.2, 1.2],
    )

    @testset "create_capacity_demand_curve(::DataFrameRow) yields season-specific curves" begin
        summer_curve = EMISAgentSimulation.create_capacity_demand_curve(
            cap_params_df[1, :], 1000.0, 1.0, 0.0, true)
        winter_curve = EMISAgentSimulation.create_capacity_demand_curve(
            cap_params_df[2, :], 1000.0, 1.0, 0.0, true)

        @test summer_curve isa CapacityMarket
        @test winter_curve isa CapacityMarket
        # Winter has higher IRM + Net CONE -> a distinct demand curve from summer.
        @test summer_curve.price_points != winter_curve.price_points
    end

    @testset "index_capacity_params_by_season keys rows by the season column" begin
        cap_params_by_season =
            EMISAgentSimulation.index_capacity_params_by_season(cap_params_df, ["summer", "winter"])
        @test Set(keys(cap_params_by_season)) == Set(["summer", "winter"])
        @test cap_params_by_season["winter"]["Net CONE per day"] == 45.0
        @test cap_params_by_season["summer"]["Net CONE per day"] == 40.0
    end

    @testset "index_capacity_params_by_season falls back to row 1 when no season column" begin
        annual_df = DataFrame(
            "introduction_year"    => [2020],
            "discontinuation_year" => [2050],
            "IRM"                  => [0.15],
            "EFORd"                => [0.08],
            "IRM perc points"      => ["0.02;0.04"],
            "Net CONE per day"     => [40.0],
            "Net CONE perc points" => ["1.0;0.75"],
            "Gross CONE per day"   => [55.0],
            "Max Clear"            => [1.2],
        )
        cap_params_by_season =
            EMISAgentSimulation.index_capacity_params_by_season(annual_df, ["annual"])
        @test collect(keys(cap_params_by_season)) == ["annual"]
        @test cap_params_by_season["annual"]["Net CONE per day"] == 40.0
    end

    @testset "resolve_capacity_seasons collapses to annual per the toggle / definition" begin
        seasonal_months = Dict("summer" => collect(5:9),
                               "winter" => vcat(collect(1:4), collect(10:12)))
        annual_months = Dict("annual" => collect(1:12))

        # Toggle OFF -> always annual, regardless of file content.
        @test EMISAgentSimulation.resolve_capacity_seasons(seasonal_months, false) == ["annual"]
        # Toggle ON with a real 2-season definition -> both seasons.
        @test Set(EMISAgentSimulation.resolve_capacity_seasons(seasonal_months, true)) ==
              Set(["summer", "winter"])
        # Toggle ON but only an annual row -> still annual.
        @test EMISAgentSimulation.resolve_capacity_seasons(annual_months, true) == ["annual"]
    end

    @testset "per-season pipeline: supply build (Phase 4a) -> clear yields one price per season" begin
        # Same project, different seasonal derating -> the season loop must clear
        # each season independently and produce a season-keyed price dict. Uses the
        # real caller helper to index params, so this exercises the production path.
        prod = make_capacity_product()
        EMISAgentSimulation.set_derating!(prod, scenario, "summer", 0.9)
        EMISAgentSimulation.set_derating!(prod, scenario, "winter", 0.5)
        project = make_thermal_project([prod])

        seasons = EMISAgentSimulation.resolve_capacity_seasons(
            Dict("summer" => collect(5:9), "winter" => vcat(collect(1:4), collect(10:12))), true)
        cap_params_by_season =
            EMISAgentSimulation.index_capacity_params_by_season(cap_params_df, seasons)

        capacity_price_dict = Dict{String, AxisArrays.AxisArray{Float64, 1}}()
        capacity_accepted_bids_dict = Dict{String, Dict{String, Float64}}()

        for season in seasons
            supply_curve = Vector{Union{String, Float64}}[]
            for p in get_products(project)
                supply_curve = EMISAgentSimulation.update_capacity_supply_curve!(
                    supply_curve, p, project, scenario, season)
            end
            sort!(supply_curve; by = x -> x[3])

            demand_curve = EMISAgentSimulation.create_capacity_demand_curve(
                cap_params_by_season[season], 1000.0, 1.0, 0.0, true)

            capacity_price_dict[season], capacity_accepted_bids_dict[season] =
                EMISAgentSimulation.capacity_market_clearing(demand_curve, supply_curve, solver)
        end

        # One price per season, and the Phase 4a seasonal derating flowed into supply [4].
        @test Set(keys(capacity_price_dict)) == Set(["summer", "winter"])
        @test length(capacity_price_dict["summer"]) == 1
        @test haskey(capacity_accepted_bids_dict["summer"], get_name(project))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 4d — realized_profits_calculator.jl: calculate_realized_profit sums over seasons
# ─────────────────────────────────────────────────────────────────────────────

@testset "Phase 4d — calculate_realized_profit(::Capacity) sums over seasons" begin
    scenario = "scenario_1"

    # Empty per-product/reserve inputs the Capacity method ignores but the signature requires.
    empty_cf   = Dict{String, Array{Float64, 2}}()
    empty_rp   = Dict{String, Dict{String, Array{Float64, 2}}}()
    empty_ip   = Dict{String, Array{Float64, 2}}()
    rec_bids   = Dict{String, Float64}()
    hour_w     = Dict{String, Dict{Int64, Vector{Float64}}}()
    rt_prods   = String[]

    function build_market_prices(price_by_season::Dict{String, Float64})
        mp = MarketPrices(nothing, nothing, nothing, nothing, nothing)
        season_prices = Dict(
            s => AxisArrays.AxisArray(reshape([v], 1), [1]) for (s, v) in price_by_season)
        EMISAgentSimulation.set_capacity_price!(mp, "realized", season_prices)
        return mp
    end

    # calculate_realized_profit(::Capacity) positional call with the Phase 4d seasonal
    # capacity_accepted_bids type (Dict{season => Dict{name => fraction}}).
    function realized_capacity_profit(project, product, mp, bids;
                                      iteration_year = 1, capacity_forward_years = 1)
        return EMISAgentSimulation.calculate_realized_profit(
            project, product, mp,
            empty_cf, empty_cf, empty_cf,
            empty_rp, empty_rp, empty_rp,
            empty_ip,
            bids, rec_bids, hour_w,
            iteration_year, capacity_forward_years,
            0.0, 60, 60, rt_prods, scenario)
    end

    @testset "profit is the sum of each season's derating x price x accepted fraction" begin
        prod = make_capacity_product()
        EMISAgentSimulation.set_derating!(prod, scenario, "summer", 0.8)
        EMISAgentSimulation.set_derating!(prod, scenario, "winter", 0.6)
        project = make_thermal_project([prod])
        name = get_name(project)
        size = get_maxcap(project)

        mp = build_market_prices(Dict("summer" => 10.0, "winter" => 20.0))
        bids = Dict("summer" => Dict(name => 1.0), "winter" => Dict(name => 1.0))

        profit, update_year = realized_capacity_profit(project, prod, mp, bids)
        # summer: size*0.8*10*1 ; winter: size*0.6*20*1
        @test profit ≈ size * (0.8 * 10.0 + 0.6 * 20.0)
        @test update_year == 1   # iteration_year + capacity_forward_years - 1
    end

    @testset "only seasons where the project cleared contribute" begin
        prod = make_capacity_product()
        EMISAgentSimulation.set_derating!(prod, scenario, "summer", 0.8)
        EMISAgentSimulation.set_derating!(prod, scenario, "winter", 0.6)
        project = make_thermal_project([prod])
        name = get_name(project)
        size = get_maxcap(project)

        mp = build_market_prices(Dict("summer" => 10.0, "winter" => 20.0))
        # Cleared in winter only (absent from summer bids).
        bids = Dict("summer" => Dict{String, Float64}(), "winter" => Dict(name => 1.0))

        profit, _ = realized_capacity_profit(project, prod, mp, bids)
        @test profit ≈ size * (0.6 * 20.0)
    end

    @testset "returns nothing when the project cleared in no season" begin
        prod = make_capacity_product()
        EMISAgentSimulation.set_derating!(prod, scenario, "summer", 0.8)
        EMISAgentSimulation.set_derating!(prod, scenario, "winter", 0.6)
        project = make_thermal_project([prod])

        mp = build_market_prices(Dict("summer" => 10.0, "winter" => 20.0))
        bids = Dict("summer" => Dict{String, Float64}(), "winter" => Dict{String, Float64}())

        profit, update_year = realized_capacity_profit(project, prod, mp, bids)
        @test isnothing(profit)
        @test update_year == 1
    end

    @testset "annual-mode single season reproduces the pre-seasonal scalar profit" begin
        prod = make_capacity_product()
        EMISAgentSimulation.set_derating!(prod, scenario, "annual", 0.75)
        project = make_thermal_project([prod])
        name = get_name(project)
        size = get_maxcap(project)

        mp = build_market_prices(Dict("annual" => 40.0))
        bids = Dict("annual" => Dict(name => 0.5))

        profit, _ = realized_capacity_profit(project, prod, mp, bids)
        @test profit ≈ size * 0.75 * 40.0 * 0.5
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 3e / 4e — HDF5 season-nesting mechanism for capacity_price
# ─────────────────────────────────────────────────────────────────────────────

@testset "capacity_price HDF5 season-nesting mechanism" begin
    # save_axisarray!/load_axisarray already exist and are generic over any single
    # AxisArray — they are NOT new code. What Phase 3e and Phase 4e add is a loop
    # that nests one extra group level (one sub-group per season) before calling
    # them. This test validates that nesting mechanism directly, since reproducing
    # the ~50-argument signatures of save_expected_market_data /
    # save_realized_market_data (and the third, previously-uncovered
    # save_realized_market_data/load_realized_market_data pair — see the
    # Pre-Implementation Review Notes in Capacity_markets_implementation.md) here
    # would test fixture-building, not the season logic. Once those functions are
    # updated, extend test_hdf5_serialization.jl's existing roundtrip-style tests
    # to call them directly for full coverage — this test only proves the
    # building block they'll be composed from already works.
    @testset "nested per-season groups round-trip an AxisArray Dict" begin
        tmpfile = tempname() * ".h5"
        try
            price_summer = AxisArrays.AxisArray(rand(3), AxisArrays.Axis{:year}(1:3))
            price_winter = AxisArrays.AxisArray(rand(3), AxisArrays.Axis{:year}(1:3))
            capacity_price = Dict("summer" => price_summer, "winter" => price_winter)

            h5open(tmpfile, "w") do f
                cap_g = HDF5.create_group(f, "capacity_price")
                for (season, prices) in capacity_price
                    EMISAgentSimulation.save_axisarray!(HDF5.create_group(cap_g, season), prices)
                end
            end

            loaded = h5open(tmpfile, "r") do f
                cap_g = f["capacity_price"]
                Dict(season => EMISAgentSimulation.load_axisarray(cap_g[season]) for season in keys(cap_g))
            end

            @test collect(loaded["summer"]) == collect(price_summer)
            @test collect(loaded["winter"]) == collect(price_winter)
        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 4e — realized-data seasonal HDF5 helpers (save/load round-trip)
# ─────────────────────────────────────────────────────────────────────────────
#
# Exercises the actual Phase 4e helpers wired into save_realized_market_data /
# load_realized_market_data, rather than re-implementing the nesting in the test.

@testset "Phase 4e — seasonal realized capacity persistence helpers" begin
    @testset "save/load_capacity_price_seasonal round-trips a season->AxisArray Dict" begin
        tmpfile = tempname() * ".h5"
        try
            price_summer = AxisArrays.AxisArray(reshape([12.0], 1), [1])
            price_winter = AxisArrays.AxisArray(reshape([25.0], 1), [1])
            capacity_price = Dict("summer" => price_summer, "winter" => price_winter)

            h5open(tmpfile, "w") do f
                EMISAgentSimulation.save_capacity_price_seasonal!(
                    HDF5.create_group(f, "capacity_price"), capacity_price)
            end
            loaded = h5open(tmpfile, "r") do f
                EMISAgentSimulation.load_capacity_price_seasonal(f["capacity_price"])
            end

            @test Set(keys(loaded)) == Set(["summer", "winter"])
            @test loaded["summer"][1] == 12.0
            @test loaded["winter"][1] == 25.0
        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end

    @testset "save/load_capacity_bids_seasonal round-trips a season->name->fraction Dict" begin
        tmpfile = tempname() * ".h5"
        try
            bids = Dict(
                "summer" => Dict("Gen1" => 1.0, "Gen2" => 0.5),
                "winter" => Dict("Gen1" => 0.25),
            )

            h5open(tmpfile, "w") do f
                EMISAgentSimulation.save_capacity_bids_seasonal!(
                    HDF5.create_group(f, "capacity_accepted_bids"), bids)
            end
            loaded = h5open(tmpfile, "r") do f
                EMISAgentSimulation.load_capacity_bids_seasonal(f["capacity_accepted_bids"])
            end

            @test Set(keys(loaded)) == Set(["summer", "winter"])
            @test loaded["summer"]["Gen1"] == 1.0
            @test loaded["summer"]["Gen2"] == 0.5
            @test loaded["winter"]["Gen1"] == 0.25
            @test !haskey(loaded["winter"], "Gen2")
        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end

    @testset "annual mode round-trips a single-season Dict" begin
        tmpfile = tempname() * ".h5"
        try
            capacity_price = Dict("annual" => AxisArrays.AxisArray(reshape([30.0], 1), [1]))
            bids = Dict("annual" => Dict("Gen1" => 0.8))

            h5open(tmpfile, "w") do f
                EMISAgentSimulation.save_capacity_price_seasonal!(
                    HDF5.create_group(f, "capacity_price"), capacity_price)
                EMISAgentSimulation.save_capacity_bids_seasonal!(
                    HDF5.create_group(f, "capacity_accepted_bids"), bids)
            end
            price_loaded, bids_loaded = h5open(tmpfile, "r") do f
                EMISAgentSimulation.load_capacity_price_seasonal(f["capacity_price"]),
                EMISAgentSimulation.load_capacity_bids_seasonal(f["capacity_accepted_bids"])
            end

            @test collect(keys(price_loaded)) == ["annual"]
            @test price_loaded["annual"][1] == 30.0
            @test bids_loaded["annual"]["Gen1"] == 0.8
        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end
end
