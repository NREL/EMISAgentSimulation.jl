# Tests for update_derating_factor! overloads.
# src: src/derating_factor_updates/derating_factor_calculator.jl
# Run from REPL after loading the package:
#   using EMISAgentSimulation
#   include("test/test_update_derating_factor.jl")

using Test
using DataFrames
using CSV
using AxisArrays
using EMISAgentSimulation
include(joinpath(@__DIR__, "helpers.jl"))

@testset "update_derating_factor!" begin
    @testset "duplicate component names" begin
        @test EMISAgentSimulation.resolve_unique_component_name(
            ["gen-543", "gen-543_1"],
            "gen-543",
        ) == "gen-543_2"
        @test EMISAgentSimulation.resolve_unique_component_name(["gen-543"], "new_name") ==
              "new_name"
    end

    # ── ThermalGenEMIS / HydroGenEMIS ────────────────────────────────────────

    @testset "ThermalGenEMIS / HydroGenEMIS" begin
        @testset "T4: Thermal — scalar applied, result not capped at 1.0" begin
            mktempdir() do tmp
                p = ThermalGenEMIS{Existing}("th1", _thermal_tech(), 2020, 2021, 2040, 2045,
                    Product[_cap()], _dummy_finance())
                _setup_derating_dict(tmp, "scenario_1")
                _write_cc_scalar(tmp, "scenario_1", "CT" => 2.0)
                update_derating_factor!(p, tmp, "scenario_1", false)
                @test _get_derating(p, "scenario_1") ≈ 0.85 * 2.0  # = 1.7
            end
        end

        @testset "T5: Thermal — no cc_scalar.csv → scalar defaults to 1.0" begin
            mktempdir() do tmp
                p = ThermalGenEMIS{Existing}("th1", _thermal_tech(), 2020, 2021, 2040, 2045,
                    Product[_cap()], _dummy_finance())
                _setup_derating_dict(tmp, "scenario_1")
                update_derating_factor!(p, tmp, "scenario_1", false)
                @test _get_derating(p, "scenario_1") ≈ 0.85
            end
        end

        @testset "T4b: HydroGenEMIS — scalar applied" begin
            mktempdir() do tmp
                p = HydroGenEMIS{Existing}("hy1", _hydro_tech(), 2020, 2021, 2040, 2045,
                    Product[_cap()], _dummy_finance())
                _setup_derating_dict(tmp, "scenario_1")
                _write_cc_scalar(tmp, "scenario_1", "HY" => 1.5)
                update_derating_factor!(p, tmp, "scenario_1", false)
                @test _get_derating(p, "scenario_1") ≈ 0.6 * 1.5
            end
        end
    end

    # ── RenewableGenEMIS ──────────────────────────────────────────────────────

    @testset "RenewableGenEMIS" begin
        @testset "T6: Existing — scalar applied exactly once" begin
            mktempdir() do tmp
                p = RenewableGenEMIS{Existing}("wind1", _renewable_tech("Wind"), 2020, 2021,
                    2040, 2045, Product[_cap()], _dummy_finance())
                _setup_derating_dict(tmp, "scenario_1")
                _write_cc_scalar(tmp, "scenario_1", "Wind" => 1.2)
                update_derating_factor!(p, tmp, "scenario_1", false)
                @test _get_derating(p, "scenario_1") ≈ 0.5 * 1.2
            end
        end

        @testset "T7: Existing — scalar = 1.0 → value unchanged" begin
            mktempdir() do tmp
                p = RenewableGenEMIS{Existing}("wind1", _renewable_tech("Wind"), 2020, 2021,
                    2040, 2045, Product[_cap()], _dummy_finance())
                _setup_derating_dict(tmp, "scenario_1")
                _write_cc_scalar(tmp, "scenario_1", "Wind" => 1.0)
                update_derating_factor!(p, tmp, "scenario_1", false)
                @test _get_derating(p, "scenario_1") ≈ 0.5
            end
        end

        @testset "T8: Existing — no cc_scalar.csv → value unchanged" begin
            mktempdir() do tmp
                p = RenewableGenEMIS{Existing}("wind1", _renewable_tech("Wind"), 2020, 2021,
                    2040, 2045, Product[_cap()], _dummy_finance())
                _setup_derating_dict(tmp, "scenario_1")
                update_derating_factor!(p, tmp, "scenario_1", false)
                @test _get_derating(p, "scenario_1") ≈ 0.5
            end
        end

        @testset "T9: Existing TopNetLoad — raw value × scalar" begin
            mktempdir() do tmp
                p = RenewableGenEMIS{Existing}("wind1", _renewable_tech("Wind"), 2020, 2021,
                    2040, 2045, Product[_cap()], _dummy_finance())
                _write_derating_dict(tmp, "scenario_1", "existing_Wind_Z1" => 0.4)
                _write_cc_scalar(tmp, "scenario_1", "Wind" => 1.2)
                update_derating_factor!(p, tmp, "scenario_1", false)
                @test _get_derating(p, "scenario_1") ≈ 0.4 * 1.2
            end
        end

        @testset "T10: Existing — scalar modifies, does not replace raw value" begin
            mktempdir() do tmp
                p = RenewableGenEMIS{Existing}("wind1", _renewable_tech("Wind"), 2020, 2021,
                    2040, 2045, Product[_cap()], _dummy_finance())
                raw, scalar = 0.4, 0.5
                _write_derating_dict(tmp, "scenario_1", "existing_Wind_Z1" => raw)
                _write_cc_scalar(tmp, "scenario_1", "Wind" => scalar)
                update_derating_factor!(p, tmp, "scenario_1", false)
                @test _get_derating(p, "scenario_1") ≈ raw * scalar
                @test _get_derating(p, "scenario_1") ≉ scalar
            end
        end

        @testset "T11: Option marginal_cc=true — uses new_ column × scalar" begin
            mktempdir() do tmp
                p = RenewableGenEMIS{Option}("wind_opt", _renewable_tech("Wind"), 2020,
                    2021,
                    2040, 2045, Product[_cap()], _dummy_finance())
                _setup_derating_dict(tmp, "scenario_1")
                _write_cc_scalar(tmp, "scenario_1", "Wind" => 1.5)
                update_derating_factor!(p, tmp, "scenario_1", true)
                @test _get_derating(p, "scenario_1") ≈ 0.35 * 1.5
            end
        end

        @testset "T12: Option TopNetLoad marginal_cc=true — scalar applied" begin
            mktempdir() do tmp
                p = RenewableGenEMIS{Option}("solar_opt", _renewable_tech("Solar"), 2020,
                    2021,
                    2040, 2045, Product[_cap()], _dummy_finance())
                _setup_derating_dict(tmp, "scenario_1")  # fixture: new_Solar_Z1=0.22
                _write_cc_scalar(tmp, "scenario_1", "Solar" => 0.9)
                update_derating_factor!(p, tmp, "scenario_1", true)
                @test _get_derating(p, "scenario_1") ≈ 0.22 * 0.9
            end
        end

        @testset "T13: Option marginal_cc=false — uses existing_ column × scalar" begin
            mktempdir() do tmp
                p = RenewableGenEMIS{Option}("wind_opt", _renewable_tech("Wind"), 2020,
                    2021,
                    2040, 2045, Product[_cap()], _dummy_finance())
                _setup_derating_dict(tmp, "scenario_1")
                _write_cc_scalar(tmp, "scenario_1", "Wind" => 1.5)
                update_derating_factor!(p, tmp, "scenario_1", false)
                @test _get_derating(p, "scenario_1") ≈ 0.5 * 1.5
            end
        end
    end

    # ── BatteryEMIS ───────────────────────────────────────────────────────────

    @testset "BatteryEMIS" begin
        @testset "T14: Existing — STOR_2 scalar applied once" begin
            mktempdir() do tmp
                p = BatteryEMIS{Existing}("batt1", _battery_tech(), 2020, 2021, 2040, 2045,
                    Product[_cap()], _dummy_finance())
                _setup_derating_dict(tmp, "scenario_1")
                _write_cc_scalar(tmp, "scenario_1", "STOR_2" => 1.2)
                update_derating_factor!(p, tmp, "scenario_1", false)
                @test _get_derating(p, "scenario_1") ≈ 0.7 * 1.2
            end
        end

        @testset "T15: Existing — anti-regression: scalar applied exactly once (not twice)" begin
            mktempdir() do tmp
                p = BatteryEMIS{Existing}("batt1", _battery_tech(), 2020, 2021, 2040, 2045,
                    Product[_cap()], _dummy_finance())
                raw, scalar = 0.7, 0.8
                _setup_derating_dict(tmp, "scenario_1")  # fixture: existing_STOR_2=0.7
                _write_cc_scalar(tmp, "scenario_1", "STOR_2" => scalar)
                update_derating_factor!(p, tmp, "scenario_1", false)
                @test _get_derating(p, "scenario_1") ≈ raw * scalar           # 0.56
                @test _get_derating(p, "scenario_1") ≉ raw * scalar * scalar  # NOT 0.448
            end
        end

        @testset "T16: Existing TopNetLoad — same update path, scalar applied once" begin
            mktempdir() do tmp
                p = BatteryEMIS{Existing}("batt1", _battery_tech(), 2020, 2021, 2040, 2045,
                    Product[_cap()], _dummy_finance())
                _write_derating_dict(tmp, "scenario_1", "existing_STOR_2" => 0.65)
                _write_cc_scalar(tmp, "scenario_1", "STOR_2" => 1.1)
                update_derating_factor!(p, tmp, "scenario_1", false)
                @test _get_derating(p, "scenario_1") ≈ 0.65 * 1.1
            end
        end

        @testset "T17: Option marginal_cc=true — uses new_STOR_N × scalar" begin
            mktempdir() do tmp
                p = BatteryEMIS{Option}("batt_opt", _battery_tech(), 2020, 2021, 2040, 2045,
                    Product[_cap()], _dummy_finance())
                _setup_derating_dict(tmp, "scenario_1")
                _write_cc_scalar(tmp, "scenario_1", "STOR_2" => 1.1)
                update_derating_factor!(p, tmp, "scenario_1", true)
                @test _get_derating(p, "scenario_1") ≈ 0.65 * 1.1
            end
        end

        @testset "T18: Option marginal_cc=false — uses existing_STOR_N × scalar" begin
            mktempdir() do tmp
                p = BatteryEMIS{Option}("batt_opt", _battery_tech(), 2020, 2021, 2040, 2045,
                    Product[_cap()], _dummy_finance())
                _setup_derating_dict(tmp, "scenario_1")
                _write_cc_scalar(tmp, "scenario_1", "STOR_2" => 1.1)
                update_derating_factor!(p, tmp, "scenario_1", false)
                @test _get_derating(p, "scenario_1") ≈ 0.7 * 1.1
            end
        end
    end

    # ── No cap at 1.0 ─────────────────────────────────────────────────────────

    @testset "No cap at 1.0" begin
        @testset "T19: Thermal scalar > 1.0 → result > 1.0, not capped" begin
            mktempdir() do tmp
                p = ThermalGenEMIS{Existing}("th1", _thermal_tech(), 2020, 2021, 2040, 2045,
                    Product[_cap()], _dummy_finance())
                _write_derating_dict(tmp, "scenario_1", "CT" => 0.9)
                _write_cc_scalar(tmp, "scenario_1", "CT" => 2.0)
                update_derating_factor!(p, tmp, "scenario_1", false)
                result = _get_derating(p, "scenario_1")
                @test result ≈ 1.8
                @test result > 1.0
            end
        end

        @testset "T20: Renewable ELCC × 1.5 = 0.525" begin
            mktempdir() do tmp
                p = RenewableGenEMIS{Existing}("solar1", _renewable_tech("Solar"), 2020,
                    2021,
                    2040, 2045, Product[_cap()], _dummy_finance())
                _setup_derating_dict(tmp, "scenario_1")  # fixture: existing_Solar_Z1=0.35
                _write_cc_scalar(tmp, "scenario_1", "Solar" => 1.5)
                update_derating_factor!(p, tmp, "scenario_1", false)
                @test _get_derating(p, "scenario_1") ≈ 0.35 * 1.5
            end
        end
    end

    # ── Regression ────────────────────────────────────────────────────────────

    @testset "T23: No cc_scalar.csv → final CC equals raw derating_dict value (all types)" begin
        mktempdir() do tmp
            _write_derating_dict(tmp, "scenario_1",
                "CT" => 0.75,
                "HY" => 0.6,
                "existing_Wind_Z1" => 0.45,
                "existing_STOR_2" => 0.55,
            )

            p_th = ThermalGenEMIS{Existing}("th1", _thermal_tech("CT"), 2020, 2021, 2040,
                2045, Product[_cap()], _dummy_finance())
            update_derating_factor!(p_th, tmp, "scenario_1", false)
            @test _get_derating(p_th, "scenario_1") ≈ 0.75

            p_hy = HydroGenEMIS{Existing}("hy1", _hydro_tech("HY"), 2020, 2021, 2040,
                2045, Product[_cap()], _dummy_finance())
            update_derating_factor!(p_hy, tmp, "scenario_1", false)
            @test _get_derating(p_hy, "scenario_1") ≈ 0.6

            p_re = RenewableGenEMIS{Existing}("wind1", _renewable_tech("Wind"), 2020, 2021,
                2040, 2045, Product[_cap()], _dummy_finance())
            update_derating_factor!(p_re, tmp, "scenario_1", false)
            @test _get_derating(p_re, "scenario_1") ≈ 0.45

            p_ba = BatteryEMIS{Existing}("batt1", _battery_tech(), 2020, 2021, 2040,
                2045, Product[_cap()], _dummy_finance())
            update_derating_factor!(p_ba, tmp, "scenario_1", false)
            @test _get_derating(p_ba, "scenario_1") ≈ 0.55
        end
    end
end
