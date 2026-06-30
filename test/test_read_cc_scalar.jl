# Tests for read_cc_scalar.
# src: src/derating_factor_updates/derating_factor_calculator.jl
# Run from REPL after loading the package:
#   using EMISAgentSimulation
#   include("test/test_read_cc_scalar.jl")

using Test
using DataFrames
using CSV
using EMISAgentSimulation
include(joinpath(@__DIR__, "helpers.jl"))

@testset "read_cc_scalar" begin
    @testset "T1: cc_scalar.csv present, column present → returns value" begin
        mktempdir() do tmp
            _write_cc_scalar(tmp, "scenario_1", "CT" => 2.0)
            @test read_cc_scalar(tmp, "scenario_1", "CT") ≈ 2.0
        end
    end

    @testset "T2: cc_scalar.csv present, column absent → returns 1.0" begin
        mktempdir() do tmp
            _write_cc_scalar(tmp, "scenario_1", "CT" => 2.0)
            @test read_cc_scalar(tmp, "scenario_1", "Wind") ≈ 1.0
        end
    end

    @testset "T3: cc_scalar.csv missing entirely → returns 1.0, no error" begin
        mktempdir() do tmp
            mkpath(joinpath(tmp, "markets_data", "derating_data", "scenario_1"))
            @test read_cc_scalar(tmp, "scenario_1", "CT") ≈ 1.0
        end
    end

    @testset "T4: cc_scalar.csv present but empty (no rows) → returns 1.0, no error" begin
        mktempdir() do tmp
            d = joinpath(tmp, "markets_data", "derating_data", "scenario_1")
            mkpath(d)
            # Write header-only CSV (no data rows)
            write(joinpath(d, "cc_scalar.csv"), "CT\n")
            @test read_cc_scalar(tmp, "scenario_1", "CT") ≈ 1.0
        end
    end

    @testset "T5: cc_scalar.csv present, column present but value is missing → returns 1.0" begin
        mktempdir() do tmp
            d = joinpath(tmp, "markets_data", "derating_data", "scenario_1")
            mkpath(d)
            # Write CSV with a missing value in the CT column
            write(joinpath(d, "cc_scalar.csv"), "CT\n\n")
            @test read_cc_scalar(tmp, "scenario_1", "CT") ≈ 1.0
        end
    end
end
