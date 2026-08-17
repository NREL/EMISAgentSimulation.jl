# Shared test helpers — safe to include multiple times (const re-declaration
# with the same value is a no-op warning in Julia, not an error).

const TEST_DATA_DIR = joinpath(@__DIR__, "test_data")

"""Minimal Finance instance. update_derating_factor! never reads finance_data."""
function _dummy_finance()
    empty_aa = AxisArray(zeros(Float64, 0, 0), Axis{:product}(Symbol[]), Axis{:year}(Int[]))
    return Finance(
        Float64[], 0.0, Float64[],
        0, 0, 0, 0.0, Float64[],
        Dict{String, Array{Float64, 2}}(),
        Dict{String, Vector{AxisArray{Float64, 2}}}(),
        empty_aa,
        0.0,
        Dict{String, Vector{Float64}}(), Float64[],
        Dict{String, Vector{Float64}}(), Float64[],
        Float64[], "",
    )
end

"""Empty Capacity product."""
_cap() = Capacity(:Capacity, Dict{String, Float64}(), Dict{String, Vector{Float64}}(), 0.0)

"""
Copy the standard derating_dict.csv fixture from test_data/ into the expected
directory structure under `dir`. Covers most tests; use `_write_derating_dict`
when a test needs values that differ from the fixture.

Fixture values (test_data/derating_dict.csv):
  CT=0.85, HY=0.6, existing_Wind_Z1=0.5, new_Wind_Z1=0.35,
  existing_Solar_Z1=0.35, new_Solar_Z1=0.22, existing_STOR_2=0.7, new_STOR_2=0.65
"""
function _setup_derating_dict(dir, scenario)
    d = joinpath(dir, "markets_data", "derating_data", scenario)
    mkpath(d)
    cp(joinpath(TEST_DATA_DIR, "derating_dict.csv"), joinpath(d, "derating_dict.csv"))
end

"""Write a single-row derating_dict.csv with custom column => value pairs (for tests
that need values other than those in the standard fixture)."""
function _write_derating_dict(dir, scenario, cols::Pair{String, <:Real}...)
    d = joinpath(dir, "markets_data", "derating_data", scenario)
    mkpath(d)
    CSV.write(joinpath(d, "derating_dict.csv"),
        DataFrame(Dict(k => [Float64(v)] for (k, v) in cols)))
end

"""Write a single-row cc_scalar.csv under markets_data/derating_data/{scenario}/."""
function _write_cc_scalar(dir, scenario, cols::Pair{String, <:Real}...)
    d = joinpath(dir, "markets_data", "derating_data", scenario)
    mkpath(d)
    CSV.write(joinpath(d, "cc_scalar.csv"),
        DataFrame(Dict(k => [Float64(v)] for (k, v) in cols)))
end

"""Return the derating stored on the first Capacity product of a project."""
function _get_derating(project, scenario)
    for prod in get_products(project)
        d = get_derating(prod)
        isnothing(d) || return d[scenario]
    end
    error("no Capacity product found on $(get_name(project))")
end

# ── Lightweight tech constructors ────────────────────────────────────────────

_thermal_tech(type = "CT") = ThermalTech(
    type, "GAS", (min = 0.0, max = 100.0), nothing, nothing, nothing,
    0.0, [(0.0, 0.0)], "bus1", "Z1", 0.05, 48,
)

_hydro_tech(type = "HY") = HydroTech(
    type, (min = 0.0, max = 100.0), nothing, nothing, nothing, "bus1", "Z1", 0.02, 48,
)

_renewable_tech(type = "Wind") = RenewableTech(
    type, (min = 0.0, max = 100.0), nothing, nothing, "bus1", "Z1", 0.0, 0,
)

# 100 MW output, 200 MWh storage → duration = 2 hours
_battery_tech() = BatteryTech(
    "BATT",
    (min = 0.0, max = 100.0), (min = 0.0, max = 100.0),
    nothing,
    (min = 0.0, max = 200.0), (min = 0.0, max = 1.0),
    0.5, 100.0, 0.5, (in = 0.9, out = 0.9),
    "bus1", "Z1", 0.05, 48, 100.0,
)
