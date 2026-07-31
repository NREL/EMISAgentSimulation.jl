# Test for HDF5 save/load roundtrip of AgentSimulation.
# Run from REPL after loading the package:
#   include("test/test_hdf5_serialization.jl")

using Test
using HDF5
using AxisArrays
using DataFrames
import PowerSystems as PSY

# ─────────────────────────────────────────────────────────────────────────────
# Recursive value comparator — returns true and is silent on success;
# prints a failure message and returns false on first mismatch.
# ─────────────────────────────────────────────────────────────────────────────

function vals_equal(a, b, path::String = "root")::Bool
    if typeof(a) != typeof(b)
        println("FAIL [$path]: type $(typeof(a)) ≠ $(typeof(b))")
        return false
    end
    if a isa AxisArrays.AxisArray
        if collect(a) != collect(b)
            println("FAIL [$path]: AxisArray data differs")
            return false
        end
        for (i, (ax_a, ax_b)) in enumerate(zip(AxisArrays.axes(a), AxisArrays.axes(b)))
            if AxisArrays.axisvalues(ax_a) != AxisArrays.axisvalues(ax_b)
                println("FAIL [$path]: AxisArray axis $i values differ")
                return false
            end
        end
        return true
    end
    if a isa Dict
        if keys(a) != keys(b)
            println("FAIL [$path]: Dict keys $(keys(a)) ≠ $(keys(b))")
            return false
        end
        return all(vals_equal(a[k], b[k], "$path.$k") for k in keys(a))
    end
    if a isa AbstractArray
        if size(a) != size(b)
            println("FAIL [$path]: array size $(size(a)) ≠ $(size(b))")
            return false
        end
        return a == b
    end
    if a isa NamedTuple
        for k in keys(a)
            vals_equal(a[k], b[k], "$path.$k") || return false
        end
        return true
    end
    ok = (a == b)
    ok || println("FAIL [$path]: $a ≠ $b")
    return ok
end

# ─────────────────────────────────────────────────────────────────────────────
# Builders
# ─────────────────────────────────────────────────────────────────────────────

SCENS = ["S1", "S2"]
N_YRS = 3
N_HRS = 4   # hours per year in tiny test arrays
PARAMS = [:cap_factor, :fuel_price]

function make_axis_2d()
    AxisArrays.AxisArray(
        rand(N_YRS, N_HRS),
        AxisArrays.Axis{:year}(1:N_YRS),
        AxisArrays.Axis{:hour}(1:N_HRS),
    )
end

function make_axis_1d()
    AxisArrays.AxisArray(rand(N_YRS), AxisArrays.Axis{:year}(1:N_YRS))
end

function make_finance()
    realized_profit = AxisArrays.AxisArray(
        rand(N_YRS, 2),
        AxisArrays.Axis{:year}(1:N_YRS),
        AxisArrays.Axis{:product}([:energy, :capacity]),
    )
    scenario_profit = Dict(s => [make_axis_2d() for _ in 1:N_YRS] for s in SCENS)
    scenario_total_utilization = Dict(s => rand(N_YRS, N_HRS) for s in SCENS)
    scenario_npv = Dict(s => rand(N_YRS) for s in SCENS)
    scenario_utility = Dict(s => rand(N_YRS) for s in SCENS)

    return Finance(
        rand(N_YRS),         # investment_cost
        1.5e6,               # effective_investment_cost
        fill(1.0, N_YRS),   # preference_multiplier
        2,                   # lag_time
        25,                  # life_time
        20,                  # capex_years
        50_000.0,            # fixed_OM_cost
        rand(N_YRS),         # queue_cost
        scenario_total_utilization,
        scenario_profit,
        realized_profit,
        0.08,                # discount_rate
        scenario_npv,
        rand(N_YRS),         # expected_npv
        scenario_utility,
        rand(N_YRS),         # expected_utility
        rand(N_YRS),         # annual_cashflow
        "TestInvestor",      # ownedby
    )
end

function make_products()
    return Product[
        Energy(
            :energy,
            Dict(s => rand(N_YRS, N_HRS) for s in SCENS),
            3.5,   # marginal_cost
            100.0, # expected_production
        ),
        Capacity(
            :capacity,
            Dict(s => 0.9 for s in SCENS),           # derating
            Dict(s => rand(N_YRS) for s in SCENS),    # accepted_perc
            200.0, # capacity_bid
        ),
        OperatingReserve{ReserveUpEMIS}(:reg_up, 0.1, 1.0),
        OperatingReserve{ReserveDownEMIS}(:reg_dn, 0.1, 1.0),
        REC(:rec, 50.0, rand(N_YRS), 15.0),
        Inertia(:inertia, true, 5.0, 0.5),
        CarbonTax(:co2, 0.05, 8.0, 3.5, rand(N_YRS)),
    ]
end

function make_thermal_tech()
    return ThermalTech(
        "CCGT",
        "Gas",
        (min = 50.0, max = 400.0),
        (up = 100.0, down = 100.0),
        (up = 2.0, down = 4.0),
        PSY.ThermalGenerationCost(;
            variable = PSY.CostCurve(
                PSY.PiecewiseLinearData([(0.0, 0.0), (400.0, 1200.0)]),
            ),
            fixed = 5000.0,
            start_up = (hot = 500.0, warm = 1000.0, cold = 2000.0),
            shut_down = 200.0,
        ),
        4.5,   # fuel_cost
        [(0.0, 8.0), (400.0, 9.5)],  # heat_rate_curve
        "BusA",
        "ZoneA",
        0.05,  # FOR
        24,    # MTTR
    )
end

function make_renewable_tech()
    return RenewableTech(
        "WT",
        (min = 0.0, max = 100.0),
        nothing,
        PSY.RenewableGenerationCost(;
            variable = PSY.CostCurve(PSY.PiecewiseLinearData([(0.0, 0.0), (100.0, 0.0)])),
            curtailment_cost = PSY.CostCurve(
                PSY.PiecewiseLinearData([(0.0, 0.0), (100.0, 0.0)]),
            ),
            fixed = 0.0,
        ),
        "BusB",
        "ZoneB",
        0.01,
        8,
    )
end

function make_hydro_tech()
    return HydroTech(
        "Hydro",
        (min = 0.0, max = 200.0),
        (up = 50.0, down = 50.0),
        nothing,
        PSY.HydroGenerationCost(;
            variable = PSY.CostCurve(PSY.PiecewiseLinearData([(0.0, 0.0), (200.0, 200.0)])),
            fixed = 1000.0,
        ),
        "BusC",
        "ZoneA",
        0.02,
        12,
    )
end

function make_battery_tech()
    return BatteryTech(
        "Li-ion",
        (min = 0.0, max = 50.0),
        (min = 0.0, max = 50.0),
        nothing,
        (min = 0.0, max = 200.0),
        (min = 0.2, max = 0.9),
        0.5,
        200.0,
        0.5,
        (in = 0.92, out = 0.92),
        "BusD",
        "ZoneB",
        0.01,
        4,
        50.0,
    )
end

function make_scenario()
    pv = [make_axis_2d() for _ in 1:N_YRS]
    return Scenario("S1", 0.6, Dict("cap_factor" => 1.05, "fuel_price" => 0.95), pv)
end

function make_scenario_no_multipliers()
    pv = [make_axis_2d() for _ in 1:N_YRS]
    return Scenario("S2", 0.4, nothing, pv)
end

function make_market_prices()
    ep = Dict(
        s => AxisArrays.AxisArray(
            rand(N_YRS, N_HRS, 2),
            AxisArrays.Axis{:year}(1:N_YRS),
            AxisArrays.Axis{:hour}(1:N_HRS),
            AxisArrays.Axis{:zone}(["ZoneA", "ZoneB"]),
        ) for s in SCENS
    )
    rp = Dict("reg_up" => Dict(s => rand(N_YRS, N_HRS) for s in SCENS))
    cp = Dict(s => make_axis_1d() for s in SCENS)
    rcp = Dict(s => make_axis_1d() for s in SCENS)
    inp = Dict(s => make_axis_2d() for s in SCENS)
    return MarketPrices(ep, rp, cp, rcp, inp)
end

function make_investor(name::String)
    projects = Project{<:BuildPhase}[
        ThermalGenEMIS{Existing}(
            "CC_plant",
            make_thermal_tech(),
            1,
            1,
            30,
            30,
            make_products(),
            make_finance(),
        ),
        RenewableGenEMIS{Option}(
            "wind_opt",
            make_renewable_tech(),
            3,
            5,
            99,
            99,
            make_products(),
            make_finance(),
        ),
        HydroGenEMIS{Planned}(
            "hydro_pl",
            make_hydro_tech(),
            2,
            4,
            50,
            50,
            make_products(),
            make_finance(),
        ),
        BatteryEMIS{Queue}(
            "batt_q",
            make_battery_tech(),
            2,
            3,
            25,
            25,
            make_products(),
            make_finance(),
        ),
    ]
    return Investor(
        name,
        "/data/$name",
        projects,
        [:energy, :capacity, :rec],
        rand(N_YRS),                   # carbon_tax
        make_market_prices(),
        24,                             # rep_period_interval
        Dict(s => Dict(y => rand(N_HRS) for y in 1:N_YRS) for s in SCENS),  # rep_hour_weight
        4,                              # avg_block_size
        true,                           # fixed_block_size
        Dict(s => Dict(y => rand(Int, N_HRS, N_HRS) for y in 1:N_YRS) for s in SCENS),  # chron_weights
        Perfect([make_scenario(), make_scenario_no_multipliers()]),
        1.0,                            # cap_cost_multiplier
        (min = 0.8, max = 1.2),             # preference_multiplier_range
        Dict(("wind_opt", "CC_plant") => rand(N_YRS)),  # portfolio_preference_multipliers
        5,                              # max_annual_projects
        RiskNeutral(),
        3,                              # retirement_lookback
    )
end

function make_simulation()
    zones = ["ZoneA", "ZoneB"]
    lines = [ZonalLine("A_to_B", "ZoneA", "ZoneB", 500.0)]
    rep_periods =
        Dict(s => Dict(y => Dict(h => h for h in 1:N_HRS) for y in 1:N_YRS) for s in SCENS)
    hour_weight = Dict(s => Dict(y => rand(N_HRS) for y in 1:N_YRS) for s in SCENS)
    peak_load = Dict(s => Dict(y => 1000.0 * y for y in 1:N_YRS) for s in SCENS)
    markets = Dict(:energy => true, :capacity => true, :rec => false)
    derating_df = DataFrame(; technology = ["CCGT", "WT"], derate = [0.9, 1.0])
    derating_data = Dict("thermal" => derating_df)
    ra = ResourceAdequacy(
        Dict("ZoneA" => 0.1, "ZoneB" => 0.1),
        rand(N_YRS),
        [Dict("eue" => rand(), "lolp" => rand()) for _ in 1:N_YRS],
    )

    return AgentSimulation(
        CaseDefinition(),
        "/results/test",
        2,
        nothing, nothing, nothing,
        Dict{String, PSY.System}(),
        zones,
        lines,
        rep_periods,
        24,
        hour_weight,
        peak_load,
        markets,
        rand(N_YRS),    # carbon_tax
        rand(N_YRS),    # rec_requirement
        [make_investor("inv1"), make_investor("inv2")],
        derating_data,
        Dict("S1" => ra),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Field-by-field comparison helpers
# ─────────────────────────────────────────────────────────────────────────────

function compare_products(orig::Vector, loaded::Vector, path::String)
    @test length(orig) == length(loaded)
    for (i, (po, pl)) in enumerate(zip(orig, loaded))
        @test typeof(po) == typeof(pl)
        for f in fieldnames(typeof(po))
            @test vals_equal(getfield(po, f), getfield(pl, f), "$path.product[$i].$f")
        end
    end
end

function compare_tech(orig, loaded, path::String)
    @test typeof(orig) == typeof(loaded)
    for f in fieldnames(typeof(orig))
        @test vals_equal(getfield(orig, f), getfield(loaded, f), "$path.tech.$f")
    end
end

function compare_finance(orig::Finance, loaded::Finance, path::String)
    for f in fieldnames(Finance)
        @test vals_equal(getfield(orig, f), getfield(loaded, f), "$path.finance.$f")
    end
end

function compare_project(orig, loaded, path::String)
    @test typeof(orig) == typeof(loaded)
    compare_tech(orig.tech, loaded.tech, "$path.$(orig.name)")
    compare_finance(orig.finance_data, loaded.finance_data, "$path.$(orig.name)")
    compare_products(orig.products, loaded.products, "$path.$(orig.name)")
    for f in (:name, :decision_year, :construction_year, :retirement_year, :end_life_year)
        @test vals_equal(getfield(orig, f), getfield(loaded, f), "$path.$(orig.name).$f")
    end
end

function compare_market_prices(orig::MarketPrices, loaded::MarketPrices, path::String)
    for f in fieldnames(MarketPrices)
        @test vals_equal(getfield(orig, f), getfield(loaded, f), "$path.market_prices.$f")
    end
end

function compare_forecast(orig, loaded, path::String)
    @test typeof(orig) == typeof(loaded)
    @test length(orig.scenario_data) == length(loaded.scenario_data)
    for (i, (so, sl)) in enumerate(zip(orig.scenario_data, loaded.scenario_data))
        @test so.name == sl.name
        @test so.probability ≈ sl.probability
        @test vals_equal(
            so.parameter_multipliers,
            sl.parameter_multipliers,
            "$path.scenario[$i].multipliers",
        )
        @test length(so.parameter_values) == length(sl.parameter_values)
        for (j, (pvo, pvl)) in enumerate(zip(so.parameter_values, sl.parameter_values))
            @test vals_equal(pvo, pvl, "$path.scenario[$i].parameter_values[$j]")
        end
    end
end

function compare_investor(orig::Investor, loaded::Investor, path::String)
    for f in (:name, :data_dir, :avg_block_size, :fixed_block_size,
        :rep_period_interval, :cap_cost_multiplier, :max_annual_projects,
        :retirement_lookback, :carbon_tax, :markets)
        @test vals_equal(getfield(orig, f), getfield(loaded, f), "$path.$f")
    end
    @test vals_equal(
        orig.preference_multiplier_range,
        loaded.preference_multiplier_range,
        "$path.preference_multiplier_range",
    )
    @test typeof(orig.risk_preference) == typeof(loaded.risk_preference)
    compare_market_prices(orig.market_prices, loaded.market_prices, path)
    compare_forecast(orig.forecast, loaded.forecast, path)
    @test vals_equal(orig.rep_hour_weight, loaded.rep_hour_weight, "$path.rep_hour_weight")
    @test vals_equal(orig.chron_weights, loaded.chron_weights, "$path.chron_weights")
    @test vals_equal(
        orig.portfolio_preference_multipliers,
        loaded.portfolio_preference_multipliers,
        "$path.portfolio_pref_mult",
    )

    @test length(orig.projects) == length(loaded.projects)
    orig_sorted = sort(orig.projects; by = p -> p.name)
    loaded_sorted = sort(loaded.projects; by = p -> p.name)
    for (po, pl) in zip(orig_sorted, loaded_sorted)
        compare_project(po, pl, path)
    end
end

function compare_simulation(orig::AgentSimulation, loaded::AgentSimulation)
    @testset "Top-level scalars" begin
        @test orig.results_dir == loaded.results_dir
        @test orig.iteration_year == loaded.iteration_year
        @test orig.rep_period_interval == loaded.rep_period_interval
        @test orig.carbon_tax == loaded.carbon_tax
        @test orig.rec_requirement == loaded.rec_requirement
        @test orig.zones == loaded.zones
        @test orig.markets == loaded.markets
    end

    @testset "Lines" begin
        @test length(orig.lines) == length(loaded.lines)
        for (lo, ll) in zip(orig.lines, loaded.lines)
            for f in fieldnames(ZonalLine)
                @test vals_equal(getfield(lo, f), getfield(ll, f), "line.$(lo.name).$f")
            end
        end
    end

    @testset "rep_periods / hour_weight / peak_load" begin
        @test vals_equal(orig.rep_periods, loaded.rep_periods, "rep_periods")
        @test vals_equal(orig.hour_weight, loaded.hour_weight, "hour_weight")
        @test vals_equal(orig.peak_load, loaded.peak_load, "peak_load")
    end

    @testset "derating_data" begin
        @test keys(orig.derating_data) == keys(loaded.derating_data)
        for k in keys(orig.derating_data)
            @test orig.derating_data[k] == loaded.derating_data[k]
        end
    end

    @testset "resource_adequacy" begin
        @test keys(orig.resource_adequacy) == keys(loaded.resource_adequacy)
        for k in keys(orig.resource_adequacy)
            rao, ral = orig.resource_adequacy[k], loaded.resource_adequacy[k]
            @test vals_equal(rao.targets, ral.targets, "ra.$k.targets")
            @test vals_equal(rao.delta_irm, ral.delta_irm, "ra.$k.delta_irm")
            @test vals_equal(rao.metrics, ral.metrics, "ra.$k.metrics")
        end
    end

    @testset "Investors" begin
        @test length(orig.investors) == length(loaded.investors)
        orig_sorted = sort(orig.investors, by = i -> i.name)
        loaded_sorted = sort(loaded.investors, by = i -> i.name)
        for (io, il) in zip(orig_sorted, loaded_sorted)
            @testset "Investor $(io.name)" begin
                compare_investor(io, il, "investor.$(io.name)")
            end
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Run
# ─────────────────────────────────────────────────────────────────────────────

@testset "HDF5 AgentSimulation roundtrip" begin
    tmpfile = tempname() * ".h5"
    try
        sim_orig = make_simulation()

        @testset "save_simulation" begin
            @test_nowarn save_simulation(tmpfile, sim_orig)
            @test isfile(tmpfile)
        end

        sim_loaded = nothing
        @testset "load_simulation" begin
            @test_nowarn (sim_loaded = load_simulation(tmpfile))
            @test sim_loaded isa AgentSimulation
        end

        @testset "roundtrip comparison" begin
            compare_simulation(sim_orig, sim_loaded)
        end
    finally
        isfile(tmpfile) && rm(tmpfile)
    end
end
