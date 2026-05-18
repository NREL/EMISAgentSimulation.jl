using Test
using EMISAgentSimulation
const EAS = EMISAgentSimulation

import PowerSystems
const PSY = PowerSystems
import PowerSimulations
const PSI = PowerSimulations
import StorageSystemsSimulations
const SSI = StorageSystemsSimulations
import DataFrames

# ─────────────────────────────────────────────────────────────────────────────
# Helper: build a minimal ThermalStandard device with a custom name and size.
# ─────────────────────────────────────────────────────────────────────────────
function make_thermal_device(name::String, max_pu::Float64)
    dev = PSY.ThermalStandard(nothing)
    PSY.set_name!(dev, name)
    PSY.set_active_power_limits!(dev, (min = 0.0, max = max_pu))
    return dev
end

# ─────────────────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────────────────
@testset "ORDC Static Reserve unit tests" begin

    # ── 1. process_ordc_data_for_siip ────────────────────────────────────────
    @testset "process_ordc_data_for_siip – flat (single-segment) ORDC" begin
        # Two-point curve: origin (0 MW, $4000) + MRR endpoint (200 MW, $4000)
        # String format from actual timeseries CSV:  [(qty, price), (qty, price)]
        flat_str = ["[(0.0, 4000.0), (200.0, 4000.0)]"]
        result = EAS.process_ordc_data_for_siip(flat_str)

        @test length(result) == 1            # one timestep
        @test length(result[1]) == 1         # one processed segment (origin skipped)
        # Second element of the tuple is MRR in MW
        @test result[1][1][2] ≈ 200.0
        # First element is cumulative cost = price × qty
        @test result[1][1][1] ≈ 4000.0 * 200.0
    end

    @testset "process_ordc_data_for_siip – curved (multi-segment) ORDC" begin
        # Four points: origin + 3 data segments (typical curved ORDC)
        curved_str = ["[(0.0, 4000.0), (1354.3, 4000.0), (2726.82, 3162.58), (4099.34, 2582.83)]"]
        result = EAS.process_ordc_data_for_siip(curved_str)

        @test length(result) == 1
        @test length(result[1]) == 3         # 3 segments after skipping origin
        # MRR is the qty (second element) of the first processed tuple
        @test result[1][1][2] ≈ 1354.3
        # Cumulative cost must be strictly increasing across segments
        @test result[1][2][1] > result[1][1][1]
        @test result[1][3][1] > result[1][2][1]
    end

    @testset "process_ordc_data_for_siip – multiple timesteps, MRR varies" begin
        mrr1, mrr2 = 200.0, 350.0
        ts = [
            "[(0.0, 4000.0), ($(mrr1), 4000.0)]",
            "[(0.0, 4000.0), ($(mrr2), 4000.0)]",
        ]
        result = EAS.process_ordc_data_for_siip(ts)

        @test length(result) == 2
        @test result[1][1][2] ≈ mrr1
        @test result[2][1][2] ≈ mrr2
    end

    # ── 2. Template service model selection ──────────────────────────────────
    # PSI template.services is Dict{Tuple{String,Symbol}, ServiceModel}.
    # Key format: (service_name, Symbol(ServiceType)).
    # PSI.get_component_type and PSI.get_formulation inspect the parameterised
    # ServiceModel{D, B} type.

    @testset "create_md_template – ordc_curved=false → StaticReserve / RangeReserve" begin
        template = EAS.create_md_template([], false)
        svc = template.services

        # Both "Synchronous" and "Primary" should be keyed as StaticReserve
        @test haskey(svc, ("Synchronous", Symbol(PSY.StaticReserve{PSY.ReserveUp})))
        @test haskey(svc, ("Primary",     Symbol(PSY.StaticReserve{PSY.ReserveUp})))
        # Must NOT be keyed as ReserveDemandCurve
        @test !haskey(svc, ("Synchronous", Symbol(PSY.ReserveDemandCurve{PSY.ReserveUp})))
        @test !haskey(svc, ("Primary",     Symbol(PSY.ReserveDemandCurve{PSY.ReserveUp})))

        synch = svc[("Synchronous", Symbol(PSY.StaticReserve{PSY.ReserveUp}))]
        prim  = svc[("Primary",     Symbol(PSY.StaticReserve{PSY.ReserveUp}))]
        @test PSI.get_component_type(synch) == PSY.StaticReserve{PSY.ReserveUp}
        @test PSI.get_formulation(synch)    == PSI.RangeReserve
        @test PSI.get_component_type(prim)  == PSY.StaticReserve{PSY.ReserveUp}
        @test PSI.get_formulation(prim)     == PSI.RangeReserve
    end

    @testset "create_md_template – ordc_curved=true → ReserveDemandCurve / StepwiseCostReserve" begin
        template = EAS.create_md_template([], true)
        svc = template.services

        @test haskey(svc, ("Synchronous", Symbol(PSY.ReserveDemandCurve{PSY.ReserveUp})))
        @test haskey(svc, ("Primary",     Symbol(PSY.ReserveDemandCurve{PSY.ReserveUp})))
        @test !haskey(svc, ("Synchronous", Symbol(PSY.StaticReserve{PSY.ReserveUp})))
        @test !haskey(svc, ("Primary",     Symbol(PSY.StaticReserve{PSY.ReserveUp})))

        synch = svc[("Synchronous", Symbol(PSY.ReserveDemandCurve{PSY.ReserveUp}))]
        prim  = svc[("Primary",     Symbol(PSY.ReserveDemandCurve{PSY.ReserveUp}))]
        @test PSI.get_component_type(synch) == PSY.ReserveDemandCurve{PSY.ReserveUp}
        @test PSI.get_formulation(synch)    == PSI.StepwiseCostReserve
        @test PSI.get_component_type(prim)  == PSY.ReserveDemandCurve{PSY.ReserveUp}
        @test PSI.get_formulation(prim)     == PSI.StepwiseCostReserve
    end

    @testset "create_uc_template – ordc_curved=false → StaticReserve / RangeReserve" begin
        template = EAS.create_uc_template([], false)
        svc = template.services

        @test haskey(svc, ("Synchronous", Symbol(PSY.StaticReserve{PSY.ReserveUp})))
        @test haskey(svc, ("Primary",     Symbol(PSY.StaticReserve{PSY.ReserveUp})))
        @test !haskey(svc, ("Synchronous", Symbol(PSY.ReserveDemandCurve{PSY.ReserveUp})))

        m = svc[("Synchronous", Symbol(PSY.StaticReserve{PSY.ReserveUp}))]
        @test PSI.get_component_type(m) == PSY.StaticReserve{PSY.ReserveUp}
        @test PSI.get_formulation(m)    == PSI.RangeReserve
    end

    @testset "create_uc_template – ordc_curved=true → ReserveDemandCurve / StepwiseCostReserve" begin
        template = EAS.create_uc_template([], true)
        svc = template.services

        @test haskey(svc, ("Synchronous", Symbol(PSY.ReserveDemandCurve{PSY.ReserveUp})))
        @test haskey(svc, ("Primary",     Symbol(PSY.ReserveDemandCurve{PSY.ReserveUp})))
        @test !haskey(svc, ("Synchronous", Symbol(PSY.StaticReserve{PSY.ReserveUp})))

        m = svc[("Synchronous", Symbol(PSY.ReserveDemandCurve{PSY.ReserveUp}))]
        @test PSI.get_component_type(m) == PSY.ReserveDemandCurve{PSY.ReserveUp}
        @test PSI.get_formulation(m)    == PSI.StepwiseCostReserve
    end

    @testset "create_ed_template – ordc_curved=false → StaticReserve / RangeReserve" begin
        template = EAS.create_ed_template([], false)
        svc = template.services

        @test haskey(svc, ("Synchronous", Symbol(PSY.StaticReserve{PSY.ReserveUp})))
        @test haskey(svc, ("Primary",     Symbol(PSY.StaticReserve{PSY.ReserveUp})))
        @test !haskey(svc, ("Synchronous", Symbol(PSY.ReserveDemandCurve{PSY.ReserveUp})))

        m = svc[("Synchronous", Symbol(PSY.StaticReserve{PSY.ReserveUp}))]
        @test PSI.get_component_type(m) == PSY.StaticReserve{PSY.ReserveUp}
        @test PSI.get_formulation(m)    == PSI.RangeReserve
    end

    @testset "create_ed_template – ordc_curved=true → ReserveDemandCurve / StepwiseCostReserve" begin
        template = EAS.create_ed_template([], true)
        svc = template.services

        @test haskey(svc, ("Synchronous", Symbol(PSY.ReserveDemandCurve{PSY.ReserveUp})))
        @test haskey(svc, ("Primary",     Symbol(PSY.ReserveDemandCurve{PSY.ReserveUp})))
        @test !haskey(svc, ("Synchronous", Symbol(PSY.StaticReserve{PSY.ReserveUp})))

        m = svc[("Synchronous", Symbol(PSY.ReserveDemandCurve{PSY.ReserveUp}))]
        @test PSI.get_component_type(m) == PSY.ReserveDemandCurve{PSY.ReserveUp}
        @test PSI.get_formulation(m)    == PSI.StepwiseCostReserve
    end

    # ── 3. feedforward_dict ED key has both feedforward types ─────────────────
    # Regression test: the pre-fix code had a duplicate "ED" key in the Dict
    # literal, silently dropping EnergyTargetFeedforward.  Verify the corrected
    # structure (matching the no-MD, no-inertia branch in create_simulation).
    @testset "feedforward_dict – ED contains both EnergyTargetFeedforward and SemiContinuousFeedforward" begin
        feedforward_dict = Dict(
            "ED" => [
                SSI.EnergyTargetFeedforward(
                    component_type  = PSY.GenericBattery,
                    source          = PSI.EnergyVariable,
                    affected_values = [PSI.EnergyVariable],
                    target_period   = 2,
                    penalty_cost    = 5000.0,
                ),
                PSI.SemiContinuousFeedforward(
                    component_type  = PSY.ThermalStandard,
                    source          = PSI.OnVariable,
                    affected_values = [PSI.ActivePowerVariable],
                ),
            ],
        )

        @test haskey(feedforward_dict, "ED")
        ed_ffs = feedforward_dict["ED"]
        @test length(ed_ffs) == 2
        @test any(ff -> isa(ff, SSI.EnergyTargetFeedforward),    ed_ffs)
        @test any(ff -> isa(ff, PSI.SemiContinuousFeedforward), ed_ffs)
    end

    # ── 4. update_realized_reserve_perc! – StaticReserve dispatch ────────────
    # Verify that the StaticReserve{ReserveUp} method reads results from the
    # "ActivePowerReserveVariable__StaticReserve__ReserveUp__<name>" keys and
    # correctly normalises by device size and base_power.
    @testset "update_realized_reserve_perc! – StaticReserve reads correct result keys" begin
        service_name = "Synchronous"
        device_name  = "Gen1"
        device_max   = 1.0   # p.u. (active_power_limits max)
        base_power   = 100.0 # MW

        provision_ed = [10.0, 20.0, 30.0, 40.0]
        provision_uc = [12.0, 22.0, 32.0, 42.0]
        n_ed = length(provision_ed)
        n_uc = length(provision_uc)

        ed_key = "ActivePowerReserveVariable__StaticReserve__ReserveUp__$(service_name)"
        uc_key = "ActivePowerReserveVariable__StaticReserve__ReserveUp__$(service_name)"

        results_ed = Dict{String, DataFrames.DataFrame}(
            ed_key => DataFrames.DataFrame(Symbol(device_name) => provision_ed),
        )
        results_uc = Dict{String, DataFrames.DataFrame}(
            uc_key => DataFrames.DataFrame(Symbol(device_name) => provision_uc),
        )

        # PSY.StaticReserve{ReserveUp}(name, available, time_frame, requirement)
        reserve = PSY.StaticReserve{PSY.ReserveUp}(service_name, true, 600.0, 2.0)

        device = make_thermal_device(device_name, device_max)

        reserve_perc_ed = Dict(device_name => Dict(service_name => zeros(Float64, 1, n_ed)))
        reserve_perc_uc = Dict(device_name => Dict(service_name => zeros(Float64, 1, n_uc)))
        reserve_perc_md = Dict{String, Dict{String, Array{Float64, 2}}}()
        inertia_perc    = Dict{String, Array{Float64, 2}}()
        rt_products     = SubString{String}[]
        da_products     = SubString{String}[]
        md_products     = SubString{String}[]

        # Should not throw – reads from __StaticReserve__ keys
        EAS.update_realized_reserve_perc!(
            device, reserve,
            results_ed, results_uc, nothing,
            reserve_perc_md, reserve_perc_uc, reserve_perc_ed,
            inertia_perc, rt_products, da_products, md_products,
            base_power, false,
        )

        expected_ed = provision_ed ./ device_max ./ base_power
        expected_uc = provision_uc ./ device_max ./ base_power
        @test reserve_perc_ed[device_name][service_name][1, :] ≈ expected_ed
        @test reserve_perc_uc[device_name][service_name][1, :] ≈ expected_uc
    end

    @testset "update_realized_reserve_perc! – StaticReserve with md_market_bool=true" begin
        service_name = "Primary"
        device_name  = "Gen2"
        device_max   = 2.0
        base_power   = 100.0

        provision_md = [5.0, 10.0]
        n = length(provision_md)

        ed_key = "ActivePowerReserveVariable__StaticReserve__ReserveUp__$(service_name)"
        uc_key = "ActivePowerReserveVariable__StaticReserve__ReserveUp__$(service_name)"
        md_key = "ActivePowerReserveVariable__StaticReserve__ReserveUp__$(service_name)"

        results_ed = Dict{String, DataFrames.DataFrame}(
            ed_key => DataFrames.DataFrame(Symbol(device_name) => zeros(Float64, n)),
        )
        results_uc = Dict{String, DataFrames.DataFrame}(
            uc_key => DataFrames.DataFrame(Symbol(device_name) => zeros(Float64, n)),
        )
        results_md = Dict{String, DataFrames.DataFrame}(
            md_key => DataFrames.DataFrame(Symbol(device_name) => provision_md),
        )

        reserve = PSY.StaticReserve{PSY.ReserveUp}(service_name, true, 1800.0, 3.0)
        device  = make_thermal_device(device_name, device_max)

        reserve_perc_ed = Dict(device_name => Dict(service_name => zeros(Float64, 1, n)))
        reserve_perc_uc = Dict(device_name => Dict(service_name => zeros(Float64, 1, n)))
        reserve_perc_md = Dict(device_name => Dict(service_name => zeros(Float64, 1, n)))
        inertia_perc    = Dict{String, Array{Float64, 2}}()

        EAS.update_realized_reserve_perc!(
            device, reserve,
            results_ed, results_uc, results_md,
            reserve_perc_md, reserve_perc_uc, reserve_perc_ed,
            inertia_perc, SubString{String}[], SubString{String}[], SubString{String}[],
            base_power, true,
        )

        expected_md = provision_md ./ device_max ./ base_power
        @test reserve_perc_md[device_name][service_name][1, :] ≈ expected_md
    end

end # @testset "ORDC Static Reserve unit tests"
