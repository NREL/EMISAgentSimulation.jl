# PSI/PSY API compatibility tests.
# Guards against regressions from the PSI 0.31→0.37 / PSY 4.x→5.11 upgrade.
# Run from REPL:
#   include("test/test_PSI.jl")

include(joinpath(@__DIR__, "includes.jl"))

@testset "PSI/PSY API compatibility (PSI ~0.37 / PSY ~5.11)" begin

    # ── Package version checks ──────────────────────────────────────────────
    @testset "Package versions are in expected range" begin
        @test pkgversion(PowerSimulations) >= v"0.37.0"
        @test pkgversion(PowerSimulations) < v"0.38.0"
        @test pkgversion(PowerSystems) >= v"5.0.0"
    end

    # ── Regression: serialize kwarg removed from PSI.build! in 0.37 ────────
    @testset "PSI.build! has no serialize kwarg" begin
        # Old call was PSI.build!(sim; serialize=false); 0.37 dropped the kwarg.
        @test !any(m -> :serialize ∈ Base.kwarg_decl(m), methods(PSI.build!))
    end

    # ── Regression: PSY 5.x renamed to plural form ─────────────────────────
    @testset "PSY.get_time_series_resolutions (plural form) exists" begin
        @test hasmethod(PSY.get_time_series_resolutions, Tuple{PSY.System})
    end

    # ── Key PSI network / simulation types ─────────────────────────────────
    @testset "Key PSI model and simulation types exist" begin
        for sym in [
            :ProblemTemplate, :NetworkModel, :AreaBalancePowerModel,
            :CopperPlateBalanceConstraint, :DecisionModel, :SimulationModels,
            :SimulationSequence, :Simulation, :SimulationResults,
            :ServiceModel, :RequirementConstraint, :InterProblemChronology,
            :OptimizationContainer, :OperationModel, :SimulationState,
        ]
            @test isdefined(PSI, sym) broken=false
        end
    end

    # ── Device formulations used in UC/ED/MD templates ─────────────────────
    @testset "PSI device formulations used in templates exist" begin
        for sym in [
            :ThermalBasicUnitCommitment, :ThermalBasicDispatch,
            :RenewableFullDispatch, :FixedOutput, :StaticPowerLoad,
            :StaticBranch, :HVDCTwoTerminalLossless, :RangeReserve,
        ]
            @test isdefined(PSI, sym) broken=false
        end
    end

    # ── Variable / constraint types ─────────────────────────────────────────
    @testset "PSI variable and constraint types exist" begin
        for sym in [
            :EnergyVariable, :OnVariable, :ActivePowerVariable,
            :SystemBalanceSlackUp, :SystemBalanceSlackDown,
            :ReserveRequirementSlack, :VariableKey,
            :ActivePowerReserveVariable, :DeviceStatus,
            :InitialEnergyLevel, :InitialCondition,
        ]
            @test isdefined(PSI, sym) broken=false
        end
    end

    # ── SimulationState APIs used in update_initial_conditions! overrides ──
    @testset "PSI simulation state APIs exist" begin
        for sym in [
            :get_system_state_value, :get_decision_state_data,
            :get_dataset_value, :get_update_timestamp, :get_system_state_data,
            :set_ic_quantity!, :get_value,
        ]
            @test isdefined(PSI, sym) broken=false
        end
    end

    # ── Internal PSI APIs used in add_constraint_dual! / adjust_reserve_voll! ─
    @testset "PSI internal APIs used in custom dispatch overrides exist" begin
        for sym in [
            :get_duals, :get_service_name, :get_available_components,
            :get_time_steps, :add_dual_container!, :assign_dual_variable!,
            :should_write_resulting_value, :get_optimization_container,
            :get_variables, :add_to_objective_variant_expression!,
            :get_component, :get_component_type, :get_component_name,
        ]
            @test isdefined(PSI, sym) broken=false
        end
    end

    # ── Feedforward types ───────────────────────────────────────────────────
    @testset "PSI feedforward types exist" begin
        @test isdefined(PSI, :SemiContinuousFeedforward)
    end

    # ── HydroPowerSimulations formulations ─────────────────────────────────
    @testset "HydroPowerSimulations formulations exist" begin
        @test isdefined(HSI, :HydroCommitmentRunOfRiver)
        @test isdefined(HSI, :HydroDispatchRunOfRiver)
    end

    # ── StorageSystemsSimulations types ────────────────────────────────────
    @testset "StorageSystemsSimulations types exist" begin
        @test isdefined(SSI, :StorageDispatchWithReserves)
        @test isdefined(SSI, :StorageEnergyShortageVariable)
        @test isdefined(SSI, :EnergyTargetFeedforward)
    end

    # ── Simulation results API ──────────────────────────────────────────────
    @testset "PSI simulation results API exists" begin
        @test hasmethod(PSI.get_decision_problem_results, Tuple{PSI.SimulationResults, String})
        @test isdefined(PSI, :read_realized_duals)
        @test isdefined(PSI, :read_realized_variables)
    end

    # ── Template construction (no PSY system needed) ───────────────────────
    @testset "UC template construction without error" begin
        template = EMISAgentSimulation.create_uc_template()
        @test template isa PSI.ProblemTemplate
    end

    @testset "ED template construction without error (no inertia)" begin
        template = EMISAgentSimulation.create_ed_template([])
        @test template isa PSI.ProblemTemplate
    end

    @testset "MD template construction without error" begin
        template = EMISAgentSimulation.create_md_template()
        @test template isa PSI.ProblemTemplate
    end

    # ── Load types: confirm StaticPowerLoad is registered for the load type(s)
    # actually present in the bundled test system, and surface any other load
    # types present so they don't silently fall out of the optimization model
    # (PSI only @warns, doesn't error, on an unmodeled component type — see the
    # PowerLoad/StandardLoad mismatch fixed in PSI_definitions.jl). ──────────
    @testset "Load device models cover the load types present in the test system" begin
        sys = PSY.System(
            joinpath(@__DIR__, "test_systems", "sys_MD_year1.json");
            runchecks = false,
        )

        function all_concrete_subtypes(T)
            out = DataType[]
            for S in InteractiveUtils.subtypes(T)
                if isconcretetype(S)
                    push!(out, S)
                else
                    append!(out, all_concrete_subtypes(S))
                end
            end
            return out
        end

        # Count every concrete PSY.StaticLoad subtype present in the system,
        # not just the ones we expect (PowerLoad/StandardLoad) — this is what
        # would have caught the original PowerLoad-vs-StandardLoad mismatch.
        load_counts = Dict{DataType, Int}()
        for T in all_concrete_subtypes(PSY.StaticLoad)
            n = length(collect(PSY.get_components(T, sys)))
            n > 0 && (load_counts[T] = n)
        end

        @info "Load component types present in MD test system" load_counts

        # PSI_definitions.jl registers PSI.StaticPowerLoad for exactly these
        # two load types across UC/MD/ED. Every load type actually present in
        # the system must be one of them, or it's being silently dropped.
        covered_load_types = Set([PSY.PowerLoad, PSY.StandardLoad])
        uncovered = setdiff(keys(load_counts), covered_load_types)
        @test isempty(uncovered)

        # Sanity check this isn't vacuously passing against a system with no
        # load components at all.
        @test sum(values(load_counts); init = 0) > 0

        # Today's test system is expected to be PowerLoad-only. This isn't a
        # correctness requirement (StandardLoad would also be fine, it's
        # covered above) — it's a canary: if this flips, the "PowerLoad-only"
        # assumption elsewhere (e.g. this test's own framing) is stale.
        @test get(load_counts, PSY.PowerLoad, 0) > 0
        @test get(load_counts, PSY.StandardLoad, 0) == 0
    end

    # ── Custom dispatch overrides ───────────────────────────────────────────
    @testset "Custom should_write_resulting_value override suppresses StorageEnergyShortageVariable" begin
        @test PSI.should_write_resulting_value(SSI.StorageEnergyShortageVariable) == false
    end

    @testset "Custom add_constraint_dual! overrides are defined" begin
        # ServiceModel dispatch (fixes dual registration for named services)
        @test hasmethod(
            PSI.add_constraint_dual!,
            Tuple{
                PSI.OptimizationContainer,
                PSY.System,
                PSI.ServiceModel{<:PSY.Service, <:PSI.AbstractServiceFormulation},
            },
        )
        # AreaBalancePowerModel dispatch (fixes dual keyed on Area vs ACBus)
        @test hasmethod(
            PSI.add_constraint_dual!,
            Tuple{
                PSI.OptimizationContainer,
                PSY.System,
                PSI.NetworkModel{PSI.AreaBalancePowerModel},
            },
        )
    end

end
