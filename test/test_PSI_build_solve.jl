# Build + solve tests for the UC/MD/ED PSI templates (src/PSI_definitions.jl,
# markets_simulation/siip_simulation_definition.jl) against the bundled
# test/test_systems fixtures.
#
# Goal: catch template/formulation regressions (wrong component type, missing
# device model, broken service model config) by actually constructing and
# solving a PSI DecisionModel, without running a full multi-year EMIS case.
#
# These are heavier than test_PSI.jl: each @testset loads a full-scale
# serialized PSY.System (thousands of components) and solves a UC/MD/ED
# DecisionModel with Xpress. Expect real wall-clock time (minutes, not
# seconds) — this is an integration-style check, not a fast unit test.
#
# Run from REPL:
#   include("test/test_PSI_build_solve.jl")

include(joinpath(@__DIR__, "includes.jl"))

const TEST_SYSTEMS_DIR = joinpath(@__DIR__, "test_systems")

# Same Xpress settings as production (run_gs_15_years.jl). NOTE: MAXTIME is
# 7200s (2hr) there — kept identical for parity, but that means a stalled
# solve here can run that long before Xpress itself gives up.
function _test_solver()
    return JuMP.optimizer_with_attributes(
        Xpress.Optimizer,
        "MIPRELSTOP" => 1e-2,
        "BARGAPSTOP" => 1e-4,
        "BARDUALSTOP" => 1e-4,
        "BARPRIMALSTOP" => 1e-4,
        "MATRIXTOL" => 1e-8,
        "BARORDER" => 1,
        "CHOLESKYTOL" => 1e-20,
        "CHOLESKYALG" => 1,
        "BARTHREADS" => 8,
        "SCALING" => 2,
        "BARSTEPSTOP" => 1e-10,
        "OUTPUTLOG" => 0,
        "PRESOLVE" => 1,
        "MAXTIME" => 7200,
        "NUMERICALEMPHASIS" => 1,
    )
end

function _load_test_system(json_name::String)
    return PSY.System(joinpath(TEST_SYSTEMS_DIR, json_name); runchecks = false)
end

# Matches the kwargs create_problem() uses in production (siip_simulation_definition.jl)
# for every problem type. `initialize_model = false` in particular is required: without
# it, PSI.DecisionModel defaults to running an internal warm-start sub-solve at build
# time (build_initial_conditions! + initialize!) that these standalone (non-Simulation)
# models aren't set up for and that fails with "Optimizer returned NO_SOLUTION".
function _test_decision_model(template::PSI.ProblemTemplate, sys::PSY.System, name::String)
    return PSI.DecisionModel(
        template,
        sys;
        optimizer = _test_solver(),
        name = name,
        optimizer_solve_log_print = false,
        warm_start = true,
        calculate_conflict = true,
        store_variable_names = true,
        export_pwl_vars = true,
        initialize_model = false,
    )
end

# PSI only @warns (doesn't error) when a system component type has no
# DeviceModel/ServiceModel in the template — those components are silently
# dropped from the optimization instead of failing the build. This makes that
# failure mode a hard test failure.
function _unmodeled_component_types(template::PSI.ProblemTemplate, sys::PSY.System)
    modeled_types = Set(PSI.get_component_types(template))
    return [t for t in PSY.get_existing_component_types(sys) if !(t in modeled_types)]
end

function _build_and_solve!(model::PSI.DecisionModel)
    @testset "model builds" begin
        build_status = PSI.build!(model; output_dir = mktempdir())
        @test build_status == PSI.ModelBuildStatus.BUILT
    end

    @testset "model solves" begin
        PSI.solve!(model; output_dir = mktempdir())
        @test PSI.get_run_status(model) == PSI.RunStatus.SUCCESSFULLY_FINALIZED
    end
end

@testset "PSI template build+solve against bundled test systems" begin

    @testset "UC template + DA test system" begin
        sys = _load_test_system("sys_UC_year1.json")
        template = EMISAgentSimulation.create_uc_template()

        @testset "no unmodeled component types" begin
            @test isempty(_unmodeled_component_types(template, sys))
        end

        model = _test_decision_model(template, sys, "UC_test")
        _build_and_solve!(model)
    end

    @testset "MD template + MD test system" begin
        sys = _load_test_system("sys_MD_year1.json")
        template = EMISAgentSimulation.create_md_template()

        @testset "no unmodeled component types" begin
            @test isempty(_unmodeled_component_types(template, sys))
        end

        model = _test_decision_model(template, sys, "MD_test")
        _build_and_solve!(model)
    end

    @testset "ED template (no inertia) + RT test system" begin
        sys = _load_test_system("sys_ED_year1.json")
        inertia_product = collect(PSY.get_components_by_name(PSY.Service, sys, "Inertia"))
        template = EMISAgentSimulation.create_ed_template(inertia_product)

        @testset "system has no Inertia service (base ED path)" begin
            # If this ever fails, the test system gained an Inertia service and
            # the ED template would take the ORDC branch — add a matching
            # "ED template (inertia) + RT test system" testset alongside this one.
            @test isempty(inertia_product)
        end

        @testset "no unmodeled component types" begin
            @test isempty(_unmodeled_component_types(template, sys))
        end

        model = _test_decision_model(template, sys, "ED_test")
        _build_and_solve!(model)
    end

end
