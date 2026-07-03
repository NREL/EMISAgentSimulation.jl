include("includes.jl")

const LOG_FILE = "emis-agent-simulation-test.log"

# Can generate with: ls -1 test | grep "test_.*\.jl"
const DISABLED_TEST_FILES = [
    "test_PSI_build_solve.jl",  # slow: real Xpress MIP solve against full-scale test systems
    "test_hdf5_serialization.jl",  # broken: make_simulation() calls CaseDefinition() with no
    # matching zero-arg constructor (real one takes 7+ positional args, ~44 validated kwargs,
    # and has a filesystem side effect via make_case_data_dir). Pre-existing, not yet fixed.
]

# Self-contained: run without any external setup. `@includetests` defaults to
# this list when called with no ARGS (see below) — it does NOT scan the
# directory, unlike the analogous macro in PowerSimulations.jl, because this
# directory also has context-dependent files (see CONTEXT_DEPENDENT_TEST_FILES)
# that error without caller-provided variables and must never be auto-discovered.
const SELF_CONTAINED_TEST_FILES = [
    "test_PSI.jl",
    "test_PSI_build_solve.jl",
    "test_hdf5_serialization.jl",
]

# Require the caller to have defined, before including this file:
#   case, simulation_years, test_system_dir  (for test_base_checks.jl)
#   projectdata, investor_dir                (for test_system_inputs.jl)
# Always run (not affected by ARGS filtering) to preserve the existing driver-
# script workflow: define those variables, then `include("test/runtests.jl")`.
const CONTEXT_DEPENDENT_TEST_FILES = [
    "test_base_checks.jl",
    "test_system_inputs.jl",
]

const LOG_LEVELS = Dict(
    "Debug" => Logging.Debug,
    "Info" => Logging.Info,
    "Warn" => Logging.Warn,
    "Error" => Logging.Error,
)

function get_logging_level(env_name::String, default)
    level = get(ENV, env_name, default)
    log_level = get(LOG_LEVELS, level, nothing)
    if isnothing(log_level)
        error("Invalid log level $level: Supported levels: $(values(LOG_LEVELS))")
    end

    return log_level
end

"""
Includes the given test files, given as a list without their ".jl" extensions
(e.g. `julia test/runtests.jl test_PSI test_hdf5_serialization` to run just
those two). If none are given, includes all of SELF_CONTAINED_TEST_FILES.
Entries in DISABLED_TEST_FILES are skipped either way.
"""
macro includetests(testarg...)
    if length(testarg) == 0
        tests = []
    elseif length(testarg) == 1
        tests = testarg[1]
    else
        error("@includetests takes zero or one argument")
    end

    quote
        tests = $tests
        if length(tests) == 0
            tests = copy(SELF_CONTAINED_TEST_FILES)
        else
            tests = map(f -> string(f, ".jl"), tests)
        end
        println()
        if !isempty(DISABLED_TEST_FILES)
            @warn("Some tests are disabled $DISABLED_TEST_FILES")
        end
        for test in tests
            test ∈ DISABLED_TEST_FILES && continue
            print(splitext(test)[1], ": ")
            include(test)
            println()
        end
    end
end

function get_logging_level_from_env(env_name::String, default)
    level = get(ENV, env_name, default)
    return IS.get_logging_level(level)
end

function run_tests()
    logging_config_filename = get(ENV, "EMIS_LOGGING_CONFIG", nothing)
    if !isnothing(logging_config_filename)
        config = IS.LoggingConfiguration(logging_config_filename)
    else
        config = IS.LoggingConfiguration(;
            filename = LOG_FILE,
            file_level = Logging.Info,
            console_level = Logging.Error,
        )
    end
    console_logger = ConsoleLogger(config.console_stream, config.console_level)

    IS.open_file_logger(LOG_FILE, config.file_level) do file_logger
        levels = (Logging.Info, Logging.Warn, Logging.Error)
        multi_logger =
            IS.MultiLogger([console_logger, file_logger], IS.LogEventTracker(levels))
        global_logger(multi_logger)

        if !isempty(config.group_levels)
            IS.set_group_levels!(multi_logger, config.group_levels)
        end

        @time @testset "Begin EMISAgentSimulation tests" begin
            @includetests ARGS

            for test in CONTEXT_DEPENDENT_TEST_FILES
                test ∈ DISABLED_TEST_FILES && continue
                if test == "test_base_checks.jl" &&
                   !(@isdefined(case) && @isdefined(simulation_years) && @isdefined(test_system_dir))
                    @info "Skipping $test: requires `case`, `simulation_years`, `test_system_dir` to be defined by the caller before include(\"test/runtests.jl\")"
                    continue
                end
                if test == "test_system_inputs.jl" &&
                   !(@isdefined(projectdata) && @isdefined(investor_dir))
                    @info "Skipping $test: requires `projectdata`, `investor_dir` to be defined by the caller before include(\"test/runtests.jl\")"
                    continue
                end
                print(splitext(test)[1], ": ")
                include(test)
                println()
            end
        end

        @test length(IS.get_log_events(multi_logger.tracker, Logging.Error)) == 0

        @info IS.report_log_summary(multi_logger)
    end
end

logger = global_logger()

try
    run_tests()
finally
    # Guarantee that the global logger is reset.
    global_logger(logger)
    nothing
end
