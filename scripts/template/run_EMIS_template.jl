using Xpress
using EMISAgentSimulation
using Dates
using Distributed
using Logging
using JLD2
using TimerOutputs

const EAS = EMISAgentSimulation
const PROJECT_ROOT = raw"{{PROJECT_ROOT}}"
const BASE_DIR = raw"{{BASE_DIR}}"
const TEST_SYSTEM_DIR = raw"{{TEST_SYSTEM_DIR}}"
const RUNS_DIR = raw"{{RUNS_DIR}}"
const SCRATCH_DIR = raw"{{SCRATCH_DIR}}"
const OUTAGE_FILEPATH = raw"{{OUTAGE_FILEPATH}}"
const SYSTEM_FILEPATH = raw"{{SYSTEM_FILEPATH}}"
const CASE_NAME = "{{CASE_NAME}}"

run_checkpoint = "--run_checkpoint" in ARGS
checkpoint_index = findfirst(==( "--run_checkpoint"), ARGS)
run_name_index = findfirst(==( "--run_name"), ARGS)
run_name = run_name_index === nothing ?
    "$(CASE_NAME)_$(Dates.format(Dates.now(), \"yyyymmdd_HHMM\"))" :
    ARGS[run_name_index + 1]
@info "Run name: $(run_name)"
run_checkpoint && checkpoint_index == length(ARGS) &&
    error("Missing simulation year after --run_checkpoint")
hpc = "--hpc" in ARGS
pras_threads_index = findfirst(==( "--pras-threads"), ARGS)
pras_threads = pras_threads_index === nothing ? 16 :
    parse(Int, ARGS[pras_threads_index + 1])

run_dir = joinpath(RUNS_DIR, CASE_NAME)
result_path = joinpath(run_dir, "Results", run_name)
project_timeseries_data_dir = joinpath(PROJECT_ROOT, "timeseries")
timeseries_data_dir = joinpath(result_path, "timeseries_data_files")
mkpath(result_path)
cd(run_dir)
ENV["SIIP_TIME_SERIES_DIRECTORY"] = SCRATCH_DIR
active_project = Base.active_project()
isnothing(active_project) && error("Run this script with the EMIS package project activated")
emis_project_dir = dirname(active_project)

if !run_checkpoint
    @info "Copying time series data from $(project_timeseries_data_dir) to $(timeseries_data_dir)"
    @info "This should be called only on first run, Checkpoints will use existing timeseries data"
    mkdir(timeseries_data_dir)
    cp(project_timeseries_data_dir, timeseries_data_dir; force=true)
end

const Xpress_optimizer = optimizer_with_attributes(
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

simulation_settings = Dict(
    row["SETTING"] => String(row["VALUE"])
    for row in EAS.read_data(joinpath(@__DIR__, "simulation_settings.csv"))
)
markets_included = Dict(
    Symbol(row["MARKET"]) => row["INCLUDED"]
    for row in EAS.read_data(joinpath(@__DIR__, "markets_included.csv"))
)

current_siip_sim = Any[1]
siip_system = Any[]
completed_sim_years = [
    parse(Int, match.captures[1])
    for filename in readdir(result_path)
    for match in [match(r"simulation_data_year_(\d+)\.h5", filename)]
    if !isnothing(match)
]
last_year = isempty(completed_sim_years) ? 0 : maximum(completed_sim_years)

if run_checkpoint
    last_year = parse(Int, ARGS[checkpoint_index + 1])
end

simulation = nothing
if run_checkpoint && last_year > 0
    @info "Restoring simulation checkpoint from year $(last_year)"
    simulation = EAS.load_simulation(result_path, last_year)
    case = EAS.get_case(simulation)
    EAS.set_solver!(case, Xpress_optimizer)
    EAS.load_sienna_systems!(simulation, result_path, last_year, simulation_settings, SCRATCH_DIR)
    current_year = last_year + 1
else
    case = EAS.CaseDefinition(
        CASE_NAME,
        BASE_DIR,
        TEST_SYSTEM_DIR,
        SCRATCH_DIR,
        OUTAGE_FILEPATH,
        simulation_settings["time_series_data_dir"],
        markets_included,
        Xpress_optimizer;
        simulation_settings=simulation_settings,
    )
    simulation = EAS.create_agent_simulation(case; results_dir=result_path)
    current_year = 1
end

EAS.create_parallel_workers(case, hpc)
pras_worker = EAS.create_pras_worker(hpc; n_threads=pras_threads)

@everywhere begin
    using Pkg
    Pkg.activate($emis_project_dir)
    using Xpress
    using EMISAgentSimulation
    const EAS = EMISAgentSimulation
    using Logging
    using Dates
    using JLD2
    using TimerOutputs
    global_logger(ConsoleLogger(stderr, Logging.Info))
end

@everywhere [pras_worker] ENV["SIIP_TIME_SERIES_DIRECTORY"] = $SCRATCH_DIR

EAS.run_agent_simulation(simulation, current_siip_sim, siip_system, current_year)