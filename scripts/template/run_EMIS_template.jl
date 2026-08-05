using Xpress
using Xpress
using EMISAgentSimulation
using AxisArrays
using DataFrames
using CSV
using Statistics
using Dates
using StatsBase
using LinearAlgebra
using Clustering
using PowerSystems
using TimeSeries
using FileIO
using PowerSimulations
using TimeZones
using Distributed
using HydroPowerSimulations
using StorageSystemsSimulations
using SiennaPRASInterface
import JLD2
using TimerOutputs
using HDF5

const PSY = PowerSystems
const TS = TimeSeries
const PSI = PowerSimulations
const HSI = HydroPowerSimulations
const SSI = StorageSystemsSimulations
const EAS = EMISAgentSimulation
const PRAS = SiennaPRASInterface

############################################################
## Update the following paths
const GS_PATH = gs_path # "/projects/gmlcmarkets/Phase2_EMIS_Analysis/GS_AAYAD"
scratch_dir = user_scratch_dir # "/kfs3/scratch/username/"
############################################################


# Parse named CLI arguments: --run_checkpoint <year> and --run_name <name>
run_checkpoint = "--run_checkpoint" in ARGS
_cp_idx = findfirst(==("--run_checkpoint"), ARGS)
_name_idx = findfirst(==("--run_name"), ARGS)

const HPC_ANALYSIS_PATH = joinpath(GS_PATH, "HPC_Analysis_Runs")
const CASE_NAME = basename(@__DIR__)

# Unique run name: auto-generated with timestamp for fresh runs so each run
# gets its own Results subfolder (prevents stale JSON contamination across runs).
# For checkpoint restarts, pass --run_name <name> matching the original run.
const run_name = _name_idx !== nothing ?
    ARGS[_name_idx + 1] :
    "$(CASE_NAME)_$(Dates.format(Dates.now(), "yyyymmdd_HHMM"))"
@info "Run name: $(run_name)"

NTP_TIMESERIES_DATA_DIR = "/projects/gmlcmarkets/Phase2_EMIS_Analysis/NTP_TimeSeries_Data"
@info "NTP Time series data directory: $(NTP_TIMESERIES_DATA_DIR)"

ENV["SIIP_TIME_SERIES_DIRECTORY"] = scratch_dir
outage_dir = "/kfs2/projects/gmlcmarkets/Phase2_EMIS_Analysis/Correlated Outages/"
outage_filename = "ThermalFOR_2011.csv"
outage_filepath = joinpath(outage_dir, outage_filename)

temp_dir = joinpath(HPC_ANALYSIS_PATH, CASE_NAME, "temp_data")
result_path = joinpath(HPC_ANALYSIS_PATH, CASE_NAME, "Results", run_name)
mkpath(result_path)
cd(joinpath(HPC_ANALYSIS_PATH, CASE_NAME))

const Xpress_optimizer = optimizer_with_attributes(Xpress.Optimizer,
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
                                                       
simulation_settings = Dict(r["SETTING"] => String(r["VALUE"]) for r in eachrow(EAS.read_data(joinpath(@__DIR__, "simulation_settings.csv"))))
markets_included = Dict(Symbol(r["MARKET"]) => r["INCLUDED"] for r in eachrow(EAS.read_data(joinpath(@__DIR__, "markets_included.csv"))))
vre_reserves_bool = EAS.parsebool(simulation_settings["vre_reserves"])

if vre_reserves_bool
    base_dir = "../../EMIS_RTS_Analysis_GS"
    test_system_dir = "../../RTS-GMLC_GS"
else
    base_dir = "../EMIS_RTS_Analysis_No_VRE_Reserves"
    test_system_dir = "../../RTS-GMLC_No_VRE_Reserves/RTS-GMLC"
 end
 

EMIS_RTS_timseries_data_dir = vre_reserves_bool ?
    joinpath(GS_PATH, "EMIS_RTS_Analysis_GS", CASE_NAME, "timeseries_data_files") :
    joinpath(GS_PATH, "EMIS_RTS_Analysis", CASE_NAME, "timeseries_data_files")
timeseries_data_dir = joinpath(result_path, "timeseries_data_files")

if !run_checkpoint
    @info "Copying time series data from $(EMIS_RTS_timseries_data_dir) to $(timeseries_data_dir)"
    @info "This should be called only on first run, Checkpoints will use existing timeseries data"
    mkdir(timeseries_data_dir)    
    cp(EMIS_RTS_timseries_data_dir, timeseries_data_dir, force=true)
end

current_sim = Any[1]
siip_system = Any[]
completed_sim_years = [parse(Int, m.captures[1]) for f in readdir(result_path)
            for m in [match(r"simulation_data_year_(\d+)\.h5", f)] if !isnothing(m)]

last_year = isempty(completed_sim_years) ? 0 : maximum(completed_sim_years)
@info "Found completed simulation years from checkpoint: $(completed_sim_years). We can start from year $(last_year + 1)."

if run_checkpoint
    @info "Overriding last year from checkpoint with command line argument --run_checkpoint"
    last_year = parse(Int, ARGS[_cp_idx + 1])
end
@info "Starting checkpoint restoration with last year: $(last_year)"

simulation = nothing
if run_checkpoint && last_year > 0
    @info "Restoring simulation from checkpoint at year $(last_year) from path\n: $(result_path)"
    simulation = EAS.load_simulation(result_path, last_year)
    case = EAS.get_case(simulation)
    EAS.set_solver!(case, Xpress_optimizer)
    EAS.load_sienna_systems!(simulation, result_path, last_year, simulation_settings, scratch_dir)
    current_year = last_year + 1
end

if isnothing(simulation)
    name = CASE_NAME
    case = EAS.CaseDefinition(name,
                          base_dir,
                          test_system_dir,
                          scratch_dir,
                          outage_filepath,
                          NTP_TIMESERIES_DATA_DIR,
                          Xpress_optimizer, #Gurobi_optimizer, #Xpress_optimizer, #HiGHS_optimizer
                          siip_market_clearing = EAS.parsebool(simulation_settings["siip_market_clearing"]),
                          pcm_scenario = simulation_settings["pcm_scenario"],
                          start_year = EAS.parseint(simulation_settings["start_year"]),
                          total_horizon = EAS.parseint(simulation_settings["total_horizon"]),
                          rolling_horizon = EAS.parseint(simulation_settings["rolling_horizon"]),
                          simulation_years = EAS.parseint(simulation_settings["simulation_years"]),
                          rep_period_interval = EAS.parseint(simulation_settings["rep_period_interval"]),
                          num_rep_periods = EAS.parseint(simulation_settings["num_rep_periods"]),
                          avg_block_size = EAS.parseint(simulation_settings["avg_block_size"]),
                          fixed_block_size = EAS.parsebool(simulation_settings["fixed_block_size"]),
                          rep_chronology_checkpoint = EAS.parseint(simulation_settings["rep_chronology_checkpoint"]),
                          da_resolution = EAS.parseint(simulation_settings["da_resolution"]),
                          rt_resolution = EAS.parseint(simulation_settings["rt_resolution"]),
                          rps_target = simulation_settings["rps_target"],
                          markets = markets_included,
                          ordc_curved = EAS.parsebool(simulation_settings["ordc_curved"]),
                          ordc_unavailability_method = simulation_settings["ordc_unavailability_method"],
                          reserve_penalty = simulation_settings["reserve_penalty"],
                          static_capacity_market = EAS.parsebool(simulation_settings["static_capacity_market"]),
                          irm_scalar = EAS.parsefloat(simulation_settings["irm_scalar"]),
                          accreditation_methodology = simulation_settings["accreditation_methodology"],
                          accreditation_metric = simulation_settings["accreditation_metric"],
                          marginal_cc_switch = EAS.parsebool(simulation_settings["marginal_cc_switch"]),
                          derating_scale = EAS.parsefloat(simulation_settings["derating_scale"]),
                          mopr = EAS.parsebool(simulation_settings["mopr"]),
                          battery_cap_mkt = EAS.parsebool(simulation_settings["battery_cap_mkt"]),
                          vre_reserves = vre_reserves_bool,
                          heterogeneity = EAS.parsebool(simulation_settings["heterogeneity"]),
                          forecast_type = simulation_settings["forecast_type"],
                          max_carbon_tax_increase =  EAS.parsefloat(simulation_settings["max_carbon_tax_increase"]),
                          info_symmetry = EAS.parsebool(simulation_settings["info_symmetry"]),
                          belief_update =EAS.parsebool(simulation_settings["belief_update"]),
                          uncertainty = EAS.parsebool(simulation_settings["uncertainty"]),
                          risk_aversion = EAS.parsebool(simulation_settings["risk_aversion"]),
                          parallel_investors = EAS.parsebool(simulation_settings["parallel_investors"]),
                          parallel_scenarios = EAS.parsebool(simulation_settings["parallel_scenarios"]),
                          md_horizon = EAS.parseint(simulation_settings["md_horizon"]),
                          md_interval = EAS.parseint(simulation_settings["md_interval"]),
                          uc_horizon = EAS.parseint(simulation_settings["uc_horizon"]),
                          uc_interval = EAS.parseint(simulation_settings["uc_interval"]),
                          ed_horizon = EAS.parseint(simulation_settings["ed_horizon"]),
                          ed_interval = EAS.parseint(simulation_settings["ed_interval"]),
                          md_market = EAS.parsebool(simulation_settings["multi_day_market"]),
                        )    
    @info "Starting simulation from scratch."
    simulation = EAS.create_agent_simulation(case);
    current_year = 1
end

hpc = false
EAS.create_parallel_workers(case, hpc)
pras_worker = EAS.create_pras_worker(hpc, n_threads = 16)

@everywhere begin
    using Pkg; Pkg.activate(".");
    using Xpress
    using EMISAgentSimulation
    const EAS = EMISAgentSimulation
    using Logging
    using Dates
    using JLD2
    using TimerOutputs
    global_logger(ConsoleLogger(stderr, Logging.Info))
end

@everywhere [pras_worker] ENV["SIIP_TIME_SERIES_DIRECTORY"] = $scratch_dir

EAS.run_agent_simulation(simulation,
                         current_sim,
                         siip_system,
                         current_year)
