const PSY_LOADS = Union{PSY.StandardLoad, PSY.PowerLoad}
const PSY_THERMAL_GENERATORS = Union{PSY.ThermalStandard, PSY.ThermalMultiStart}

const DEFAULT_TIME_RESOLUTION = 60
const DEFAULT_HOURS_PER_YEAR = 8760

const DEFAULT_LOAD_YEAR = 2020
const DEFAULT_RTS_LOAD = 75.0 # GW

const BASE_POWER = 100.0 # MW

# PSI constants
const BALANCE_SLACK_COST = 1e6
const SERVICES_SLACK_COST = 1e5

# Simulations constants
const PENALTY_COST = 5000.0

# PRAS constants
const PRAS_N_SAMPLES = 1000
const PRAS_WORKER = Ref{Union{Nothing, Int}}(nothing)

const SIM_START_DATE = Dates.DateTime("2018-01-01T00:00:00")
const SIM_END_DATE = Dates.DateTime("2019-01-01T00:00:00")

# File pointers
const TIMESERIES_DATA_DIR = "/projects/gmlcmarkets/Phase2_EMIS_Analysis/NTP_TimeSeries_Data"
const POINTER_FILE = Dict(
    :NTPS_TS_DATA_DIR => joinpath(TIMESERIES_DATA_DIR, "input_processing"),
   
)


const OBJ_SCALE = 1 # Scale factor to convert objective function values from $ to millions of dollars

# Reading and writing data
const SCHEMA_VERSION = 1