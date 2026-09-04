using Test
using Logging
using EMISAgentSimulation
using DataFrames
using CSV
using AxisArrays

# Self-contained tests — run without any external setup.
include("test_PSI.jl")
include("test_hdf5_serialization.jl")
include("test_read_cc_scalar.jl")
include("test_update_derating_factor.jl")
include("test_seasonal_capacity_markets.jl")
include("test_system_config.jl")
include("test_project_input_validation.jl")

# Context-dependent tests — require the caller to have defined:
#   case, simulation_years, test_system_dir  (for test_base_checks.jl)
#   projectdata, investor_dir                (for test_system_inputs.jl)
include("test_base_checks.jl")
include("test_system_inputs.jl")
