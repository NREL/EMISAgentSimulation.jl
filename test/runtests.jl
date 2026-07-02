include("includes.jl")

# Self-contained tests — run without any external setup.
include("test_PSI.jl")
include("test_PSI_build_solve.jl")
include("test_hdf5_serialization.jl")

# Context-dependent tests — require the caller to have defined:
#   case, simulation_years, test_system_dir  (for test_base_checks.jl)
#   projectdata, investor_dir                (for test_system_inputs.jl)
include("test_base_checks.jl")
include("test_system_inputs.jl")
