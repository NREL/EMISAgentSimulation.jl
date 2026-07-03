# Single place for the package imports shared across test/*.jl files.
#
# Each test file starts with:
#   include(joinpath(@__DIR__, "includes.jl"))
# so it stays runnable standalone from the REPL (not just via runtests.jl),
# without every file repeating its own `using` list.

using Test
using Logging
using InteractiveUtils

using PowerSimulations
using PowerSystems
using StorageSystemsSimulations
using HydroPowerSimulations
using InfrastructureSystems
using JuMP
using Xpress
using HDF5
using AxisArrays
using DataFrames
using CSV

using EMISAgentSimulation

# Not `const`: this file is `include`d once per test file, and redefining a
# `const` at top level warns on every re-inclusion.
PSI = PowerSimulations
PSY = PowerSystems
SSI = StorageSystemsSimulations
HSI = HydroPowerSimulations
IS = InfrastructureSystems
