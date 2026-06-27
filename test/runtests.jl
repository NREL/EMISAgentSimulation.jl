using Test
using EMISAgentSimulation
using DataFrames
using CSV
using AxisArrays

include("helpers.jl")

@testset "EMISAgentSimulation" begin
    include("test_read_cc_scalar.jl")
    include("test_update_derating_factor.jl")
end
