using Test
using EMISAgentSimulation
using CSV
using DataFrames

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))

@testset "Project init validation" begin
    spec_dir = joinpath(PROJECT_ROOT, "config", "project_templates", "project_spec")
    init_dir = joinpath(mktempdir(), "project_init_validation")
    initialize_emis_project(spec_dir; output_dir=init_dir, reference_case_dir=nothing)

    @test validate_emis_project(init_dir) == init_dir
end
