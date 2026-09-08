using Test
using EMISAgentSimulation
using PowerSystems
using DataFrames
using CSV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const TECHNOLOGIES_FILE = joinpath(
    PROJECT_ROOT, "config", "ercot_est", "system_config", "technologies.csv"
)
const MAPPING_FILE = joinpath(
    PROJECT_ROOT, "config", "psy5", "psy_classification_mapping.csv"
)
const DEFAULTS_FILE = joinpath(PROJECT_ROOT, "config", "project_defaults.csv")
const OPTIONS_TEMPLATE = joinpath(
    PROJECT_ROOT, "config", "project_templates", "project_spec", "projectoptions.csv"
)
const EXISTING_TEMPLATE = joinpath(
    PROJECT_ROOT, "config", "project_templates", "project_spec", "projectexisting.csv"
)

@testset "Project input validation" begin
    mapping = load_psy_classification_mapping(MAPPING_FILE)
    technologies = DataFrame(CSV.File(TECHNOLOGIES_FILE))

    @test validate_psy_classification_mapping(mapping, technologies)
    @test nrow(load_project_defaults(DEFAULTS_FILE)) > 0
    @test nrow(validate_project_input_template(
        OPTIONS_TEMPLATE;
        kind=:options,
        investors=["investor_1"],
    )) == 1
    @test nrow(validate_project_input_template(
        EXISTING_TEMPLATE;
        kind=:existing,
        investors=["investor_1"],
    )) == 1

    invalid_mapping = copy(mapping)
    invalid_mapping.unit_type[1] = "UNKNOWN"
    @test_throws ErrorException validate_psy_classification_mapping(
        invalid_mapping, technologies
    )

    invalid_options = DataFrame(
        [["investor_1"], ["new_WT_2"], ["WT"], ["L"], [missing], [missing]],
        [:Investor, :GEN_UID, Symbol("Unit Type"), :Size, :Zone, Symbol("Bus ID")],
    )
    invalid_options_file = joinpath(mktempdir(), "projectoptions.csv")
    CSV.write(invalid_options_file, invalid_options)
    @test_throws ErrorException validate_project_input_template(
        invalid_options_file;
        kind=:options,
        investors=["investor_1"],
    )

    empty_system = PowerSystems.System(100.0; runchecks=false)
    @test_throws ErrorException extract_zones(empty_system)
    @test isempty(extract_branches(empty_system).ac)
    @test isempty(extract_branches(empty_system).dc)
    @test isempty(extract_reserves(empty_system))
    reserve_defaults = DataFrame(
        Symbol.(EMISAgentSimulation.RESERVE_COLUMNS) .=> [
            ["Curve"], [60], [12.0], ["(1)"], ["(Generator)"], [""], ["Up"]
        ],
    )
    demand_curve = PowerSystems.ReserveDemandCurve{PowerSystems.ReserveUp}(nothing)
    @test ismissing(EMISAgentSimulation._reserve_requirement(
        demand_curve,
        nothing,
        "Curve",
    ))
    @test EMISAgentSimulation._reserve_requirement(
        demand_curve,
        reserve_defaults,
        "Curve",
    ) == 12.0
    ownership = DataFrame(Investor=String[], GEN_UID=String[])
    technologies = DataFrame(CSV.File(TECHNOLOGIES_FILE))
    @test nrow(extract_fleet(
        empty_system,
        mapping,
        technologies,
        ownership;
        defaults=load_project_defaults(DEFAULTS_FILE),
    ).projects) == 0
    output_dir = mktempdir()
    @test_throws ErrorException write_system_inputs(
        empty_system,
        output_dir;
        mapping=mapping,
        technologies=technologies,
    )

    spec_dir = joinpath(PROJECT_ROOT, "config", "project_templates", "project_spec")
    init_dir = joinpath(mktempdir(), "project_init_test")
    result = initialize_emis_project(spec_dir; output_dir=init_dir, reference_case_dir=nothing)
    @test isdir(joinpath(init_dir, "EMIS_RTS_Analysis", "Heterogeneous", "system_config"))
    @test isdir(joinpath(init_dir, "EMIS_RTS_Analysis", "Heterogeneous", "markets_data"))
    @test isdir(joinpath(init_dir, "EMIS_RTS_Analysis", "Heterogeneous", "investors", "investor_1", "markets_data"))
    @test haskey(result, :base_dir)
    @test isdir(joinpath(init_dir, "case_templates"))

    stale_root = mktempdir()
    stale_base = joinpath(stale_root, "EMIS_RTS_Analysis")
    mkpath(joinpath(stale_base, "case_1"))
    @test_throws ErrorException initialize_emis_project(
        spec_dir;
        output_dir=stale_root,
        reference_case_dir=nothing,
    )
end