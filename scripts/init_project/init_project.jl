#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using EMISAgentSimulation

function usage()
    println("""
    Usage:
      julia scripts/init_project/init_project.jl init <spec_dir> --output-dir <dir> [--reference-case <dir>]
      julia scripts/init_project/init_project.jl case <project_dir> <case_name> [--override key=value ...]

    Commands:
      init      Create a reusable EMIS project tree from a project spec directory.
      case      Stamp a new case directory under the project runs tree from the generated templates.
    """)
end

function parse_override_pairs(raw_args)
    overrides = Dict{String, String}()
    i = 1
    while i <= length(raw_args)
        arg = raw_args[i]
        if arg == "--override"
            i >= length(raw_args) && error("Missing value after --override")
            pair = raw_args[i + 1]
            i += 2
        elseif startswith(arg, "--override=")
            pair = replace(arg, "--override=" => ""; count = 1)
            i += 1
        else
            i += 1
            continue
        end

        parts = split(pair, "="; limit = 2)
        length(parts) == 2 || error("Override must use key=value format: $(pair)")
        key = strip(parts[1])
        value = strip(parts[2])
        isempty(key) && error("Override key cannot be empty")
        overrides[key] = value
    end
    return overrides
end

function parse_init_args(raw_args)
    spec_dir = ""
    output_dir = ""
    reference_case = nothing
    i = 1
    while i <= length(raw_args)
        arg = raw_args[i]
        if arg == "--help" || arg == "-h"
            usage(); exit(0)
        elseif arg == "--output-dir" || arg == "-o"
            i >= length(raw_args) && error("Missing value after $(arg)")
            output_dir = raw_args[i + 1]
            i += 2
        elseif arg == "--reference-case"
            i >= length(raw_args) && error("Missing value after $(arg)")
            reference_case = raw_args[i + 1]
            i += 2
        elseif isempty(spec_dir)
            spec_dir = arg
            i += 1
        else
            error("Unknown argument: $(arg)")
        end
    end

    isempty(spec_dir) && error("A project spec directory is required")
    isempty(output_dir) && error("An output directory is required; use --output-dir")
    return (spec_dir, output_dir, reference_case)
end

function parse_case_args(raw_args)
    project_dir = ""
    case_name = ""
    overrides = parse_override_pairs(raw_args)

    i = 1
    while i <= length(raw_args)
        arg = raw_args[i]
        if arg == "--help" || arg == "-h"
            usage(); exit(0)
        elseif arg == "--override"
            i += 2
            continue
        elseif startswith(arg, "--override=")
            i += 1
            continue
        elseif isempty(project_dir)
            project_dir = arg
        elseif isempty(case_name)
            case_name = arg
        else
            error("Unknown argument: $(arg)")
        end
        i += 1
    end

    isempty(project_dir) && error("A project directory is required")
    isempty(case_name) && error("A case name is required")
    return (project_dir, case_name, overrides)
end

function main()
    if isempty(ARGS) || ARGS[1] == "--help" || ARGS[1] == "-h"
        usage()
        return
    end

    command = lowercase(ARGS[1])
    rest = ARGS[2:end]

    if command == "init"
        spec_dir, output_dir, reference_case = parse_init_args(rest)
        result = initialize_emis_project(spec_dir; output_dir = output_dir, reference_case_dir = reference_case)
        println("Initialized EMIS project at $(output_dir)")
        println("  base_dir: $(result[:base_dir])")
        println("  runs_dir: $(result[:runs_dir])")
        println("  test_system_dir: $(result[:test_system_dir])")
        println("  case_templates_dir: $(result[:case_templates_dir])")
        return
    end

    if command == "case"
        project_dir, case_name, overrides = parse_case_args(rest)
        result = new_emis_case(project_dir, case_name; overrides = overrides)
        println("Created case $(case_name) under $(result[:run_dir])")
        println("  case_data_dir: $(result[:case_data_dir])")
        return
    end

    usage()
    error("Unknown command: $(command)")
end

main()
