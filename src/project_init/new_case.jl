function _read_project_metadata(project_dir::AbstractString)
    path = joinpath(project_dir, "project_metadata.csv")
    isfile(path) || error("Project metadata does not exist: $(path)")
    metadata = Dict{String, String}()
    for row in CSV.File(path; stringtype=String)
        metadata[strip(string(row.key))] = strip(string(row.value))
    end
    return metadata
end

function _render_case_tree(src_dir::AbstractString, dst_dir::AbstractString, replacements::Dict{String, String})
    mkpath(dst_dir)
    for entry in readdir(src_dir)
        src = joinpath(src_dir, entry)
        dst = joinpath(dst_dir, entry)
        if isdir(src)
            _render_case_tree(src, dst, replacements)
        else
            mkpath(dirname(dst))
            write(dst, _render_template_text(read(src, String), replacements))
        end
    end
    return dst_dir
end

function _apply_case_overrides(settings_path::AbstractString, overrides)
    isempty(overrides) && return
    rows = collect(CSV.File(settings_path; stringtype=String))
    values = Dict(string(row["SETTING"]) => string(row["VALUE"]) for row in rows)
    for (key, value) in overrides
        key = String(key)
        haskey(values, key) || error("Unknown case setting override: $(key)")
        values[key] = String(value)
    end
    open(settings_path, "w") do io
        write(io, "SETTING,VALUE,,Comments\n")
        for row in rows
            key = string(row["SETTING"])
            comments = string(row["Comments"])
            write(io, "$(key),$(values[key]),,$(comments)\n")
        end
    end
end

function new_emis_case(project_dir::AbstractString, case_name::AbstractString; overrides=Dict())
    case_name = strip(String(case_name))
    isempty(case_name) && error("case_name must be non-empty")
    isdir(project_dir) || error("Project directory does not exist: $(project_dir)")

    project_template_dir = joinpath(project_dir, "case_templates")
    isdir(project_template_dir) || error("Project templates do not exist: $(project_template_dir)")
    metadata = _read_project_metadata(project_dir)
    runs_dir = get(metadata, "runs_dir", "")
    isempty(runs_dir) && error("Project metadata does not define runs_dir")
    run_dir = joinpath(runs_dir, case_name)
    isdir(run_dir) && error("Case run directory already exists: $(run_dir)")

    replacements = Dict("CASE_NAME" => case_name)
    _render_case_tree(project_template_dir, run_dir, replacements)
    _apply_case_overrides(joinpath(run_dir, "simulation_settings.csv"), overrides)

    return Dict(
        :project_dir => project_dir,
        :case_name => case_name,
        :case_dir => run_dir,
        :run_dir => run_dir,
        :case_data_dir => joinpath(metadata["base_dir"], case_name),
    )
end
