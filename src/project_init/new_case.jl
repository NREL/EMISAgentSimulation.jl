function new_emis_case(project_dir::AbstractString, case_name::AbstractString; overrides=Dict())
    case_name = strip(String(case_name))
    isempty(case_name) && error("case_name must be non-empty")
    isdir(project_dir) || error("Project directory does not exist: $(project_dir)")

    base_dir = joinpath(project_dir, "EMIS_RTS_Analysis")
    case_dir = joinpath(base_dir, case_name)
    isdir(case_dir) && error("Case directory already exists: $(case_dir)")

    project_template_dir = joinpath(project_dir, "case_templates")
    mkpath(case_dir)
    if isdir(project_template_dir)
        for entry in readdir(project_template_dir)
            src = joinpath(project_template_dir, entry)
            dst = joinpath(case_dir, entry)
            if isdir(src)
                cp(src, dst; force=true)
            else
                cp(src, dst; force=true)
            end
        end
    end

    for (key, value) in overrides
        key = String(key)
        value = String(value)
        if occursin("{{", key)
            continue
        end
    end

    return Dict(
        :project_dir => project_dir,
        :case_name => case_name,
        :case_dir => case_dir,
    )
end
