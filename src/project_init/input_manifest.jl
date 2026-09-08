const PROJECT_INIT_MANIFEST_COLUMNS = [:key, :scope, :path, :required, :source]

"""A minimal declarative manifest for the project initializer."""
function build_project_init_manifest()
    DataFrames.DataFrame(
        key = [
            "project_spec",
            "zones",
            "technologies",
            "scenarios",
            "system_config",
            "projectoptions",
            "projectexisting",
            "outages",
            "devices_to_remove",
        ],
        scope = [
            "system",
            "system",
            "system",
            "system",
            "system",
            "investor",
            "investor",
            "system",
            "system",
        ],
        path = [
            "project_spec.csv",
            "zones.csv",
            "technologies.csv",
            "scenarios.csv",
            "system_config.csv",
            "projectoptions.csv",
            "projectexisting.csv",
            "outages.csv",
            "devices_to_remove.csv",
        ],
        required = [
            true,
            true,
            true,
            true,
            true,
            true,
            false,
            false,
            false,
        ],
        source = [
            "project_spec",
            "project_spec",
            "project_spec",
            "project_spec",
            "project_spec",
            "project_spec",
            "project_spec",
            "project_spec",
            "project_spec",
        ],
    )
end

function validate_project_init_manifest(manifest::DataFrames.DataFrame)
    required_columns = Set(string.(PROJECT_INIT_MANIFEST_COLUMNS))
    missing = setdiff(required_columns, Set(string.(names(manifest))))
    isempty(missing) || error(
        "Project init manifest is missing columns: $(join(sort!(collect(missing)), ", "))"
    )
    return manifest
end
