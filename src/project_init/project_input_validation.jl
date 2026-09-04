const DEFAULT_PSY_CLASSIFICATION_MAPPING = joinpath(
    @__DIR__, "..", "..", "config", "psy5", "psy_classification_mapping.csv"
)
const DEFAULT_PROJECT_DEFAULTS = joinpath(
    @__DIR__, "..", "..", "config", "project_defaults.csv"
)

const PSY_CLASSIFICATION_COLUMNS = [
    :component_type,
    :prime_mover,
    :fuel,
    :unit_type,
]

const PROJECT_OPTION_INPUT_COLUMNS = [
    :Investor,
    :GEN_UID,
    Symbol("Unit Type"),
    :Size,
    :Zone,
    Symbol("Bus ID"),
]

const PROJECT_EXISTING_INPUT_COLUMNS = [:Investor, :GEN_UID]

function _require_columns(df::DataFrames.DataFrame, columns, path::AbstractString)
    missing_columns = setdiff(columns, Symbol.(names(df)))
    isempty(missing_columns) || error(
        "$(path) is missing required columns: $(join(string.(missing_columns), ", "))"
    )
    return nothing
end

function _require_nonempty_column(df::DataFrames.DataFrame, column::Symbol, path::AbstractString)
    invalid_rows = findall(row -> ismissing(row[column]) || isempty(strip(string(row[column]))), eachrow(df))
    isempty(invalid_rows) || error(
        "$(path) has empty $(column) values at rows $(join(invalid_rows, ", "))"
    )
    return nothing
end

"""
    load_psy_classification_mapping(path=DEFAULT_PSY_CLASSIFICATION_MAPPING)

Load the user-maintained PSY component classification mapping. The mapping is
deliberately not version-tagged; users update it when their installed PSY types
or enum values change.
"""
function load_psy_classification_mapping(
    path::AbstractString=DEFAULT_PSY_CLASSIFICATION_MAPPING,
)
    isfile(path) || error("PSY classification mapping does not exist: $(path)")
    mapping = DataFrames.DataFrame(CSV.File(path))
    _require_columns(mapping, PSY_CLASSIFICATION_COLUMNS, path)

    for column in (:component_type, :unit_type)
        _require_nonempty_column(mapping, column, path)
    end

    keys = Tuple.(eachrow(mapping[:, PSY_CLASSIFICATION_COLUMNS[1:3]]))
    length(unique(keys)) == length(keys) || error(
        "$(path) contains duplicate component_type/prime_mover/fuel mappings"
    )
    return mapping
end

"""
    validate_psy_classification_mapping(mapping, technologies)

Validate mapping targets against the phase 1 technology configuration.
"""
function validate_psy_classification_mapping(
    mapping::DataFrames.DataFrame,
    technologies::DataFrames.DataFrame,
)
    _require_columns(mapping, PSY_CLASSIFICATION_COLUMNS, "PSY classification mapping")
    _require_columns(technologies, [:unit_type], "technologies.csv")
    valid_unit_types = Set(string.(technologies.unit_type))
    invalid = sort!(unique([
        unit_type for unit_type in string.(mapping.unit_type) if
        !in(unit_type, valid_unit_types)
    ]))
    isempty(invalid) || error(
        "PSY classification mapping references unknown unit_type values: $(join(invalid, ", "))"
    )
    return true
end

function load_project_defaults(path::AbstractString=DEFAULT_PROJECT_DEFAULTS)
    isfile(path) || error("Project defaults do not exist: $(path)")
    defaults = DataFrames.DataFrame(CSV.File(path))
    _require_columns(defaults, [:unit_type, :field, :value], path)
    for column in [:unit_type, :field]
        _require_nonempty_column(defaults, column, path)
    end
    keys = Tuple.(eachrow(defaults[:, [:unit_type, :field]]))
    length(unique(keys)) == length(keys) || error(
        "$(path) contains duplicate unit_type/field defaults"
    )
    return defaults
end

"""Validate one of the minimal user-facing project input templates."""
function validate_project_input_template(
    path::AbstractString;
    kind::Symbol,
    investors::Union{Nothing, AbstractVector}=nothing,
)
    isfile(path) || error("Project input template does not exist: $(path)")
    input = DataFrames.DataFrame(CSV.File(path; stringtype=String))
    required = kind === :options ? PROJECT_OPTION_INPUT_COLUMNS :
               kind === :existing ? PROJECT_EXISTING_INPUT_COLUMNS :
               error("kind must be :options or :existing")
    _require_columns(input, required, path)

    for column in required[1:2]
        _require_nonempty_column(input, column, path)
    end
    if kind === :options
        _require_nonempty_column(input, Symbol("Unit Type"), path)
        _require_nonempty_column(input, :Size, path)
        for row in eachrow(input)
            has_zone = !(ismissing(row.Zone) || isempty(strip(string(row.Zone))))
            has_bus = !(ismissing(row[Symbol("Bus ID")]) || isempty(strip(string(row[Symbol("Bus ID")]))))
            has_zone || has_bus || error(
                "$(path) requires Zone or Bus ID for GEN_UID $(row.GEN_UID)"
            )
        end
    end

    if investors !== nothing
        valid_investors = Set(string.(investors))
        invalid = sort!(unique([
            investor for investor in string.(input.Investor) if
            !in(investor, valid_investors)
        ]))
        isempty(invalid) || error(
            "$(path) references investors not listed in project_spec.csv: $(join(invalid, ", "))"
        )
    end
    return input
end