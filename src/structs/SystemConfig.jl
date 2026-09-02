# Minimal system-level configuration for generic project initialization.
# This keeps RTS-specific defaults working for legacy cases while allowing a generated
# project folder to supply its own zones/scenario metadata.

"""
    SystemConfig

Container for the system-level config data that is presently hard-coded in the RTS/
ERCOT workflow. The config is read from <data_dir>/system_config/*.csv when present,
otherwise the legacy defaults are used.
"""
struct SystemConfig
    data_dir::String
    system_name::String
    zone_names::Vector{String}
    scenario_names::Vector{String}
    default_rts_load::Float64
    base_load_year::Int
    zone_lookup::Dict{String, String}
    scenario_lookup::Dict{String, Dict{String, Any}}
    technology_lookup::Dict{String, Dict{String, Any}}
    devices_to_remove::Dict{String, Vector{String}}
    zone_column_start::Int
    has_custom_config::Bool
end

function _safe_read_config_csv(path::String)
    if isfile(path)
        return read_data(path)
    end
    return DataFrames.DataFrame()
end

function _resolve_system_config_dir(data_dir::AbstractString)
    candidates = String[]
    push!(candidates, joinpath(data_dir, "system_config"))
    push!(candidates, joinpath(data_dir, "Heterogeneous", "system_config"))
    push!(candidates, joinpath(data_dir, "Homogeneous", "system_config"))

    for candidate in candidates
        if isdir(candidate)
            return candidate
        end
    end

    return joinpath(data_dir, "system_config")
end

function _coerce_string(value)
    if value === missing
        return ""
    end
    return string(value)
end

function _parse_bool(value)
    if value === missing
        return false
    end
    if value isa Bool
        return value
    end
    return lowercase(String(value)) in ("true", "t", "1", "yes", "y")
end

function _build_zone_lookup(zone_names::Vector{String}, zones_df::DataFrames.DataFrame)
    lookup = Dict{String, String}()
    for (idx, zone_name) in enumerate(zone_names)
        lookup[string(idx)] = zone_name
        lookup[zone_name] = zone_name
    end
    if hasproperty(zones_df, :zone_name)
        for row in eachrow(zones_df)
            zone_name = _coerce_string(row[:zone_name])
            for alias_column in (:zone_id, :area_name, :load_column)
                if hasproperty(zones_df, alias_column)
                    alias = _coerce_string(row[alias_column])
                    isempty(alias) || (lookup[alias] = zone_name)
                end
            end
        end
    end
    return lookup
end

function _build_technology_lookup(df::DataFrames.DataFrame)
    lookup = Dict{String, Dict{String, Any}}()

    for row in eachrow(df)
        unit_type = _coerce_string(row[:unit_type])
        isempty(unit_type) && continue
        lookup[unit_type] = Dict(
            "class" => _coerce_string(get(row, :class, "thermal")),
            "category" => _coerce_string(get(row, :category, unit_type)),
            "mopr_exempt" => _parse_bool(get(row, :mopr_exempt, false)),
        )
    end
    return lookup
end

function _build_devices_to_remove(df::DataFrames.DataFrame)
    devices = Dict{String, Vector{String}}()
    if nrow(df) == 0
        return devices
    end

    if hasproperty(df, :device_type) && hasproperty(df, :device_name)
        for row in eachrow(df)
            device_type = _coerce_string(row[:device_type])
            device_name = _coerce_string(row[:device_name])
            if isempty(device_type) || isempty(device_name)
                continue
            end
            push!(get!(devices, device_type, String[]), device_name)
        end
    elseif hasproperty(df, :component_type) && hasproperty(df, :device_names)
        for row in eachrow(df)
            device_type = _coerce_string(row[:component_type])
            device_names = split(_coerce_string(row[:device_names]), ';'; keepempty = false)
            if isempty(device_type) || isempty(device_names)
                continue
            end
            devices[device_type] = strip.(device_names)
        end
    end
    return devices
end

function _get_config_setting(df::DataFrames.DataFrame, key::String, default)
    if nrow(df) == 0
        return default
    end
    key_column = hasproperty(df, :key) ? :key : (hasproperty(df, :SETTING) ? :SETTING : nothing)
    value_column = hasproperty(df, :value) ? :value : (hasproperty(df, :VALUE) ? :VALUE : nothing)
    if isnothing(key_column) || isnothing(value_column)
        return default
    end
    row_index = findfirst(==("$(key)"), string.(df[:, key_column]))
    return isnothing(row_index) ? default : df[row_index, value_column]
end

function _build_scenario_lookup(df::DataFrames.DataFrame)
    lookup = Dict{String, Dict{String, Any}}()
    for row in eachrow(df)
        name = get(row, :scenario_name, get(row, :Scenario, ""))
        if name === missing || String(name) == ""
            continue
        end
        scenario_name = String(name)
        probability = tryparse(Float64, string(get(row, :probability, 1.0)))
        if probability === nothing
            probability = 1.0
        end
        pcm_label = String(get(row, :pcm_label, "baseline"))
        weather_year = tryparse(Int, string(get(row, :weather_year, "")))
        if weather_year === nothing
            error("Invalid weather_year for scenario $(scenario_name)")
        end
        lookup[scenario_name] = Dict(
            "probability" => probability,
            "pcm_label" => pcm_label,
            "weather_year" => weather_year,
        )
    end

    return lookup
end

"""
    load_system_config(data_dir::AbstractString)

Load the system configuration from the project's required `system_config/` directory.
The repository's `config/legacy_ercot/system_config/` bundle can be copied into legacy
case data directories to preserve ERCOT behavior.
"""
function load_system_config(data_dir::AbstractString)
    config_dir = _resolve_system_config_dir(data_dir)
    required_files = ["system_config.csv", "zones.csv", "scenarios.csv", "technologies.csv", "devices_to_remove.csv"]
    missing_files = filter(file -> !isfile(joinpath(config_dir, file)), required_files)
    isempty(missing_files) || error("Missing system configuration in $(config_dir): $(join(missing_files, ", ")). Copy config/legacy_ercot/system_config for legacy ERCOT cases.")
    settings_df = _safe_read_config_csv(joinpath(config_dir, "system_config.csv"))
    zones_df = _safe_read_config_csv(joinpath(config_dir, "zones.csv"))
    scenarios_df = _safe_read_config_csv(joinpath(config_dir, "scenarios.csv"))
    technologies_df = _safe_read_config_csv(joinpath(config_dir, "technologies.csv"))
    devices_df = _safe_read_config_csv(joinpath(config_dir, "devices_to_remove.csv"))

    hasproperty(zones_df, :zone_name) || error("zones.csv must include a zone_name column")
    zone_names = String.(zones_df[:, :zone_name])
    isempty(zone_names) && error("zones.csv must define at least one zone")

    scenario_lookup = _build_scenario_lookup(scenarios_df)
    hasproperty(scenarios_df, :scenario_name) || error("scenarios.csv must include a scenario_name column")
    scenario_names = String.(scenarios_df[:, :scenario_name])
    isempty(scenario_names) && error("scenarios.csv must define at least one scenario")

    lookup = _build_zone_lookup(zone_names, zones_df)
    technology_lookup = _build_technology_lookup(technologies_df)
    devices_to_remove = _build_devices_to_remove(devices_df)
    default_rts_load = something(
        tryparse(Float64, string(_get_config_setting(settings_df, "default_rts_load", DEFAULT_RTS_LOAD))),
        DEFAULT_RTS_LOAD,
    )
    base_load_year_value = _get_config_setting(settings_df, "base_load_year", DEFAULT_LOAD_YEAR)
    base_load_year = base_load_year_value isa Number ?
        Int(base_load_year_value) :
        something(tryparse(Int, string(base_load_year_value)), DEFAULT_LOAD_YEAR)
    zone_column_start_value = _get_config_setting(settings_df, "zone_column_start", 5)
    zone_column_start = if zone_column_start_value isa Number
        Int(zone_column_start_value)
    else
        something(tryparse(Int, string(zone_column_start_value)), 5)
    end
    zone_column_start > 0 || (zone_column_start = 5)
    system_name = splitdir(abspath(data_dir))[end]
    if system_name == "." || isempty(system_name)
        system_name = "EMISProject"
    end

    has_custom_config = isdir(config_dir) && (nrow(settings_df) > 0 || nrow(zones_df) > 0 || nrow(scenarios_df) > 0 || nrow(technologies_df) > 0 || nrow(devices_df) > 0)

    return SystemConfig(
        data_dir,
        system_name,
        zone_names,
        scenario_names,
        default_rts_load,
        base_load_year,
        lookup,
        scenario_lookup,
        technology_lookup,
        devices_to_remove,
        zone_column_start,
        has_custom_config,
    )
end

function get_zone_name(cfg::SystemConfig, zone_id::Union{Int, String})
    if zone_id isa Int
        idx = zone_id
        if 1 <= idx <= length(cfg.zone_names)
            return cfg.zone_names[idx]
        end
        error("Unknown zone index: $(idx)")
    end
    key = String(zone_id)
    haskey(cfg.zone_lookup, key) || error("Unknown zone identifier: $(key)")
    return cfg.zone_lookup[key]
end

function get_default_scenario(cfg::SystemConfig)
    return cfg.scenario_names[1]
end

function get_scenario_pcm_label(cfg::SystemConfig, scenario_name::String)
    scenario = get(cfg.scenario_lookup, scenario_name, nothing)
    isnothing(scenario) && error("Unknown scenario: $(scenario_name)")
    return String(scenario["pcm_label"])
end

function get_scenario_probability(cfg::SystemConfig, scenario_name::String)
    scenario = get(cfg.scenario_lookup, scenario_name, nothing)
    isnothing(scenario) && error("Unknown scenario: $(scenario_name)")
    return Float64(scenario["probability"])
end

function get_technology_class(cfg::SystemConfig, unit_type::AbstractString)
    technology = get(cfg.technology_lookup, String(unit_type), nothing)
    isnothing(technology) && error("Unknown technology unit type: $(unit_type)")
    return String(technology["class"])
end

function get_technology_category(cfg::SystemConfig, unit_type::AbstractString)
    technology = get(cfg.technology_lookup, String(unit_type), nothing)
    isnothing(technology) && error("Unknown technology unit type: $(unit_type)")
    return String(technology["category"])
end

function is_mopr_exempt(cfg::SystemConfig, unit_type::AbstractString)
    technology = get(cfg.technology_lookup, String(unit_type), nothing)
    isnothing(technology) && error("Unknown technology unit type: $(unit_type)")
    return Bool(technology["mopr_exempt"])
end
