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
    zone_lookup::Dict{String, String}
    scenario_lookup::Dict{String, Dict{String, Any}}
    has_custom_config::Bool
end

const DEFAULT_SYSTEM_SCENARIO_NAMES = ["scenario_1", "scenario_2", "scenario_3"]

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

function _default_zone_names()
    # Defaults match the ERCOT/RTS-GMLC system (8 zones)
    default_names = String[]
    for i in 1:8
        push!(default_names, "zone_$(i)")
    end
    return default_names
end

function _default_scenario_names()
    return copy(DEFAULT_SYSTEM_SCENARIO_NAMES)
end

function _coerce_string(value)
    if value === missing
        return ""
    end
    return String(value)
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

function _build_zone_lookup(zone_names::Vector{String})
    lookup = Dict{String, String}()
    for (idx, zone_name) in enumerate(zone_names)
        lookup[string(idx)] = zone_name
        lookup[zone_name] = zone_name
    end
    return lookup
end

function _build_scenario_lookup(df::DataFrames.DataFrame)
    lookup = Dict{String, Dict{String, Any}}()
    if nrow(df) == 0
        # ERCOT defaults: scenario_1→baseline, scenario_2→central, scenario_3→ira
        lookup["scenario_1"] = Dict("probability" => 1.0, "pcm_label" => "baseline", "weather_year" => DEFAULT_LOAD_YEAR)
        lookup["scenario_2"] = Dict("probability" => 1.0, "pcm_label" => "central", "weather_year" => DEFAULT_LOAD_YEAR)
        lookup["scenario_3"] = Dict("probability" => 1.0, "pcm_label" => "ira", "weather_year" => DEFAULT_LOAD_YEAR)
        return lookup
    end

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
        weather_year = tryparse(Int, string(get(row, :weather_year, DEFAULT_LOAD_YEAR)))
        if weather_year === nothing
            weather_year = DEFAULT_LOAD_YEAR
        end
        lookup[scenario_name] = Dict(
            "probability" => probability,
            "pcm_label" => pcm_label,
            "weather_year" => weather_year,
        )
    end

    if isempty(lookup)
        # ERCOT defaults: scenario_1→baseline, scenario_2→central, scenario_3→ira
        lookup["scenario_1"] = Dict("probability" => 1.0, "pcm_label" => "baseline", "weather_year" => DEFAULT_LOAD_YEAR)
        lookup["scenario_2"] = Dict("probability" => 1.0, "pcm_label" => "central", "weather_year" => DEFAULT_LOAD_YEAR)
        lookup["scenario_3"] = Dict("probability" => 1.0, "pcm_label" => "ira", "weather_year" => DEFAULT_LOAD_YEAR)
    end
    return lookup
end

"""
    load_system_config(data_dir::AbstractString)

Load the system configuration for a project. On a generated project, this reads the CSV
files under the project's `system_config/` directory. If the directory is absent, a
legacy RTS-style default is returned so existing workflows continue to work.
"""
function load_system_config(data_dir::AbstractString)
    config_dir = _resolve_system_config_dir(data_dir)
    zones_df = _safe_read_config_csv(joinpath(config_dir, "zones.csv"))
    scenarios_df = _safe_read_config_csv(joinpath(config_dir, "scenarios.csv"))
    technologies_df = _safe_read_config_csv(joinpath(config_dir, "technologies.csv"))
    devices_df = _safe_read_config_csv(joinpath(config_dir, "devices_to_remove.csv"))

    zone_names = _default_zone_names()
    if nrow(zones_df) > 0
        if hasproperty(zones_df, :zone_name)
            zone_names = String.(zones_df[:, :zone_name])
        elseif hasproperty(zones_df, :zone_id)
            zone_names = ["zone_$(x)" for x in zones_df[:, :zone_id]]
        end
    end

    scenario_names = _default_scenario_names()
    scenario_lookup = _build_scenario_lookup(scenarios_df)
    if nrow(scenarios_df) > 0 && hasproperty(scenarios_df, :scenario_name)
        scenario_names = String.(scenarios_df[:, :scenario_name])
    end

    lookup = _build_zone_lookup(zone_names)
    system_name = splitdir(abspath(data_dir))[end]
    if system_name == "." || isempty(system_name)
        system_name = "EMISProject"
    end

    has_custom_config = isdir(config_dir) && (nrow(zones_df) > 0 || nrow(scenarios_df) > 0 || nrow(technologies_df) > 0 || nrow(devices_df) > 0)

    return SystemConfig(
        data_dir,
        system_name,
        zone_names,
        scenario_names,
        DEFAULT_RTS_LOAD,
        lookup,
        scenario_lookup,
        has_custom_config,
    )
end

function get_zone_name(cfg::SystemConfig, zone_id::Union{Int, String})
    if zone_id isa Int
        idx = zone_id
        if 1 <= idx <= length(cfg.zone_names)
            return cfg.zone_names[idx]
        end
        return "zone_$(idx)"
    end
    key = String(zone_id)
    return get(cfg.zone_lookup, key, key)
end

function get_default_scenario(cfg::SystemConfig)
    if isempty(cfg.scenario_names)
        return "scenario_1"
    end
    return cfg.scenario_names[1]
end

function get_scenario_pcm_label(cfg::SystemConfig, scenario_name::String)
    scenario = get(cfg.scenario_lookup, scenario_name, Dict("pcm_label" => "baseline"))
    return String(scenario["pcm_label"])
end

function get_scenario_probability(cfg::SystemConfig, scenario_name::String)
    scenario = get(cfg.scenario_lookup, scenario_name, Dict("probability" => 1.0))
    return Float64(scenario["probability"])
end
