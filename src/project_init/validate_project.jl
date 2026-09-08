function _project_csv_file(project_root::AbstractString, filename::AbstractString)
    path = joinpath(project_root, filename)
    isfile(path) || error("Missing project input: $(path)")
    return path
end

function _read_csv_df(path::AbstractString)
    isfile(path) || error("CSV missing: $(path)")
    return DataFrames.DataFrame(CSV.File(path; stringtype=String))
end

function _as_string(value)
    return strip(string(value))
end

function _must_have_column(df::DataFrames.DataFrame, column::Symbol, path::AbstractString)
    column in Symbol.(names(df)) || error("$(path) is missing required column $(column)")
    return nothing
end

function _base_project_paths(project_root::AbstractString)
    metadata_path = joinpath(project_root, "project_metadata.csv")
    if isfile(metadata_path)
        metadata = Dict{String, String}()
        for row in CSV.File(metadata_path; stringtype=String)
            metadata[strip(string(row.key))] = strip(string(row.value))
        end
        base_dir = get(metadata, "base_dir", project_root)
        heterogeneity = get(metadata, "heterogeneity", "Heterogeneous")
        config_root = joinpath(base_dir, heterogeneity, "system_config")
        investor_root = joinpath(base_dir, heterogeneity, "investors")
        ts_root = joinpath(project_root, "timeseries")
        return config_root, investor_root, ts_root
    end
    return project_root, joinpath(project_root, "investors"), joinpath(project_root, "timeseries")
end

function _project_zone_names(config_root::AbstractString)
    zones_path = joinpath(config_root, "zones.csv")
    isfile(zones_path) || error("Missing zones.csv in $(config_root)")
    zones = _read_csv_df(zones_path)
    _must_have_column(zones, :zone_name, zones_path)
    return Set(_as_string.(zones.zone_name))
end

function _project_scenario_names(config_root::AbstractString)
    scenarios_path = joinpath(config_root, "scenarios.csv")
    isfile(scenarios_path) || error("Missing scenarios.csv in $(config_root)")
    scenarios = _read_csv_df(scenarios_path)
    _must_have_column(scenarios, :scenario_name, scenarios_path)
    return Set(_as_string.(scenarios.scenario_name))
end

function _project_technology_names(config_root::AbstractString)
    tech_path = joinpath(config_root, "technologies.csv")
    isfile(tech_path) || error("Missing technologies.csv in $(config_root)")
    techs = _read_csv_df(tech_path)
    _must_have_column(techs, :unit_type, tech_path)
    return Set(_as_string.(techs.unit_type))
end

function _validate_probability_sum(config_root::AbstractString)
    scenarios_path = joinpath(config_root, "scenarios.csv")
    scenarios = _read_csv_df(scenarios_path)
    _must_have_column(scenarios, :probability, scenarios_path)
    prob_total = 0.0
    for row in eachrow(scenarios)
        val = _as_string(row.probability)
        isempty(val) && error("Scenario probability is empty in $(scenarios_path)")
        prob_total += parse(Float64, val)
    end
    abs(prob_total - 1.0) < 1e-8 || error("Scenario probabilities sum to $(prob_total), expected 1.0")
    return nothing
end

function _validate_load_zone_columns(config_root::AbstractString, ts_root::AbstractString)
    zones_path = joinpath(config_root, "zones.csv")
    zones = _read_csv_df(zones_path)
    _must_have_column(zones, :zone_name, zones_path)
    zone_names = Set(_as_string.(zones.zone_name))
    isdir(ts_root) || return nothing

    valid_scenarios = _project_scenario_names(config_root)
    for scenario in readdir(ts_root)
        scenario_path = joinpath(ts_root, scenario)
        isdir(scenario_path) || continue
        if !(scenario in valid_scenarios)
            error("Timeseries scenario folder $(scenario) is not present in scenarios.csv")
        end
        for sim_year in readdir(scenario_path)
            year_path = joinpath(scenario_path, sim_year)
            isdir(year_path) || continue
            load_dir = joinpath(year_path, "Load")
            isdir(load_dir) || continue
            load_file = joinpath(load_dir, "DAY_AHEAD_regional_Load.csv")
            isfile(load_file) || continue
            load_df = _read_csv_df(load_file)
            if :Year in names(load_df) && :Month in names(load_df) && :Day in names(load_df) && :Period in names(load_df)
                for zone in zone_names
                    if !(String(zone) in String.(names(load_df)))
                        error("Load timeseries for $(scenario)/$(sim_year) is missing zone column $(zone)")
                    end
                end
            end
        end
    end
    return nothing
end

function _validate_option_techs(config_root::AbstractString, investor_root::AbstractString)
    valid_techs = _project_technology_names(config_root)
    isdir(investor_root) || return nothing
    for investor in readdir(investor_root)
        investor_path = joinpath(investor_root, investor)
        isdir(investor_path) || continue
        options_path = joinpath(investor_path, "projectoptions.csv")
        if !isfile(options_path)
            continue
        end
        options = _read_csv_df(options_path)
        if Symbol("Unit Type") in names(options)
            for row in eachrow(options)
                tech = _as_string(row[Symbol("Unit Type")])
                isempty(tech) && continue
                if !(tech in valid_techs)
                    error("Investor $(investor) projectoptions.csv references unknown technology $(tech)")
                end
            end
        end
    end
    return nothing
end

function _validate_psy_generator_names(config_root::AbstractString, ts_root::AbstractString)
    isdir(ts_root) || return nothing
    zone_names = _project_zone_names(config_root)
    for scenario in readdir(ts_root)
        scenario_path = joinpath(ts_root, scenario)
        isdir(scenario_path) || continue
        for sim_year in readdir(scenario_path)
            year_path = joinpath(scenario_path, sim_year)
            isdir(year_path) || continue
            for kind in ("WIND", "PV", "RTPV", "Hydro")
                kind_dir = joinpath(year_path, kind)
                isdir(kind_dir) || continue
                for file in readdir(kind_dir)
                    path = joinpath(kind_dir, file)
                    isfile(path) || continue
                    if endswith(lowercase(file), ".csv")
                        df = _read_csv_df(path)
                        for name in String.(names(df))
                            if !(name in zone_names || name in ("Year", "Month", "Day", "Period"))
                                continue
                            end
                        end
                    end
                end
            end
        end
    end
    return nothing
end

function validate_emis_project(project_root::AbstractString)
    isdir(project_root) || error("Project root does not exist: $(project_root)")

    config_root, investor_root, ts_root = _base_project_paths(project_root)

    for filename in ("zones.csv", "technologies.csv", "scenarios.csv")
        _project_csv_file(config_root, filename)
    end
    if isfile(joinpath(project_root, "system_config.csv"))
        _project_csv_file(project_root, "system_config.csv")
    elseif isfile(joinpath(config_root, "system_config.csv"))
        _project_csv_file(config_root, "system_config.csv")
    end

    _validate_probability_sum(config_root)
    _validate_load_zone_columns(config_root, ts_root)
    _validate_option_techs(config_root, investor_root)
    _validate_psy_generator_names(config_root, ts_root)

    return project_root
end
