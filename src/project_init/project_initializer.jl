function _read_project_spec(spec_dir::AbstractString)
    path = joinpath(spec_dir, "project_spec.csv")
    isfile(path) || error("Project spec does not exist: $(path)")
    spec = DataFrames.DataFrame(CSV.File(path; stringtype=String))
    isempty(spec) && error("Project spec is empty: $(path)")
    if !("key" in names(spec) && "value" in names(spec))
        error("$(path) must contain key/value columns")
    end
    spec_lookup = Dict{String, String}()
    for row in eachrow(spec)
        key = strip(string(row.key))
        isempty(key) && continue
        spec_lookup[key] = strip(string(row.value))
    end
    return spec_lookup
end

function _project_spec_value(spec::Dict{String, String}, key::AbstractString; default::Union{Nothing, AbstractString}=nothing)
    haskey(spec, key) && !isempty(strip(spec[key])) && return spec[key]
    default === nothing && error("Required project_spec key missing: $(key)")
    return default
end

function _split_investors(spec::Dict{String, String})
    values = _project_spec_value(spec, "investors")
    return [strip(part) for part in split(values, ";") if !isempty(strip(part))]
end

function _copy_template_file(src::AbstractString, dst::AbstractString)
    isfile(src) || error("Template file does not exist: $(src)")
    mkpath(dirname(dst))
    cp(src, dst; force=true)
    return dst
end

function _resolve_path(base_dir::AbstractString, path_value::AbstractString)
    isempty(strip(path_value)) && return ""
    if isabspath(path_value)
        return normpath(path_value)
    end
    return normpath(joinpath(base_dir, path_value))
end

function _write_system_config_bundle(base_dir::AbstractString, spec_dir::AbstractString, spec::Dict{String, String})
    source_dir = joinpath(spec_dir, "system_config")
    heterogeneity = _project_spec_value(spec, "heterogeneity")
    target_dir = joinpath(base_dir, heterogeneity, "system_config")
    mkpath(target_dir)

    for filename in ("system_config.csv", "zones.csv", "technologies.csv", "scenarios.csv", "devices_to_remove.csv")
        src_path = joinpath(spec_dir, filename)
        dst_path = joinpath(target_dir, filename)
        if isfile(src_path)
            cp(src_path, dst_path; force=true)
        elseif isfile(joinpath(source_dir, filename))
            cp(joinpath(source_dir, filename), dst_path; force=true)
        end
    end

    if !isfile(joinpath(target_dir, "system_config.csv"))
        maybe_default = joinpath(@__DIR__, "..", "..", "config", "ercot_est", "system_config", "system_config.csv")
        isfile(maybe_default) && cp(maybe_default, joinpath(target_dir, "system_config.csv"); force=true)
    end
    if !isfile(joinpath(target_dir, "zones.csv"))
        maybe_default = joinpath(@__DIR__, "..", "..", "config", "ercot_est", "system_config", "zones.csv")
        isfile(maybe_default) && cp(maybe_default, joinpath(target_dir, "zones.csv"); force=true)
    end
    if !isfile(joinpath(target_dir, "technologies.csv"))
        maybe_default = joinpath(@__DIR__, "..", "..", "config", "ercot_est", "system_config", "technologies.csv")
        isfile(maybe_default) && cp(maybe_default, joinpath(target_dir, "technologies.csv"); force=true)
    end
    if !isfile(joinpath(target_dir, "scenarios.csv"))
        maybe_default = joinpath(@__DIR__, "..", "..", "config", "ercot_est", "system_config", "scenarios.csv")
        isfile(maybe_default) && cp(maybe_default, joinpath(target_dir, "scenarios.csv"); force=true)
    end
    if !isfile(joinpath(target_dir, "devices_to_remove.csv"))
        maybe_default = joinpath(@__DIR__, "..", "..", "config", "ercot_est", "system_config", "devices_to_remove.csv")
        isfile(maybe_default) && cp(maybe_default, joinpath(target_dir, "devices_to_remove.csv"); force=true)
    end

    return target_dir
end

function _template_source_dir()
    return normpath(joinpath(@__DIR__, "..", "..", "scripts", "template"))
end

function _render_template_text(text::AbstractString, replacements::Dict{String, String})
    rendered = String(text)
    for (key, value) in replacements
        rendered = replace(rendered, "{{$(key)}}" => value)
    end
    return rendered
end

function _copy_rendered_template(src::AbstractString, dst::AbstractString, replacements::Dict{String, String})
    mkpath(dirname(dst))
    if isfile(src)
        write(dst, _render_template_text(read(src, String), replacements))
    elseif isdir(src)
        _copy_directory_contents(src, dst)
    else
        error("Template path does not exist: $(src)")
    end
    return dst
end

function _write_case_template(project_root::AbstractString, replacements::Dict{String, String})
    template_dir = joinpath(project_root, "case_templates")
    mkpath(template_dir)
    source_dir = _template_source_dir()
    if isdir(source_dir)
        for entry in readdir(source_dir)
            _copy_rendered_template(
                joinpath(source_dir, entry),
                joinpath(template_dir, entry),
                replacements,
            )
        end
    end
    return template_dir
end

function _copy_case_template_overrides(spec_dir::AbstractString, template_dir::AbstractString, replacements::Dict{String, String})
    for filename in ("simulation_settings.csv", "markets_included.csv", "options.csv")
        src = joinpath(spec_dir, filename)
        if isfile(src)
            _copy_rendered_template(src, joinpath(template_dir, filename), replacements)
        end
    end
    return template_dir
end

function _copy_directory_contents(src_dir::AbstractString, dst_dir::AbstractString)
    isdir(src_dir) || return dst_dir
    mkpath(dst_dir)
    for entry in readdir(src_dir)
        src_path = joinpath(src_dir, entry)
        dst_path = joinpath(dst_dir, entry)
        if isdir(src_path)
            cp(src_path, dst_path; force=true)
        elseif isfile(src_path)
            cp(src_path, dst_path; force=true)
        end
    end
    return dst_dir
end

function _copy_system_file_bundle(system_path::AbstractString, target_dir::AbstractString)
    isfile(system_path) || error("System file does not exist: $(system_path)")
    mkpath(target_dir)
    copied_system_path = joinpath(target_dir, basename(system_path))
    cp(system_path, copied_system_path; force=true)

    root, ext = splitext(system_path)
    if lowercase(ext) == ".json"
        storage_path = root * "_time_series_storage.h5"
        isfile(storage_path) || error(
            "PowerSystems time-series storage file does not exist: $(storage_path)"
        )
        cp(storage_path, joinpath(target_dir, basename(storage_path)); force=true)

        for optional_path in (root * "_metadata.json", root * "_validation_descriptors.json")
            if isfile(optional_path)
                cp(optional_path, joinpath(target_dir, basename(optional_path)); force=true)
            end
        end
    end

    return copied_system_path
end

const PROJECT_OPTION_DOWNSTREAM_COLUMNS = [
    "GEN_UID", "Unit Type", "Fuel", "Category", "Size", "Min Gen pu",
    "Min Up Time Hr", "Min Down Time Hr", "Shut Down Cost", "Ramp Rate pu/Hr",
    "Input Power Rating pu", "Output Power Rating pu", "Min Storgae pu", "Duration Hr",
    "Round Trip Efficiency pu", "Fuel Price \$/MMBTU", "Start Heat Cold MBTU",
    "Output_pct_0", "Output_pct_1", "Output_pct_2", "Output_pct_3", "Output_pct_4",
    "HR_avg_0", "HR_incr_1", "HR_incr_2", "HR_incr_3", "HR_incr_4",
    "CO2_Emissions ton/MMBTU", "Inertia MJ/MW", "FOR", "MTTR Hr", "Lagtime",
    "Capex Years", "Lifetime", "Fixed OM Cost per MW", "Preference Multiplier",
    "Zone", "Bus ID", "Capacity Eligible", "REC Eligible", "Synchronous_Inertia",
]

function _string_column_names(table::DataFrames.DataFrame)
    return String.(names(table))
end

function _is_complete_projectoptions(table::DataFrames.DataFrame)
    return all(column -> column in _string_column_names(table), PROJECT_OPTION_DOWNSTREAM_COLUMNS)
end

function _technology_lookup(technologies::DataFrames.DataFrame)
    lookup = Dict{String, DataFrames.DataFrameRow}()
    for row in eachrow(technologies)
        lookup[string(row.unit_type)] = row
    end
    return lookup
end

function _default_lookup(defaults::DataFrames.DataFrame)
    lookup = Dict{String, Dict{String, String}}()
    for row in eachrow(defaults)
        unit_defaults = get!(lookup, string(row.unit_type), Dict{String, String}())
        unit_defaults[string(row.field)] = string(row.value)
    end
    return lookup
end

function _option_default(defaults, unit_type::AbstractString, field::AbstractString, fallback)
    unit_defaults = get(defaults, string(unit_type), Dict{String, String}())
    value = get(unit_defaults, string(field), fallback)
    value === missing || isempty(strip(string(value))) ? string(fallback) : string(value)
end

function _option_technology_default(technology, field::AbstractString, fallback)
    symbol = Symbol(field)
    if symbol in Symbol.(names(technology))
        value = technology[symbol]
        return value === missing || isempty(strip(string(value))) ? string(fallback) : string(value)
    end
    return string(fallback)
end

function _complete_projectoptions_for_investor(
    options::DataFrames.DataFrame,
    investor::AbstractString,
    technologies::DataFrames.DataFrame,
    defaults::DataFrames.DataFrame,
)
    if _is_complete_projectoptions(options)
        investor_options = if "Investor" in _string_column_names(options)
            filter(row -> string(row.Investor) == investor, options)
        else
            options
        end
        return investor_options[:, PROJECT_OPTION_DOWNSTREAM_COLUMNS]
    end

    _require_columns(options, PROJECT_OPTION_INPUT_COLUMNS, "projectoptions.csv")
    technology_by_unit = _technology_lookup(technologies)
    defaults_by_unit = _default_lookup(defaults)
    rows = DataFrames.DataFrame([String[] for _ in PROJECT_OPTION_DOWNSTREAM_COLUMNS], PROJECT_OPTION_DOWNSTREAM_COLUMNS)
    for input_row in eachrow(options)
        string(input_row.Investor) == investor || continue
        unit_type = string(input_row[Symbol("Unit Type")])
        haskey(technology_by_unit, unit_type) || error("projectoptions.csv references unknown technology $(unit_type)")
        technology = technology_by_unit[unit_type]
        technology_class = string(technology.class)
        capacity_default = _option_default(
            defaults_by_unit,
            unit_type,
            "Capacity Eligible",
            _option_technology_default(technology, "capacity_eligible", "true"),
        )
        rec_default = _option_default(defaults_by_unit, unit_type, "REC Eligible", "false")
        duration_default = _option_default(
            defaults_by_unit,
            unit_type,
            "Duration Hr",
            _option_technology_default(technology, "duration_hr", "n/a"),
        )
        row = Dict{String, String}(
            "GEN_UID" => string(input_row.GEN_UID),
            "Unit Type" => unit_type,
            "Fuel" => unit_type in ("WT", "WIND") ? "WIND" : unit_type in ("PVe", "PV", "RTPV") ? "SOLAR" : "N/A",
            "Category" => string(technology.category),
            "Size" => string(input_row.Size),
            "Min Gen pu" => _option_default(defaults_by_unit, unit_type, "Min Gen pu", "0"),
            "Min Up Time Hr" => _option_default(defaults_by_unit, unit_type, "Min Up Time Hr", "0"),
            "Min Down Time Hr" => _option_default(defaults_by_unit, unit_type, "Min Down Time Hr", "0"),
            "Shut Down Cost" => _option_default(defaults_by_unit, unit_type, "Shut Down Cost", "0"),
            "Ramp Rate pu/Hr" => _option_default(defaults_by_unit, unit_type, "Ramp Rate pu/Hr", "1"),
            "Input Power Rating pu" => technology_class == "storage" ? _option_default(defaults_by_unit, unit_type, "Input Power Rating pu", "1") : "n/a",
            "Output Power Rating pu" => technology_class == "storage" ? _option_default(defaults_by_unit, unit_type, "Output Power Rating pu", "1") : "n/a",
            "Min Storgae pu" => technology_class == "storage" ? _option_default(defaults_by_unit, unit_type, "Min Storgae pu", "0") : "n/a",
            "Duration Hr" => technology_class == "storage" ? duration_default : "n/a",
            "Round Trip Efficiency pu" => technology_class == "storage" ? _option_default(defaults_by_unit, unit_type, "Round Trip Efficiency pu", "0.85") : "n/a",
            "Fuel Price \$/MMBTU" => _option_default(defaults_by_unit, unit_type, "Fuel Price \$/MMBTU", "0"),
            "Start Heat Cold MBTU" => _option_default(defaults_by_unit, unit_type, "Start Heat Cold MBTU", "0"),
            "CO2_Emissions ton/MMBTU" => _option_default(defaults_by_unit, unit_type, "CO2_Emissions ton/MMBTU", "0"),
            "Inertia MJ/MW" => _option_default(defaults_by_unit, unit_type, "Inertia MJ/MW", "0"),
            "FOR" => _option_default(defaults_by_unit, unit_type, "FOR", "0"),
            "MTTR Hr" => _option_default(defaults_by_unit, unit_type, "MTTR Hr", "0"),
            "Lagtime" => _option_default(defaults_by_unit, unit_type, "Lagtime", "1"),
            "Capex Years" => _option_default(defaults_by_unit, unit_type, "Capex Years", "20"),
            "Lifetime" => _option_default(defaults_by_unit, unit_type, "Lifetime", "30"),
            "Fixed OM Cost per MW" => _option_default(defaults_by_unit, unit_type, "Fixed OM Cost per MW", "0"),
            "Preference Multiplier" => _option_default(defaults_by_unit, unit_type, "Preference Multiplier", "1"),
            "Zone" => string(input_row.Zone),
            "Bus ID" => string(input_row[Symbol("Bus ID")]),
            "Capacity Eligible" => capacity_default,
            "REC Eligible" => rec_default,
            "Synchronous_Inertia" => technology_class == "thermal" ? "true" : "false",
        )
        for field in ("Output_pct_0", "Output_pct_1", "Output_pct_2", "Output_pct_3", "Output_pct_4", "HR_avg_0", "HR_incr_1", "HR_incr_2", "HR_incr_3", "HR_incr_4")
            row[field] = _option_default(defaults_by_unit, unit_type, field, "0")
        end
        push!(rows, [row[column] for column in PROJECT_OPTION_DOWNSTREAM_COLUMNS])
    end
    return rows
end

function _write_projectoptions_for_investor(
    spec_dir::AbstractString,
    investor_dir::AbstractString,
    investor::AbstractString,
    technologies::DataFrames.DataFrame,
)
    options_path = joinpath(spec_dir, "projectoptions.csv")
    isfile(options_path) || return nothing
    options = DataFrames.DataFrame(CSV.File(options_path; stringtype=String))
    completed_options = _complete_projectoptions_for_investor(
        options,
        investor,
        technologies,
        load_project_defaults(),
    )
    CSV.write(joinpath(investor_dir, "projectoptions.csv"), completed_options)
    return completed_options
end

function _write_project_metadata(project_root::AbstractString, values::Dict{String, String})
    metadata_path = joinpath(project_root, "project_metadata.csv")
    open(metadata_path, "w") do io
        write(io, "key,value\n")
        for key in sort!(collect(keys(values)))
            write(io, "$(key),$(replace(values[key], "," => "\\,"))\n")
        end
    end
    return metadata_path
end

function _find_ref_het_dir(ref_case::AbstractString, base_dir_name::AbstractString, heterogeneity::AbstractString)
    d1 = joinpath(ref_case, base_dir_name, heterogeneity)
    isdir(d1) && return d1
    d2 = joinpath(ref_case, heterogeneity)
    isdir(d2) && return d2
    return ""
end

function _copy_reference_investor_assets(ref_case::Union{Nothing, AbstractString}, base_dir_name::AbstractString, heterogeneity::AbstractString, project_dir::AbstractString, investor::AbstractString)
    isnothing(ref_case) && return project_dir
    ref_case = normpath(String(ref_case))
    ref_het = _find_ref_het_dir(ref_case, base_dir_name, heterogeneity)
    isempty(ref_het) && return project_dir

    candidates = [
        investor,
        replace(investor, "_" => ""),
        replace(investor, "investor_" => "investor"),
        replace(investor, "investor" => "investor_"),
    ]

    ref_investor_dir = ""
    for cand in unique(candidates)
        d = joinpath(ref_het, "investors", cand)
        if isdir(d)
            ref_investor_dir = d
            break
        end
    end

    isempty(ref_investor_dir) && return project_dir
    _copy_directory_contents(ref_investor_dir, project_dir)
    return project_dir
end

function initialize_emis_project(spec_dir::AbstractString; output_dir::AbstractString, reference_case_dir=nothing)
    isdir(spec_dir) || error("Project spec directory does not exist: $(spec_dir)")
    spec = _read_project_spec(spec_dir)

    system_name = _project_spec_value(spec, "system_name")
    scratch_dir = _project_spec_value(spec, "scratch_dir")
    outage_filepath = _project_spec_value(spec, "outage_filepath")
    heterogeneity = _project_spec_value(spec, "heterogeneity")
    investors = _split_investors(spec)
    isempty(investors) && error("project_spec.csv contains no investors")

    base_dir_name = get(spec, "base_dir_name", "EMIS_RTS_Analysis")
    test_system_dir_name = get(spec, "test_system_dir_name", "RTS-GMLC")
    runs_dir_name = get(spec, "runs_dir_name", "HPC_Analysis_Runs")
    system_filepath = get(spec, "system_filepath", "")
    time_series_data_dir = get(spec, "time_series_data_dir", "")
    start_year = get(spec, "start_year", "2020")
    pcm_scenario = get(spec, "pcm_scenario", "scenario_1")

    base_dir = joinpath(output_dir, base_dir_name)
    test_system_dir = joinpath(output_dir, test_system_dir_name)
    runs_dir = joinpath(output_dir, runs_dir_name)
    project_timeseries_dir = joinpath(test_system_dir, "RTS_Data", "timeseries_data_files")

    template_replacements = Dict(
        "PROJECT_ROOT" => normpath(output_dir),
        "BASE_DIR" => normpath(base_dir),
        "TEST_SYSTEM_DIR" => normpath(test_system_dir),
        "RUNS_DIR" => normpath(runs_dir),
        "HETEROGENEITY" => heterogeneity,
        "HETEROGENEITY_BOOL" => string(lowercase(heterogeneity) == "heterogeneous"),
        "SYSTEM_NAME" => system_name,
        "SCRATCH_DIR" => _resolve_path(spec_dir, scratch_dir),
        "OUTAGE_FILEPATH" => _resolve_path(spec_dir, outage_filepath),
        "SYSTEM_FILEPATH" => _resolve_path(spec_dir, system_filepath),
        "TIME_SERIES_DATA_DIR" => normpath(project_timeseries_dir),
        "CASE_NAME" => "{{CASE_NAME}}",
        "START_YEAR" => start_year,
        "PCM_SCENARIO" => pcm_scenario,
    )

    if isdir(base_dir)
        existing_entries = readdir(base_dir)
        if !isempty(existing_entries)
            error("Project initialization target already exists and is not empty: $(base_dir)")
        end
    end
    if isdir(test_system_dir) || isdir(runs_dir)
        error("Project initialization target already exists: $(output_dir)")
    end

    mkpath(output_dir)
    mkpath(joinpath(base_dir, heterogeneity))
    mkpath(joinpath(base_dir, heterogeneity, "markets_data"))
    template_dir = _write_case_template(output_dir, template_replacements)
    _copy_case_template_overrides(spec_dir, template_dir, template_replacements)
    metadata_path = _write_project_metadata(output_dir, Dict(
        "base_dir" => normpath(base_dir),
        "runs_dir" => normpath(runs_dir),
        "test_system_dir" => normpath(test_system_dir),
        "timeseries_data_dir" => normpath(project_timeseries_dir),
        "heterogeneity" => heterogeneity,
    ))
    mkpath(runs_dir)
    mkpath(test_system_dir)

    target_config_dir = _write_system_config_bundle(base_dir, spec_dir, spec)
    technologies = DataFrames.DataFrame(CSV.File(joinpath(target_config_dir, "technologies.csv"); stringtype=String))

    spec_markets_dir = joinpath(spec_dir, "markets_data")
    if isdir(spec_markets_dir)
        _copy_directory_contents(spec_markets_dir, joinpath(base_dir, heterogeneity, "markets_data"))
    end
    spec_queue_cost = joinpath(spec_dir, "queue_cost_data.csv")
    if isfile(spec_queue_cost)
        cp(spec_queue_cost, joinpath(base_dir, "queue_cost_data.csv"); force=true)
    end

    ref_case = isnothing(reference_case_dir) ? get(spec, "reference_case_dir", nothing) : String(reference_case_dir)
    if !isnothing(ref_case) && !isempty(strip(ref_case))
        ref_case = normpath(strip(ref_case))
        ref_het = _find_ref_het_dir(ref_case, base_dir_name, heterogeneity)
        if !isempty(ref_het)
            ref_markets = joinpath(ref_het, "markets_data")
            if isdir(ref_markets)
                _copy_directory_contents(ref_markets, joinpath(base_dir, heterogeneity, "markets_data"))
            end
            ref_queue = joinpath(ref_het, "queue_cost_data.csv")
            if !isfile(ref_queue)
                ref_queue = joinpath(ref_case, "queue_cost_data.csv")
            end
            if isfile(ref_queue)
                cp(ref_queue, joinpath(base_dir, "queue_cost_data.csv"); force=true)
            end
        end
    end

    for investor in investors
        investor_dir = joinpath(base_dir, heterogeneity, "investors", investor)
        investor_markets_dir = joinpath(investor_dir, "markets_data")
        mkpath(investor_markets_dir)

        if !isnothing(ref_case) && !isempty(strip(ref_case))
            _copy_reference_investor_assets(ref_case, base_dir_name, heterogeneity, investor_dir, investor)
        end

        spec_inv_dir = joinpath(spec_dir, "investors", investor)
        if isdir(spec_inv_dir)
            _copy_directory_contents(spec_inv_dir, investor_dir)
        end

        for fn in ("characteristics.csv", "finance_params.csv", "MACRS Schedule.csv", "project_capex.csv", "sizedict.csv")
            spec_file = joinpath(spec_dir, fn)
            if isfile(spec_file)
                cp(spec_file, joinpath(investor_dir, fn); force=true)
            end
        end

        if isfile(joinpath(spec_dir, "projectexisting.csv"))
            cp(joinpath(spec_dir, "projectexisting.csv"), joinpath(investor_dir, "projectexisting.csv"); force=true)
        end
        if isfile(joinpath(spec_dir, "projectoptions.csv"))
            _write_projectoptions_for_investor(spec_dir, investor_dir, investor, technologies)
        end
    end

    copied_system_path = ""
    if !isempty(system_filepath)
        resolved_system = _resolve_path(spec_dir, system_filepath)
        if !isempty(resolved_system) && isfile(resolved_system)
            copied_system_path = _copy_system_file_bundle(resolved_system, test_system_dir)
        end
    end

    if !isempty(time_series_data_dir)
        resolved_ts = _resolve_path(spec_dir, time_series_data_dir)
        if !isempty(resolved_ts) && isdir(resolved_ts)
            mkpath(dirname(project_timeseries_dir))
            cp(resolved_ts, project_timeseries_dir; force=true)
        end
    end

    if !isempty(copied_system_path)
        source_data_dir = joinpath(test_system_dir, "RTS_Data", "SourceData")
        sys = PSY.System(copied_system_path; runchecks=false)
        write_system_inputs(sys, source_data_dir; technologies=technologies)
    end

    if !isnothing(reference_case_dir)
        reference_case_dir = normpath(String(reference_case_dir))
        if !isempty(reference_case_dir) && !isdir(reference_case_dir)
            error("reference_case_dir does not exist: $(reference_case_dir)")
        end
    end

    manifest = build_project_init_manifest()
    validate_project_init_manifest(manifest)

    return Dict(
        :spec_dir => spec_dir,
        :output_dir => output_dir,
        :base_dir => base_dir,
        :test_system_dir => test_system_dir,
        :runs_dir => runs_dir,
        :timeseries_data_dir => project_timeseries_dir,
        :heterogeneity => heterogeneity,
        :system_name => system_name,
        :scratch_dir => scratch_dir,
        :outage_filepath => outage_filepath,
        :investors => investors,
        :system_config_dir => target_config_dir,
        :case_templates_dir => template_dir,
        :project_metadata => metadata_path,
        :template_replacements => template_replacements,
        :manifest => manifest,
    )
end
