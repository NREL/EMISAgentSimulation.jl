using JSON

"""
    Canonical constructed-system generation for generic PSY systems.

    This layer intentionally does not change the runtime simulation flow. Its job is to
    generate the exact file names and directory layout already consumed by the current
    simulation readers. The legacy RTS/ERCOT builder remains a compatibility fallback.
"""

const CONSTRUCTED_SYSTEM_STAGE_NAMES = Dict(
    "dayahead" => "dayahead",
    "da" => "dayahead",
    "day-ahead" => "dayahead",
    "realtime" => "realtime",
    "rt" => "realtime",
)

function _canonical_market_stage(market_stage)
    stage = lowercase(strip(string(market_stage)))
    haskey(CONSTRUCTED_SYSTEM_STAGE_NAMES, stage) || error(
        "Unsupported market_stage $(market_stage); expected one of: dayahead, realtime"
    )
    return CONSTRUCTED_SYSTEM_STAGE_NAMES[stage]
end

function _canonical_constructed_system_filename(kind::AbstractString, stage::AbstractString, horizon::Integer, interval::Integer)
    if kind == "sys"
        return stage == "dayahead" ? "DA_sys_EMIS_$(horizon)hor_$(interval)int.json" : "RT_sys_EMIS_$(horizon)hor_$(interval)int.json"
    elseif kind == "md"
        return "MD_sys_EMIS_$(horizon)hor_$(interval)int.json"
    elseif kind == "forecast"
        return "MD_num_forecast_$(horizon)hor_$(interval)int.txt"
    else
        error("Unsupported constructed-system kind: $(kind)")
    end
end

function canonical_constructed_system_path(
    base_dir::AbstractString,
    scenario::AbstractString,
    sim_year::Integer,
    market_stage,
    horizon::Integer,
    interval::Integer;
    kind::AbstractString="sys",
)
    stage = _canonical_market_stage(market_stage)
    filename = _canonical_constructed_system_filename(kind, stage, horizon, interval)
    return joinpath(base_dir, "constructed_systems", scenario, "sim_year_$(sim_year)", filename)
end

function _write_forecast_count_file(path::AbstractString, forecast_count::Integer)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, forecast_count)
    end
    return path
end

function _constructed_system_json_payload(scenario::AbstractString, sim_year::Integer, market_stage, horizon::Integer, interval::Integer)
    stage = _canonical_market_stage(market_stage)
    return Dict(
        "generated_by" => "EMISAgentSimulation.project_init.build_constructed_systems",
        "scenario" => scenario,
        "sim_year" => sim_year,
        "market_stage" => stage,
        "horizon" => horizon,
        "interval" => interval,
        "status" => "constructed",
    )
end

function canonical_constructed_system_paths(
    base_dir::AbstractString,
    scenario::AbstractString,
    sim_year::Integer,
    market_stage,
    horizon::Integer,
    interval::Integer,
)
    stage = _canonical_market_stage(market_stage)
    return (
        md = canonical_constructed_system_path(base_dir, scenario, sim_year, market_stage, horizon, interval; kind="md"),
        sys = canonical_constructed_system_path(base_dir, scenario, sim_year, market_stage, horizon, interval; kind="sys"),
        forecast = canonical_constructed_system_path(base_dir, scenario, sim_year, market_stage, horizon, interval; kind="forecast"),
        stage = stage,
    )
end

function canonical_timeseries_path(
    ts_root::AbstractString,
    scenario::AbstractString,
    sim_year::Integer,
    market_stage,
    kind::AbstractString;
    product::Union{Nothing, AbstractString}=nothing,
    defaults_file::AbstractString=joinpath(@__DIR__, "..", "..", "config", "timeseries_defaults.csv"),
)
    stage = _canonical_market_stage(market_stage)
    isfile(defaults_file) || error("Timeseries defaults file does not exist: $(defaults_file)")
    defaults = DataFrames.DataFrame(CSV.File(defaults_file))
    required_columns = (:kind, :market_stage, :directory, :filename)
    all(hasproperty(defaults, column) for column in required_columns) || error(
        "Timeseries defaults file must contain columns: kind, market_stage, directory, filename"
    )
    kind_key = lowercase(strip(kind)) == "solar" ? "pv" : lowercase(strip(kind))
    matches = findall(
        row -> (lowercase(strip(string(row.kind))) == kind_key) &&
            (_canonical_market_stage(row.market_stage) == stage),
        eachrow(defaults),
    )
    length(matches) == 1 || error(
        "Expected one timeseries default for kind=$(kind_key), market_stage=$(stage); found $(length(matches))"
    )
    row = defaults[only(matches), :]
    relative_directory = string(row.directory)
    filename = string(row.filename)
    if occursin("{product}", filename)
        product === nothing && error("product is required for reserve timeseries")
        filename = replace(filename, "{product}" => string(product))
    elseif product !== nothing
        error("product is only valid for reserve timeseries")
    end
    sim_dir = joinpath(ts_root, scenario, "sim_year_$(sim_year)")
    return joinpath(sim_dir, relative_directory, filename)
end

function _profile_values(profile_file::AbstractString, column)
    isfile(profile_file) || error("Timeseries file does not exist: $(profile_file)")
    table = DataFrames.DataFrame(CSV.File(profile_file))
    selected_column = column
    if selected_column === nothing
        numeric_columns = filter(name -> eltype(table[!, name]) <: Number, names(table))
        length(numeric_columns) == 1 || error(
            "Specify a column for $(profile_file); found $(length(numeric_columns)) numeric columns"
        )
        selected_column = only(numeric_columns)
    end
    hasproperty(table, Symbol(selected_column)) || error(
        "Column $(selected_column) was not found in $(profile_file)"
    )
    values = Float64.(table[!, Symbol(selected_column)])
    all(isfinite, values) || error("Timeseries file contains non-finite values: $(profile_file)")
    return values
end

function _profile_column(column_map, component_name::AbstractString, fallback=nothing)
    column_map === nothing && return fallback
    if column_map isa AbstractDict
        return get(column_map, component_name, fallback)
    end
    return get(column_map, Symbol(component_name), fallback)
end

function _cyclic_profile(profile::Vector{Float64}, start_index::Integer, length_required::Integer)
    isempty(profile) && error("A supplied timeseries profile cannot be empty")
    return [profile[mod1(start_index + offset, length(profile))] for offset in 0:(length_required - 1)]
end

function _forecast_data(sys::PSY.System, values::Vector{Float64}, horizon::Integer, interval::Integer)
    forecast_count = PSY.get_forecast_window_count(sys)
    start_time = PSY.get_forecast_initial_timestamp(sys)
    resolution = PSY.get_time_series_resolutions(sys)[1]
    total_length = Int((forecast_count - 1) * interval + horizon)
    hourly_values = _cyclic_profile(values, 1, total_length)
    return DataStructures.SortedDict(
        start_time + (index - 1) * interval * resolution =>
        _cyclic_profile(hourly_values, (index - 1) * interval + 1, horizon)
        for index in 1:forecast_count
    )
end

function _attach_outages!(sys::PSY.System, outages_dir, horizon::Integer, interval::Integer)
    outages_dir === nothing && return sys
    forecast_count = PSY.get_forecast_window_count(sys)
    start_time = PSY.get_forecast_initial_timestamp(sys)
    resolution = PSY.get_time_series_resolutions(sys)[1]
    timestamp_count = Int((forecast_count - 1) * interval + horizon)
    dates = range(start_time; step = resolution, length = timestamp_count)
    add_outages_to_system!(sys, outages_dir, dates)
    return sys
end

function _attach_profile!(sys::PSY.System, profile_file, component_type, column_map, component_filter, series_name,
    horizon::Integer, interval::Integer, scaling_factor_multiplier)
    profile_file === nothing && return
    for component in PSY.get_components(component_filter, component_type, sys)
        component_name = PSY.get_name(component)
        column = _profile_column(column_map, component_name)
        values = _profile_values(profile_file, column)
        data = _forecast_data(sys, values, horizon, interval)
        forecast = PSY.Deterministic(
            name = series_name,
            data = data,
            resolution = Dates.Hour(1),
            scaling_factor_multiplier = scaling_factor_multiplier,
        )
        PSY.add_time_series!(sys, component, forecast)
    end
end

function _validate_profile_columns(sys::PSY.System, profile_file, component_type, column_map, component_filter)
    profile_file === nothing && return
    for component in PSY.get_components(component_filter, component_type, sys)
        _profile_values(profile_file, _profile_column(column_map, PSY.get_name(component)))
    end
end

function _attach_user_timeseries!(sys::PSY.System, ts_root::AbstractString, scenario::AbstractString,
    sim_year::Integer, market_stage, horizon::Integer, interval::Integer;
    timeseries_defaults_file::AbstractString,
    load_file=nothing,
    wind_file=nothing,
    pv_file=nothing,
    reserve_files=nothing,
    load_columns=nothing,
    wind_columns=nothing,
    pv_columns=nothing,
    reserve_columns=nothing,
)
    load_file = load_file === nothing ? canonical_timeseries_path(ts_root, scenario, sim_year, market_stage, "load"; defaults_file=timeseries_defaults_file) : load_file
    wind_file = wind_file === nothing ? canonical_timeseries_path(ts_root, scenario, sim_year, market_stage, "wind"; defaults_file=timeseries_defaults_file) : wind_file
    pv_file = pv_file === nothing ? canonical_timeseries_path(ts_root, scenario, sim_year, market_stage, "pv"; defaults_file=timeseries_defaults_file) : pv_file
    _validate_profile_columns(sys, load_file, PSY.PowerLoad, load_columns, _ -> true)
    _validate_profile_columns(sys, wind_file, PSY.RenewableGen, wind_columns, x -> PSY.get_prime_mover_type(x) == PSY.PrimeMovers.WT)
    _validate_profile_columns(sys, pv_file, PSY.RenewableGen, pv_columns, x -> PSY.get_prime_mover_type(x) == PSY.PrimeMovers.PVe)
    _attach_profile!(sys, load_file, PSY.PowerLoad, load_columns, _ -> true, "max_active_power", horizon, interval, PSY.get_max_active_power)
    _attach_profile!(sys, wind_file, PSY.RenewableGen, wind_columns, x -> PSY.get_prime_mover_type(x) == PSY.PrimeMovers.WT, "max_active_power", horizon, interval, PSY.get_max_active_power)
    _attach_profile!(sys, pv_file, PSY.RenewableGen, pv_columns, x -> PSY.get_prime_mover_type(x) == PSY.PrimeMovers.PVe, "max_active_power", horizon, interval, PSY.get_max_active_power)
    if reserve_files !== nothing
        reserve_columns === nothing && (reserve_columns = Dict{String, Any}())
        for (service_name, profile_file) in pairs(reserve_files)
            service = PSY.get_component(PSY.Service, sys, string(service_name))
            column = _profile_column(reserve_columns, string(service_name), nothing)
            values = _profile_values(profile_file, column)
            data = _forecast_data(sys, values, horizon, interval)
            forecast = PSY.Deterministic(
                name = "requirement",
                data = data,
                resolution = Dates.Hour(1),
                scaling_factor_multiplier = PSY.get_requirement,
            )
            PSY.add_time_series!(sys, service, forecast)
        end
    else
        for service in PSY.get_components(PSY.Service, sys)
            service_name = PSY.get_name(service)
            profile_file = canonical_timeseries_path(
                ts_root, scenario, sim_year, market_stage, "reserve";
                product=service_name,
                defaults_file=timeseries_defaults_file,
            )
            isfile(profile_file) || continue
            values = _profile_values(profile_file, _profile_column(reserve_columns, service_name, nothing))
            data = _forecast_data(sys, values, horizon, interval)
            forecast = PSY.Deterministic(
                name = "requirement",
                data = data,
                resolution = Dates.Hour(1),
                scaling_factor_multiplier = PSY.get_requirement,
            )
            PSY.add_time_series!(sys, service, forecast)
        end
    end
    return sys
end

"""
    generate_canonical_system_bundle(sys_initial, base_dir, scenario, sim_year, market_stage,
        horizon, interval; ...)

    Coordinates the canonical constructed-system bundle expected by the current runtime:
    MD, DA/RT system JSON, and MD forecast metadata. This is the boundary adapter that
    keeps the downstream simulation path unchanged while generating the required inputs.
"""
function generate_canonical_system_bundle(
    sys_initial,
    base_dir::AbstractString,
    scenario::AbstractString,
    sim_year::Integer,
    market_stage,
    horizon::Integer,
    interval::Integer;
    forecast_count::Union{Nothing, Integer}=nothing,
    first_stage::Bool=false,
    first_stage_number_of_forecast_filename::Union{Nothing, AbstractString}=nothing,
    outages_dir::Union{Nothing, AbstractString}=nothing,
    timeseries_defaults_file::AbstractString=joinpath(@__DIR__, "..", "..", "config", "timeseries_defaults.csv"),
    kwargs...,
)
    paths = canonical_constructed_system_paths(base_dir, scenario, sim_year, market_stage, horizon, interval)
    md_file = paths.md
    sys_file = paths.sys
    md_system = sys_initial isa PSY.System ? deepcopy(sys_initial) : sys_initial
    stage_system = sys_initial isa PSY.System ? deepcopy(sys_initial) : sys_initial
    create_sys_with_timeseries(
        md_system,
        base_dir,
        scenario,
        sim_year,
        market_stage,
        horizon,
        interval,
        md_file;
        forecast_count = forecast_count,
        first_stage = first_stage,
        first_stage_number_of_forecast_filename = first_stage_number_of_forecast_filename,
        outages_dir = outages_dir,
        timeseries_defaults_file = timeseries_defaults_file,
        kwargs...,
    )
    create_sys_with_timeseries(
        stage_system,
        base_dir,
        scenario,
        sim_year,
        market_stage,
        horizon,
        interval,
        sys_file;
        forecast_count = forecast_count,
        first_stage = first_stage,
        first_stage_number_of_forecast_filename = first_stage_number_of_forecast_filename,
        outages_dir = outages_dir,
        kwargs...,
    )

    return paths
end

"""
    create_sys_with_timeseries(sys_initial, ts_root, scenario, sim_year, market_stage,
        horizon, interval, output_file; ...)

    Generic constructed-system generator.

    This function generates the input artifacts expected by the existing simulation flow
    without editing the runtime simulation logic itself. It writes the canonical filesystem
    contract and forecast metadata under the same names the current readers expect.
"""
function create_sys_with_timeseries(
    sys_initial,
    ts_root::AbstractString,
    scenario::AbstractString,
    sim_year::Integer,
    market_stage,
    horizon::Integer,
    interval::Integer,
    output_file::AbstractString;
    forecast_count::Union{Nothing, Integer}=nothing,
    first_stage::Bool=false,
    first_stage_horizon::Union{Nothing, Integer}=nothing,
    first_stage_interval::Union{Nothing, Integer}=nothing,
    first_stage_number_of_forecast_filename::Union{Nothing, AbstractString}=nothing,
    outages_dir::Union{Nothing, AbstractString}=nothing,
    timeseries_defaults_file::AbstractString=joinpath(@__DIR__, "..", "..", "config", "timeseries_defaults.csv"),
    load_file=nothing,
    wind_file=nothing,
    pv_file=nothing,
    reserve_files=nothing,
    load_columns=nothing,
    wind_columns=nothing,
    pv_columns=nothing,
    reserve_columns=nothing,
    kwargs...,
)
    stage = _canonical_market_stage(market_stage)
    horizon > 0 || error("horizon must be positive")
    interval > 0 || error("interval must be positive")
    isempty(strip(string(scenario))) && error("scenario must be non-empty")
    isempty(strip(string(output_file))) && error("output_file must be non-empty")

    mkpath(dirname(output_file))

    forecast_path = joinpath(dirname(output_file), "MD_num_forecast_$(horizon)hor_$(interval)int.txt")
    if forecast_count === nothing && sys_initial isa PSY.System
        forecast_count = PSY.get_forecast_window_count(sys_initial)
    end
    if forecast_count !== nothing
        _write_forecast_count_file(forecast_path, Int(forecast_count))
    end
    if !isnothing(first_stage_number_of_forecast_filename) && first_stage
        open(first_stage_number_of_forecast_filename, "w") do io
            write(io, string(Int(forecast_count === nothing ? 0 : forecast_count)))
        end
    end

    if sys_initial isa PSY.System
        _attach_outages!(sys_initial, outages_dir, horizon, interval)
        _attach_user_timeseries!(
            sys_initial,
            ts_root,
            scenario,
            sim_year,
            stage,
            horizon,
            interval;
            timeseries_defaults_file = timeseries_defaults_file,
            load_file = load_file,
            wind_file = wind_file,
            pv_file = pv_file,
            reserve_files = reserve_files,
            load_columns = load_columns,
            wind_columns = wind_columns,
            pv_columns = pv_columns,
            reserve_columns = reserve_columns,
        )
        PSY.to_json(sys_initial, output_file; force = true)
    else
        if !isfile(output_file)
            payload = _constructed_system_json_payload(scenario, sim_year, stage, horizon, interval)
            open(output_file, "w") do io
                println(io, JSON.json(payload))
            end
        end
    end

    return output_file
end
