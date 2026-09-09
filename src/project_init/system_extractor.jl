function _mapping_key_value(value)
    value === missing ? "" : strip(string(value))
end

function _psy_component_type(component)
    if component isa PSY.ThermalGen
        return "ThermalGen"
    elseif component isa PSY.RenewableGen
        return "RenewableGen"
    elseif component isa PSY.HydroGen
        return "HydroGen"
    elseif component isa PSY.Storage
        return "Storage"
    end
    error("Unsupported PSY component type $(typeof(component))")
end

function _psy_mapping_key(component)
    component_type = _psy_component_type(component)
    if component isa PSY.Storage
        return (component_type, "", "")
    end
    prime_mover = _mapping_key_value(PSY.get_prime_mover_type(component))
    fuel = if component isa PSY.ThermalGen
        _mapping_key_value(PSY.get_fuel(component))
    elseif component isa PSY.HydroGen
        "HYDRO"
    elseif component isa PSY.RenewableGen
        prime_mover in ("WT", "Wind") ? "WIND" : "SOLAR"
    else
        ""
    end
    return (component_type, prime_mover, fuel)
end

function _build_psy_mapping(mapping::DataFrames.DataFrame)
    _require_columns(mapping, PSY_CLASSIFICATION_COLUMNS, "PSY classification mapping")
    result = Dict{NTuple{3, String}, String}()
    for row in eachrow(mapping)
        key = (
            _mapping_key_value(row.component_type),
            _mapping_key_value(row.prime_mover),
            _mapping_key_value(row.fuel),
        )
        result[key] = _mapping_key_value(row.unit_type)
    end
    return result
end

"""Return the configured EMIS unit type for a PSY generator or storage device."""
function classify_psy_component(component, mapping::DataFrames.DataFrame)
    key = _psy_mapping_key(component)
    lookup = _build_psy_mapping(mapping)
    haskey(lookup, key) || error(
        "No PSY classification mapping for $(PSY.get_name(component)): " *
        "component_type=$(key[1]), prime_mover=$(key[2]), fuel=$(key[3])"
    )
    return lookup[key]
end

function _sorted_components(component_type, sys)
    components = [
        component for component in PSY.get_components(component_type, sys) if
        PSY.get_available(component)
    ]
    sort!(components; by=PSY.get_name)
    return components
end

"""Extract stable zone identifiers and PSY area aliases."""
function extract_zones(sys::PSY.System)
    areas = sort!(collect(PSY.get_components(PSY.Area, sys)); by=PSY.get_name)
    isempty(areas) && error("PSY system does not contain any Area components")
    names = PSY.get_name.(areas)
    return DataFrames.DataFrame(
        zone_id=["zone_$(index)" for index in eachindex(names)],
        zone_name=names,
        area_name=names,
        load_column=["load_$(name)" for name in names],
    )
end

function _component_capacity_mw(component, sys::PSY.System)
    try
        limits = if component isa PSY.Storage
            PSY.get_output_active_power_limits(component)
        else
            PSY.get_active_power_limits(component)
        end
        return limits.max * PSY.get_base_power(sys)
    catch exception
        if exception isa MethodError || exception isa ArgumentError
            return PSY.get_rating(component)
        end
        rethrow()
    end
end

function _component_min_capacity_mw(component, sys::PSY.System)
    try
        limits = component isa PSY.Storage ?
                 PSY.get_output_active_power_limits(component) :
                 PSY.get_active_power_limits(component)
        return limits.min * PSY.get_base_power(sys)
    catch exception
        if exception isa MethodError || exception isa ArgumentError
            return 0.0
        end
        rethrow()
    end
end

function _component_status(component)
    try
        return PSY.get_status(component)
    catch exception
        if exception isa MethodError
            return PSY.get_available(component)
        end
        rethrow()
    end
end

"""Extract the generator/storage rows consumed by the RTS reader."""
function extract_gen_table(
    sys::PSY.System,
    mapping::DataFrames.DataFrame;
    technologies::Union{Nothing, DataFrames.DataFrame}=nothing,
)
    technologies !== nothing && validate_psy_classification_mapping(mapping, technologies)
    components = vcat(
        _sorted_components(PSY.Generator, sys),
        _sorted_components(PSY.Storage, sys),
    )
    rows = DataFrames.DataFrame(
        "GEN UID" => String[],
        "Unit Type" => String[],
        "PMax MW" => Float64[],
        "PMin MW" => Float64[],
        "Bus" => Int[],
        "Zone" => String[],
        "Status" => Bool[],
    )
    failures = String[]
    for component in components
        try
            bus = PSY.get_bus(component)
            area = PSY.get_area(bus)
            push!(rows, (
                PSY.get_name(component),
                classify_psy_component(component, mapping),
                _component_capacity_mw(component, sys),
                _component_min_capacity_mw(component, sys),
                PSY.get_number(bus),
                PSY.get_name(area),
                _component_status(component),
            ))
        catch exception
            push!(failures, "$(PSY.get_name(component)): $(sprint(showerror, exception))")
        end
    end
    isempty(failures) || error(
        "PSY generator extraction failed for $(length(failures)) component(s):\n" *
        join(failures, "\n")
    )
    return rows
end

function _branch_endpoints(branch)
    arc = PSY.get_arc(branch)
    return PSY.get_from(arc), PSY.get_to(arc)
end

function _numeric_or_default(value, default::Float64=0.0)
    if value === nothing || value === missing
        return default
    end
    return Float64(value)
end

function _branch_value(getter, branch, default::Float64=0.0)
    try
        return _numeric_or_default(getter(branch), default)
    catch exception
        if exception isa MethodError || exception isa ArgumentError
            return default
        end
        rethrow()
    end
end

function _branch_row(branch, sys::PSY.System; dc::Bool)
    from_bus, to_bus = _branch_endpoints(branch)
    rating = if dc
        (
            _numeric_or_default(PSY.get_active_power_limits_from(branch).max),
            _numeric_or_default(PSY.get_active_power_limits_to(branch).max),
        )
    else
        rate = _branch_value(PSY.get_rating, branch)
        (rate, rate)
    end
    resistance = dc ? 0.0 : _branch_value(PSY.get_r, branch)
    reactance = dc ? 0.0 : _branch_value(PSY.get_x, branch)
    susceptance = dc ? 0.0 : _branch_value(PSY.get_b, branch)
    return (
        PSY.get_name(branch),
        string(PSY.get_number(from_bus)),
        string(PSY.get_number(to_bus)),
        resistance,
        reactance,
        susceptance,
        rating[1],
        rating[2],
    )
end

"""Extract supported AC and DC branches into reader-compatible tables."""
function extract_branches(sys::PSY.System)
    ac_rows = DataFrames.DataFrame(
        "UID" => String[],
        "From Bus" => String[],
        "To Bus" => String[],
        "R" => Float64[],
        "X" => Float64[],
        "B" => Float64[],
        "Cont Rating" => Float64[],
    )
    dc_rows = DataFrames.DataFrame(
        "UID" => String[],
        "From Bus" => String[],
        "To Bus" => String[],
        "MW Load" => Float64[],
    )
    failures = String[]
    for (component_type, destination, is_dc) in (
        (PSY.ACTransmission, ac_rows, false),
        (PSY.DCBranch, dc_rows, true),
    )
        for branch in _sorted_components(component_type, sys)
            try
                row = _branch_row(branch, sys; dc=is_dc)
                if is_dc
                    push!(destination, row[1:4])
                else
                    push!(destination, row[1:7])
                end
            catch exception
                push!(failures, "$(PSY.get_name(branch)): $(sprint(showerror, exception))")
            end
        end
    end
    isempty(failures) || error(
        "PSY branch extraction failed for $(length(failures)) component(s):\n" *
        join(failures, "\n")
    )
    return (ac=ac_rows, dc=dc_rows)
end

const RESERVE_COLUMNS = [
    "Reserve Product",
    "Timeframe (sec)",
    "Requirement (MW)",
    "Eligible Regions",
    "Eligible Device Categories",
    "Eligible Device SubCategories",
    "Direction",
]

function _reserve_services(sys::PSY.System)
    services = PSY.Service[]
    for device_type in (PSY.Generator, PSY.Storage, PSY.HydroGen)
        for device in PSY.get_components(device_type, sys)
            for service in PSY.get_services(device)
                service in services || push!(services, service)
            end
        end
    end
    sort!(services; by=PSY.get_name)
    return services
end

function _reserve_default(defaults, product, field, fallback)
    defaults === nothing && return fallback
    default_row = findfirst(row -> string(row[Symbol("Reserve Product")]) == product, eachrow(defaults))
    default_row === nothing && return fallback
    value = defaults[default_row, Symbol(field)]
    ismissing(value) || isempty(strip(string(value))) ? fallback : value
end

function _reserve_requirement(service, defaults, product)
    _reserve_requirement(service, defaults, product, Val(:default))
end

function _reserve_requirement(service::Union{
    PSY.VariableReserve,
    PSY.ConstantReserve,
    PSY.VariableReserveNonSpinning,
    PSY.ConstantReserveNonSpinning,
    PSY.ConstantReserveGroup,
}, defaults, product, ::Val{:default})
    return Float64(PSY.get_requirement(service))
end

function _reserve_requirement(
    service::PSY.ReserveDemandCurve,
    defaults,
    product,
    ::Val{:default},
)
    fallback = _reserve_default(defaults, product, "Requirement (MW)", missing)
    return fallback === missing ? missing : Float64(fallback)
end

function _reserve_requirement(service, defaults, product, ::Val{:default})
    fallback = _reserve_default(defaults, product, "Requirement (MW)", nothing)
    return fallback === nothing ? missing : Float64(fallback)
end

function _reserve_timeframe(service, defaults, product)
    return Int(_reserve_default(defaults, product, "Timeframe (sec)", 0))
end

function _reserve_timeframe(service::PSY.ReserveDemandCurve, defaults, product)
    return Int(PSY.get_time_frame(service))
end

"""Extract reserve products using the existing EMIS reserve-table schema."""
function extract_reserves(
    sys::PSY.System;
    defaults::Union{Nothing, DataFrames.DataFrame}=nothing,
)
    defaults !== nothing && _require_columns(defaults, Symbol.(RESERVE_COLUMNS), "reserve defaults")
    rows = DataFrames.DataFrame(
        Symbol.(RESERVE_COLUMNS[1]) => String[],
        Symbol.(RESERVE_COLUMNS[2]) => Int[],
        Symbol.(RESERVE_COLUMNS[3]) => Union{Missing, Float64}[],
        Symbol.(RESERVE_COLUMNS[4]) => String[],
        Symbol.(RESERVE_COLUMNS[5]) => String[],
        Symbol.(RESERVE_COLUMNS[6]) => String[],
        Symbol.(RESERVE_COLUMNS[7]) => String[],
    )
    for service in _reserve_services(sys)
        product = PSY.get_name(service)
        requirement = _reserve_requirement(service, defaults, product)
        service_type = string(typeof(service))
        direction = occursin("Down", service_type) ? "Down" : "Up"
        push!(rows, (
            product,
            _reserve_timeframe(service, defaults, product),
            requirement,
            string(_reserve_default(defaults, product, "Eligible Regions", "")),
            string(_reserve_default(defaults, product, "Eligible Device Categories", "(Generator)")),
            string(_reserve_default(defaults, product, "Eligible Device SubCategories", "")),
            string(_reserve_default(defaults, product, "Direction", direction)),
        ))
    end
    return rows
end

const PROJECT_EXISTING_COLUMNS = [
    "GEN_UID", "Unit Type", "Size", "Fuel", "Category", "Min Gen pu",
    "Ramp Rate pu/Hr", "Input Power Rating pu", "Output Power Rating pu",
    "Min Storgae pu", "Duration Hr", "Round Trip Efficiency pu",
        "Fuel Price \$/MMBTU", "Output_pct_0", "Output_pct_1", "Output_pct_2",
    "Output_pct_3", "Output_pct_4", "HR_avg_0", "HR_incr_1", "HR_incr_2",
    "HR_incr_3", "HR_incr_4", "CO2_Emissions ton/MMBTU", "Inertia MJ/MW",
    "FOR", "MTTR Hr", "Lagtime", "Online Year", "Capex Years", "Lifetime",
    "Fixed OM Cost per MW", "Zone", "Bus ID", "Capacity Eligible", "REC Eligible",
    "Synchronous_Inertia",
]

function _defaults_by_unit(defaults::DataFrames.DataFrame)
    result = Dict{String, Dict{String, Any}}()
    for row in eachrow(defaults)
        get!(result, string(row.unit_type), Dict{String, Any}())[string(row.field)] = row.value
    end
    return result
end

function _default_value(defaults, unit_type, field, fallback)
    value = get(get(defaults, unit_type, Dict{String, Any}()), field, fallback)
    value === missing || isempty(strip(string(value))) ? fallback : value
end

function _outage_for(outage, name, fallback)
    outage === nothing && return fallback
    name in string.(names(outage)) || return fallback
    values = Float64[]
    for value in outage[!, Symbol(name)]
        value === missing || push!(values, Float64(value))
    end
    isempty(values) ? fallback : Statistics.mean(values)
end

function _outage_metadata(outage, device, unit_type, size_mw)
    outage === nothing && return nothing
    required = (:PrimeMovers, :ThermalFuels, :NameplateLimit_MW, :FOR, :MTTR)
    all(column -> column in Symbol.(names(outage)), required) || return nothing
    prime_mover = string(PSY.get_prime_mover_type(device))
    fuel = device isa PSY.ThermalGen ? string(PSY.get_fuel(device)) : "NA"
    candidates = filter(eachrow(outage)) do row
        mover_matches = string(row.PrimeMovers) == prime_mover
        fuel_value = string(row.ThermalFuels)
        fuel_matches = fuel_value == fuel || fuel_value == "NA"
        mover_matches && fuel_matches
    end
    isempty(candidates) && return nothing
    ordered = sort(collect(candidates); by=row -> Float64(row.NameplateLimit_MW))
    selected = findfirst(row -> size_mw <= Float64(row.NameplateLimit_MW), ordered)
    selected === nothing && (selected = lastindex(ordered))
    row = ordered[selected]
    for_value = Float64(row.FOR)
    for_value > 1.0 && (for_value /= 100.0)
    return (FOR=for_value, MTTR=Float64(row.MTTR))
end

function _project_row_value(defaults, unit_type, field, fallback)
    value = _default_value(defaults, unit_type, field, fallback)
    value isa AbstractString && lowercase(value) in ("true", "false") &&
        return lowercase(value) == "true"
    value isa AbstractString && tryparse(Float64, value) !== nothing &&
        return parse(Float64, value)
    return value
end

"""Extract complete existing-project rows from PSY and explicit ownership input."""
function extract_fleet(
    sys::PSY.System,
    mapping::DataFrames.DataFrame,
    technologies::DataFrames.DataFrame,
    ownership::DataFrames.DataFrame;
    defaults::DataFrames.DataFrame=load_project_defaults(),
    outage=nothing,
)
    _require_columns(ownership, [:Investor, :GEN_UID], "existing project ownership")
    _require_columns(technologies, [:unit_type, :category], "technologies.csv")
    validate_psy_classification_mapping(mapping, technologies)
    default_lookup = _defaults_by_unit(defaults)
    technology_lookup = Dict(string(row.unit_type) => string(row.category) for row in eachrow(technologies))
    devices = Dict(PSY.get_name(device) => device for device in vcat(
        _sorted_components(PSY.Generator, sys), _sorted_components(PSY.Storage, sys)
    ))
    rows = DataFrames.DataFrame([Any[] for _ in PROJECT_EXISTING_COLUMNS], PROJECT_EXISTING_COLUMNS)
    report = String[]
    for owner_row in eachrow(ownership)
        name = string(owner_row.GEN_UID)
        haskey(devices, name) || error("Existing project $(name) is not available in the PSY system")
        device = devices[name]
        unit_type = classify_psy_component(device, mapping)
        haskey(technology_lookup, unit_type) || error(
            "No technology configuration for $(unit_type)"
        )
        category = technology_lookup[unit_type]
        bus = PSY.get_bus(device)
        area = PSY.get_area(bus)
        limits = device isa PSY.Storage ? PSY.get_output_active_power_limits(device) : PSY.get_active_power_limits(device)
        size = limits.max * PSY.get_base_power(sys)
        min_pu = limits.min
        row = Dict{String, Any}(column => "NA" for column in PROJECT_EXISTING_COLUMNS)
        row["GEN_UID"] = name
        row["Unit Type"] = unit_type
        row["Size"] = size
        row["Category"] = category
        row["Min Gen pu"] = min_pu
        row["Zone"] = PSY.get_name(area)
        row["Bus ID"] = PSY.get_number(bus)
        metadata = _outage_metadata(outage, device, unit_type, size)
        row["FOR"] = metadata === nothing ?
                  _outage_for(outage, name, _project_row_value(default_lookup, unit_type, "FOR", 0.0)) :
                  metadata.FOR
        row["MTTR Hr"] = metadata === nothing ?
                  _project_row_value(default_lookup, unit_type, "MTTR Hr", 24.0) :
                  metadata.MTTR
        row["Capacity Eligible"] = _project_row_value(default_lookup, unit_type, "Capacity Eligible", true)
        row["REC Eligible"] = _project_row_value(default_lookup, unit_type, "REC Eligible", false)
        row["Synchronous_Inertia"] = unit_type in ("ST", "CT", "CC", "GT", "NU_ST")
        if device isa PSY.Storage
            row["Input Power Rating pu"] = PSY.get_input_active_power_limits(device).max
            row["Output Power Rating pu"] = PSY.get_output_active_power_limits(device).max
            row["Duration Hr"] = PSY.get_storage_capacity(device) / size
            efficiency = PSY.get_efficiency(device)
            row["Round Trip Efficiency pu"] = efficiency.in * efficiency.out
        elseif device isa PSY.ThermalGen
            operation_cost = PSY.get_operation_cost(device)
            row["Fuel"] = string(PSY.get_fuel(device))
            row["Ramp Rate pu/Hr"] = something(PSY.get_ramp_limits(device), (up=1.0,)).up
            row["Min Up Time Hr"] = PSY.get_time_limits(device)
            operation_cost === nothing && push!(report, "$(name): heat-rate fallback used")
        end
        for column in PROJECT_EXISTING_COLUMNS
            row[column] == "NA" && (row[column] = _project_row_value(default_lookup, unit_type, column, "NA"))
        end
        row["Investor"] = string(owner_row.Investor)
        push!(rows, [row[column] for column in PROJECT_EXISTING_COLUMNS])
    end
    return (projects=rows, report=report)
end

"""Write the extracted system inputs used by the current simulation readers."""
function write_system_inputs(
    sys::PSY.System,
    source_data_dir::AbstractString;
    mapping::DataFrames.DataFrame=load_psy_classification_mapping(),
    technologies::DataFrames.DataFrame,
    reserve_defaults::Union{Nothing, DataFrames.DataFrame}=nothing,
    ownership::Union{Nothing, DataFrames.DataFrame}=nothing,
    investor_dir::Union{Nothing, AbstractString}=nothing,
    defaults::DataFrames.DataFrame=load_project_defaults(),
    outage=nothing,
)
    mkpath(source_data_dir)
    zones = extract_zones(sys)
    gen = extract_gen_table(sys, mapping; technologies=technologies)
    gen = gen[:, ["GEN UID", "Unit Type", "PMax MW"]]
    branches = extract_branches(sys)
    reserves = extract_reserves(sys; defaults=reserve_defaults)
    CSV.write(joinpath(source_data_dir, "zones.csv"), zones)
    CSV.write(joinpath(source_data_dir, "gen.csv"), gen)
    CSV.write(joinpath(source_data_dir, "branch.csv"), branches.ac)
    CSV.write(joinpath(source_data_dir, "dc_branch.csv"), branches.dc)
    CSV.write(joinpath(source_data_dir, "reserves.csv"), reserves)

    report = String[]
    project_files = String[]
    if ownership !== nothing
        investor_dir === nothing && error(
            "investor_dir is required when ownership is provided"
        )
        fleet = extract_fleet(
            sys,
            mapping,
            technologies,
            ownership;
            defaults=defaults,
            outage=outage,
        )
        append!(report, fleet.report)
        for investor in sort!(unique(string.(ownership.Investor)))
            owned_names = string.(ownership.GEN_UID[string.(ownership.Investor) .== investor])
            investor_rows = filter(row -> string(row.GEN_UID) in owned_names, fleet.projects)
            directory = joinpath(investor_dir, investor)
            mkpath(directory)
            path = joinpath(directory, "projectexisting.csv")
            CSV.write(path, investor_rows)
            push!(project_files, path)
        end
    end
    return (zones=zones, gen=gen, branches=branches, reserves=reserves,
            project_files=project_files, report=report)
end