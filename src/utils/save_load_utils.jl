# HDF5 field-by-field serialization for AgentSimulation and all nested EMIS structs.
#
# Design rules
# ─────────────────────────────────────────────────────────────────────────────
# • Every concrete EMIS struct gets a save_X!(group, value) and load_X(group).
# • Abstract types (Project, Tech, Forecast, RiskPreference, Product) store their
#   concrete type name in an HDF5 attribute "type" so load dispatch can reconstruct
#   the right concrete type.
# • BuildPhase is stored as attribute "build_phase" on the project group.
# • Union{Nothing, T} fields: write attribute is_nothing=true on the group and skip
#   writing child data; on load check the attribute before reading.
# • Symbol values are stored as plain HDF5 strings and converted back with Symbol().
# • NamedTuple{(:a,:b),...} → group with datasets named "a", "b".
# • AxisArray{T,N} → group with dataset "data" + "axis<i>_names" per axis.
# • DataFrame → group with dataset "col_names" + one dataset per column.
# • PSY.System / PSY.*GenerationCost – systems are excluded entirely (PSY.to_json);
#   generation costs are saved as scalar fields + breakpoint arrays.
#
# ─────────────────────────────────────────────────────────────────────────────

SCHEMA_VERSION = 1

"""
    save_simulation(path, simulation)

Write an `AgentSimulation` to an HDF5 file at `path`.
PSY systems (system_MDs, system_UCs, system_EDs, system_PRAS) are NOT written here;
they must be saved separately with PSY.to_json.
"""
function save_simulation(simulation::AgentSimulation, save_dir::String, iteration_year::Int)
    @info "Saving simulation to directory: $save_dir"
    isdir(save_dir) || mkpath(save_dir)
    h5_path = joinpath(save_dir, "simulation_data_year_$(iteration_year).h5")
    h5open(h5_path, "w") do f
        attributes(f)["schema_version"] = SCHEMA_VERSION
        save_simulation!(f, simulation)
    end

    @info "Saving Sienna systems..."
    save_Sienna_systems(simulation, save_dir, iteration_year)
end

"""
    load_simulation(save_dir, restore_year, case) -> AgentSimulation

Reconstruct an `AgentSimulation` from an HDF5 checkpoint.
`case` is the caller's `CaseDefinition`; its `solver` field is spliced into the
reconstructed case (the solver cannot be serialized to HDF5).
PSY systems are not stored in the HDF5 file; callers must reload them separately.
"""
function load_simulation(save_dir::String, restore_year::Int)
    h5_path = joinpath(save_dir, "simulation_data_year_$(restore_year).h5")
    h5open(h5_path, "r") do f
        return load_simulation(f)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# CaseDefinition  (solver field not serialized — must be provided at load time)
# ─────────────────────────────────────────────────────────────────────────────

function save_case_definition!(g::HDF5.Group, c::CaseDefinition)
    write(g, "name", get_name(c))
    write(g, "base_dir", get_base_dir(c))
    write(g, "sys_dir", get_sys_dir(c))
    write(g, "scratch_dir", get_scratch_dir(c))
    write(g, "outage_dir", get_outage_dir(c))
    write(g, "timeseries_data_dir", get_timeseries_data_dir(c))
    write(g, "siip_market_clearing", get_siip_market_clearing(c))
    # solver is not serialized — MOI.OptimizerWithAttributes contains Any-typed fields
    write(g, "pcm_scenario", get_pcm_scenario(c))
    write(g, "start_year", get_start_year(c))
    write(g, "total_horizon", get_total_horizon(c))
    write(g, "rolling_horizon", get_rolling_horizon(c))
    write(g, "simulation_years", get_simulation_years(c))
    write(g, "rep_period_interval", get_rep_period_interval(c))
    write(g, "num_rep_periods", get_num_rep_periods(c))
    write(g, "avg_block_size", get_avg_block_size(c))
    write(g, "fixed_block_size", get_fixed_block_size(c))
    write(g, "rep_chronology_checkpoint", get_rep_chronology_checkpoint(c))
    write(g, "da_resolution", get_da_resolution(c))
    write(g, "rt_resolution", get_rt_resolution(c))
    write(g, "rps_target", get_rps_target(c))
    save_dict_symbol_bool!(create_group(g, "markets"), get_markets(c))
    write(g, "ordc_curved", get_ordc_curved(c))
    write(g, "ordc_unavailability_method", get_ordc_unavailability_method(c))
    write(g, "reserve_penalty", get_reserve_penalty(c))
    write(g, "static_capacity_market", get_static_capacity_market(c))
    write(g, "irm_scalar", get_irm_scalar(c))
    write(g, "accreditation_methodology", get_accreditation_methodology(c))
    write(g, "accreditation_metric", get_accreditation_metric(c))
    write(g, "marginal_cc_switch", get_marginal_cc_switch(c))
    write(g, "derating_scale", get_derating_scale(c))
    write(g, "mopr", get_mopr(c))
    write(g, "battery_cap_mkt", get_battery_cap_mkt(c))
    write(g, "vre_reserves", get_vre_reserves(c))
    write(g, "heterogeneity", get_heterogeneity(c))
    write(g, "forecast_type", get_forecast_type(c))
    write(g, "max_carbon_tax_increase", get_max_carbon_tax_increase(c))
    write(g, "info_symmetry", get_info_symmetry(c))
    write(g, "belief_update", get_belief_update(c))
    write(g, "uncertainty", get_uncertainty(c))
    write(g, "risk_aversion", get_risk_aversion(c))
    write(g, "parallel_investors", get_parallel_investors(c))
    write(g, "parallel_scenarios", get_parallel_scenarios(c))
    write(g, "md_horizon", get_md_horizon(c))
    write(g, "md_interval", get_md_interval(c))
    write(g, "uc_horizon", get_uc_horizon(c))
    write(g, "uc_interval", get_uc_interval(c))
    write(g, "ed_horizon", get_ed_horizon(c))
    write(g, "ed_interval", get_ed_interval(c))
    write(g, "md_market", get_md_market(c))
    write(g, "single_stage", get_single_stage(c))
    write(g, "step_size", get_step_size(c))
end

function load_case_definition(g::HDF5.Group)
    return CaseDefinition(
        read(g, "name"),
        read(g, "base_dir"),
        read(g, "sys_dir"),
        read(g, "scratch_dir"),
        read(g, "outage_dir"),
        read(g, "timeseries_data_dir"),
        nothing,  # solver — not serialized; caller must re-inject after loading
        read(g, "siip_market_clearing"),
        read(g, "pcm_scenario"),
        read(g, "start_year"),
        read(g, "total_horizon"),
        read(g, "rolling_horizon"),
        read(g, "simulation_years"),
        read(g, "rep_period_interval"),
        read(g, "num_rep_periods"),
        read(g, "avg_block_size"),
        read(g, "fixed_block_size"),
        read(g, "rep_chronology_checkpoint"),
        read(g, "da_resolution"),
        read(g, "rt_resolution"),
        read(g, "rps_target"),
        load_dict_symbol_bool(g["markets"]),
        read(g, "ordc_curved"),
        read(g, "ordc_unavailability_method"),
        read(g, "reserve_penalty"),
        read(g, "static_capacity_market"),
        read(g, "irm_scalar"),
        read(g, "accreditation_methodology"),
        read(g, "accreditation_metric"),
        read(g, "marginal_cc_switch"),
        read(g, "derating_scale"),
        read(g, "mopr"),
        read(g, "battery_cap_mkt"),
        read(g, "vre_reserves"),
        read(g, "heterogeneity"),
        read(g, "forecast_type"),
        read(g, "max_carbon_tax_increase"),
        read(g, "info_symmetry"),
        read(g, "belief_update"),
        read(g, "uncertainty"),
        read(g, "risk_aversion"),
        read(g, "parallel_investors"),
        read(g, "parallel_scenarios"),
        read(g, "md_horizon"),
        read(g, "md_interval"),
        read(g, "uc_horizon"),
        read(g, "uc_interval"),
        read(g, "ed_horizon"),
        read(g, "ed_interval"),
        read(g, "md_market"),
        read(g, "single_stage"),
        read(g, "step_size"),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# AgentSimulation
# ─────────────────────────────────────────────────────────────────────────────
function save_simulation!(f::Union{HDF5.File, HDF5.Group}, sim::AgentSimulation)
    @info "Saving AgentSimulation ($(length(sim.investors)) investors)"
    save_case_definition!(create_group(f, "case"), sim.case)
    write(f, "results_dir", sim.results_dir)
    write(f, "iteration_year", sim.iteration_year)
    write(f, "rep_period_interval", sim.rep_period_interval)
    write(f, "carbon_tax", sim.carbon_tax)
    write(f, "rec_requirement", sim.rec_requirement)
    write(f, "zones", sim.zones)

    @info "  Saving simulation metadata (markets, lines, rep_periods, weights, derating, RA)"
    save_dict_symbol_bool!(create_group(f, "markets"), sim.markets)
    save_lines!(create_group(f, "lines"), sim.lines)
    save_rep_periods!(create_group(f, "rep_periods"), sim.rep_periods)
    save_nested_dict_vf!(create_group(f, "hour_weight"), sim.hour_weight)
    save_nested_dict_f!(create_group(f, "peak_load"), sim.peak_load)
    save_derating_data!(create_group(f, "derating_data"), sim.derating_data)
    save_resource_adequacy_dict!(
        create_group(f, "resource_adequacy"),
        sim.resource_adequacy,
    )

    inv_g = create_group(f, "investors")
    for investor in sim.investors
        @info "  Saving investor: $(investor.name) ($(length(investor.projects)) projects)"
        save_investor!(create_group(inv_g, investor.name), investor)
    end
    @info "Simulation save complete"
end

function load_simulation(f::Union{HDF5.File, HDF5.Group})
    @info "Loading AgentSimulation"
    case = load_case_definition(f["case"])
    results_dir = read(f, "results_dir")
    iteration_year = read(f, "iteration_year")
    rep_period_interval = read(f, "rep_period_interval")
    carbon_tax = read(f, "carbon_tax")
    rec_requirement = read(f, "rec_requirement")
    zones = read(f, "zones")

    @info "  Loading simulation metadata (markets, lines, rep_periods, weights, derating, RA)"
    markets = load_dict_symbol_bool(f["markets"])
    lines = load_lines(f["lines"])
    rep_periods = load_rep_periods(f["rep_periods"])
    hour_weight = load_nested_dict_vf(f["hour_weight"])
    peak_load = load_nested_dict_f(f["peak_load"])
    derating_data = load_derating_data(f["derating_data"])
    resource_adequacy = load_resource_adequacy_dict(f["resource_adequacy"])

    investors = Investor[]
    inv_keys = keys(f["investors"])
    @info "  Loading $(length(inv_keys)) investor(s)"
    for k in inv_keys
        @info "  Loading investor: $k"
        push!(investors, load_investor(f["investors"][k]))
    end
    @info "Simulation load complete"

    return AgentSimulation(
        case,
        results_dir,
        iteration_year,
        nothing,                    # system_MDs — reload via PSY.to_json
        nothing,                    # system_UCs
        nothing,                    # system_EDs
        Dict{String, PSY.System}(), # system_PRAS
        zones,
        lines,
        rep_periods,
        rep_period_interval,
        hour_weight,
        peak_load,
        markets,
        carbon_tax,
        rec_requirement,
        investors,
        derating_data,
        resource_adequacy,
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Lines  (Vector{ZonalLine} stored as parallel arrays)
# ─────────────────────────────────────────────────────────────────────────────

function save_lines!(g::HDF5.Group, lines::Vector{ZonalLine})
    isempty(lines) && return
    write(g, "names", [l.name for l in lines])
    write(g, "from_zone", [l.from_zone for l in lines])
    write(g, "to_zone", [l.to_zone for l in lines])
    write(g, "active_power_limit", [l.active_power_limit for l in lines])
end

function load_lines(g::HDF5.Group)
    isempty(keys(g)) && return ZonalLine[]
    names = read(g, "names")
    from = read(g, "from_zone")
    to = read(g, "to_zone")
    limits = read(g, "active_power_limit")
    return [ZonalLine(names[i], from[i], to[i], limits[i]) for i in eachindex(names)]
end

# ─────────────────────────────────────────────────────────────────────────────
# rep_periods  Dict{String, Dict{Int64, Dict{Int64,Int64}}}
# ─────────────────────────────────────────────────────────────────────────────

function save_rep_periods!(g::HDF5.Group, rp::Dict)
    for (scen, yr_dict) in rp
        sg = create_group(g, scen)
        for (yr, mapping) in yr_dict
            yg = create_group(sg, string(yr))
            ks = collect(Int64, keys(mapping))
            vs = [mapping[k] for k in ks]
            write(yg, "keys", ks)
            write(yg, "values", vs)
        end
    end
end

function load_rep_periods(g::HDF5.Group)
    out = Dict{String, Dict{Int64, Dict{Int64, Int64}}}()
    for scen in keys(g)
        out[scen] = Dict{Int64, Dict{Int64, Int64}}()
        for yr_str in keys(g[scen])
            yr = parse(Int64, yr_str)
            yg = g[scen][yr_str]
            ks = read(yg, "keys")
            vs = read(yg, "values")
            out[scen][yr] = Dict(ks[i] => vs[i] for i in eachindex(ks))
        end
    end
    return out
end

# ─────────────────────────────────────────────────────────────────────────────
# hour_weight  Dict{String, Dict{Int64, Vector{Float64}}}
# ─────────────────────────────────────────────────────────────────────────────

function save_nested_dict_vf!(g::HDF5.Group, d::Dict{String, Dict{Int64, Vector{Float64}}})
    for (scen, yr_dict) in d
        sg = create_group(g, scen)
        for (yr, v) in yr_dict
            write(sg, string(yr), v)
        end
    end
end

# JLD2 may restore rep_hour_weight as Vector{Vector{Float64}} (flattened CEM block weights).
function save_nested_dict_vf!(g::HDF5.Group, vv::Vector{<:AbstractVector{<:Real}})
    attributes(g)["format"] = "vector_of_vectors"
    write(g, "n", length(vv))
    for (i, v) in enumerate(vv)
        write(g, string(i), Vector{Float64}(v))
    end
end

function load_nested_dict_vf(g::HDF5.Group)
    if haskey(HDF5.attributes(g), "format") && read_attribute(g, "format") == "vector_of_vectors"
        n = read(g, "n")
        return [read(g, string(i)) for i in 1:n]
    end
    out = Dict{String, Dict{Int64, Vector{Float64}}}()
    for scen in keys(g)
        out[scen] = Dict{Int64, Vector{Float64}}()
        for yr_str in keys(g[scen])
            out[scen][parse(Int64, yr_str)] = read(g[scen], yr_str)
        end
    end
    return out
end

# ─────────────────────────────────────────────────────────────────────────────
# peak_load  Dict{String, Dict{Int64, Float64}}
# ─────────────────────────────────────────────────────────────────────────────

function save_nested_dict_f!(g::HDF5.Group, d::Dict{String, Dict{Int64, Float64}})
    for (scen, yr_dict) in d
        sg = create_group(g, scen)
        ks = collect(Int64, keys(yr_dict))
        vs = [yr_dict[k] for k in ks]
        write(sg, "keys", ks)
        write(sg, "values", vs)
    end
end

function load_nested_dict_f(g::HDF5.Group)
    out = Dict{String, Dict{Int64, Float64}}()
    for scen in keys(g)
        ks = read(g[scen], "keys")
        vs = read(g[scen], "values")
        out[scen] = Dict(ks[i] => vs[i] for i in eachindex(ks))
    end
    return out
end

# ─────────────────────────────────────────────────────────────────────────────
# markets  Dict{Symbol, Bool}
# ─────────────────────────────────────────────────────────────────────────────

function save_dict_symbol_bool!(g::HDF5.Group, d::Dict{Symbol, Bool})
    ks = [string(k) for k in keys(d)]
    vs = [d[Symbol(k)] for k in ks]
    write(g, "keys", ks)
    write(g, "values", vs)
end

function load_dict_symbol_bool(g::HDF5.Group)
    ks = read(g, "keys")
    vs = read(g, "values")
    return Dict(Symbol(ks[i]) => vs[i] for i in eachindex(ks))
end

# ─────────────────────────────────────────────────────────────────────────────
# derating_data  Dict{String, DataFrame}
# ─────────────────────────────────────────────────────────────────────────────

function save_derating_data!(g::HDF5.Group, d::Dict{String, DataFrames.DataFrame})
    for (tech, df) in d
        save_dataframe!(create_group(g, tech), df)
    end
end

function load_derating_data(g::HDF5.Group)
    out = Dict{String, DataFrames.DataFrame}()
    for tech in keys(g)
        out[tech] = load_dataframe(g[tech])
    end
    return out
end

# ─────────────────────────────────────────────────────────────────────────────
# ResourceAdequacy  Dict{String, ResourceAdequacy}
# ─────────────────────────────────────────────────────────────────────────────

function save_resource_adequacy_dict!(g::HDF5.Group, d::Dict{String, ResourceAdequacy})
    for (scen, ra) in d
        save_resource_adequacy!(create_group(g, scen), ra)
    end
end

function load_resource_adequacy_dict(g::HDF5.Group)
    out = Dict{String, ResourceAdequacy}()
    for scen in keys(g)
        out[scen] = load_resource_adequacy(g[scen])
    end
    return out
end

function save_resource_adequacy!(g::HDF5.Group, ra::ResourceAdequacy)
    target_keys = collect(String, keys(ra.targets))
    target_values = [ra.targets[k] for k in target_keys]
    write(g, "target_keys", target_keys)
    write(g, "target_values", target_values)
    write(g, "delta_irm", ra.delta_irm)
    mg = create_group(g, "metrics")
    for (i, m) in enumerate(ra.metrics)
        yg = create_group(mg, string(i))
        ks = collect(String, keys(m))
        vs = [m[k] for k in ks]
        write(yg, "keys", ks)
        write(yg, "values", vs)
    end
end

function load_resource_adequacy(g::HDF5.Group)
    target_keys = read(g, "target_keys")
    target_values = read(g, "target_values")
    targets = Dict(target_keys[i] => target_values[i] for i in eachindex(target_keys))
    delta_irm = read(g, "delta_irm")
    mg = g["metrics"]
    n = length(keys(mg))
    metrics = [
        begin
            yg = mg[string(i)]
            ks = read(yg, "keys")
            vs = read(yg, "values")
            Dict(ks[j] => vs[j] for j in eachindex(ks))
        end for i in 1:n
    ]
    return ResourceAdequacy(targets, delta_irm, metrics)
end

# ─────────────────────────────────────────────────────────────────────────────
# DataFrame
# ─────────────────────────────────────────────────────────────────────────────

function save_dataframe!(g::HDF5.Group, df::DataFrames.DataFrame)
    col_names = string.(names(df))
    write(g, "col_names", col_names)
    for col in col_names
        data = df[!, col]
        # HDF5 needs concrete element type; convert if needed
        if eltype(data) <: AbstractString
            write(g, col, Vector{String}(data))
        elseif eltype(data) <: Integer
            write(g, col, Vector{Int64}(data))
        elseif eltype(data) <: AbstractFloat
            write(g, col, Vector{Float64}(data))
        elseif eltype(data) <: Bool
            write(g, col, Vector{Bool}(data))
        else
            write(g, col, string.(data))   # fallback: stringify
        end
    end
end

function load_dataframe(g::HDF5.Group)
    col_names = read(g, "col_names")
    df = DataFrames.DataFrame()
    for col in col_names
        DataFrames.insertcols!(df, col => read(g, col))
    end
    return df
end

# ─────────────────────────────────────────────────────────────────────────────
# AxisArray  (generic helper)
# ─────────────────────────────────────────────────────────────────────────────

function save_axisarray!(g::HDF5.Group, a::AxisArrays.AxisArray)
    write(g, "data", collect(a))
    axes_vals = AxisArrays.axisvalues(a)
    axes_names = AxisArrays.axisnames(a)
    write(g, "axis_count", length(axes_vals))
    for (i, (ax_name, ax_vals)) in enumerate(zip(axes_names, axes_vals))
        write(g, "axis$(i)_name", string(ax_name))
        if eltype(ax_vals) <: Symbol
            write(g, "axis$(i)_values", string.(ax_vals))
            write(g, "axis$(i)_type", "Symbol")
        elseif eltype(ax_vals) <: Integer
            write(g, "axis$(i)_values", collect(Int64, ax_vals))
            write(g, "axis$(i)_type", "Int64")
        elseif eltype(ax_vals) <: AbstractFloat
            write(g, "axis$(i)_values", collect(Float64, ax_vals))
            write(g, "axis$(i)_type", "Float64")
        else
            write(g, "axis$(i)_values", string.(ax_vals))
            write(g, "axis$(i)_type", "String")
        end
    end
end

function load_axisarray(g::HDF5.Group)
    data = read(g, "data")
    axis_count = read(g, "axis_count")
    axes = []
    for i in 1:axis_count
        ax_name = Symbol(read(g, "axis$(i)_name"))
        ax_type = read(g, "axis$(i)_type")
        ax_raw = read(g, "axis$(i)_values")
        ax_vals = if ax_type == "Symbol"
            Symbol.(ax_raw)
        elseif ax_type == "Int64"
            collect(Int64, ax_raw)
        elseif ax_type == "Float64"
            collect(Float64, ax_raw)
        else
            ax_raw
        end
        push!(axes, AxisArrays.Axis{ax_name}(ax_vals))
    end
    return AxisArrays.AxisArray(data, axes...)
end

# ─────────────────────────────────────────────────────────────────────────────
# Investor
# ─────────────────────────────────────────────────────────────────────────────

function save_investor!(g::HDF5.Group, inv::Investor)
    @info "    Saving investor data for $(inv.name)"
    write(g, "name", inv.name)
    write(g, "data_dir", inv.data_dir)
    write(g, "avg_block_size", inv.avg_block_size)
    write(g, "fixed_block_size", inv.fixed_block_size)
    write(g, "rep_period_interval", inv.rep_period_interval)
    write(g, "cap_cost_multiplier", inv.cap_cost_multiplier)
    write(g, "max_annual_projects", inv.max_annual_projects)
    write(g, "retirement_lookback", inv.retirement_lookback)
    write(g, "carbon_tax", inv.carbon_tax)
    write(g, "markets", string.(inv.markets))

    # NamedTuple
    pmr = create_group(g, "preference_multiplier_range")
    write(pmr, "min", inv.preference_multiplier_range.min)
    write(pmr, "max", inv.preference_multiplier_range.max)

    save_risk_preference!(create_group(g, "risk_preference"), inv.risk_preference)
    save_forecast!(create_group(g, "forecast"), inv.forecast)
    save_nested_dict_vf!(create_group(g, "rep_hour_weight"), inv.rep_hour_weight)
    save_chron_weights!(create_group(g, "chron_weights"), inv.chron_weights)
    save_portfolio_preference_multipliers!(
        create_group(g, "portfolio_preference_multipliers"),
        inv.portfolio_preference_multipliers,
    )
    save_market_prices!(create_group(g, "market_prices"), inv.market_prices)

    proj_g = create_group(g, "projects")
    for proj in inv.projects
        @info "      Saving project: $(proj.name) ($(typeof(proj)))"
        save_project!(create_group(proj_g, proj.name), proj)
    end
end

function load_investor(g::HDF5.Group)
    @info "    Loading investor data"
    name = read(g, "name")
    data_dir = read(g, "data_dir")
    avg_block_size = read(g, "avg_block_size")
    fixed_block_size = read(g, "fixed_block_size")
    rep_period_interval = read(g, "rep_period_interval")
    cap_cost_multiplier = read(g, "cap_cost_multiplier")
    max_annual_projects = read(g, "max_annual_projects")
    retirement_lookback = read(g, "retirement_lookback")
    carbon_tax = read(g, "carbon_tax")
    markets = Symbol.(read(g, "markets"))

    pmr_g = g["preference_multiplier_range"]
    preference_multiplier_range = (min = read(pmr_g, "min"), max = read(pmr_g, "max"))

    risk_preference = load_risk_preference(g["risk_preference"])
    forecast = load_forecast(g["forecast"])
    rep_hour_weight = load_nested_dict_vf(g["rep_hour_weight"])
    chron_weights = load_chron_weights(g["chron_weights"])
    portfolio_preference_multipliers =
        load_portfolio_preference_multipliers(g["portfolio_preference_multipliers"])
    market_prices = load_market_prices(g["market_prices"])

    projects = Vector{Project{<:BuildPhase}}()
    proj_keys = keys(g["projects"])
    @info "      Loading $(length(proj_keys)) project(s) for investor: $name"
    for k in proj_keys
        @info "      Loading project: $k"
        push!(projects, load_project(g["projects"][k]))
    end

    return Investor(
        name, data_dir, projects, markets, carbon_tax,
        market_prices, rep_period_interval, rep_hour_weight,
        avg_block_size, fixed_block_size, chron_weights,
        forecast, cap_cost_multiplier, preference_multiplier_range,
        portfolio_preference_multipliers, max_annual_projects,
        risk_preference, retirement_lookback,
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# chron_weights  Dict{String, Dict{Int64, Matrix{Int64}}}
# ─────────────────────────────────────────────────────────────────────────────

function save_chron_weights!(g::HDF5.Group, d::Dict{String, Dict{Int64, Matrix{Int64}}})
    for (scen, yr_dict) in d
        sg = create_group(g, scen)
        for (yr, mat) in yr_dict
            write(sg, string(yr), mat)
        end
    end
end

function load_chron_weights(g::HDF5.Group)
    out = Dict{String, Dict{Int64, Matrix{Int64}}}()
    for scen in keys(g)
        out[scen] = Dict{Int64, Matrix{Int64}}()
        for yr_str in keys(g[scen])
            out[scen][parse(Int64, yr_str)] = read(g[scen], yr_str)
        end
    end
    return out
end

# ─────────────────────────────────────────────────────────────────────────────
# portfolio_preference_multipliers  Dict{Tuple{String,String}, Vector{Float64}}
# ─────────────────────────────────────────────────────────────────────────────

function save_portfolio_preference_multipliers!(
    g::HDF5.Group,
    d::Dict{Tuple{String, String}, Vector{Float64}},
)
    isempty(d) && return
    entries = collect(d)
    key_a = [e[1][1] for e in entries]
    key_b = [e[1][2] for e in entries]
    write(g, "key_a", key_a)
    write(g, "key_b", key_b)
    vals_g = create_group(g, "values")
    for (i, e) in enumerate(entries)
        write(vals_g, string(i), e[2])
    end
end

function load_portfolio_preference_multipliers(g::HDF5.Group)
    out = Dict{Tuple{String, String}, Vector{Float64}}()
    isempty(keys(g)) && return out
    key_a = read(g, "key_a")
    key_b = read(g, "key_b")
    vals_g = g["values"]
    for i in eachindex(key_a)
        out[(key_a[i], key_b[i])] = read(vals_g, string(i))
    end
    return out
end

# ─────────────────────────────────────────────────────────────────────────────
# RiskPreference
# ─────────────────────────────────────────────────────────────────────────────

function save_risk_preference!(g::HDF5.Group, r::RiskNeutral)
    attributes(g)["type"] = "RiskNeutral"
end

function save_risk_preference!(g::HDF5.Group, r::RiskAverse)
    attributes(g)["type"] = "RiskAverse"
    write(g, "constant", r.constant)
    write(g, "multiplier", r.multiplier)
    write(g, "risk_coefficient", r.risk_coefficient)
end

function load_risk_preference(g::HDF5.Group)
    t = read_attribute(g, "type")
    if t == "RiskNeutral"
        return RiskNeutral()
    elseif t == "RiskAverse"
        return RiskAverse(
            read(g, "constant"),
            read(g, "multiplier"),
            read(g, "risk_coefficient"),
        )
    else
        error("Unknown RiskPreference type: $t")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Forecast  (Perfect / Imperfect)
# ─────────────────────────────────────────────────────────────────────────────

function save_forecast!(g::HDF5.Group, f::Perfect)
    attributes(g)["type"] = "Perfect"
    sg = create_group(g, "scenarios")
    for s in f.scenario_data
        save_scenario!(create_group(sg, s.name), s)
    end
end

function save_forecast!(g::HDF5.Group, f::Imperfect)
    attributes(g)["type"] = "Imperfect"
    if isnothing(f.kalman_filter)
        attributes(g)["kalman_filter_is_nothing"] = true
    else
        attributes(g)["kalman_filter_is_nothing"] = false
        save_kalman_filter!(create_group(g, "kalman_filter"), f.kalman_filter)
    end
    sg = create_group(g, "scenarios")
    for s in f.scenario_data
        save_scenario!(create_group(sg, s.name), s)
    end
end

function load_forecast(g::HDF5.Group)
    t = read_attribute(g, "type")
    sg = g["scenarios"]
    scenarios = [load_scenario(sg[k]) for k in keys(sg)]
    if t == "Perfect"
        return Perfect(scenarios)
    elseif t == "Imperfect"
        kf_nothing = read_attribute(g, "kalman_filter_is_nothing")
        kf = kf_nothing ? nothing : load_kalman_filter(g["kalman_filter"])
        return Imperfect(kf, scenarios)
    else
        error("Unknown Forecast type: $t")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Scenario
# ─────────────────────────────────────────────────────────────────────────────

function save_scenario!(g::HDF5.Group, s::Scenario)
    write(g, "name", s.name)
    write(g, "probability", s.probability)

    pm_g = create_group(g, "parameter_multipliers")
    if isnothing(s.parameter_multipliers)
        attributes(pm_g)["is_nothing"] = true
    else
        attributes(pm_g)["is_nothing"] = false
        ks = collect(String, keys(s.parameter_multipliers))
        vs = [s.parameter_multipliers[k] for k in ks]
        write(pm_g, "keys", ks)
        write(pm_g, "values", vs)
    end

    pv_g = create_group(g, "parameter_values")
    for (i, ax) in enumerate(s.parameter_values)
        save_axisarray!(create_group(pv_g, string(i)), ax)
    end
end

function load_scenario(g::HDF5.Group)
    name = read(g, "name")
    probability = read(g, "probability")

    pm_g = g["parameter_multipliers"]
    pm_none = read_attribute(pm_g, "is_nothing")
    parameter_multipliers = if pm_none
        nothing
    else
        ks = read(pm_g, "keys")
        vs = read(pm_g, "values")
        Dict(ks[i] => vs[i] for i in eachindex(ks))
    end

    pv_g = g["parameter_values"]
    n = length(keys(pv_g))
    parameter_values = [load_axisarray(pv_g[string(i)]) for i in 1:n]

    return Scenario(name, probability, parameter_multipliers, parameter_values)
end

# ─────────────────────────────────────────────────────────────────────────────
# KalmanFilter / InvestorBelief
# ─────────────────────────────────────────────────────────────────────────────

function save_kalman_filter!(g::HDF5.Group, kf::KalmanFilter)
    ib_g = create_group(g, "investor_belief")
    save_axisarray!(
        create_group(ib_g, "process_covariance"),
        kf.investor_data.process_covariance,
    )
    save_axisarray!(
        create_group(ib_g, "measurement_covariance"),
        kf.investor_data.measurement_covariance,
    )
    write(ib_g, "state_transition_matrix", kf.investor_data.state_transition_matrix)
    write(ib_g, "state_measurement_matrix", kf.investor_data.state_measurement_matrix)

    save_axisarray!(create_group(g, "state_estimate"), kf.state_estimate)
    save_axisarray!(
        create_group(g, "error_covariance_estimate"),
        kf.error_covariance_estimate,
    )
end

function load_kalman_filter(g::HDF5.Group)
    ib_g = g["investor_belief"]
    investor_data = InvestorBelief(
        load_axisarray(ib_g["process_covariance"]),
        load_axisarray(ib_g["measurement_covariance"]),
        read(ib_g, "state_transition_matrix"),
        read(ib_g, "state_measurement_matrix"),
    )
    state_estimate = load_axisarray(g["state_estimate"])
    error_covariance_estimate = load_axisarray(g["error_covariance_estimate"])
    return KalmanFilter(investor_data, state_estimate, error_covariance_estimate)
end

# ─────────────────────────────────────────────────────────────────────────────
# MarketPrices
# ─────────────────────────────────────────────────────────────────────────────

function save_market_prices!(g::HDF5.Group, mp::MarketPrices)
    _save_optional_axis_dict!(create_group(g, "energy_price"), mp.energy_price)
    _save_optional_reserve_price!(create_group(g, "reserve_price"), mp.reserve_price)
    _save_optional_axis_dict!(create_group(g, "capacity_price"), mp.capacity_price)
    _save_optional_axis_dict!(create_group(g, "rec_price"), mp.rec_price)
    _save_optional_axis_dict!(create_group(g, "inertia_price"), mp.inertia_price)
end

# Dict{String, AxisArray} — used by energy_price, capacity_price, rec_price, inertia_price
function _save_optional_axis_dict!(g::HDF5.Group, d)
    if isnothing(d)
        attributes(g)["is_nothing"] = true
        return
    end
    attributes(g)["is_nothing"] = false
    for (scen, ax) in d
        save_axisarray!(create_group(g, scen), ax)
    end
end

function _load_optional_axis_dict(g::HDF5.Group)
    read_attribute(g, "is_nothing") && return nothing
    return Dict(k => load_axisarray(g[k]) for k in keys(g) if k != "is_nothing")
end

# reserve_price: Dict{String, Dict{String, Array{Float64,2}}}
function _save_optional_reserve_price!(g::HDF5.Group, d)
    if isnothing(d)
        attributes(g)["is_nothing"] = true
        return
    end
    attributes(g)["is_nothing"] = false
    for (prod, scen_dict) in d
        pg = create_group(g, prod)
        for (scen, arr) in scen_dict
            write(pg, scen, arr)
        end
    end
end

function _load_optional_reserve_price(g::HDF5.Group)
    read_attribute(g, "is_nothing") && return nothing
    out = Dict{String, Dict{String, Array{Float64, 2}}}()
    for prod in keys(g)
        prod == "is_nothing" && continue
        out[prod] = Dict(scen => read(g[prod], scen) for scen in keys(g[prod]))
    end
    return out
end

function load_market_prices(g::HDF5.Group)
    return MarketPrices(
        _load_optional_axis_dict(g["energy_price"]),
        _load_optional_reserve_price(g["reserve_price"]),
        _load_optional_axis_dict(g["capacity_price"]),
        _load_optional_axis_dict(g["rec_price"]),
        _load_optional_axis_dict(g["inertia_price"]),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Project{BuildPhase}  (dispatches on concrete type)
# ─────────────────────────────────────────────────────────────────────────────

# Map type name string → BuildPhase abstract type for parameterisation
const _BUILDPHASE_MAP = Dict(
    "Existing" => Existing, "Planned" => Planned, "Queue" => Queue,
    "Option" => Option, "Retired" => Retired,
)

function save_project!(g::HDF5.Group, p::ThermalGenEMIS{T}) where {T <: BuildPhase}
    attributes(g)["project_type"] = "ThermalGenEMIS"
    attributes(g)["build_phase"] = string(T)
    _save_project_common!(g, p)
    save_tech!(create_group(g, "tech"), p.tech)
end

function save_project!(g::HDF5.Group, p::RenewableGenEMIS{T}) where {T <: BuildPhase}
    attributes(g)["project_type"] = "RenewableGenEMIS"
    attributes(g)["build_phase"] = string(T)
    _save_project_common!(g, p)
    save_tech!(create_group(g, "tech"), p.tech)
end

function save_project!(g::HDF5.Group, p::HydroGenEMIS{T}) where {T <: BuildPhase}
    attributes(g)["project_type"] = "HydroGenEMIS"
    attributes(g)["build_phase"] = string(T)
    _save_project_common!(g, p)
    save_tech!(create_group(g, "tech"), p.tech)
end

function save_project!(g::HDF5.Group, p::BatteryEMIS{T}) where {T <: BuildPhase}
    attributes(g)["project_type"] = "BatteryEMIS"
    attributes(g)["build_phase"] = string(T)
    _save_project_common!(g, p)
    save_tech!(create_group(g, "tech"), p.tech)
end

function _save_project_common!(g::HDF5.Group, p)
    write(g, "name", p.name)
    write(g, "decision_year", p.decision_year)
    write(g, "construction_year", p.construction_year)
    write(g, "retirement_year", p.retirement_year)
    write(g, "end_life_year", p.end_life_year)
    save_finance!(create_group(g, "finance_data"), p.finance_data)
    prod_g = create_group(g, "products")
    for prod in p.products
        save_product!(create_group(prod_g, string(prod.name)), prod)
    end
end

function load_project(g::HDF5.Group)
    pt = read_attribute(g, "project_type")
    bp = read_attribute(g, "build_phase")
    # Strip module path if present (e.g. "EMISAgentSimulation.Existing" → "Existing")
    bp_short = split(bp, ".")[end]
    phase = _BUILDPHASE_MAP[bp_short]

    name = read(g, "name")
    decision_year = read(g, "decision_year")
    construction_year = read(g, "construction_year")
    retirement_year = read(g, "retirement_year")
    end_life_year = read(g, "end_life_year")
    finance_data = load_finance(g["finance_data"])
    products = [load_product(g["products"][k]) for k in keys(g["products"])]
    tech = load_tech(g["tech"])

    if pt == "ThermalGenEMIS"
        return ThermalGenEMIS{phase}(
            name,
            tech,
            decision_year,
            construction_year,
            retirement_year,
            end_life_year,
            products,
            finance_data,
        )
    elseif pt == "RenewableGenEMIS"
        return RenewableGenEMIS{phase}(
            name,
            tech,
            decision_year,
            construction_year,
            retirement_year,
            end_life_year,
            products,
            finance_data,
        )
    elseif pt == "HydroGenEMIS"
        return HydroGenEMIS{phase}(
            name,
            tech,
            decision_year,
            construction_year,
            retirement_year,
            end_life_year,
            products,
            finance_data,
        )
    elseif pt == "BatteryEMIS"
        return BatteryEMIS{phase}(
            name,
            tech,
            decision_year,
            construction_year,
            retirement_year,
            end_life_year,
            products,
            finance_data,
        )
    else
        error("Unknown project_type: $pt")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Tech  (ThermalTech / RenewableTech / HydroTech / BatteryTech)
# ─────────────────────────────────────────────────────────────────────────────

function save_tech!(g::HDF5.Group, t::ThermalTech)
    attributes(g)["tech_type"] = "ThermalTech"
    _save_tech_common!(g, t)
    write(g, "fuel", t.fuel)
    write(g, "fuel_cost", t.fuel_cost)
    _save_optional_named_tuple_ud!(create_group(g, "time_limits"), t.time_limits)
    _save_heat_rate_curve!(create_group(g, "heat_rate_curve"), t.heat_rate_curve)
    save_operation_cost!(create_group(g, "operation_cost"), t.operation_cost)
end

function save_tech!(g::HDF5.Group, t::RenewableTech)
    attributes(g)["tech_type"] = "RenewableTech"
    _save_tech_common!(g, t)
    save_operation_cost!(create_group(g, "operation_cost"), t.operation_cost)
end

function save_tech!(g::HDF5.Group, t::HydroTech)
    attributes(g)["tech_type"] = "HydroTech"
    _save_tech_common!(g, t)
    _save_optional_named_tuple_ud!(create_group(g, "time_limits"), t.time_limits)
    save_operation_cost!(create_group(g, "operation_cost"), t.operation_cost)
end

function save_tech!(g::HDF5.Group, t::BatteryTech)
    attributes(g)["tech_type"] = "BatteryTech"
    # BatteryTech has no active_power_limits — write shared fields directly
    write(g, "type", t.type)
    write(g, "bus", t.bus)
    write(g, "zone", t.zone)
    write(g, "FOR", t.FOR)
    write(g, "MTTR", t.MTTR)
    _save_optional_named_tuple_ud!(create_group(g, "ramp_limits"), t.ramp_limits)
    _save_named_tuple_minmax!(
        create_group(g, "input_active_power_limits"),
        t.input_active_power_limits,
    )
    _save_named_tuple_minmax!(
        create_group(g, "output_active_power_limits"),
        t.output_active_power_limits,
    )
    _save_named_tuple_minmax!(create_group(g, "storage_capacity"), t.storage_capacity)
    _save_named_tuple_minmax!(
        create_group(g, "storage_level_limits"),
        t.storage_level_limits,
    )
    write(g, "initial_storage_capacity_level", t.initial_storage_capacity_level)
    write(g, "rating", t.rating)
    write(g, "soc", t.soc)
    write(g, "base_power", t.base_power)
    eff_g = create_group(g, "efficiency")
    write(eff_g, "in", t.efficiency.in)
    write(eff_g, "out", t.efficiency.out)
end

function _save_tech_common!(g::HDF5.Group, t)
    write(g, "type", t.type)
    write(g, "bus", t.bus)
    write(g, "zone", t.zone)
    write(g, "FOR", t.FOR)
    write(g, "MTTR", t.MTTR)
    _save_named_tuple_minmax!(create_group(g, "active_power_limits"), t.active_power_limits)
    _save_optional_named_tuple_ud!(create_group(g, "ramp_limits"), t.ramp_limits)
end

function _save_named_tuple_minmax!(g::HDF5.Group, nt::NamedTuple)
    write(g, "min", nt.min)
    write(g, "max", nt.max)
end

function _save_optional_named_tuple_ud!(g::HDF5.Group, nt::Nothing)
    attributes(g)["is_nothing"] = true
end
function _save_optional_named_tuple_ud!(g::HDF5.Group, nt::NamedTuple)
    attributes(g)["is_nothing"] = false
    write(g, "up", nt.up)
    write(g, "down", nt.down)
end

function _save_heat_rate_curve!(g::HDF5.Group, curve::Vector{Tuple{Float64, Float64}})
    write(g, "x", [c[1] for c in curve])
    write(g, "y", [c[2] for c in curve])
end

function load_tech(g::HDF5.Group)
    tt = read_attribute(g, "tech_type")
    type = read(g, "type")
    bus = read(g, "bus")
    zone = read(g, "zone")
    FOR = read(g, "FOR")
    MTTR = read(g, "MTTR")
    rl = _load_optional_named_tuple_ud(g["ramp_limits"])

    if tt == "ThermalTech"
        apl = (
            min = read(g["active_power_limits"], "min"),
            max = read(g["active_power_limits"], "max"),
        )
        fuel = read(g, "fuel")
        fuel_cost = read(g, "fuel_cost")
        tl = _load_optional_named_tuple_ud(g["time_limits"])
        hrc = _load_heat_rate_curve(g["heat_rate_curve"])
        op_cost = load_operation_cost(g["operation_cost"])
        return ThermalTech(
            type,
            fuel,
            apl,
            rl,
            tl,
            op_cost,
            fuel_cost,
            hrc,
            bus,
            zone,
            FOR,
            MTTR,
        )

    elseif tt == "RenewableTech"
        apl = (
            min = read(g["active_power_limits"], "min"),
            max = read(g["active_power_limits"], "max"),
        )
        op_cost = load_operation_cost(g["operation_cost"])
        return RenewableTech(type, apl, rl, op_cost, bus, zone, FOR, MTTR)

    elseif tt == "HydroTech"
        apl = (
            min = read(g["active_power_limits"], "min"),
            max = read(g["active_power_limits"], "max"),
        )
        tl = _load_optional_named_tuple_ud(g["time_limits"])
        op_cost = load_operation_cost(g["operation_cost"])
        return HydroTech(type, apl, rl, tl, op_cost, bus, zone, FOR, MTTR)

    elseif tt == "BatteryTech"
        iapl = (
            min = read(g["input_active_power_limits"], "min"),
            max = read(g["input_active_power_limits"], "max"),
        )
        oapl = (
            min = read(g["output_active_power_limits"], "min"),
            max = read(g["output_active_power_limits"], "max"),
        )
        sc = (
            min = read(g["storage_capacity"], "min"),
            max = read(g["storage_capacity"], "max"),
        )
        sll = (
            min = read(g["storage_level_limits"], "min"),
            max = read(g["storage_level_limits"], "max"),
        )
        iscl = read(g, "initial_storage_capacity_level")
        rating = read(g, "rating")
        soc = read(g, "soc")
        base_power = read(g, "base_power")
        eff = (in = read(g["efficiency"], "in"), out = read(g["efficiency"], "out"))
        return BatteryTech(
            type,
            iapl,
            oapl,
            rl,
            sc,
            sll,
            iscl,
            rating,
            soc,
            eff,
            bus,
            zone,
            FOR,
            MTTR,
            base_power,
        )

    else
        error("Unknown tech_type: $tt")
    end
end

function _load_optional_named_tuple_ud(g::HDF5.Group)
    read_attribute(g, "is_nothing") && return nothing
    return (up = read(g, "up"), down = read(g, "down"))
end

function _load_heat_rate_curve(g::HDF5.Group)
    xs = read(g, "x")
    ys = read(g, "y")
    return [(xs[i], ys[i]) for i in eachindex(xs)]
end

# ─────────────────────────────────────────────────────────────────────────────
# operation_cost  (PSY *GenerationCost — save essential numeric fields only)
# ─────────────────────────────────────────────────────────────────────────────

function save_operation_cost!(g::HDF5.Group, oc::Nothing)
    attributes(g)["is_nothing"] = true
end

function save_operation_cost!(g::HDF5.Group, oc::PSY.ThermalGenerationCost)
    attributes(g)["is_nothing"] = false
    attributes(g)["cost_type"] = "ThermalGenerationCost"
    write(g, "fixed", oc.fixed)
    write(g, "shut_down", oc.shut_down)
    _save_psy_cost_curve!(create_group(g, "variable"), oc.variable)
    su_g = create_group(g, "start_up")
    if oc.start_up isa NamedTuple
        attributes(su_g)["is_named_tuple"] = true
        write(su_g, "hot", oc.start_up.hot)
        write(su_g, "warm", oc.start_up.warm)
        write(su_g, "cold", oc.start_up.cold)
    else
        attributes(su_g)["is_named_tuple"] = false
        write(su_g, "value", Float64(oc.start_up))
    end
end

function save_operation_cost!(g::HDF5.Group, oc::PSY.RenewableGenerationCost)
    attributes(g)["is_nothing"] = false
    attributes(g)["cost_type"] = "RenewableGenerationCost"
    write(g, "fixed", oc.fixed)
    _save_psy_cost_curve!(create_group(g, "variable"), oc.variable)
    _save_psy_cost_curve!(create_group(g, "curtailment_cost"), oc.curtailment_cost)
end

function save_operation_cost!(g::HDF5.Group, oc::PSY.HydroGenerationCost)
    attributes(g)["is_nothing"] = false
    attributes(g)["cost_type"] = "HydroGenerationCost"
    write(g, "fixed", oc.fixed)
    _save_psy_cost_curve!(create_group(g, "variable"), oc.variable)
end

function _save_psy_cost_curve!(g::HDF5.Group, curve)
    if isnothing(curve)
        attributes(g)["is_nothing"] = true
        return
    end
    attributes(g)["is_nothing"] = false
    vc = PSY.get_value_curve(curve)
    if vc isa PSY.LinearCurve
        attributes(g)["curve_kind"] = "Linear"
        write(g, "proportional_term", PSY.get_proportional_term(vc))
        write(g, "constant_term", PSY.get_constant_term(vc))
    else
        # Try to extract piecewise breakpoints; if we can't get ≥2, fall back to
        # zero-cost LinearCurve so the load path stays unambiguous.
        pts = nothing
        try
            extracted = PSY.get_breakpoint_vars(curve)
            length(extracted) >= 2 && (pts = extracted)
        catch
        end
        if !isnothing(pts)
            attributes(g)["curve_kind"] = "Piecewise"
            write(g, "x", Float64[p[1] for p in pts])
            write(g, "y", Float64[p[2] for p in pts])
        else
            attributes(g)["curve_kind"] = "Linear"
            write(g, "proportional_term", 0.0)
            write(g, "constant_term", 0.0)
            @warn "Could not extract ≥2 breakpoints from PSY cost curve; saved as zero-cost LinearCurve"
        end
    end
end

function load_operation_cost(g::HDF5.Group)
    read_attribute(g, "is_nothing") && return nothing
    ct = read_attribute(g, "cost_type")
    if ct == "ThermalGenerationCost"
        fixed = read(g, "fixed")
        shut_down = read(g, "shut_down")
        variable = _load_psy_cost_curve(g["variable"])
        su_g = g["start_up"]
        start_up = if read_attribute(su_g, "is_named_tuple")
            (hot = read(su_g, "hot"), warm = read(su_g, "warm"), cold = read(su_g, "cold"))
        else
            read(su_g, "value")
        end
        return PSY.ThermalGenerationCost(;
            variable = variable,
            fixed = fixed,
            start_up = start_up,
            shut_down = shut_down,
        )
    elseif ct == "RenewableGenerationCost"
        fixed = read(g, "fixed")
        variable = _load_psy_cost_curve(g["variable"])
        curtailment_cost = _load_psy_cost_curve(g["curtailment_cost"])
        return PSY.RenewableGenerationCost(;
            variable = variable,
            curtailment_cost = curtailment_cost,
            fixed = fixed,
        )
    elseif ct == "HydroGenerationCost"
        fixed = read(g, "fixed")
        variable = _load_psy_cost_curve(g["variable"])
        return PSY.HydroGenerationCost(; variable = variable, fixed = fixed)
    else
        error("Unknown cost_type: $ct")
    end
end

function _load_psy_cost_curve(g::HDF5.Group)
    read_attribute(g, "is_nothing") && return nothing
    # curve_kind absent in files saved before this fix → assume Piecewise
    curve_kind = if haskey(HDF5.attributes(g), "curve_kind")
        read_attribute(g, "curve_kind")
    else
        "Piecewise"
    end
    if curve_kind == "Linear"
        prop = read(g, "proportional_term")
        const_term = read(g, "constant_term")
        return PSY.CostCurve(PSY.LinearCurve(prop, const_term))
    else
        xs = read(g, "x")
        ys = read(g, "y")
        # Old files may have stored only 1 point → fall back to zero-cost LinearCurve
        length(xs) < 2 && return PSY.CostCurve(PSY.LinearCurve(0.0, 0.0))
        pts = [(xs[i], ys[i]) for i in eachindex(xs)]
        # PiecewiseLinearData is function data; wrap in InputOutputCurve to get a ValueCurve
        return PSY.CostCurve(PSY.InputOutputCurve(PSY.PiecewiseLinearData(pts)))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Finance
# ─────────────────────────────────────────────────────────────────────────────

function save_finance!(g::HDF5.Group, f::Finance)
    write(g, "investment_cost", f.investment_cost)
    write(g, "effective_investment_cost", f.effective_investment_cost)
    write(g, "preference_multiplier", f.preference_multiplier)
    write(g, "lag_time", f.lag_time)
    write(g, "life_time", f.life_time)
    write(g, "capex_years", f.capex_years)
    write(g, "fixed_OM_cost", f.fixed_OM_cost)
    write(g, "queue_cost", f.queue_cost)
    write(g, "discount_rate", f.discount_rate)
    write(g, "expected_npv", f.expected_npv)
    write(g, "expected_utility", f.expected_utility)
    write(g, "annual_cashflow", f.annual_cashflow)
    write(g, "ownedby", f.ownedby)

    stu_g = create_group(g, "scenario_total_utilization")
    for (scen, arr) in f.scenario_total_utilization
        write(stu_g, scen, arr)
    end

    snpv_g = create_group(g, "scenario_npv")
    for (scen, v) in f.scenario_npv
        write(snpv_g, scen, v)
    end

    su_g = create_group(g, "scenario_utility")
    for (scen, v) in f.scenario_utility
        write(su_g, scen, v)
    end

    sp_g = create_group(g, "scenario_profit")
    for (scen, yr_vec) in f.scenario_profit
        sg = create_group(sp_g, scen)
        for (i, ax) in enumerate(yr_vec)
            save_axisarray!(create_group(sg, string(i)), ax)
        end
    end

    save_axisarray!(create_group(g, "realized_profit"), f.realized_profit)
end

function load_finance(g::HDF5.Group)
    investment_cost = read(g, "investment_cost")
    effective_investment_cost = read(g, "effective_investment_cost")
    preference_multiplier = read(g, "preference_multiplier")
    lag_time = read(g, "lag_time")
    life_time = read(g, "life_time")
    capex_years = read(g, "capex_years")
    fixed_OM_cost = read(g, "fixed_OM_cost")
    queue_cost = read(g, "queue_cost")
    discount_rate = read(g, "discount_rate")
    expected_npv = read(g, "expected_npv")
    expected_utility = read(g, "expected_utility")
    annual_cashflow = read(g, "annual_cashflow")
    ownedby = read(g, "ownedby")

    stu_g = g["scenario_total_utilization"]
    scenario_total_utilization = Dict(k => read(stu_g, k) for k in keys(stu_g))

    snpv_g = g["scenario_npv"]
    scenario_npv = Dict(k => read(snpv_g, k) for k in keys(snpv_g))

    su_g = g["scenario_utility"]
    scenario_utility = Dict(k => read(su_g, k) for k in keys(su_g))

    sp_g = g["scenario_profit"]
    scenario_profit = Dict{String, Vector{AxisArrays.AxisArray{Float64, 2}}}()
    for scen in keys(sp_g)
        n = length(keys(sp_g[scen]))
        scenario_profit[scen] = [load_axisarray(sp_g[scen][string(i)]) for i in 1:n]
    end

    realized_profit = load_axisarray(g["realized_profit"])

    return Finance(
        investment_cost, effective_investment_cost, preference_multiplier,
        lag_time, life_time, capex_years, fixed_OM_cost, queue_cost,
        scenario_total_utilization, scenario_profit, realized_profit,
        discount_rate, scenario_npv, expected_npv, scenario_utility,
        expected_utility, annual_cashflow, ownedby,
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Products  (Energy / Capacity / OperatingReserve / REC / Inertia)
# ─────────────────────────────────────────────────────────────────────────────

function save_product!(g::HDF5.Group, p::Energy)
    attributes(g)["product_type"] = "Energy"
    write(g, "name", string(p.name))
    write(g, "marginal_cost", p.marginal_cost)
    write(g, "expected_production", p.expected_production)
    cf_g = create_group(g, "capacity_factors")
    for (scen, arr) in p.capacity_factors
        write(cf_g, scen, arr)
    end
end

function save_product!(g::HDF5.Group, p::Capacity)
    attributes(g)["product_type"] = "Capacity"
    write(g, "name", string(p.name))
    write(g, "capacity_bid", p.capacity_bid)
    der_g = create_group(g, "derating")
    ks = collect(String, keys(p.derating))
    vs = [p.derating[k] for k in ks]
    write(der_g, "keys", ks)
    write(der_g, "values", vs)
    ap_g = create_group(g, "accepted_perc")
    for (scen, v) in p.accepted_perc
        write(ap_g, scen, v)
    end
end

function save_product!(g::HDF5.Group, p::OperatingReserve{T}) where {T}
    attributes(g)["product_type"] = "OperatingReserve{$(T)}"
    write(g, "name", string(p.name))
    write(g, "max_limit", p.max_limit)
    write(g, "marginal_cost", p.marginal_cost)
end

function save_product!(g::HDF5.Group, p::CarbonTax)
    attributes(g)["product_type"] = "CarbonTax"
    write(g, "name", string(p.name))
    write(g, "emission_intensity", p.emission_intensity)
    write(g, "avg_heat_rate", p.avg_heat_rate)
    write(g, "fuel_cost", p.fuel_cost)
    write(g, "total_emission", p.total_emission)
end

function save_product!(g::HDF5.Group, p::REC)
    attributes(g)["product_type"] = "REC"
    write(g, "name", string(p.name))
    write(g, "expected_certificates", p.expected_certificates)
    write(g, "correction_factor", p.correction_factor)
    write(g, "rec_bid", p.rec_bid)
end

function save_product!(g::HDF5.Group, p::Inertia)
    attributes(g)["product_type"] = "Inertia"
    write(g, "name", string(p.name))
    write(g, "synchronous", p.synchronous)
    write(g, "h_constant", p.h_constant)
    write(g, "marginal_cost", p.marginal_cost)
end

function load_product(g::HDF5.Group)
    pt = read_attribute(g, "product_type")

    if pt == "Energy"
        name = Symbol(read(g, "name"))
        mc = read(g, "marginal_cost")
        ep = read(g, "expected_production")
        cf_g = g["capacity_factors"]
        cf = Dict(k => read(cf_g, k) for k in keys(cf_g))
        return Energy(name, cf, mc, ep)

    elseif pt == "Capacity"
        name = Symbol(read(g, "name"))
        bid = read(g, "capacity_bid")
        der_g = g["derating"]
        ks = read(der_g, "keys")
        vs = read(der_g, "values")
        der = Dict(ks[i] => vs[i] for i in eachindex(ks))
        ap_g = g["accepted_perc"]
        ap = Dict(k => read(ap_g, k) for k in keys(ap_g))
        return Capacity(name, der, ap, bid)

    elseif startswith(pt, "OperatingReserve{")
        name = Symbol(read(g, "name"))
        ml = read(g, "max_limit")
        mc = read(g, "marginal_cost")
        if occursin("ReserveUpEMIS", pt)
            return OperatingReserve{ReserveUpEMIS}(name, ml, mc)
        else
            return OperatingReserve{ReserveDownEMIS}(name, ml, mc)
        end

    elseif pt == "REC"
        name = Symbol(read(g, "name"))
        ec = read(g, "expected_certificates")
        cf = read(g, "correction_factor")
        rb = read(g, "rec_bid")
        return REC(name, ec, cf, rb)

    elseif pt == "Inertia"
        name = Symbol(read(g, "name"))
        sync = read(g, "synchronous")
        hc = read(g, "h_constant")
        mc = read(g, "marginal_cost")
        return Inertia(name, sync, hc, mc)

    elseif pt == "CarbonTax"
        name = Symbol(read(g, "name"))
        emission_intensity = read(g, "emission_intensity")
        avg_heat_rate = read(g, "avg_heat_rate")
        fuel_cost = read(g, "fuel_cost")
        total_emission = read(g, "total_emission")
        return CarbonTax(name, emission_intensity, avg_heat_rate, fuel_cost, total_emission)

    else
        error("Unknown product_type: $pt")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Save Sienna Systems (MD, UC, ED, PRAS) to files or checkpoints
# ─────────────────────────────────────────────────────────────────────────────
function save_Sienna_systems(simulation::AgentSimulation, result_path::String, iteration_year::Int)

    scenario_names = String.(get_all_scenario_names(get_data_dir(get_case(simulation))))
    # Save from iteration_year onward only. Past years were already saved at their own
    # checkpoints and are immutable (finish_construction! only modifies construction_year:end,
    # never retroactively changes past-year systems).
    all_years = length(simulation.system_MDs)
    for year in iteration_year:all_years
        @info "Saving systems for year $(year) to checkpoint."
        PSY.to_json(simulation.system_MDs[year], joinpath(result_path, "sys_MD_year$(year).json"), force = true)
        PSY.to_json(simulation.system_UCs[year], joinpath(result_path, "sys_UC_year$(year).json"), force = true)
        PSY.to_json(simulation.system_EDs[year], joinpath(result_path, "sys_ED_year$(year).json"), force = true)
    end
    # PRAS is a single object per scenario (not year-indexed). Only save for iteration_year;
    # load_sienna_systems! is updated to load only restore_year instead of 1:restore_year.
    for scenario in scenario_names
        PSY.to_json(simulation.system_PRAS[scenario],
        joinpath(get_results_dir(simulation), "sys_PRAS_$(scenario)_year$(iteration_year).json"), force = true)
    end
end


# ─────────────────────────────────────────────────────────────────────────────
# Load Sienna Systems (MD, UC, ED, PRAS) from files or checkpoints
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# Derating factors  (DataFrame → HDF5)
# ─────────────────────────────────────────────────────────────────────────────
function save_derating_factors(path::String, df::DataFrames.DataFrame)
    isdir(dirname(path)) || mkpath(dirname(path))
    h5open(path, "w") do f
        attributes(f)["schema_version"] = SCHEMA_VERSION
        save_dataframe!(create_group(f, "derating_factors"), df)
    end
end

function load_derating_factors(path::String)
    h5open(path, "r") do f
        return load_dataframe(f["derating_factors"])
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Clean energy percentage vector  (Vector{Float64} → HDF5)
# ─────────────────────────────────────────────────────────────────────────────
function save_clean_energy_percentage(path::String, v::Vector{Float64})
    isdir(dirname(path)) || mkpath(dirname(path))
    h5open(path, "w") do f
        attributes(f)["schema_version"] = SCHEMA_VERSION
        write(f, "clean_energy_percentage", v)
    end
end

function load_clean_energy_percentage(path::String)
    h5open(path, "r") do f
        return read(f, "clean_energy_percentage")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Shortfall data  (PRAS ShortfallResult → HDF5, write-only)
# Saves all Float64 fields from a PRAS ShortfallResult.
# Timestamps are stored as ISO8601 strings; no round-trip load is needed.
# ─────────────────────────────────────────────────────────────────────────────
function save_shortfall_data(path::String, sf)
    isdir(dirname(path)) || mkpath(dirname(path))
    h5open(path, "w") do f
        attributes(f)["schema_version"] = SCHEMA_VERSION
        g = create_group(f, "shortfall")
        if isnothing(sf.nsamples)
            attributes(g)["nsamples_is_nothing"] = true
        else
            attributes(g)["nsamples_is_nothing"] = false
            write(g, "nsamples", sf.nsamples)
        end
        write(g, "region_names", collect(String, sf.regions.names))
        write(g, "timestamps", string.(collect(sf.timestamps)))
        write(g, "eventperiod_mean", sf.eventperiod_mean)
        write(g, "eventperiod_std", sf.eventperiod_std)
        write(g, "shortfall_std", sf.shortfall_std)
        write(g, "eventperiod_region_mean", sf.eventperiod_region_mean)
        write(g, "eventperiod_region_std", sf.eventperiod_region_std)
        write(g, "eventperiod_period_mean", sf.eventperiod_period_mean)
        write(g, "eventperiod_period_std", sf.eventperiod_period_std)
        write(g, "shortfall_region_std", sf.shortfall_region_std)
        write(g, "shortfall_period_std", sf.shortfall_period_std)
        write(g, "eventperiod_regionperiod_mean", sf.eventperiod_regionperiod_mean)
        write(g, "eventperiod_regionperiod_std", sf.eventperiod_regionperiod_std)
        write(g, "shortfall_mean", sf.shortfall_mean)
        write(g, "shortfall_regionperiod_std", sf.shortfall_regionperiod_std)
    end
end


function load_sienna_systems!(simulation::AgentSimulation, result_path::String, restore_year::Int, simulation_settings::Dict, scratch_dir::String)

    simulation.system_MDs = Vector{PSY.System}()
    simulation.system_UCs = Vector{PSY.System}()
    simulation.system_EDs = Vector{PSY.System}()
    simulation.system_PRAS = Dict{String, PSY.System}()
    case = get_case(simulation)
    data_dir = get_data_dir(case);
    scenarios = string.(get_all_scenario_names(data_dir));

    for year in 1:restore_year
        @info "Restoring systems for year $(year) from checkpoint."
        sys_MD = PSY.System(joinpath(result_path, "sys_MD_year$(year).json"), runchecks = false)
        push!(simulation.system_MDs, sys_MD)

        sys_UC = PSY.System(joinpath(result_path, "sys_UC_year$(year).json"), runchecks = false)
        push!(simulation.system_UCs, sys_UC)

        sys_ED = PSY.System(joinpath(result_path, "sys_ED_year$(year).json"), runchecks = false)
        push!(simulation.system_EDs, sys_ED)
    end
    # PRAS is not year-indexed; only the restore_year file exists and is needed.
    for scenario in scenarios
        @info "Restoring PRAS system for scenario $(scenario) from year $(restore_year) checkpoint."
        sys_PRAS = PSY.System(joinpath(result_path, "sys_PRAS_$(scenario)_year$(restore_year).json"), runchecks = false)
        simulation.system_PRAS[scenario] = sys_PRAS
    end

    simulation_years = get_total_horizon(case)
    timeseries_data_dir = get_timeseries_data_dir(case)
    rts_dir = get_sys_dir(case)
    ntp_ts_data_dir = joinpath(timeseries_data_dir, "input_processing")
    runchecks = false
    MD_horizon = get_md_horizon(case)
    MD_interval = get_md_interval(case)
    UC_horizon = get_uc_horizon(case)
    UC_interval = get_uc_interval(case)
    ED_horizon = get_ed_horizon(case)
    ED_interval = get_ed_interval(case)
    outage_dir = get_outage_dir(case)

    initial_Sienna_system_name = "DA_sys_zonal_with_storage_capacities.json"
    pcm_scenario = simulation_settings["pcm_scenario"]

    if pcm_scenario == "scenario_1"
        supercc_scenario = "baseline"
    elseif pcm_scenario == "scenario_2"
        supercc_scenario = "central"
    elseif pcm_scenario == "scenario_3"
        supercc_scenario = "ira"
    else
        "Not a pre-defined scenario."
    end

    for sim_year in restore_year+1:simulation_years
        md_json = joinpath(result_path, "sys_MD_year$(sim_year).json")
        uc_json = joinpath(result_path, "sys_UC_year$(sim_year).json")
        ed_json = joinpath(result_path, "sys_ED_year$(sim_year).json")

        if isfile(md_json) && isfile(uc_json) && isfile(ed_json)
            # Load from checkpoint: these systems already have investor-constructed devices
            # from finish_construction! calls in years 1..restore_year, plus all post-load
            # modifications (fix_multistart_cost_curves!, component removals, renaming, units).
            @info "Restoring systems for year $(sim_year) from checkpoint."
            sys_MD = PSY.System(md_json, time_series_directory = scratch_dir, runchecks = runchecks)
            sys_UC = PSY.System(uc_json, time_series_directory = scratch_dir, runchecks = runchecks)
            sys_ED = PSY.System(ed_json, time_series_directory = scratch_dir, runchecks = runchecks)
        else
            # Fallback: no checkpoint for this year — load base system and re-apply modifications.
            @info "No checkpoint found for year $(sim_year), reading from constructed_systems."
            MD_sys_filename = joinpath(rts_dir, "constructed_systems", pcm_scenario,
            "sim_year_$(sim_year)", "MD_sys_EMIS_$(MD_horizon)hor_$(MD_interval)int.json")
            UC_filename = joinpath(rts_dir, "constructed_systems", pcm_scenario,
            "sim_year_$(sim_year)",
            "DA_sys_EMIS_$(UC_horizon)hor_$(UC_interval)int_$(MD_horizon)mdhor_$(MD_interval)mdint.json")
            ED_filename = joinpath(rts_dir, "constructed_systems", pcm_scenario,
            "sim_year_$(sim_year)", "RT_sys_EMIS_$(ED_horizon)hor_$(ED_interval)int_$(MD_horizon)mdhor_$(MD_interval)mdint.json")

            sys_MD = PSY.System(MD_sys_filename, time_series_directory = scratch_dir, runchecks = runchecks)
            sys_UC = PSY.System(UC_filename, time_series_directory = scratch_dir, runchecks = runchecks)
            sys_ED = PSY.System(ED_filename, time_series_directory = scratch_dir, runchecks = runchecks)

            # fix_multistart_cost_curves!(sys_MD)
            # fix_multistart_cost_curves!(sys_UC)
            # fix_multistart_cost_curves!(sys_ED)

            removegen_name = ["AUSTIN_1","AUSTIN_2"]
            for sys in [sys_MD, sys_UC, sys_ED]
                for d in PSY.get_components(PSY.Generator, sys)
                    d.name in removegen_name && PSY.remove_component!(sys, d)
                end
                PSY.remove_component!(sys, PSY.get_component(PSY.VariableReserve, sys, "SPIN"))
                PSY.remove_component!(sys, PSY.get_component(PSY.VariableReserveNonSpinning, sys, "NONSPIN"))
                PSY.set_name!(sys, PSY.get_component(PSY.VariableReserve, sys, "REG_DN"), "Reg_Down")
                PSY.set_name!(sys, PSY.get_component(PSY.VariableReserve, sys, "REG_UP"), "Reg_Up")
                PSY.set_units_base_system!(sys, PSY.IS.UnitSystem.DEVICE_BASE)
            end
        end

        push!(simulation.system_MDs, sys_MD)
        push!(simulation.system_UCs, sys_UC)
        push!(simulation.system_EDs, sys_ED)
    end

    # Propagate ThermalMultiStart end-of-year-restore_year states into year restore_year+1 fresh systems
    if restore_year + 1 <= length(simulation.system_UCs)
        for gen_prev in PSY.get_components(PSY.ThermalMultiStart, simulation.system_UCs[restore_year])
            name = PSY.get_name(gen_prev)
            gen_next = PSY.get_component(PSY.ThermalMultiStart, simulation.system_UCs[restore_year + 1], name)
            if gen_next !== nothing
                PSY.set_status!(gen_next, PSY.get_status(gen_prev))
                PSY.set_time_at_status!(gen_next, PSY.get_time_at_status(gen_prev))
            end
        end
        for gen_prev in PSY.get_components(PSY.ThermalMultiStart, simulation.system_EDs[restore_year])
            name = PSY.get_name(gen_prev)
            gen_next = PSY.get_component(PSY.ThermalMultiStart, simulation.system_EDs[restore_year + 1], name)
            if gen_next !== nothing
                PSY.set_status!(gen_next, PSY.get_status(gen_prev))
                PSY.set_time_at_status!(gen_next, PSY.get_time_at_status(gen_prev))
            end
        end
        @info "Propagated ThermalMultiStart initial conditions from year $(restore_year) to year $(restore_year + 1)"
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# expected_market_data  (standalone .h5 file per scenario × year)
#
# Private helpers for types not covered by the existing save/load utilities:
#   _save_dict_str_float! / _load_dict_str_float   Dict{String, Float64}
#   _save_dict_str_int!   / _load_dict_str_int     Dict{String, Int64}
#   _save_dict_str_matrix!/ _load_dict_str_matrix  Dict{String, Matrix{Float64}}
#   _save_dict_str_vec_float!/_load_dict_str_vec_float Dict{String, Vector{Float64}}
#   _save_dict_str_vec_str! / _load_dict_str_vec_str   Dict{String, Vector{String}}
# ─────────────────────────────────────────────────────────────────────────────

function _save_dict_str_float!(g::HDF5.Group, d::AbstractDict)
    ks = collect(String, keys(d))
    vs = [Float64(d[k]) for k in ks]
    write(g, "keys", ks)
    write(g, "values", vs)
end

# JLD2 may restore parameter-vector AxisArrays (String axis, Float64 data).
function _save_dict_str_float!(g::HDF5.Group, a::AxisArrays.AxisArray{<:Real, 1})
    ks = string.(collect(AxisArrays.axisvalues(a)[1]))
    vs = Float64.(collect(a))
    write(g, "keys", ks)
    write(g, "values", vs)
end

function _load_dict_str_float(g::HDF5.Group)
    ks = read(g, "keys")
    vs = read(g, "values")
    return Dict(ks[i] => vs[i] for i in eachindex(ks))
end

function _save_dict_str_int!(g::HDF5.Group, d::AbstractDict)
    ks = collect(String, keys(d))
    vs = [Int64(d[k]) for k in ks]
    write(g, "keys", ks)
    write(g, "values", vs)
end

function _save_dict_str_int!(g::HDF5.Group, a::AxisArrays.AxisArray{<:Real, 1})
    ks = string.(collect(AxisArrays.axisvalues(a)[1]))
    vs = Int64.(collect(a))
    write(g, "keys", ks)
    write(g, "values", vs)
end

function _load_dict_str_int(g::HDF5.Group)
    ks = read(g, "keys")
    vs = read(g, "values")
    return Dict(ks[i] => vs[i] for i in eachindex(ks))
end

# Each key maps to a 2-D Float64 matrix stored as a named dataset inside g.
function _save_dict_str_matrix!(g::HDF5.Group, d)
    ks = collect(String, keys(d))
    write(g, "keys", ks)
    for k in ks
        write(g, k, Matrix{Float64}(d[k]))
    end
end

function _load_dict_str_matrix(g::HDF5.Group)
    ks = read(g, "keys")
    return Dict(k => read(g, k) for k in ks)
end

# Each key maps to a 1-D Float64 vector stored as a named dataset inside g.
function _save_dict_str_vec_float!(g::HDF5.Group, d)
    ks = collect(String, keys(d))
    write(g, "keys", ks)
    for k in ks
        write(g, k, Vector{Float64}(d[k]))
    end
end

function _load_dict_str_vec_float(g::HDF5.Group)
    ks = read(g, "keys")
    return Dict(k => read(g, k) for k in ks)
end

# Each key maps to a Vector{String}; empty vectors written as a 0-element string array.
function _save_dict_str_vec_str!(g::HDF5.Group, d)
    ks = collect(String, keys(d))
    write(g, "keys", ks)
    for k in ks
        v = Vector{String}(d[k])
        if isempty(v)
            write(g, k, String[])
        else
            write(g, k, v)
        end
    end
end

function _load_dict_str_vec_str(g::HDF5.Group)
    ks = read(g, "keys")
    return Dict(k => read(g, k) for k in ks)
end

# AxisArray helper with an is_empty guard for 0-length first-axis arrays
# (e.g. p_in_ru_detail when storage_projects is empty).
function _save_axisarray_guarded!(g::HDF5.Group, a::AxisArrays.AxisArray)
    if isempty(a)
        attributes(g)["is_empty"] = true
        # Store axis metadata so load can reconstruct the shape.
        axes_names = AxisArrays.axisnames(a)
        axes_vals  = AxisArrays.axisvalues(a)
        write(g, "axis_count", length(axes_vals))
        for (i, (ax_name, ax_vals)) in enumerate(zip(axes_names, axes_vals))
            write(g, "axis$(i)_name", string(ax_name))
            if eltype(ax_vals) <: Symbol
                write(g, "axis$(i)_values", string.(collect(ax_vals)))
                write(g, "axis$(i)_type", "Symbol")
            elseif eltype(ax_vals) <: Integer
                write(g, "axis$(i)_values", collect(Int64, ax_vals))
                write(g, "axis$(i)_type", "Int64")
            else
                write(g, "axis$(i)_values", string.(collect(ax_vals)))
                write(g, "axis$(i)_type", "String")
            end
        end
        return
    end
    attributes(g)["is_empty"] = false
    save_axisarray!(g, a)
end

function _load_axisarray_guarded(g::HDF5.Group)
    if read_attribute(g, "is_empty")
        axis_count = read(g, "axis_count")
        axes = []
        sizes = Int[]
        for i in 1:axis_count
            ax_name = Symbol(read(g, "axis$(i)_name"))
            ax_type = read(g, "axis$(i)_type")
            ax_raw  = read(g, "axis$(i)_values")
            ax_vals = if ax_type == "Symbol"
                Symbol.(ax_raw)
            elseif ax_type == "Int64"
                collect(Int64, ax_raw)
            else
                ax_raw
            end
            push!(axes, AxisArrays.Axis{ax_name}(ax_vals))
            push!(sizes, length(ax_vals))
        end
        return AxisArrays.AxisArray(zeros(Float64, sizes...), axes...)
    end
    return load_axisarray(g)
end

# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

"""
    save_expected_market_data(path, <all 49 CEM output variables>)

Write the outputs of one CEM solve to a standalone HDF5 file at `path`.
The file is created (or overwritten) by this function; it is separate from
simulation_data_year_N.h5.
"""
function save_expected_market_data(path::String,
    capacity_price,
    energy_price,
    reserve_price,
    rec_price,
    inertia_price,
    capacity_factors,
    total_utilization,
    capacity_accepted_perc,
    new_options,
    new_options_by_type,
    REC_slack,
    REC_supply,
    REC_demand,
    rps_compliant_projects,
    rec_requirement,
    investment,
    retirement,
    max_new_options,
    max_gen,
    demand_e,
    rep_hour_weight,
    p_e_print,
    in_flow_print,
    out_flow_print,
    v_e_print,
    p_in_print,
    linepowerlimit,
    rec_correction,
    p_e_storage_print,
    zone_storage,
    zone_projects,
    lines_to_zone,
    lines_from_zone,
    init_storage,
    p_e_detail,
    remaining_buildtime,
    projects,
    generator_projects,
    storage_projects,
    reserve_up_products,
    ordc_products,
    reserve_down_products,
    p_in_ru_detail,
    p_in_ordc_detail,
    p_in_inertia_detail,
    p_in_rd_detail,
    p_out_ru_detail,
    p_out_ordc_detail,
    p_out_inertia_detail,
    p_out_rd_detail,
)
    isdir(dirname(path)) || mkpath(dirname(path))
    h5open(path, "w") do f
        attributes(f)["schema_version"] = SCHEMA_VERSION

        save_axisarray!(create_group(f, "capacity_price"),   capacity_price)
        save_axisarray!(create_group(f, "energy_price"),     energy_price)
        _save_dict_str_matrix!(create_group(f, "reserve_price"),         reserve_price)
        save_axisarray!(create_group(f, "rec_price"),        rec_price)
        save_axisarray!(create_group(f, "inertia_price"),    inertia_price)

        _save_dict_str_matrix!(create_group(f, "capacity_factors"),      capacity_factors)
        _save_dict_str_matrix!(create_group(f, "total_utilization"),     total_utilization)
        _save_dict_str_vec_float!(create_group(f, "capacity_accepted_perc"), capacity_accepted_perc)

        _save_dict_str_int!(create_group(f, "new_options"),              new_options)
        _save_dict_str_float!(create_group(f, "new_options_by_type"),    new_options_by_type)

        save_axisarray!(create_group(f, "REC_slack"),  REC_slack)
        save_axisarray!(create_group(f, "REC_supply"), REC_supply)
        save_axisarray!(create_group(f, "REC_demand"), REC_demand)

        write(f, "rps_compliant_projects", Vector{String}(rps_compliant_projects))
        save_axisarray!(create_group(f, "rec_requirement"), rec_requirement)

        save_axisarray!(create_group(f, "investment"), investment)
        save_axisarray!(create_group(f, "retirement"), retirement)

        _save_dict_str_float!(create_group(f, "max_new_options"),  max_new_options)
        _save_dict_str_float!(create_group(f, "max_gen"),          max_gen)

        save_axisarray!(create_group(f, "demand_e"), demand_e)
        save_nested_dict_vf!(create_group(f, "rep_hour_weight"), rep_hour_weight)

        save_axisarray!(create_group(f, "p_e_print"),         p_e_print)
        save_axisarray!(create_group(f, "in_flow_print"),     in_flow_print)
        save_axisarray!(create_group(f, "out_flow_print"),    out_flow_print)
        save_axisarray!(create_group(f, "v_e_print"),         v_e_print)
        save_axisarray!(create_group(f, "p_in_print"),        p_in_print)
        save_axisarray!(create_group(f, "p_e_storage_print"), p_e_storage_print)

        _save_dict_str_float!(create_group(f, "linepowerlimit"),  linepowerlimit)
        _save_dict_str_float!(create_group(f, "rec_correction"),  rec_correction)

        _save_dict_str_vec_str!(create_group(f, "zone_storage"),    zone_storage)
        _save_dict_str_vec_str!(create_group(f, "zone_projects"),   zone_projects)
        _save_dict_str_vec_str!(create_group(f, "lines_to_zone"),   lines_to_zone)
        _save_dict_str_vec_str!(create_group(f, "lines_from_zone"), lines_from_zone)

        _save_dict_str_float!(create_group(f, "init_storage"),      init_storage)
        save_axisarray!(create_group(f, "p_e_detail"), p_e_detail)
        _save_dict_str_float!(create_group(f, "remaining_buildtime"), remaining_buildtime)

        write(f, "projects",             Vector{String}(projects))
        write(f, "generator_projects",   Vector{String}(generator_projects))
        write(f, "storage_projects",     Vector{String}(storage_projects))
        write(f, "reserve_up_products",  Vector{String}(collect(reserve_up_products)))
        write(f, "ordc_products",        Vector{String}(collect(ordc_products)))
        write(f, "reserve_down_products",Vector{String}(collect(reserve_down_products)))

        _save_axisarray_guarded!(create_group(f, "p_in_ru_detail"),      p_in_ru_detail)
        _save_axisarray_guarded!(create_group(f, "p_in_ordc_detail"),    p_in_ordc_detail)
        _save_axisarray_guarded!(create_group(f, "p_in_inertia_detail"), p_in_inertia_detail)
        _save_axisarray_guarded!(create_group(f, "p_in_rd_detail"),      p_in_rd_detail)
        _save_axisarray_guarded!(create_group(f, "p_out_ru_detail"),     p_out_ru_detail)
        _save_axisarray_guarded!(create_group(f, "p_out_ordc_detail"),   p_out_ordc_detail)
        _save_axisarray_guarded!(create_group(f, "p_out_inertia_detail"),p_out_inertia_detail)
        _save_axisarray_guarded!(create_group(f, "p_out_rd_detail"),     p_out_rd_detail)
    end
end

"""
    load_expected_market_data(path) -> Dict{String, Any}

Read a CEM output file written by `save_expected_market_data`.
Returns a `Dict{String, Any}` with the same keys as the old JLD2 format so
that call sites in investor_iteration.jl need no changes.
"""
function load_expected_market_data(path::String)
    h5open(path, "r") do f
        return Dict{String, Any}(
            "capacity_price"          => load_axisarray(f["capacity_price"]),
            "energy_price"            => load_axisarray(f["energy_price"]),
            "reserve_price"           => _load_dict_str_matrix(f["reserve_price"]),
            "rec_price"               => load_axisarray(f["rec_price"]),
            "inertia_price"           => load_axisarray(f["inertia_price"]),
            "capacity_factors"        => _load_dict_str_matrix(f["capacity_factors"]),
            "total_utilization"       => _load_dict_str_matrix(f["total_utilization"]),
            "capacity_accepted_perc"  => _load_dict_str_vec_float(f["capacity_accepted_perc"]),
            "new_options"             => _load_dict_str_int(f["new_options"]),
            "new_options_by_type"     => _load_dict_str_float(f["new_options_by_type"]),
            "REC_slack"               => load_axisarray(f["REC_slack"]),
            "REC_supply"              => load_axisarray(f["REC_supply"]),
            "REC_demand"              => load_axisarray(f["REC_demand"]),
            "rps_compliant_projects"  => read(f, "rps_compliant_projects"),
            "rec_requirement"         => load_axisarray(f["rec_requirement"]),
            "investment"              => load_axisarray(f["investment"]),
            "retirement"              => load_axisarray(f["retirement"]),
            "max_new_options"         => _load_dict_str_float(f["max_new_options"]),
            "max_gen"                 => _load_dict_str_float(f["max_gen"]),
            "demand_e"                => load_axisarray(f["demand_e"]),
            "rep_hour_weight"         => load_nested_dict_vf(f["rep_hour_weight"]),
            "p_e_print"               => load_axisarray(f["p_e_print"]),
            "in_flow_print"           => load_axisarray(f["in_flow_print"]),
            "out_flow_print"          => load_axisarray(f["out_flow_print"]),
            "v_e_print"               => load_axisarray(f["v_e_print"]),
            "p_in_print"              => load_axisarray(f["p_in_print"]),
            "linepowerlimit"          => _load_dict_str_float(f["linepowerlimit"]),
            "rec_correction"          => _load_dict_str_float(f["rec_correction"]),
            "p_e_storage_print"       => load_axisarray(f["p_e_storage_print"]),
            "zone_storage"            => _load_dict_str_vec_str(f["zone_storage"]),
            "zone_projects"           => _load_dict_str_vec_str(f["zone_projects"]),
            "lines_to_zone"           => _load_dict_str_vec_str(f["lines_to_zone"]),
            "lines_from_zone"         => _load_dict_str_vec_str(f["lines_from_zone"]),
            "init_storage"            => _load_dict_str_float(f["init_storage"]),
            "p_e_detail"              => load_axisarray(f["p_e_detail"]),
            "remaining_buildtime"     => _load_dict_str_float(f["remaining_buildtime"]),
            "projects"                => read(f, "projects"),
            "generator_projects"      => read(f, "generator_projects"),
            "storage_projects"        => read(f, "storage_projects"),
            "reserve_up_products"     => read(f, "reserve_up_products"),
            "ordc_products"           => read(f, "ordc_products"),
            "reserve_down_products"   => read(f, "reserve_down_products"),
            "p_in_ru_detail"          => _load_axisarray_guarded(f["p_in_ru_detail"]),
            "p_in_ordc_detail"        => _load_axisarray_guarded(f["p_in_ordc_detail"]),
            "p_in_inertia_detail"     => _load_axisarray_guarded(f["p_in_inertia_detail"]),
            "p_in_rd_detail"          => _load_axisarray_guarded(f["p_in_rd_detail"]),
            "p_out_ru_detail"         => _load_axisarray_guarded(f["p_out_ru_detail"]),
            "p_out_ordc_detail"       => _load_axisarray_guarded(f["p_out_ordc_detail"]),
            "p_out_inertia_detail"    => _load_axisarray_guarded(f["p_out_inertia_detail"]),
            "p_out_rd_detail"         => _load_axisarray_guarded(f["p_out_rd_detail"]),
        )
    end
end




# investor_iteration.jl calls save_expected_market_data
# This code snippet transforms the jld2 to new h5 format
# Note: this is a one-time utility function to convert existing JLD2 files to the new HDF5 format.
# base_path = "/projects/gmlcmarkets/Phase2_EMIS_Analysis/GS_AAYAD/EMIS_RTS_Analysis_GS/20250310_no_sdes_High_RECT_Static_ORDC_RA_Cap_wo_md_storff_High_RPS/investors"

function transform_jld2_to_h5(base_path::String)

    for investor_id in 1:4
        investor_name = "investor$(investor_id)"
        for scenario_name in ["scenario_1", "scenario_2", "scenario_3"]
            for iteration_year in 1:15
                @info "Loading OLD format expected market data for $(investor_name) iteration year $(iteration_year)"
                data_path = joinpath(base_path, investor_name, "expected_market_data", "$(scenario_name)_year_$(iteration_year).jld2")
                expected_data = FileIO.load(data_path)
                
                capacity_price = expected_data["capacity_price"]
                energy_price = expected_data["energy_price"]
                reserve_price = expected_data["reserve_price"]
                rec_price = expected_data["rec_price"]
                inertia_price = expected_data["inertia_price"]
                capacity_factors = expected_data["capacity_factors"]
                total_utilization = expected_data["total_utilization"]
                capacity_accepted_perc = expected_data["capacity_accepted_perc"]
                new_options = expected_data["new_options"]
                new_options_by_type = expected_data["new_options_by_type"]
                REC_slack = expected_data["REC_slack"]
                REC_supply = expected_data["REC_supply"]
                REC_demand = expected_data["REC_demand"]
                rps_compliant_projects = expected_data["rps_compliant_projects"]
                rec_requirement = expected_data["rec_requirement"]
                investment = expected_data["investment"]
                retirement = expected_data["retirement"]
                max_new_options = expected_data["max_new_options"]
                max_gen = expected_data["max_gen"]
                demand_e = expected_data["demand_e"]
                rep_hour_weight = expected_data["rep_hour_weight"]
                p_e_print = expected_data["p_e_print"]
                in_flow_print = expected_data["in_flow_print"]
                out_flow_print = expected_data["out_flow_print"]
                v_e_print = expected_data["v_e_print"]
                p_in_print = expected_data["p_in_print"]
                linepowerlimit = expected_data["linepowerlimit"]
                rec_correction = expected_data["rec_correction"]    
                p_e_storage_print = expected_data["p_e_storage_print"]
                zone_storage = expected_data["zone_storage"]
                zone_projects = expected_data["zone_projects"]
                lines_to_zone = expected_data["lines_to_zone"]
                lines_from_zone = expected_data["lines_from_zone"]
                init_storage = expected_data["init_storage"]
                p_e_detail = expected_data["p_e_detail"]
                remaining_buildtime = expected_data["remaining_buildtime"]
                projects = expected_data["projects"]
                generator_projects = expected_data["generator_projects"]
                storage_projects = expected_data["storage_projects"]
                reserve_up_products = expected_data["reserve_up_products"]
                ordc_products = expected_data["ordc_products"]
                reserve_down_products = expected_data["reserve_down_products"]
                p_in_ru_detail = expected_data["p_in_ru_detail"]
                p_in_ordc_detail = expected_data["p_in_ordc_detail"]
                p_in_inertia_detail = expected_data["p_in_inertia_detail"]
                p_in_rd_detail = expected_data["p_in_rd_detail"]
                p_out_ru_detail = expected_data["p_out_ru_detail"]
                p_out_ordc_detail = expected_data["p_out_ordc_detail"]
                p_out_inertia_detail = expected_data["p_out_inertia_detail"]
                p_out_rd_detail = expected_data["p_out_rd_detail"]

                @info "Saving NEW format expected market data for $(investor_name) iteration year $(iteration_year)"
                h5_data_path = joinpath(base_path, investor_name, "expected_market_data", "$(scenario_name)_year_$(iteration_year).h5")
                save_expected_market_data(h5_data_path,
                capacity_price, energy_price, reserve_price, rec_price, inertia_price,
                capacity_factors, total_utilization, capacity_accepted_perc,
                new_options, new_options_by_type,
                REC_slack, REC_supply, REC_demand,
                rps_compliant_projects, rec_requirement,
                investment, retirement,
                max_new_options, max_gen, demand_e, rep_hour_weight,
                p_e_print, in_flow_print, out_flow_print, v_e_print, p_in_print,
                linepowerlimit, rec_correction, p_e_storage_print,
                zone_storage, zone_projects, lines_to_zone, lines_from_zone,
                init_storage, p_e_detail, remaining_buildtime,
                projects, generator_projects, storage_projects,
                reserve_up_products, ordc_products, reserve_down_products,
                p_in_ru_detail, p_in_ordc_detail, p_in_inertia_detail, p_in_rd_detail,
                p_out_ru_detail, p_out_ordc_detail, p_out_inertia_detail, p_out_rd_detail,
            )
            end
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# realized_market_data  (standalone .h5 file per simulation year)
#
# Additional private helpers:
#   _save_nested_dict_str_matrix! / _load_nested_dict_str_matrix
#       Dict{String, Dict{String, Matrix{Float64}}}  — reserve_perc_{md/uc/ed}
#   _save_supply_curve! / _load_supply_curve
#       Vector{Vector{Union{String,Float64}}} with rows [name, val1, val2]
# ─────────────────────────────────────────────────────────────────────────────

# Outer key → inner key → 2-D Float64 matrix.
function _save_nested_dict_str_matrix!(g::HDF5.Group, d)
    outer_keys = collect(String, keys(d))
    write(g, "outer_keys", outer_keys)
    for ok in outer_keys
        og = create_group(g, ok)
        inner = d[ok]
        inner_keys = collect(String, keys(inner))
        write(og, "inner_keys", inner_keys)
        for ik in inner_keys
            write(og, ik, Matrix{Float64}(inner[ik]))
        end
    end
end

function _load_nested_dict_str_matrix(g::HDF5.Group)
    outer_keys = read(g, "outer_keys")
    return Dict(
        ok => begin
            og = g[ok]
            inner_keys = read(og, "inner_keys")
            Dict(ik => read(og, ik) for ik in inner_keys)
        end
        for ok in outer_keys
    )
end

# rec_supply_curve rows: [name::String, val1::Float64, val2::Float64]
function _save_supply_curve!(g::HDF5.Group, sc::Vector)
    n = length(sc)
    write(g, "n", n)
    if n == 0
        return
    end
    write(g, "names",  String[row[1] for row in sc])
    ncols = length(sc[1]) - 1   # number of float columns
    write(g, "ncols", ncols)
    for c in 1:ncols
        write(g, "col$(c)", Float64[row[c + 1] for row in sc])
    end
end

function _load_supply_curve(g::HDF5.Group)
    n = read(g, "n")
    n == 0 && return Vector{Union{String, Float64}}[]
    names = read(g, "names")
    ncols = read(g, "ncols")
    cols = [read(g, "col$(c)") for c in 1:ncols]
    return [Union{String, Float64}[names[i], [cols[c][i] for c in 1:ncols]...] for i in 1:n]
end

# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

"""
    save_realized_market_data(path, <all 30 actual-market variables>)

Write the outputs of one year's actual market clearing to a standalone HDF5
file at `path`.  The file is separate from simulation_data_year_N.h5.
"""
function save_realized_market_data(path::String,
    capacity_price,
    energy_price_ed,
    energy_price_uc,
    energy_price_md,
    reserve_price_ed,
    reserve_price_uc,
    reserve_price_md,
    rec_price,
    inertia_price,
    capacity_factors_md,
    capacity_factors_uc,
    capacity_factors_ed,
    reserve_perc_md,
    reserve_perc_uc,
    reserve_perc_ed,
    capacity_accepted_bids,
    rec_accepted_bids,
    inertia_perc,
    start_up_costs,
    shut_down_costs,
    energy_voll,
    energy_voll_uc,
    energy_voll_md,
    reserve_voll,
    reserve_voll_uc,
    reserve_voll_md,
    inertia_voll,
    rec_supply_curve,
    rec_energy_requirment,
    cet_achieved_ratio,
)
    isdir(dirname(path)) || mkpath(dirname(path))
    h5open(path, "w") do f
        attributes(f)["schema_version"] = SCHEMA_VERSION

        save_axisarray!(create_group(f, "capacity_price"), capacity_price)
        save_axisarray!(create_group(f, "energy_price_ed"), energy_price_ed)
        save_axisarray!(create_group(f, "energy_price_uc"), energy_price_uc)
        save_axisarray!(create_group(f, "energy_price_md"), energy_price_md)

        _save_dict_str_matrix!(create_group(f, "reserve_price_ed"), reserve_price_ed)
        _save_dict_str_matrix!(create_group(f, "reserve_price_uc"), reserve_price_uc)
        _save_dict_str_matrix!(create_group(f, "reserve_price_md"), reserve_price_md)

        save_axisarray!(create_group(f, "rec_price"),      rec_price)
        save_axisarray!(create_group(f, "inertia_price"),  inertia_price)

        _save_dict_str_matrix!(create_group(f, "capacity_factors_md"), capacity_factors_md)
        _save_dict_str_matrix!(create_group(f, "capacity_factors_uc"), capacity_factors_uc)
        _save_dict_str_matrix!(create_group(f, "capacity_factors_ed"), capacity_factors_ed)

        _save_nested_dict_str_matrix!(create_group(f, "reserve_perc_md"), reserve_perc_md)
        _save_nested_dict_str_matrix!(create_group(f, "reserve_perc_uc"), reserve_perc_uc)
        _save_nested_dict_str_matrix!(create_group(f, "reserve_perc_ed"), reserve_perc_ed)

        _save_dict_str_float!(create_group(f, "capacity_accepted_bids"), capacity_accepted_bids)
        _save_dict_str_float!(create_group(f, "rec_accepted_bids"),      rec_accepted_bids)

        _save_dict_str_matrix!(create_group(f, "inertia_perc"),   inertia_perc)
        _save_dict_str_matrix!(create_group(f, "start_up_costs"), start_up_costs)
        _save_dict_str_matrix!(create_group(f, "shut_down_costs"), shut_down_costs)

        save_axisarray!(create_group(f, "energy_voll"),    energy_voll)
        save_axisarray!(create_group(f, "energy_voll_uc"), energy_voll_uc)
        save_axisarray!(create_group(f, "energy_voll_md"), energy_voll_md)

        _save_dict_str_matrix!(create_group(f, "reserve_voll"),    reserve_voll)
        _save_dict_str_matrix!(create_group(f, "reserve_voll_uc"), reserve_voll_uc)
        _save_dict_str_matrix!(create_group(f, "reserve_voll_md"), reserve_voll_md)

        save_axisarray!(create_group(f, "inertia_voll"), inertia_voll)

        _save_supply_curve!(create_group(f, "rec_supply_curve"), rec_supply_curve)
        write(f, "rec_energy_requirment", Float64(rec_energy_requirment))
        write(f, "cet_achieved_ratio",    Float64(cet_achieved_ratio))
    end
end

"""
    load_realized_market_data(path) -> Dict{String, Any}

Read a realized market data file written by `save_realized_market_data`.
Returns a `Dict{String, Any}` with the same keys as the old JLD2 format.
"""
function load_realized_market_data(path::String)
    h5open(path, "r") do f
        return Dict{String, Any}(
            "capacity_price"        => load_axisarray(f["capacity_price"]),
            "energy_price_ed"       => load_axisarray(f["energy_price_ed"]),
            "energy_price_uc"       => load_axisarray(f["energy_price_uc"]),
            "energy_price_md"       => load_axisarray(f["energy_price_md"]),
            "reserve_price_ed"      => _load_dict_str_matrix(f["reserve_price_ed"]),
            "reserve_price_uc"      => _load_dict_str_matrix(f["reserve_price_uc"]),
            "reserve_price_md"      => _load_dict_str_matrix(f["reserve_price_md"]),
            "rec_price"             => load_axisarray(f["rec_price"]),
            "inertia_price"         => load_axisarray(f["inertia_price"]),
            "capacity_factors_md"   => _load_dict_str_matrix(f["capacity_factors_md"]),
            "capacity_factors_uc"   => _load_dict_str_matrix(f["capacity_factors_uc"]),
            "capacity_factors_ed"   => _load_dict_str_matrix(f["capacity_factors_ed"]),
            "reserve_perc_md"       => _load_nested_dict_str_matrix(f["reserve_perc_md"]),
            "reserve_perc_uc"       => _load_nested_dict_str_matrix(f["reserve_perc_uc"]),
            "reserve_perc_ed"       => _load_nested_dict_str_matrix(f["reserve_perc_ed"]),
            "capacity_accepted_bids"=> _load_dict_str_float(f["capacity_accepted_bids"]),
            "rec_accepted_bids"     => _load_dict_str_float(f["rec_accepted_bids"]),
            "inertia_perc"          => _load_dict_str_matrix(f["inertia_perc"]),
            "start_up_costs"        => _load_dict_str_matrix(f["start_up_costs"]),
            "shut_down_costs"       => _load_dict_str_matrix(f["shut_down_costs"]),
            "energy_voll"           => load_axisarray(f["energy_voll"]),
            "energy_voll_uc"        => load_axisarray(f["energy_voll_uc"]),
            "energy_voll_md"        => load_axisarray(f["energy_voll_md"]),
            "reserve_voll"          => _load_dict_str_matrix(f["reserve_voll"]),
            "reserve_voll_uc"       => _load_dict_str_matrix(f["reserve_voll_uc"]),
            "reserve_voll_md"       => _load_dict_str_matrix(f["reserve_voll_md"]),
            "inertia_voll"          => load_axisarray(f["inertia_voll"]),
            "rec_supply_curve"      => _load_supply_curve(f["rec_supply_curve"]),
            "rec_energy_requirment" => read(f, "rec_energy_requirment"),
            "cet_achieved_ratio"    => read(f, "cet_achieved_ratio"),
        )
    end
end