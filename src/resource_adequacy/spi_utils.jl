using DataFrames
using CSV
using Dates
using TimeSeries

"""
    attach_outage_data_from_csv!(sys, outage_csv_file; mttr_hours=24)

Read a forced-outage-rate (FOR) CSV with one column per generator name and one row per
timestep - the same format consumed by `PSY2PRAS.make_pras_system`'s
`outage_csv_location`/`outage_ts_flag` kwargs - and attach a
`PSY.GeometricDistributionForcedOutage` supplemental attribute to each matching generator
in `sys`, with time-varying `"outage_probability"`/`"recovery_probability"` time series
attached to it in the format `SiennaPRASInterface.generate_pras_system` expects (see
`SiennaPRASInterface`'s "How do I add outage data to Sienna?" how-to guide).

`mttr_hours` is the fixed mean time to recovery (in hours) assumed for every generator,
matching the assumption already made in `PSY2PRAS.make_pras_system`.

If `outage_csv_file === nothing`, this is a no-op: `sys` is returned unchanged, so any
`GeometricDistributionForcedOutage` attributes already on `sys` are used as-is, and
`SiennaPRASInterface.generate_pras_system` falls back to its own nominal outage data for
any generator that still has none.

Any pre-existing `GeometricDistributionForcedOutage` attribute on a matched generator is
removed before the new one is attached, so this function can be called more than once on
the same system without duplicating attributes.
"""
function attach_outage_data_from_csv!(
    sys::PSY.System,
    outage_csv_file::Union{Nothing, String};
    mttr_hours::Real = 24,
)
    if outage_csv_file === nothing
        return sys
    end

    outage_df = read_data(outage_csv_file)
    gen_names = DataFrames.names(outage_df)
    n_rows = DataFrames.nrow(outage_df)

    all_ts = PSY.get_time_series_multiple(sys)
    isempty(all_ts) &&
        error(
            "Cannot attach outage time series: system has no existing time series to infer a start time and resolution from.",
        )
    start_datetime = PSY.IS.get_initial_timestamp(first(all_ts))
    resolution = PSY.get_time_series_resolutions(sys)[1]
    timestamps = range(start_datetime; step = resolution, length = n_rows)

    for gen in PSY.get_components(PSY.Generator, sys)
        gname = PSY.get_name(gen)
        gname in gen_names || continue

        for_values = Float64.(outage_df[1:n_rows, gname])
        rates = SPI.rate_to_probability.(for_values, Int(mttr_hours))
        λ_values = getfield.(rates, :λ)
        μ_values = getfield.(rates, :μ)

        for existing_attr in
            PSY.get_supplemental_attributes(PSY.GeometricDistributionForcedOutage, gen)
            PSY.remove_supplemental_attribute!(sys, gen, existing_attr)
        end

        transition_data = PSY.GeometricDistributionForcedOutage(;
            mean_time_to_recovery = mttr_hours,
            outage_transition_probability = λ_values[1],
        )
        PSY.add_supplemental_attribute!(sys, gen, transition_data)

        outage_ts = PSY.SingleTimeSeries(;
            name = "outage_probability",
            data = TimeArray(timestamps, λ_values),
        )
        recovery_ts = PSY.SingleTimeSeries(;
            name = "recovery_probability",
            data = TimeArray(timestamps, μ_values),
        )
        PSY.add_time_series!(sys, transition_data, outage_ts)
        PSY.add_time_series!(sys, transition_data, recovery_ts)
    end

    return sys
end

"""
    attach_outage_data_from_ext!(sys)

Convert each `PSY.Generator`/`PSY.Storage` component's existing `ext["outage_probability"]`
/`ext["recovery_probability"]` (set by `add_outage_info!` from the component's own Tech
FOR/MTTR data at creation time - see `ra_utils.jl`) into a `PSY.GeometricDistributionForcedOutage`
supplemental attribute, so `SiennaPRASInterface.generate_pras_system` uses the same outage
assumptions `PSY2PRAS.make_pras_system`'s non-CSV fallback path already does.

Components without both ext keys are left untouched (no attribute attached), so
`SiennaPRASInterface`'s own nominal/default outage data applies to them instead.
"""
function attach_outage_data_from_ext!(sys::PSY.System)
    for comp in Iterators.flatten((
        PSY.get_components(PSY.Generator, sys),
        PSY.get_components(PSY.Storage, sys),
    ))
        ext = PSY.get_ext(comp)
        (haskey(ext, "outage_probability") && haskey(ext, "recovery_probability")) ||
            continue

        λ = ext["outage_probability"]
        μ = ext["recovery_probability"]
        mttr_hours = iszero(μ) ? 0.0 : 1 / μ

        for existing_attr in
            PSY.get_supplemental_attributes(PSY.GeometricDistributionForcedOutage, comp)
            PSY.remove_supplemental_attribute!(sys, comp, existing_attr)
        end

        transition_data = PSY.GeometricDistributionForcedOutage(;
            mean_time_to_recovery = mttr_hours,
            outage_transition_probability = λ,
        )
        PSY.add_supplemental_attribute!(sys, comp, transition_data)
    end

    return sys
end

"""
    set_line_capacities!(sys, copper_plate; copper_plate_capacity_mw=99999.0)

Set the forward/backward transfer capacity of every `PSY.Line`, `PSY.MonitoredLine`, and
`PSY.TwoTerminalGenericHVDCLine` in `sys`. If `copper_plate` is `true`, every line is set
to `copper_plate_capacity_mw` so transmission never binds in `SiennaPRASInterface`'s
`LinePRAS` formulation - replicating `PSY2PRAS.make_pras_system`'s `"Single-Node"`
behavior while still going through SPI's normal zonal/regions+lines path. If `false`,
each line is (re-)assigned its own current rating, i.e. real transmission limits are used.
"""
function set_line_capacities!(
    sys::PSY.System,
    copper_plate::Bool;
    copper_plate_capacity_mw::Real = 99999.0,
)
    for line in PSY.get_components(PSY.Branch, sys)
        if line isa PSY.TwoTerminalGenericHVDCLine
            cap_from = if copper_plate
                copper_plate_capacity_mw
            else
                PSY.get_active_power_limits_from(line).max
            end
            cap_to = if copper_plate
                copper_plate_capacity_mw
            else
                PSY.get_active_power_limits_to(line).max
            end
            PSY.set_active_power_limits_from!(line, (min = -cap_from, max = cap_from))
            PSY.set_active_power_limits_to!(line, (min = -cap_to, max = cap_to))
        elseif line isa Union{PSY.Line, PSY.MonitoredLine}
            cap = copper_plate ? copper_plate_capacity_mw : PSY.get_rating(line)
            PSY.set_rating!(line, cap)
        end
    end

    return sys
end

"""
    make_pras_system_spi(sys, aggregation=PSY.Area; outage_csv_file=nothing, mttr_hours=24,
                          copper_plate=true, copper_plate_capacity_mw=99999.0,
                          lump_region_renewable_gens=false, export_location=nothing)

`SiennaPRASInterface`-based equivalent of `PSY2PRAS.make_pras_system`, for testing SPI as
a drop-in replacement before rewiring the live call sites in `derating_factor_calculator.jl`.
Operates on a deep copy of `sys`, so the caller's system is left untouched.

- `outage_csv_file`: forwarded to `attach_outage_data_from_csv!`, run after
  `attach_outage_data_from_ext!` populates static outage data from each component's
  existing `ext["outage_probability"]`/`ext["recovery_probability"]`. Pass `nothing` to
  skip CSV-based (time-varying) outage attachment and rely solely on the ext-based
  static values (or SPI's own nominal defaults for anything without either).
- `copper_plate`: when `true` (default), line capacities are relaxed to
  `copper_plate_capacity_mw` so the assessment is effectively single-node/copper-plate;
  when `false`, each line's real transmission rating is used.
"""
function make_pras_system_spi(
    sys::PSY.System,
    aggregation::Type{<:PSY.AggregationTopology} = PSY.Area;
    outage_csv_file::Union{Nothing, String} = nothing,
    mttr_hours::Real = 24,
    copper_plate::Bool = true,
    copper_plate_capacity_mw::Real = 99999.0,
    lump_region_renewable_gens::Bool = false,
    export_location::Union{Nothing, String} = nothing,
)
    sys = deepcopy(sys)
    attach_outage_data_from_ext!(sys)
    attach_outage_data_from_csv!(sys, outage_csv_file; mttr_hours = mttr_hours)
    set_line_capacities!(
        sys,
        copper_plate;
        copper_plate_capacity_mw = copper_plate_capacity_mw,
    )

    return SPI.generate_pras_system(
        sys,
        aggregation,
        lump_region_renewable_gens,
        export_location,
    )
end
