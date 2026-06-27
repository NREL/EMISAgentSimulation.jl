"""
This function constructs the ORDCs for spinning and primary reserve products.
"""
function construct_ordc(sys::PSY.System,
    simulation_dir::String,
    scenario::String,
    sim_year::Int64,
    investors::Vector{Investor},
    iteration_year::Int64,
    representative_periods::Union{
        Dict{Int64, Int64},
        OrderedCollections.OrderedDict{Int64, Int64},
    },
    rep_period_interval::Int64,
    ordc_curved::Bool,
    ordc_unavailability_method::String,
    reserve_penalty::String,
    timeseries_data_dir::String)

    products = split(
        read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[
            1,
            "ordc_products",
        ],
        "; ",
    )
    @info "Creating ORDC for products: $(products)"
    existing_gens = leaftypes(GeneratorEMIS{Existing})
    generators =
        filter(p -> in(typeof(p), existing_gens), vcat(get_existing.(investors)...))

    ##TODO AA: This should be moved from here
    zonal = false

    smc_unavailability_timeseries =
        construct_smc_unavailabilities(sys, ordc_unavailability_method)
    conv_unavail_mean, conv_unavail_std = construct_conv_unavailabilities(
        simulation_dir,
        scenario,
        sim_year,
        generators,
        zonal,
        ordc_unavailability_method,
        timeseries_data_dir
    )

    load_n_vg_df = read_data(
        joinpath(
            timeseries_data_dir,
            scenario,
            "sim_year_$(sim_year)",
            "Net Load Data",
            "load_n_vg_data.csv",
        ),
    )
    load_n_vg_df_rt = read_data(
        joinpath(
            timeseries_data_dir,
            scenario,
            "sim_year_$(sim_year)",
            "Net Load Data",
            "load_n_vg_data_rt.csv",
        ),
    )
    zones = chop.(filter(n -> occursin("load", n), names(load_n_vg_df)); head = 5, tail = 0)

    num_rt_intervals =
        round(Int, DataFrames.nrow(load_n_vg_df_rt)/DataFrames.nrow(load_n_vg_df))

    ## TODO AA: Check this fix. 
    # Previous error: ArgumentError: column name "investor1_PVe_zone_5_year_1_1" not found in the data frame; existing most similar names are: "investor1_PVe_zone_5_year_5_1", "investor1_PVe_zone_8_year_5_1", "investor3_PVe_zone_5_year_5_1" and "investor4_PVe_zone_5_year_5_1"
    # The load_n_vg_df_rt file contains sim_year data, but generators 
    renewable_generators = filter(
        g -> typeof(g) == RenewableGenEMIS{Existing} && get_name(g) in names(load_n_vg_df_rt),
        generators)

    for product in products
        product_data = read_data(
            joinpath(
                simulation_dir,
                "markets_data",
                "$(reserve_penalty)_reserve_penalty",
                "$(product).csv",
            ),
        )

        MRR_scale = product_data[1, "MRR_scale"]
        penalty = product_data[1, "penalty"]
        step_size = product_data[1, "stepsize (MW)"]

        MRR = calculate_min_reserve_req(
            simulation_dir,
            scenario,
            sim_year,
            generators,
            MRR_scale,
            zonal,
            timeseries_data_dir
        )

        #println("MRR_$(product): ", MRR)
        #println("Penalty_$(product): ", penalty)

        seasons = chop.(split(product_data[1, "seasons"], ";"); head = 1, tail = 1)
        timeblocks = chop.(split(product_data[1, "timeblocks"], ";"); head = 1, tail = 1)

        season_months = Dict{String, Vector{Int64}}()
        timeblock_hours = Dict{String, Vector{Int64}}()

        for season in seasons
            months = split(season, "-")
            start_month = month_lookup(months[1])
            end_month = month_lookup(months[2])

            if start_month <= end_month
                season_months[season] = collect(start_month:end_month)
            else
                season_months[season] = collect(1:end_month)
                append!(season_months[season], collect(start_month:12))
            end
        end

        for timeblock in timeblocks
            hours = parse.(Int64, split(timeblock, "-"))/100
            start_hour = hours[1]
            end_hour = hours[2]

            if start_hour <= end_hour
                timeblock_hours[timeblock] = collect(start_hour:end_hour)
            else
                timeblock_hours[timeblock] = collect(1:end_hour)
                append!(timeblock_hours[timeblock], collect(start_hour:24))
            end
        end

        ordc_df = load_n_vg_df[:, 1:4]
        ordc_df_rt = load_n_vg_df_rt[:, 1:4]
        #Need code fix to work with zonal ORDC markets
        if zonal
            for zone in zones
                ordc_df[:, "$(product)_$(zone)"] =
                    Vector{Vector{Tuple{Float64, Float64}}}(undef, DataFrames.nrow(ordc_df))
                ordc_df_rt[:, "$(product)_$(zone)"] =
                    Vector{Vector{Tuple{Float64, Float64}}}(
                        undef,
                        DataFrames.nrow(ordc_df_rt),
                    )
            end
        else
            ordc_df[:, "$(product)"] =
                Vector{Vector{Tuple{Float64, Float64}}}(undef, DataFrames.nrow(ordc_df))
            ordc_df_rt[:, "$(product)"] =
                Vector{Vector{Tuple{Float64, Float64}}}(undef, DataFrames.nrow(ordc_df_rt))
        end

        for months_key in keys(season_months)
            months = season_months[months_key]
            for hours_key in keys(timeblock_hours)
                hours = timeblock_hours[hours_key]
                rt_periods = find_rt_periods(hours, num_rt_intervals)

                error_mean, error_var, meanload =
                    construct_net_load_forecast_error_distribution(
                        simulation_dir,
                        scenario,
                        sim_year,
                        renewable_generators,
                        months,
                        hours,
                        zonal,
                        timeseries_data_dir,
                    )
                #println("$(product), $(months_key): error_mean is $(error_mean) & error_var is $(error_var) & meanload is $(meanload)")
                unavail_mean, unavail_std = construct_gen_unavail_distribution(
                    simulation_dir,
                    scenario,
                    sim_year,
                    smc_unavailability_timeseries,
                    conv_unavail_mean,
                    conv_unavail_std,
                    months,
                    hours,
                    timeseries_data_dir,
                )
                #println("$(product), $(months_key): unavail_mean is $(unavail_mean) & unavail_std is $(unavail_std)")

                if zonal
                    aggregate_distribution = Dict{String, Distributions.Normal}()

                    for zone in zones
                        aggregate_distribution_mean = error_mean[zone] + unavail_mean[zone]
                        aggregate_distribution_std =
                            sqrt(error_var[zone]^2 + unavail_std[zone]^2)
                        aggregate_distribution[zone] = Distributions.Normal(
                            aggregate_distribution_mean,
                            aggregate_distribution_std,
                        )

                        three_std =
                            aggregate_distribution_mean + 3 * aggregate_distribution_std
                        maximum_error = round(ceil(three_std/step_size))*step_size ## Round to the next step_size
                        if maximum_error < min(aggregate_distribution_std, 0.1*meanload)
                            maximum_error = min(aggregate_distribution_std, 0.1*meanload) #this step is used when the mean net load error+3std is negative; then we use a small positive number to represent the x-interval (above MRR)
                        end

                        initial_points =
                            [(0.0, Float64(penalty)), (MRR[zone], Float64(penalty))]

                        ordc_points = [
                            ((step * step_size) + MRR[zone],
                                (
                                    1 - Distributions.cdf(
                                        aggregate_distribution[zone],
                                        step * step_size,
                                    )
                                ) * penalty)
                            for step in 1:Int64(maximum_error / step_size)
                        ]

                        if ordc_curved
                            ordc_points = append!(initial_points, ordc_points)
                        else
                            ordc_points = initial_points
                        end
                    end

                    for month in months
                        for hour in hours
                            rows = findall(
                                (ordc_df.Month .== month) .& (ordc_df.Period .== hour),
                            )
                            for row in rows
                                ordc_df[row, "$(product)_$(zone)"] = ordc_points
                            end
                        end

                        for rt_period in rt_periods
                            rows = findall(
                                (ordc_df_rt.Month .== month) .&
                                (ordc_df_rt.Period .== rt_period),
                            )
                            for row in rows
                                ordc_df_rt[row, "$(product)_$(zone)"] = ordc_points
                            end
                        end
                    end

                else
                    aggregate_distribution_mean = error_mean + unavail_mean
                    aggregate_distribution_std = sqrt(error_var^2 + unavail_std^2)

                    aggregate_distribution = Distributions.Normal(
                        aggregate_distribution_mean,
                        aggregate_distribution_std,
                    )

                    three_std = aggregate_distribution_mean + 3*aggregate_distribution_std
                    maximum_error = round(ceil(three_std/step_size))*step_size ## Round to the next step_size MW
                    if maximum_error < min(aggregate_distribution_std, 0.1*meanload)
                        maximum_error = min(aggregate_distribution_std, 0.1*meanload) #this step is used when the mean net load error+3std is negative; then we use a small positive number to represent the x-interval (above MRR)
                    end
                    #println("$(product), $(months_key): maximum_error is $(maximum_error)")

                    initial_points = [(0.0, Float64(penalty)), (MRR, Float64(penalty))]

                    n_steps = 10

                    step_size = maximum_error/n_steps

                    ordc_points = [
                        ((step * step_size) + MRR,
                            (
                                1 - Distributions.cdf(
                                    aggregate_distribution,
                                    step * step_size,
                                )
                            ) * penalty)
                        for step in 1:n_steps
                    ]
                    if ordc_curved
                        ordc_points = append!(initial_points, ordc_points)
                    else
                        ordc_points = initial_points
                    end

                    for month in months
                        for hour in hours
                            rows = findall(
                                (ordc_df.Month .== month) .& (ordc_df.Period .== hour),
                            )
                            for row in rows
                                ordc_df[row, "$(product)"] = ordc_points
                            end
                        end

                        for rt_period in rt_periods
                            rows = findall(
                                (ordc_df_rt.Month .== month) .&
                                (ordc_df_rt.Period .== rt_period),
                            )
                            for row in rows
                                ordc_df_rt[row, "$(product)"] = ordc_points
                            end
                        end
                    end
                end
            end
        end

        write_data(
            joinpath(
                timeseries_data_dir,
                scenario,
                "sim_year_$(sim_year)",
                "Reserves",
            ),
            "$(product).csv",
            ordc_df,
        )
        write_data(
            joinpath(
                timeseries_data_dir,
                scenario,
                "sim_year_$(sim_year)",
                "Reserves",
            ),
            "$(product)_REAL_TIME.csv",
            ordc_df_rt,
        )

        ordc_df[!, "Period_Number"] = 1:size(ordc_df, 1)
        ordc_df[!, "Representative_Period"] =
            add_representative_period.(ordc_df[:, "Period_Number"], rep_period_interval)

        rep_ordc_df = filter(
            row -> in(row["Representative_Period"], keys(representative_periods)),
            ordc_df,
        )

        write_data(
            joinpath(
                timeseries_data_dir,
                scenario,
                "sim_year_$(sim_year)",
                "Reserves",
            ),
            "rep_$(product).csv",
            rep_ordc_df,
        )
    end
    return
end

function process_ordc_data_for_siip(
    raw_data::Union{
        Vector{String},
        PooledArrays.PooledVector{String, UInt32, Vector{UInt32}},
    },
)
    T = length(raw_data)

    product_da_ts = Vector{PSY.PiecewiseStepData}(undef, T)

    # for t = 1:T
    #     tuples = split.(chop.(split(chop(raw_data[t], head = 1, tail = 2), "), "), head = 1, tail = 0), ", ")
    #     product_da_ts[t] = [(parse.(Float64, tuple)[2], parse.(Float64, tuple)[1]) for tuple in tuples]
    #     l = length(product_da_ts[t])
    #     if l > 2
    #         for s in 3:2:(l * 2) - 2
    #             temp_1 = copy(product_da_ts[t][1:(s-1)])
    #             temp_2 = copy(product_da_ts[t][s:end])
    #             push!(temp_1, (product_da_ts[t][s][1], product_da_ts[t][s-1][2]))
    #             new = vcat(temp_1, temp_2)
    #             product_da_ts[t] = copy(new)
    #         end
    #      end
    # end

    # reconstruct the ORDC curve to be a PiecewiseStepCurve of n power points and n -1 price points
    for t in 1:T
        tuples = split.(
            chop.(split(chop(raw_data[t]; head = 1, tail = 2), "), "); head = 1, tail = 0),
            ", ",
        )
        # (price, quantity) pairs
        product_da =
            [(parse.(Float64, tuple)[2], parse.(Float64, tuple)[1]) for tuple in tuples]
        l = length(product_da)
        prices = zeros(Float64, l - 1)
        quantities = zeros(Float64, l)
        prices = [point[1] for point in product_da][2:end]
        quantities = [point[2] for point in product_da]
        # reconstruct as (quantity, price) pairs for PiecewiseStepCurve
        product_da_ts[t] = PSY.PiecewiseStepData(quantities, prices)
        # l = length(product_da_ts[t])
        # for i in 1:l
        #     if i == 1
        #         product_da_ts[t][i] = (product_da_ts[t][i][1] * product_da_ts[t][i][2], product_da_ts[t][i][2])
        #     else
        #         product_da_ts[t][i] = (product_da_ts[t][i-1][1] + product_da_ts[t][i][1] * (product_da_ts[t][i][2] - product_da_ts[t][i-1][2]), product_da_ts[t][i][2])
        #     end
        # end
    end

    return product_da_ts
end

function add_psy_ordc!(simulation_dir::String,
    markets_dict::Dict{Symbol, Bool},
    sys::Nothing,
    type::String,
    scenario::String,
    iteration_year::Int64,
    da_resolution::Int64,
    rt_resolution::Int64,
    reserve_penalty::String,
    timeseries_data_dir::String)
    return
end

function add_psy_ordc!(simulation_dir::String,
    markets_dict::Dict{Symbol, Bool},
    sys::PSY.System,
    type::String,
    scenario::String,
    iteration_year::Int64,
    da_resolution::Int64,
    rt_resolution::Int64,
    reserve_penalty::String,
    timeseries_data_dir::String
)
    products = split(
        read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[
            1,
            "ordc_products",
        ],
        "; ",
    )
    da_products = split(
        read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[
            1,
            "da_products",
        ],
        "; ",
    )
    rt_products = split(
        read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[
            1,
            "rt_products",
        ],
        "; ",
    )

    sys_interval = PSY.get_forecast_interval(sys)
    sys_horizon = PSY.get_forecast_horizon(sys)
    sys_horizon = Dates.Hour(sys_horizon).value
    forecast_count = PSY.get_forecast_window_count(sys)
    sys_resolution = PSY.get_time_series_resolutions(sys)[1]
    start_datetime = PSY.get_forecast_initial_timestamp(sys)
    finish_datetime =
        start_datetime + Dates.Hour((
            forecast_count * sys_interval/sys_resolution +
            (sys_horizon - sys_interval/sys_resolution) - 1
        ))
    # start_datetime = Dates.DateTime("2019-01-01T00:00:00")
    # finish_datetime = start_datetime + Dates.Hour(8759)
    time_stamps = StepRange(start_datetime, Dates.Hour(1), finish_datetime);

    additional_timestep = length(time_stamps) - DEFAULT_HOURS_PER_YEAR

    for product in products
        if markets_dict[Symbol(product)]
            product_data = read_data(
                joinpath(
                    simulation_dir,
                    "markets_data",
                    "$(reserve_penalty)_reserve_penalty",
                    "$(product).csv",
                ),
            )
            eligible_categories = product_data[1, "eligible categories"]

            ####### Adding ORDC reserve
            reserve = PSY.ReserveDemandCurve{PSY.ReserveUp}(
                nothing,    #InfrastructureSystems.TimeSeriesKey
                product,
                true,
                product_data[1, "timescale (min)"] * DEFAULT_TIME_RESOLUTION,
            )

            PSY.add_service!(sys, reserve, PSY.get_components(PSY.ThermalStandard, sys))

            # for component in PSY.get_components(PSYE.ThermalCleanEnergy, sys)
            #     PSY.add_service!(component, reserve, sys)
            # end

            for component in PSY.get_components(ThermalFastStartSIIP, sys)
                PSY.add_service!(component, reserve, sys)
            end

            for component in PSY.get_components(PSY.ThermalMultiStart, sys)
                PSY.add_service!(component, reserve, sys)
            end

            if occursin("Hydro", eligible_categories)
                for component in PSY.get_components(PSY.HydroDispatch, sys)
                    PSY.add_service!(component, reserve, sys)
                end
                for component in PSY.get_components(PSY.HydroTurbine, sys)
                    PSY.add_service!(component, reserve, sys)
                end
            end
            if occursin("Wind", eligible_categories)
                for component in PSY.get_components(PSY.RenewableDispatch, sys)
                    PSY.add_service!(component, reserve, sys)
                end
            end
            if occursin("Battery", eligible_categories)
                for component in PSY.get_components(PSY.EnergyReservoirStorage, sys)
                    PSY.add_service!(component, reserve, sys)
                end
            end

            # time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
            #     PSY.SingleTimeSeries,
            #     first(PSY.get_components(PSY.ElectricLoad, sys)),
            #     "max_active_power"
            #     )))

            # time_stamps = StepRange(start_datetime, Dates.Hour(1), finish_datetime);
            if type == "ED"
                product_ts_raw = read_data(
                    joinpath(
                        timeseries_data_dir,
                        scenario,
                        "sim_year_$(iteration_year)",
                        "Reserves",
                        "$(product)_REAL_TIME.csv",
                    ),
                )[
                    :,
                    product,
                ]
            elseif type in ["UC", "MD"]
                product_ts_raw = read_data(
                    joinpath(
                        timeseries_data_dir,
                        scenario,
                        "sim_year_$(iteration_year)",
                        "Reserves",
                        "$(product).csv",
                    ),
                )[
                    :,
                    product,
                ]
            end

            if type in ["MD", "UC", "ED"]
                @info "Processing ORDC timeseries for $(product) product - $(type) model"
                product_data_ts = process_ordc_data_for_siip(product_ts_raw)
                product_data_ts = [product_data_ts; product_data_ts[1:additional_timestep]]

                # Use Deterministic (ForecastKey) instead of SingleTimeSeries (StaticTimeSeriesKey)
                # PSI's _get_reserve_pwl_data only handles ForecastKey, not StaticTimeSeriesKey
                interval_steps = round(Int, sys_interval / sys_resolution)
                data_dict = Dict{Dates.DateTime, Vector{PSY.PiecewiseStepData}}()
                for w in 1:forecast_count
                    t_start = start_datetime + (w - 1) * sys_interval
                    idx_start = (w - 1) * interval_steps + 1
                    data_dict[t_start] =
                        product_data_ts[idx_start:(idx_start + sys_horizon - 1)]
                end
                forecast = PSY.Deterministic("variable_cost", data_dict, sys_resolution)
                PSY.set_variable_cost!(sys, reserve, forecast)
            end
        end
    end

    return
end

#=
system = SystemModel("../PLEXOS2PRAS/test/rts/rts_interfaces.pras")
nsamples = 100

unavailable = unavailabilities(system, nsamples) # samples x timesteps
mus = vec(mean(unavailable, dims=1))
sigmas = vec(std(unavailable, dims=1))
=#
