"""
This function constructs the ORDCs for spinning and primary reserve products.
"""
function construct_ordc(simulation_dir::String,
                        investors::Vector{Investor},
                        iteration_year::Int64,
                        rep_days::Dict{Dates.Date,Int64})

    products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"ordc_products"], "; ")
    println("Creating ORDC for products: $(products)")
    generators = filter(p -> in(typeof(p), leaftypes(GeneratorEMIS{Existing})), vcat(get_existing.(investors)...))
    zonal = false

    load_n_vg_df = read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))
    load_n_vg_df_rt = read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data_rt.csv"))
    zones = chop.(filter(n -> occursin("load", n), names(load_n_vg_df)), head = 5, tail = 0)

    num_rt_intervals = round(Int, DataFrames.nrow(load_n_vg_df_rt)/DataFrames.nrow(load_n_vg_df))

    renewable_generators = filter(g -> typeof(g) == RenewableGenEMIS{Existing}, generators)
    unavail_mean, unavail_std = construct_gen_unavail_distribution(simulation_dir, generators, zonal)

    for product in products

        product_data = read_data(joinpath(simulation_dir, "markets_data", "$(product).csv"))

        MRR_scale = product_data[1,"MRR_scale"]
        penalty = product_data[1,"penalty"]
        step_size = product_data[1, "stepsize (MW)"]
        curve = product_data[1, "curve"]

        MRR = calculate_min_reserve_req(simulation_dir, generators, MRR_scale, zonal)

        println("MRR_$(product): ", MRR)
        println("Penalty_$(product): ", penalty)

        seasons = chop.(split(product_data[1, "seasons"], ";"), head = 1, tail = 1)
        timeblocks = chop.(split(product_data[1, "timeblocks"], ";"), head = 1, tail = 1)

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

        if zonal
            for zone in zones
                ordc_df[:,"$(product)_$(zone)"] = Vector{Vector{Tuple{Float64, Float64}}}(undef, DataFrames.nrow(ordc_df))
                ordc_df_rt[:,"$(product)_$(zone)"] = Vector{Vector{Tuple{Float64, Float64}}}(undef, DataFrames.nrow(ordc_df_rt))
            end
        else
            ordc_df[:,"$(product)"] = Vector{Vector{Tuple{Float64, Float64}}}(undef, DataFrames.nrow(ordc_df))
            ordc_df_rt[:,"$(product)"] = Vector{Vector{Tuple{Float64, Float64}}}(undef, DataFrames.nrow(ordc_df_rt))
        end

        for months_key in keys(season_months)
            months = season_months[months_key]
            for hours_key in keys(timeblock_hours)
                hours = timeblock_hours[hours_key]
                rt_periods = find_rt_periods(hours, num_rt_intervals)

                error_mean, error_var = construct_net_load_forecast_error_distribution(simulation_dir, renewable_generators, months, hours, zonal)
                if zonal
                    aggregate_distribution = Dict{String, Distributions.Normal}()

                    for zone in zones
                        aggregate_distribution_mean = error_mean[zone] + unavail_mean[zone]
                        aggregate_distribution_std = sqrt(error_mean[zone]^2 + unavail_mean[zone]^2)
                        aggregate_distribution[zone] = Distributions.Normal(aggregate_distribution_mean, aggregate_distribution_std)

                        three_std = aggregate_distribution_mean + 3 * aggregate_distribution_std
                        maximum_error = round(ceil(three_std/step_size))*step_size ## Round to the next step_size

                        initial_points = [(0.0, Float64(penalty)), (MRR[zone], Float64(penalty))]

                        ordc_points = [((step * step_size) + MRR[zone],
                                    (1 - Distributions.cdf(aggregate_distribution[zone], step * step_size)) * penalty)
                                    for step in 1:Int64(maximum_error/step_size)]

                        ordc_points = append!(initial_points, ordc_points)

                    end

                    for month in months
                        for hour in hours
                            rows = findall((ordc_df.Month .== month) .& (ordc_df.Period .== hour))
                            for row in rows
                                ordc_df[row, "$(product)_$(zone)"] = ordc_points
                            end
                        end

                        for rt_period in rt_periods
                            rows = findall((ordc_df_rt.Month .== month) .& (ordc_df_rt.Period .== rt_period))
                            for row in rows
                                ordc_df_rt[row, "$(product)_$(zone)"] = ordc_points
                            end
                        end
                    end

                else
                    aggregate_distribution_mean = error_mean + unavail_mean
                    aggregate_distribution_std = sqrt(error_mean^2 + unavail_mean^2)

                    aggregate_distribution = Distributions.Normal(aggregate_distribution_mean, aggregate_distribution_std)

                    three_std = aggregate_distribution_mean + 3*aggregate_distribution_std
                    maximum_error = round(ceil(three_std/step_size))*step_size ## Round to the next step_size MW

                    initial_points = [(0.0, Float64(penalty)), (MRR, Float64(penalty))]

                    n_steps = 10

                    step_size = maximum_error/n_steps

                    ordc_points = [((step * step_size) + MRR,
                                (1 - Distributions.cdf(aggregate_distribution, step * step_size)) * penalty)
                                for step in 1:n_steps]
                    if curve
                        ordc_points = append!(initial_points, ordc_points)
                    end

                    for month in months
                        for hour in hours
                            rows = findall((ordc_df.Month .== month) .& (ordc_df.Period .== hour))
                            for row in rows
                                ordc_df[row, "$(product)"] = ordc_points
                            end
                        end

                        for rt_period in rt_periods
                            rows = findall((ordc_df_rt.Month .== month) .& (ordc_df_rt.Period .== rt_period))
                            for row in rows
                                ordc_df_rt[row, "$(product)"] = ordc_points
                            end
                        end
                    end

                end
            end
        end

        write_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves"), "$(product)_$(iteration_year).csv", ordc_df)
        write_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves"), "$(product)_REAL_TIME_$(iteration_year).csv", ordc_df_rt)

        rep_ordc_df = filter(row -> in(Dates.Date(row[:Year], row[:Month], row[:Day]), keys(rep_days)), ordc_df)

        write_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves"), "rep_$(product)_$(iteration_year).csv", rep_ordc_df)
    end
    return
end

function add_psy_ordc!(simulation_dir::String,
             markets_dict::Dict{Symbol, Bool},
             sys::Nothing,
             type::String)
    return
end

function add_psy_ordc!(simulation_dir::String,
             markets_dict::Dict{Symbol, Bool},
             sys::PSY.System,
             type::String
            )

    products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"ordc_products"], "; ")
    da_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"da_products"], "; ")
    rt_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"rt_products"], "; ")

    for product in products
        if markets_dict[Symbol(product)]
            product_data = read_data(joinpath(simulation_dir, "markets_data", "$(product).csv"))
            eligible_categories = product_data[1, "eligible categories"]

            ####### Adding ORDC reserve
            reserve = PSY.ReserveDemandCurve{PSY.ReserveUp}(
                nothing,
                product,
                true,
                product_data[1, "timescale (min)"],
            )

                PSY.add_service!(sys, reserve, PSY.get_components(PSY.ThermalStandard, sys))

                if occursin("Hydro", eligible_categories)
                    for component in PSY.get_components(PSY.HydroDispatch, sys)
                        PSY.add_service!(component, reserve, sys)
                    end
                    for component in PSY.get_components(PSY.HydroEnergyReservoir, sys)
                        PSY.add_service!(component, reserve, sys)
                    end
                end
                if occursin("Wind", eligible_categories)
                    for component in PSY.get_components(PSY.RenewableDispatch, sys)
                        PSY.add_service!(component, reserve, sys)
                    end
                end
                if occursin("Battery", eligible_categories)
                    for component in PSY.get_components(PSY.GenericBattery, sys)
                        PSY.add_service!(component, reserve_UC, sys)
                    end
                end


            time_stamps = TS.timestamp(PSY.get_data(PSY.get_time_series(
                                                    PSY.SingleTimeSeries,
                                                    first(PSY.get_components(PSY.ElectricLoad, sys)),
                                                    "max_active_power"
                                                    )))

            if type == "UC"
                product_da_ts_raw = read_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "$(product)_0.csv"))[:, product]
            elseif type == "ED"
                product_da_ts_raw = read_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "$(product)_REAL_TIME_0.csv"))[:, product]
            else
                error("Type should be either DA or ED")
            end

            T = length(product_da_ts_raw)

            product_da_ts = Vector{Vector{Tuple{Float64, Float64}}}(undef, T)

            for t = 1:T
                tuples = split.(chop.(split(chop(product_da_ts_raw[t], head = 1, tail = 2), "), "), head = 1, tail = 0), ", ")
                product_da_ts[t] = [(parse.(Float64, tuple)[2], parse.(Float64, tuple)[1]) for tuple in tuples]
            end

            ts = TS.TimeArray(time_stamps, product_da_ts);
            forecast = PSY.SingleTimeSeries("variable_cost", ts)
            PSY.set_variable_cost!(sys, reserve, forecast)
        end
    end

    return
end
