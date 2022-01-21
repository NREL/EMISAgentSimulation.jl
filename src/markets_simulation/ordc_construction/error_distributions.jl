"""
This function constructs the net load forecast error distribution based on Day-Ahead and Real-Time data.
"""
function construct_net_load_forecast_error_distribution(simulation_dir::String,
                                                        renewable_generators::Vector{Project},
                                                        months::Vector{Int64},
                                                        hours::Vector{Int64},
                                                        zonal::Bool)
    load_n_vg_df = read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))
    load_n_vg_df_rt = read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data_rt.csv"))
    zones = chop.(filter(n -> occursin("load", n), names(load_n_vg_df)), head = 5, tail = 0)

    num_rt_intervals = round(Int, DataFrames.nrow(load_n_vg_df_rt)/DataFrames.nrow(load_n_vg_df))

    rt_periods = find_rt_periods(hours, num_rt_intervals)
    filter!(row -> in(row["Month"], months) && in(row["Period"], hours), load_n_vg_df)

    filter!(row -> in(row["Month"], months) && in(row["Period"], rt_periods), load_n_vg_df_rt)

    #half hourly forecasts
    filter!(row -> row["Period"] % (num_rt_intervals/2) == 0, load_n_vg_df_rt)

    load_n_vg_df_extrap = repeat(load_n_vg_df; inner = 2, outer = 1)
    load_n_vg_df_extrap[:, "Period"] = load_n_vg_df_rt[:, "Period"]

    net_load_forecast_error = load_n_vg_df_extrap[:, 1:4]

    if zonal

        zonal_gens = Dict{String, Vector{Project}}()
        error_mean = Dict{String, Float64}()
        error_std = Dict{String, Float64}()

        for zone in zones
            zonal_gens[zone] = filter(p -> get_zone(get_tech(p)) == zone, renewable_generators)

            total_gen_rt = zeros(DataFrames.nrow(load_n_vg_df_rt))
            total_gen_extrap = zeros(DataFrames.nrow(load_n_vg_df_extrap))
            for gen in get_name.(zonal_gens[zone])
                total_gen_rt += load_n_vg_df_rt[:, gen]
                total_gen_extrap += load_n_vg_df_extrap[:, gen]
            end

            load_n_vg_df_rt[:, "Total Generation_$(zone)"] = total_gen_rt
            load_n_vg_df_extrap[:, "Total Generation_$(zone)"] = total_gen_extrap

            load_n_vg_df_rt[:, "Net Load_$(zone)"] = load_n_vg_df_rt[:, "load_$(zone)"] - total_gen_rt
            load_n_vg_df_extrap[:, "Net Load_$(zone)"] = load_n_vg_df_extrap[:, "load_$(zone)"] - total_gen_extrap

            net_load_forecast_error[:, zone] = load_n_vg_df_rt[:, "Net Load_$(zone)"] - load_n_vg_df_extrap[:, "Net Load_$(zone)"]

            error_mean[zone] = Statistics.mean(net_load_forecast_error[:, zone])
            error_std[zone] = Statistics.std(net_load_forecast_error[:, zone])
        end
    else
        total_gen_rt = zeros(DataFrames.nrow(load_n_vg_df_rt))
        total_gen_extrap = zeros(DataFrames.nrow(load_n_vg_df_extrap))
        total_load_rt = zeros(DataFrames.nrow(load_n_vg_df_extrap))
        total_load_extrap = zeros(DataFrames.nrow(load_n_vg_df_extrap))

        for gen in get_name.(renewable_generators)
            total_gen_rt += load_n_vg_df_rt[:, gen]
            total_gen_extrap += load_n_vg_df_extrap[:, gen] * 2/3 + load_n_vg_df_rt[:, gen] * 1/3 # blended forecasts
        end

        for zone in zones
            total_load_rt += load_n_vg_df_rt[:, "load_$(zone)"]
            total_load_extrap += load_n_vg_df_extrap[:, "load_$(zone)"] * 2/3 + load_n_vg_df_rt[:, "load_$(zone)"] * 1/3 # blended forecasts
        end

        load_n_vg_df_rt[:, "Total Generation"] = total_gen_rt
        load_n_vg_df_extrap[:, "Total Generation"] = total_gen_extrap

        load_n_vg_df_rt[:, "Net Load"] = total_load_rt - total_gen_rt
        load_n_vg_df_extrap[:, "Net Load"] = total_load_extrap - total_gen_extrap


        net_load_forecast_error[:, "Total"] = load_n_vg_df_rt[:, "Net Load"] - load_n_vg_df_extrap[:, "Net Load"]

        error_mean = Statistics.mean(net_load_forecast_error[:, "Total"])
        error_std = Statistics.std(net_load_forecast_error[:, "Total"])
    end

    return error_mean, error_std

end


"""
This function constructs the generator unavailability distribution based on Forced Outage Rates using colvolution.
"""
function construct_conv_unavailabilities(simulation_dir::String,
                                            generators::Vector{Project},
                                            zonal::Bool,
                                            ordc_unavailability_method::String)

    if ordc_unavailability_method == "Convolution"
        load_n_vg_df = read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))
        zones = chop.(filter(n -> occursin("load", n), names(load_n_vg_df)), head = 5, tail = 0)

        if zonal
            unavail_mean = Dict{String, Float64}()
            unavail_std = Dict{String, Float64}()
            for zone in zones
                gens = filter(p -> get_zone(get_tech(p)) == zone, generators)
                capacity_vector = round.(Int, get_maxcap.(gens)) #PRAS needs integer capacities
                FOR_vector = get_FOR.(get_tech.(gens))

                distribution = PRAS.ResourceAdequacy.spconv(capacity_vector, FOR_vector)

                unavail_mean[zone] = Distributions.mean(distribution)
                unavail_std[zone] = sqrt(Distributions.var(distribution))
            end
        else
            capacity_vector = round.(Int, get_maxcap.(generators)) #PRAS needs integer capacities
            FOR_vector = get_FOR.(get_tech.(generators))

            distribution = spconv(capacity_vector, FOR_vector)

            unavail_mean = Distributions.mean(distribution)
            unavail_std = sqrt(Distributions.var(distribution))
        end
    else
        unavail_mean = nothing
        unavail_std = nothing
    end

    return unavail_mean, unavail_std
end

"""
This function constructs the generator unavailability distribution using PRAS unavailibility timeseries.
"""
function construct_gen_unavail_distribution(simulation_dir::String,
                                            smc_unavailability_timeseries::Array{Float64, 2},
                                            conv_unavailability_mean::Nothing,
                                            conv_unavailability_std::Nothing,
                                            months::Vector{Int64},
                                            hours::Vector{Int64})

    gen_unavail_df = read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))[:, 1:4]
    @assert DataFrames.nrow(gen_unavail_df) == size(smc_unavailability_timeseries, 2)
    for i in 1:size(smc_unavailability_timeseries, 1)
        gen_unavail_df[!, "$(i)"] = smc_unavailability_timeseries[1, :]
    end
    filter!(row -> in(row["Month"], months) && in(row["Period"], hours), gen_unavail_df)

    unavail_mean = Statistics.mean(Matrix(gen_unavail_df[:, 5:end]))
    unavail_std = Statistics.std(Matrix(gen_unavail_df[:, 5:end]))
    return unavail_mean, unavail_std
end

function construct_gen_unavail_distribution(simulation_dir::String,
                                            smc_unavailability_timeseries::Nothing,
                                            conv_unavailability_mean::Float64,
                                            conv_unavailability_std::Float64,
                                            months::Vector{Int64},
                                            hours::Vector{Int64})

    unavail_mean = conv_unavailability_mean
    unavail_std = conv_unavailability_std
    return unavail_mean, unavail_std
end


function calculate_min_reserve_req(simulation_dir::String,
                                   generators::Vector{Project},
                                   MRR_scale::Real,
                                   zonal::Bool)

    load_n_vg_df = read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))
    zones = chop.(filter(n -> occursin("load", n), names(load_n_vg_df)), head = 5, tail = 0)

    thermal_generators = filter(p -> typeof(p) == ThermalGenEMIS{Existing}, generators)

    if zonal
        req = Dict{String, Float64}()

        for zone in zones
            largest_gen = maximum(get_maxcap.(filter(p -> get_zone(get_tech(p)) == zone, thermal_generators)))
            req[zone] = largest_gen * MRR_scale
        end
    else
        largest_gen = maximum(get_maxcap.(thermal_generators))
        req = largest_gen * MRR_scale
    end

    return req
end
