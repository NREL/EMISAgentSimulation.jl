
function cluster(period_labels::AbstractVector,
                 load::AbstractVector{Float64},
                 wind::AbstractVector{Float64},
                 pv::AbstractVector{Float64},
                 n_clusters::Int)

    n_labels = length(period_labels)
    n_load   = length(load)
    n_wind   = length(wind)
    n_pv     = length(pv)

    (n_labels == n_load == n_wind == n_pv) ||
        error("Time series inputs must have matching length")

    freqs = StatsBase.countmap(period_labels)
    periods = collect(keys(freqs))
    counts = values(freqs)
    period_ntimesteps = first(counts)
    all(x -> x == period_ntimesteps, counts) ||
        error("Each period must contain the same number of timesteps")

    n_periods = length(periods)
    period_key = Dict(zip(keys(freqs), 1:n_periods))

    # Create feature vectors for each period
    features = Matrix{Float64}(undef, 3 * period_ntimesteps, n_periods)
    idxs = 1:period_ntimesteps
    for p in keys(freqs)
        p_idx = period_key[p]
        include_timesteps = findall(isequal(p), period_labels)
        features[idxs, p_idx] = load[include_timesteps]
        features[period_ntimesteps .+ idxs, p_idx] = wind[include_timesteps]
        features[2*period_ntimesteps .+ idxs, p_idx] = pv[include_timesteps]
    end

    # Normalize feature vectors
    features .-= Statistics.mean(features, dims=2)
    stds = Statistics.std(features, dims=2)
    for i in 1:size(features, 1)
        stds[i] > 0 && (features[i, :] ./= stds[i])
    end

    # Create period distance matrix
    feature_dists = [LinearAlgebra.norm(ci .- cj, 2)
                     for ci in eachcol(features), cj in eachcol(features)]

    result = Clustering.kmedoids(feature_dists, n_clusters)
    medoid_periods = getindex.(Ref(periods), result.medoids)
    cluster_nperiods = sort(Dict(zip(medoid_periods, result.counts)))
    period_clusters = sort(Dict(periods[i] => medoid_periods[result.assignments[i]] for i in 1:length(periods)))

    return cluster_nperiods, period_clusters

end


function add_representative_period(period, interval)
    return Int(floor((period - 1) / interval) + 1)
end

function prune_extra_rows!(df, interval)
    last_row = findlast(x -> x == Int(floor(nrow(df) / interval)), df[:, "Representative_Period"])
    return df = df[1:last_row, :]
end

function find_representative_periods(simulation_dir::String,
                                  test_system_dir::String,
                                  base_dir::String,
                                  interval::Int64,
                                  n_clusters::Int64)

    representative_periods_file = joinpath(base_dir, "representative_periods_$(interval)h_interval_$(n_clusters)_clusters.csv")
    cluster_indices_file = joinpath(base_dir, "cluster_indices_$(interval)h_interval_$(n_clusters)_clusters.csv")

    # Use the same representative periods if file already exists, otherwise create new clusters
    if isfile(representative_periods_file) && isfile(cluster_indices_file)
        representative_periods_data = read_data(representative_periods_file)
        clusters = sort(Dict(representative_periods_data[i, "Period"] => representative_periods_data[i, "Weight"] for i in 1:DataFrames.nrow(representative_periods_data)))
        
        cluster_indices_data = read_data(cluster_indices_file)
        allperiods = sort(Dict(cluster_indices_data[i, "Period"] => cluster_indices_data[i, "Rep_Period"] for i in 1:DataFrames.nrow(cluster_indices_data)))
        
    else
        net_load_data = read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))
        net_load_data[!, "Period_Number"] = 1:size(net_load_data, 1)
        net_load_data[!, "Representative_Period"] = add_representative_period.(net_load_data[:, "Period_Number"], interval)
        net_load_data = prune_extra_rows!(net_load_data, interval)
        
        labels = net_load_data[!, "Representative_Period"]

        load = vec(sum(Matrix(net_load_data[:, r"load"]), dims=2))
        wind = zeros(length(load))
        pv = zeros(length(load))

        existing_generator_data = DataFrames.DataFrame(CSV.File(joinpath(test_system_dir, "RTS_Data", "SourceData", "gen.csv")))

        for i in names(net_load_data)
            if length(filter( row -> row."GEN UID" in [i], existing_generator_data)."Unit Type") > 0
                if occursin("WIND", filter( row -> row."GEN UID" in [i], existing_generator_data)."Unit Type"[1]) #occursin("WIND", i) || occursin("WT", i)
                    wind += net_load_data[:, i]
                elseif occursin("PV", filter( row -> row."GEN UID" in [i], existing_generator_data)."Unit Type"[1])#occursin("PV", i) || occursin("PVe", i)
                    pv += net_load_data[:, i]
                end
            end
        end

        clusters, allperiods = cluster(labels, load, wind, pv, n_clusters)

        representative_periods_data = DataFrames.DataFrame(Period= Int64[], Weight = Int64[])
        cluster_indices_data = DataFrames.DataFrame(Period = Int64[], Rep_Period= Int64[])
        for i in keys(clusters)
            push!(representative_periods_data, (i, clusters[i]))
        end
        for i in keys(allperiods)
            push!(cluster_indices_data, [i, allperiods[i]])
        end

        CSV.write(representative_periods_file, representative_periods_data)
        CSV.write(cluster_indices_file, cluster_indices_data)
    end
    
    return clusters, allperiods
end
