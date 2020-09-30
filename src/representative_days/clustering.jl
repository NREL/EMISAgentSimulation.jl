
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
    cluster_nperiods = Dict(zip(medoid_periods, result.counts))
    period_clusters = Dict(zip(periods, result.assignments))

    return cluster_nperiods, period_clusters

end

function find_representative_days(simulation_dir::String,
                                  base_dir::String,
                                  n_clusters::Int64)

    representative_days_file = joinpath(base_dir, "representative_days_data_$(n_clusters).csv")

    # Use the same representative days if file already exists, otherwise create new clusters
    if isfile(representative_days_file)
        representative_days_data = read_data(representative_days_file)
        clusters = Dict(representative_days_data[i, "Date"] => representative_days_data[i, "Weight"] for i in 1:DataFrames.nrow(representative_days_data))
    else
        net_load_data = read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))
        labels = [Dates.Date(net_load_data[i, :Year], net_load_data[i, :Month], net_load_data[i, :Day]) for i in 1:DataFrames.nrow(net_load_data)]
        load = vec(sum(Matrix(net_load_data[:, r"load"]), dims=2))
        wind = zeros(length(load))
        pv = zeros(length(load))

        for i in names(net_load_data)
            if occursin("WIND", i) || occursin("WT", i)
                wind += net_load_data[:, i]
            elseif occursin("PV", i) || occursin("PVe", i)
                pv += net_load_data[:, i]
            end
        end
        clusters, allperiods = cluster(labels, load, wind, pv, n_clusters)
        representative_days_data = DataFrames.DataFrame(Date = Dates.Date[], Weight = Int64[])
        for i in keys(clusters)
            push!(representative_days_data, (i, clusters[i]))
        end
        CSV.write(representative_days_file, representative_days_data)
    end
    return clusters
end
