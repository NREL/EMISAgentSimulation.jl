using Clustering
using CSV
using Dates
using DataFrames
using Statistics

# Define a custom function to extract the upper triangular part of a matrix.
function custom_triu(A::AbstractMatrix)
    m, n = size(A)
    upper_triangular_part = similar(A, m, n)
    
    for i in 1:m, j in i:n
        upper_triangular_part[i, j] = A[i, j]
    end
    
    return upper_triangular_part
end

function normalize_vector(vector::Vector{Float64})
    # Calculate mean and standard deviation
    mean_value = mean(vector)
    std_value = std(vector)
    
    # Normalize the vector
    normalized_vector = (vector .- mean_value) / std_value
    
    return normalized_vector
end

function calculate_centroid(data)
    num_points = size(data, 1)  # Get the number of data points.

    # Calculate the centroid by taking the sum of data points in each dimension within the cluster and dividing by the number of data points.
    centroid = sum(data, dims=1) / num_points

    return vec(centroid)  # Convert the centroid to a 1D vector.
end

function calculate_dissimilarity(data, cluster1, cluster2)
    # Extract data points belonging to cluster1 and cluster2
    data_cluster1 = data[cluster1, :]
    data_cluster2 = data[cluster2, :]

    # Calculate the squared Euclidean distance between the centroids of cluster1 and cluster2
    squared_distance = sum((mean(data_cluster1, dims=1) .- mean(data_cluster2, dims=1)).^2)

    # Calculate the dissimilarity based on the provided formula
    dissimilarity = 2 * length(cluster1) * length(cluster2)  * squared_distance / (length(cluster1) + length(cluster2))

    return dissimilarity
end

# Define a function to check if two time periods are adjacent.
function is_adjacent(hours1, hours2)
    for i in hours1
        for j in hours2
            if abs(i - j) == 1
                return true  # At least one pair of consecutive elements found.
            end
        end
    end
    return false
end


function chronological_clustering(xi, K)
    # Define the number of hours in a time period (e.g., 24 for daily periods).
    τt = 1


    # Step 1: Set the initial number of clusters n to the total number of hours N.
    N = size(xi, 1)
    n = N

    # Initialize cluster assignments.
    clusters = collect(1:N)
    cluster_hours = Dict(x => [x] for x in clusters)
        
    while n > K
        cluster_centroids = Dict(x => calculate_centroid(xi[cluster_hours[x], :]) for x in clusters)

        adjacent_clusters = Tuple[]
        dissimilarity_matrix = Float64[]

        for I in clusters           
            for J in clusters
                if I > J && is_adjacent(cluster_hours[I], cluster_hours[J])
                    push!(adjacent_clusters, (I, J))
                    push!(dissimilarity_matrix, calculate_dissimilarity(xi, cluster_hours[I], cluster_hours[J]))
                end
            end
        end

        similar_clusters = adjacent_clusters[argmin(dissimilarity_matrix)]

        cluster_hours[similar_clusters[1]] = vcat(cluster_hours[similar_clusters[1]], cluster_hours[similar_clusters[2]])

        filter!(x -> x != similar_clusters[2], clusters)
        delete!(cluster_hours, similar_clusters[2])

        n = length(clusters)

    end

    # Step 7: Determine the set of representative periods as the clusters' centroids xI.
    unique_clusters = unique(clusters)
    centroids = Dict(x => calculate_centroid(xi[cluster_hours[x], :]) for x in unique_clusters)
    representative_periods = [centroids[i] for i in unique_clusters]

    # Step 8: The number of hours belonging to each cluster corresponds to the value of the time-period duration τt.
    time_period_durations = [size(clusters[clusters .== cluster], 1) for cluster in unique_clusters]

    # Print the representative periods and their corresponding time-period durations.
    #println("Representative Periods:", representative_periods)
    #println("Time-Period Durations:", time_period_durations)

    return cluster_hours
end


#= 
net_load_data = CSV.read(joinpath("data", "load_n_vg_data.csv"), DataFrame)

labels = [Dates.Date(net_load_data[i, :Year], net_load_data[i, :Month], net_load_data[i, :Day]) for i in 1:DataFrames.nrow(net_load_data)]
global load = vec(sum(Matrix(net_load_data[:, r"load"]), dims=2))
global wind = zeros(length(load))
global pv = zeros(length(load))

existing_generator_data = DataFrames.DataFrame(CSV.File(joinpath("data", "gen.csv")))

for i in names(net_load_data)
    if length(filter( row -> row."GEN UID" in [i], existing_generator_data)."Unit Type") > 0
        if occursin("WIND", filter( row -> row."GEN UID" in [i], existing_generator_data)."Unit Type"[1]) #occursin("WIND", i) || occursin("WT", i)
            global wind += net_load_data[:, i]
        elseif occursin("PV", filter( row -> row."GEN UID" in [i], existing_generator_data)."Unit Type"[1])#occursin("PV", i) || occursin("PVe", i)
            global pv += net_load_data[:, i]
        end
    end
end

normalized_load = normalize_vector(load)
normalized_wind = normalize_vector(wind)
normalized_pv = normalize_vector(pv)

features = hcat(load, wind, pv)
xi = hcat(normalized_load, normalized_wind, normalized_pv)

# Set the constant for the reduced number of time periods.
K = 240

@time clusters, cluster_hours = chronological_clustering(xi, K)
sorted_dict = sort(cluster_hours)

df = DataFrame(cluster = string.(keys(sorted_dict)), hours = string.(values(sorted_dict)))
CSV.write("cluster.csv", df)
 =#
