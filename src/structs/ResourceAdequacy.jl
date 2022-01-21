mutable struct ResourceAdequacy
    targets :: Dict{String, Float64}
    delta_irm :: Vector{Float64}
    metrics :: Vector{Dict{String, Float64}}
end

get_targets(ra::ResourceAdequacy) = ra.targets
get_delta_irm(ra::ResourceAdequacy) = ra.delta_irm
get_delta_irm(ra::ResourceAdequacy, iteration_year::Int64) = ra.delta_irm[iteration_year]
get_metrics(ra::ResourceAdequacy) = ra.metrics
get_metrics(ra::ResourceAdequacy, iteration_year::Int64) = ra.metrics[iteration_year]


function set_delta_irm!(ra::ResourceAdequacy, iteration_year::Int64, value::Float64)
    ra.delta_irm[iteration_year] = value
end

function set_metrics!(ra::ResourceAdequacy, iteration_year::Int64, value::Dict{String, Float64})
    ra.metrics[iteration_year] = value
end
