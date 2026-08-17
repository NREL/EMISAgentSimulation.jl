@testset "Test output point fields and heat rate fields" begin
    unit_type = String(projectdata["Unit Type"])
    size_raw = projectdata["Size"]
    if typeof(size_raw) !== Float64
        size_raw = String(size_raw)
    end

    size = Float64(size_in_MW(investor_dir,
        unit_type,
        size_raw))

    output_point_fields = String[]
    heat_rate_fields = String[]
    fields = names(projectdata)
    project_scale = 1e3

    for field in fields
        if occursin("Output_pct_", field)
            push!(output_point_fields, field)
        elseif occursin("HR_", field)
            push!(heat_rate_fields, field)
        end
    end

    @assert length(output_point_fields) > 0
end
