"""
This function returns a vector of the products provided by the project
based on the data provided in the input files.
"""
function create_products(simulation_data::AgentSimulationData,
                         projectdata::DataFrames.DataFrameRow)
    products = Product[]

    markets = get_markets(simulation_data)

    reserve_up_mkt_df = read_data(joinpath(get_data_dir(get_case(simulation_data)), "markets_data", "reserve_up_mkt_param.csv"))
    reserve_down_mkt_df = read_data(joinpath(get_data_dir(get_case(simulation_data)), "markets_data", "reserve_down_mkt_param.csv"))

    reserve_up_eligible_categories = reserve_up_mkt_df[findfirst(x-> x == projectdata["Zone"], reserve_up_mkt_df[:, "zones"]), "eligible categories"]
    reserve_down_eligible_categories = reserve_down_mkt_df[findfirst(x-> x == projectdata["Zone"], reserve_down_mkt_df[:, "zones"]), "eligible categories"]

    # Energy market product
    variable_cost = projectdata["Fuel Price \$/MMBTU"] * projectdata["HR_avg_0"] / 1000
        push!(products, Energy(:Energy,
                               Dict{String, Array{Float64, 2}}(),
                               variable_cost))

    if markets[:Reserves] && occursin(projectdata["Category"], reserve_up_eligible_categories) > 0
        push!(products, OperatingReserve{ReserveUpEMIS}(:ReserveUp,
                                        projectdata["Ramp Rate pu/Hr"],
                                        0.0))
    end

    if markets[:Reserves] && occursin(projectdata["Category"], reserve_down_eligible_categories) > 0
        push!(products, OperatingReserve{ReserveDownEMIS}(:ReserveDown,
                                        projectdata["Ramp Rate pu/Hr"],
                                        0.0))
    end

    if markets[:Capacity] && projectdata["Capacity Eligible"]
        push!(products, Capacity(:Capacity,
                                 0.0,
                                 Dict{String, Array{Float64, 1}}(),
                                 0.0))
    end

    if markets[:REC] && projectdata["REC Eligible"]
        push!(products, REC(:REC,
                            0.0,
                            0.0))
    end

    return products
end


"""
This function returns a vector of the products provided by the project
based on the data provided in the input files and PSY Device data.
"""
function create_products(simulation_data::AgentSimulationData,
                         device::T,
                         projectdata::DataFrames.DataFrameRow) where T <: Union{PSY.Generator, PSY.Storage}

    products = Product[]

    markets = get_markets(simulation_data)

    service_types = typeof.(PSY.get_services(device))

    # Energy market product
    variable_cost = projectdata["Fuel Price \$/MMBTU"] * projectdata["HR_avg_0"] / 1000.0

    push!(products, Energy(:Energy,
                               Dict{String, Array{Float64, 2}}(),
                               variable_cost))

    fields = collect(fieldnames(typeof(device)))

    reserve_up_lim = 1.0
    reserve_down_lim = 1.0

    if in(:ramp_limits, fields)
        ramp_limits = PSY.get_ramp_limits(device)
        rating = get_device_size(device)
        if !isnothing(ramp_limits)
            reserve_up_lim = ramp_limits[:up] * 60 / rating
            reserve_down_lim = ramp_limits[:down] * 60 / rating
        end
    end



    if in(PSY.VariableReserve{PSY.ReserveUp}, service_types)
        push!(products, OperatingReserve{ReserveUpEMIS}(:ReserveUp,
                                        reserve_up_lim,
                                        0.0))
    end



    if in(PSY.VariableReserve{PSY.ReserveDown}, service_types)
        push!(products, OperatingReserve{ReserveDownEMIS}(:ReserveDown,
                                        reserve_down_lim,
                                        0.0))
    end

    if markets[:Capacity] && projectdata["Capacity Eligible"]
        push!(products, Capacity(:Capacity,
                                 0.0,
                                 Dict{String, Array{Float64, 1}}(),
                                 0.0))
    end

    if markets[:REC] && projectdata["REC Eligible"]
        push!(products, REC(:REC,
                            0.0,
                            0.0))
    end

    return products
end
