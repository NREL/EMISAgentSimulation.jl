"""
This function returns a vector of the products provided by the project
based on the data provided in the input files.
"""
function create_products(simulation_data::AgentSimulationData,
                         projectdata::DataFrames.DataFrameRow)
    products = Product[]

    markets = get_markets(simulation_data)

    reserve_penalty = get_reserve_penalty(get_case(simulation_data))
    mopr = get_mopr(get_case(simulation_data))
    bat_cap = get_battery_cap_mkt(get_case(simulation_data))

    simulation_horizon = get_simulation_years(get_case(simulation_data))

    # Energy market product
    variable_cost = projectdata["Fuel Price \$/MMBTU"] * projectdata["HR_avg_0"] / 1000
    push!(products, Energy(:Energy, Dict{String, Array{Float64, 2}}(), variable_cost, 0.0))

    data_dir = get_data_dir(get_case(simulation_data))
    reserve_definition = read_data(joinpath(data_dir, "markets_data", "reserve_products.csv"))

    system_cfg = load_system_config(data_dir)
    reserve_products = split(reserve_definition[1, "all_products"], "; ")

    for product in reserve_products
        product_data = read_data(joinpath(data_dir, "markets_data", "$(reserve_penalty)_reserve_penalty", "$(product).csv"))
        eligible_zones = [get_zone_name(system_cfg, n) for n in split(product_data[1, "eligible_zones"], ";")]

        if markets[Symbol(product)] && (projectdata["Zone"] in eligible_zones) && occursin(projectdata["Category"], product_data[1, "eligible categories"])
            time_scale = product_data[1, "timescale (min)"]
            direction = product_data[1, "direction"]
            if lowercase(direction) == "up"
                push!(products, OperatingReserve{ReserveUpEMIS}(Symbol(product),
                                        projectdata["Ramp Rate pu/Hr"] * time_scale / 60,
                                        0.0))
            elseif lowercase(direction) == "down"
                push!(products, OperatingReserve{ReserveDownEMIS}(Symbol(product),
                                        projectdata["Ramp Rate pu/Hr"] * time_scale / 60,
                                        0.0))
            end
        end
    end

    unit_type = String(projectdata["Unit Type"])
    technology_class = get_technology_class(system_cfg, unit_type)
    capacity_eligible = projectdata["Capacity Eligible"]
    if mopr && is_mopr_exempt(system_cfg, unit_type)
        capacity_eligible = false
    end

    if !bat_cap && technology_class == "storage"
        capacity_eligible = false
    end

    if markets[:Capacity] && capacity_eligible
        push!(products, Capacity(:Capacity,
                                 Dict{String, Dict{String, Float64}}(),
                                 Dict{String, Dict{String, Vector{Float64}}}(),
                                 0.0))
    end

    if markets[:REC] && projectdata["REC Eligible"]
        push!(products, REC(:REC,
                            0.0,
                            ones(simulation_horizon),
                            0.0))
    end

    push!(products, CarbonTax(:CarbonTax,
                              projectdata["CO2_Emissions ton/MMBTU"],
                              projectdata["HR_avg_0"] / 1000,
                              projectdata["Fuel Price \$/MMBTU"],
                              zeros(simulation_horizon)))

    if markets[:Inertia]
        push!(products, Inertia(:Inertia,
                                projectdata["Synchronous_Inertia"],
                                projectdata["Inertia MJ/MW"],
                                0.0))
    end
    return products
end

