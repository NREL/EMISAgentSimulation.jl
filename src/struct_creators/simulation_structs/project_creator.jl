# This file contains functions which help creating a Project.

"""
This function returns a dataframe containing project data if the data file exists,
otherwise it returns nothing.
"""
function extract_projectdata(inv_dir::String,
                        file_name::String)

    file_path = joinpath(inv_dir, file_name)
    if isfile(file_path)
        projectdata = read_data(file_path)
        return projectdata
    else
        return nothing
    end

end


"""
This function returns the queue cost data for a project
based on its size.
"""
function create_project_queuecost(size::Float64,
                                queue_cost_data::DataFrames.DataFrame,
                                lag_bool::Bool)
    if lag_bool
        D1_row = findfirst(x -> x == "D1", queue_cost_data[:, "Ref"])
        D2_row = findfirst(x -> x == "D2", queue_cost_data[:, "Ref"])
        M2_row = findfirst(x -> x == "M2 (per MW)", queue_cost_data[:, "Ref"])

        if size < 6
            D1 = queue_cost_data[D1_row, "less than 6 MW"]
            D2 = queue_cost_data[D2_row, "less than 6 MW"]
            M2 = queue_cost_data[M2_row, "less than 6 MW"] * size
        elseif 6 <= size <= 20
            D1 = queue_cost_data[D1_row, "6 to 20 MW"]
            D2 = queue_cost_data[D2_row, "6 to 20 MW"]
            M2 = queue_cost_data[M2_row, "6 to 20 MW"] * size
        elseif 20 < size <= 50
            D1 = queue_cost_data[D1_row, "20 to 50 MW"]
            D2 = queue_cost_data[D2_row, "20 to 50 MW"]
            M2 = queue_cost_data[M2_row, "20 to 50 MW"] * size
        elseif 50 < size <= 100
            D1 = queue_cost_data[D1_row, "50 to 100 MW"]
            D2 = queue_cost_data[D2_row, "50 to 100 MW"]
            M2 = queue_cost_data[M2_row, "50 to 100 MW"] * size
        elseif 100 < size <= 200
            D1 = queue_cost_data[D1_row, "100 to 200 MW"]
            D2 = queue_cost_data[D2_row, "100 to 200 MW"]
            M2 = queue_cost_data[M2_row, "100 to 200 MW"] * size
        elseif 200 < size <= 500
            D1 = queue_cost_data[D1_row, "200 to 500 MW"]
            D2 = queue_cost_data[D2_row, "200 to 500 MW"]
            M2 = queue_cost_data[M2_row, "200 to 500 MW"] * size
        elseif 500 < size < 1000
            D1 = queue_cost_data[D1_row, "500 to 1000 MW"]
            D2 = queue_cost_data[D2_row, "500 to 1000 MW"]
            M2 = queue_cost_data[M2_row, "500 to 1000 MW"] * size
        elseif size >= 1000
            D1 = queue_cost_data[D1_row, "greater than 1000 MW"]
            D2 = queue_cost_data[D2_row, "greater than 1000 MW"]
            M2 = queue_cost_data[M2_row, "greater than 1000 MW"] * size
        end
    else
        D1 = 0; D2 = 0; M2 = 0
    end
    queue_cost = [D1 + D2, M2]
    return queue_cost
end

"""
This function writes the project's hourly availability data
to the "DAY_AHEAD_availability.csv" file stored in
simulation_dir/timeseries_data_files/Availability
"""
function add_investor_project_availability!(simulation_dir::String,
                                              projects::Vector{Project},
                                              sys_UC::Union{Nothing, PSY.System})
    system_availability_data = DataFrames.DataFrame(CSV.File(joinpath(simulation_dir, "timeseries_data_files", "Availability", "DAY_AHEAD_availability.csv")))
    system_availability_data_rt = DataFrames.DataFrame(CSV.File(joinpath(simulation_dir, "timeseries_data_files", "Availability", "REAL_TIME_availability.csv")))
    gennames = names(system_availability_data)[5:length(names(system_availability_data))] #################

    for project in projects
        # println("this is $(get_name(project))")
        project_name = get_name(project)
        tech = get_tech(project)
        type_zone_id = "$(get_type(tech))_$(get_zone(tech))"

        if (!in(project_name, names(system_availability_data)) && !in(type_zone_id, names(system_availability_data)))
            if typeof(project) == ThermalGenEMIS{Existing} || typeof(project) == BatteryEMIS{Existing}
                system_availability_data[:, project_name] = ones(DataFrames.nrow(system_availability_data))
                system_availability_data_rt[:, project_name] = ones(DataFrames.nrow(system_availability_data_rt))
            elseif typeof(project) == RenewableGenEMIS{Existing}
                error("Timeseries for existing renewable generator not provided")
            elseif typeof(project) == ThermalGenEMIS{Option} || typeof(project) == BatteryEMIS{Option}
                system_availability_data[:, type_zone_id] = ones(DataFrames.nrow(system_availability_data))
                system_availability_data_rt[:, type_zone_id] = ones(DataFrames.nrow(system_availability_data_rt))
            elseif typeof(project) == RenewableGenEMIS{Option}
                availability = zeros(DataFrames.nrow(system_availability_data))
                availability_rt = zeros(DataFrames.nrow(system_availability_data_rt))
                bus = get_bus(tech)
                if get_type(tech) == "WT"
                    # gens_in_zone = filter(g -> (first("$(bus)", 1) == first(g, 1)) && (occursin("wind", lowercase(g)) || occursin("wt", lowercase(g))),
                    #                                 names(system_availability_data))
                    gens_in_zone = filter(g -> (first("$(bus)", 1) == first(string(PSY.get_number(PSY.get_bus(PSY.get_component(PSY.StaticInjection, sys_UC,g)))), 1)) &&
                                                occursin("WT", string(PSY.get_prime_mover(PSY.get_component(PSY.StaticInjection, sys_UC,g)))),
                                                    gennames)
                elseif get_type(tech) == "PVe"
                    # gens_in_zone = filter(g -> (first("$(bus)", 1)== first(g, 1)) && occursin("pv", lowercase(g)),
                    #                                 names(system_availability_data))
                    gens_in_zone = filter(g -> (first("$(bus)", 1) == first(string(PSY.get_number(PSY.get_bus(PSY.get_component(PSY.StaticInjection, sys_UC,g)))), 1)) &&
                                                occursin("PVe", string(PSY.get_prime_mover(PSY.get_component(PSY.StaticInjection, sys_UC,g)))),
                                                    gennames)
                elseif get_type(tech) == "HY"
                    # gens_in_zone = filter(g -> (first("$(bus)", 1) == first(g, 1)) && occursin("hy", lowercase(g)),
                    #                                 names(system_availability_data))
                    gens_in_zone = filter(g -> (first("$(bus)", 1) == first(string(PSY.get_number(PSY.get_bus(PSY.get_component(PSY.StaticInjection, sys_UC,g)))), 1)) &&
                                                occursin("HY", string(PSY.get_prime_mover(PSY.get_component(PSY.StaticInjection, sys_UC,g)))),
                                                    gennames)
                end

                if length(gens_in_zone) < 1
                    error("Similar renewable generator in zone not found")
                else
                    for g in gens_in_zone
                        availability += system_availability_data[:, g]
                        availability_rt += system_availability_data_rt[:, g]
                    end
                    system_availability_data[:, type_zone_id] = availability / length(gens_in_zone)
                    system_availability_data_rt[:, type_zone_id] = availability_rt / length(gens_in_zone)
                end
            end
        end
    end

    write_data(joinpath(simulation_dir, "timeseries_data_files", "Availability"), "DAY_AHEAD_availability.csv", system_availability_data)
    write_data(joinpath(simulation_dir, "timeseries_data_files", "Availability"), "REAL_TIME_availability.csv", system_availability_data_rt)

    return
end



"""
This function creates a piece-wise linear operation cost for the project according to the structure used in PSY.
This allows smooth conversion of local Project structs to PSY Device structs.
Returns PSY.TwoPartCost or PSY.ThreePartCost based on the type of the Project.
"""
function create_operation_cost(projectdata::DataFrames.DataFrameRow, size::Float64, base_power::Float64, sys_base_power::Float64, products::Vector{Product}, carbon_tax::Vector{Float64})

    output_point_fields = String[]
    heat_rate_fields = String[]

    fields = names(projectdata)

    for field in fields
        if occursin("Output_pct_", field)
            push!(output_point_fields, field)
        elseif occursin("HR_", field)
            push!(heat_rate_fields, field)
        end
    end

    @assert length(output_point_fields) > 0
    cost_colnames = zip(heat_rate_fields, output_point_fields)

    fuel_cost =  projectdata["Fuel Price \$/MMBTU"] / 1000.0
    var_cost = [(projectdata[hr], projectdata[mw]) for (hr, mw) in cost_colnames]

    var_cost = unique([
        (tryparse(Float64, string(c[1])), tryparse(Float64, string(c[2])))
        for c in var_cost if !in("NA", c)
    ])

    if length(var_cost) > 1
        var_cost[2:end] = [
            (
                var_cost[i][1] *
                (var_cost[i][2] - var_cost[i - 1][2]) *
                fuel_cost,
                var_cost[i][2] / base_power,
            ) .* size for i in 2:length(var_cost)
        ]
        var_cost[1] =
            (var_cost[1][1] * var_cost[1][2] * fuel_cost, var_cost[1][2] / base_power) .* size

        fixed = max(
            0.0,
            var_cost[1][1] -
            (var_cost[2][1] / (var_cost[2][2] - var_cost[1][2]) * var_cost[1][2]),
        )
        var_cost[1] = (var_cost[1][1] - fixed, var_cost[1][2] * size)
        for i in 2:length(var_cost)
            var_cost[i] = (var_cost[i - 1][1] + var_cost[i][1], var_cost[i][2] * size)
        end
    elseif length(var_cost) == 1
        # if there is only one point, use it to determine the constant $/MW cost
        var_cost = var_cost[1][1] * var_cost[1][2] * fuel_cost * base_power
        fixed = 0.0
    else
        var_cost = 0.0
        fixed = 0.0
    end

    start_up_cost = fuel_cost * projectdata["Start Heat Cold MBTU"] * 1000.0
    shut_down_cost = 0.0

    if occursin("CC", projectdata["Unit Type"]) || occursin("CT", projectdata["Unit Type"]) || occursin("GT", projectdata["Unit Type"]) || occursin("ST", projectdata["Unit Type"]) || occursin("NUCLEAR", projectdata["Unit Type"])
        operation_cost = PSY.ThreePartCost(var_cost, fixed, start_up_cost, shut_down_cost)
    else
        operation_cost = PSY.TwoPartCost(var_cost, fixed)
    end

    return operation_cost
end

"""
This function populates project investment parameters and creates Finance struct.
"""
function create_investment_data(size::Float64,
                              projectdata::DataFrames.DataFrameRow,
                              investor_dir::String,
                              simulation_data::AgentSimulationData,
                              products::Vector{<: Product},
                              investor_name::String,
                              scenario_names::Vector{String},
                              lag_bool::Bool
                              )
    start_year = get_start_year(get_case(simulation_data))
    max_horizon = get_total_horizon(get_case(simulation_data))

    decision_year = 1
    queue_cost_data = get_queue_cost_data(simulation_data)
    queue_cost = create_project_queuecost(size, queue_cost_data, lag_bool)
    queue_time = length(queue_cost)
    construction_year = decision_year + (queue_time + projectdata["Lagtime"]) * lag_bool

    capex_years = projectdata["Capex Years"]

    capex_data = read_data(joinpath(investor_dir, "project_capex.csv"))

    category = String(projectdata["Category"])
    capex_row = findfirst(x-> x == category, capex_data[:, "Category"])

    if lag_bool
        online_year = start_year + decision_year + queue_time - 1
        age = 0
        remaining_lifetime = projectdata["Lifetime"]
        project_capex = [capex_data[capex_row, "$(start_year + y - 1)"] for y in 1:max_horizon]
        preference_multiplier = projectdata["Preference Multiplier"] * ones(max_horizon)
    else
        online_year = projectdata["Online Year"]
        age = start_year - online_year
        remaining_lifetime = projectdata["Lifetime"] - age
        project_capex = zeros(max_horizon)

        for y = 1:min(max_horizon, (capex_years - age))
            project_capex[y] = capex_data[capex_row, "pre $(start_year)"]
        end

        preference_multiplier = ones(max_horizon)
    end

    investment_cost, discount_rate = calculate_capcost_and_wacc(project_capex, category, investor_dir, online_year)

    investment_cost = investment_cost * size * 1000.0

    effective_investment_cost = investment_cost[queue_time + 1]
    fixedOM_cost = projectdata["Fixed OM Cost per MW"] * size

    end_life_year = (queue_time + projectdata["Lagtime"]) * lag_bool + remaining_lifetime

    retirement_year = (queue_time + projectdata["Lagtime"]) * lag_bool + remaining_lifetime

    num_profit_years = max(max_horizon, end_life_year)

    scenario_total_utilization = Dict{String, Array{Float64, 2}}()

    yearly_profit = [AxisArrays.AxisArray(zeros(length(products), num_profit_years),
                                                 AxisArrays.Axis{:prod}(get_name.(products)),
                                                 AxisArrays.Axis{:year}(1:1:num_profit_years))
                                                 for y in 1:num_profit_years]

    scenario_profit = Dict(scenario_name => [AxisArrays.AxisArray(zeros(length(products), num_profit_years),
                                                 AxisArrays.Axis{:prod}(get_name.(products)),
                                                 AxisArrays.Axis{:year}(1:1:num_profit_years))
                                                 for y in 1:num_profit_years]  for scenario_name in scenario_names)

    realized_profit = AxisArrays.AxisArray(zeros(length(products), num_profit_years),
                                                 AxisArrays.Axis{:prod}(get_name.(products)),
                                                 AxisArrays.Axis{:year}(1:1:num_profit_years))

    scenario_npv = Dict(scenario_name => zeros(num_profit_years) for scenario_name in scenario_names)
    expected_npv = zeros(num_profit_years)
    scenario_utility = Dict(scenario_name => zeros(num_profit_years) for scenario_name in scenario_names)
    expected_utility = zeros(num_profit_years)

    annual_cashflow = zeros(num_profit_years)

    finance_data = Finance(investment_cost,
                           effective_investment_cost,
                           preference_multiplier,
                           projectdata["Lagtime"],
                           remaining_lifetime,
                           capex_years,
                           fixedOM_cost,
                           queue_cost,
                           scenario_total_utilization,
                           scenario_profit,
                           realized_profit,
                           discount_rate,
                           scenario_npv,
                           expected_npv,
                           scenario_utility,
                           expected_utility,
                           annual_cashflow,
                           investor_name)

    return decision_year, construction_year, retirement_year, end_life_year, finance_data
end

"""
This function creates the Tech struct for projects
based only on the CSV data.
"""
function create_tech_type(name::String,
                          projectdata::DataFrames.DataFrameRow,
                          size::Float64,
                          base_power::Float64,
                          decision_year::Int64,
                          construction_year::Int64,
                          retirement_year::Int64,
                          end_life_year::Int64,
                          products::Vector{Product},
                          finance_data::Finance,
                          sys_UC::Union{PSY.System, Nothing},
                          carbon_tax::Vector{Float64}
                        )

    type = projectdata["Unit Type"]
    min_cap = projectdata["Min Gen pu"] * size

    active_power_limits = (min = min_cap, max = size)
    ramp_limits = (up = projectdata["Ramp Rate pu/Hr"]  * size, down = projectdata["Ramp Rate pu/Hr"] * size)
    fuel_cost = projectdata["Fuel Price \$/MMBTU"]
    zone = projectdata["Zone"]
    bus = projectdata["Bus ID"]
    FOR = projectdata["FOR"]
    MTTR = projectdata["MTTR Hr"]

    if !isnothing(sys_UC)
        sys_base_power = PSY.get_base_power(sys_UC)
        operation_cost = create_operation_cost(projectdata, size, base_power, sys_base_power, products, carbon_tax)
        up_down_time = (up = projectdata["Min Up Time Hr"], down = projectdata["Min Down Time Hr"])
    else
        operation_cost = nothing
        up_down_time = nothing
    end

    if type == "ST" || type == "CT" || type == "CC" || type == "NU_ST" || type == "GT" || type == "RE_CT"

        output_point_fields = String[]
        heat_rate_fields = String[]
        fields = names(projectdata)

        for field in fields
            if occursin("Output_pct_", field)
                push!(output_point_fields, field)
            elseif occursin("HR_", field)
                push!(heat_rate_fields, field)
            end
        end

        @assert length(output_point_fields) > 0
        cost_colnames = zip(heat_rate_fields, output_point_fields)

        heat_rate = [(projectdata[hr], projectdata[mw]) for (hr, mw) in cost_colnames]

        heat_rate = unique([
            (tryparse(Float64, string(c[1])), tryparse(Float64, string(c[2])))
            for c in heat_rate if !in("NA", c)
        ])

        if length(heat_rate) > 1
            heat_rate[2:end] = [
                (
                    heat_rate[i][1] *
                    (heat_rate[i][2] - heat_rate[i - 1][2]),
                    heat_rate[i][2]
                ) for i in 2:length(heat_rate)
            ]
            heat_rate[1] =
                (heat_rate[1][1] * heat_rate[1][2], heat_rate[1][2])

                fixed = max(
                    0.0,
                    heat_rate[1][1] -
                    (heat_rate[2][1] / (heat_rate[2][2] - heat_rate[1][2]) * heat_rate[1][2]),
                )
                heat_rate[1] = (heat_rate[1][1] / 1000.0, heat_rate[1][2] * size)
                for i in 2:length(heat_rate)
                    heat_rate[i] = (heat_rate[i - 1][1] + heat_rate[i][1] / 1000.0, heat_rate[i][2] * size)
                end
                pushfirst!(heat_rate, (fixed / 1000.0, 0.0))
        elseif length(heat_rate) == 1
            # if there is only one point, use it to determine the constant $/MW cost
            heat_rate = (heat_rate[1][1] * heat_rate[1][2] / 1000.0, heat_rate[1][2] * size)
            fixed = 0.0
        else
            heat_rate = [(0.0, 0.0)]
            fixed = 0.0
        end
        tech = ThermalTech(type,
                       projectdata["Fuel"],
                       active_power_limits,
                       ramp_limits,
                       up_down_time,
                       operation_cost,
                       fuel_cost,
                       heat_rate,
                       bus,
                       zone,
                       FOR,
                       MTTR)

        return  ThermalGenEMIS{BuildPhase}(name,
                                 tech,
                                 decision_year,
                                 construction_year,
                                 retirement_year,
                                 end_life_year,
                                 products,
                                 finance_data)
    elseif type == "WT" || type == "PVe"
        tech = RenewableTech(type,
                        active_power_limits,
                        ramp_limits,
                        operation_cost,
                        bus,
                        zone,
                        FOR,
                        MTTR)

        return  RenewableGenEMIS{BuildPhase}(name,
                                 tech,
                                 decision_year,
                                 construction_year,
                                 retirement_year,
                                 end_life_year,
                                 products,
                                 finance_data)

    elseif type == "HY"
        tech = HydroTech(type,
                       active_power_limits,
                       ramp_limits,
                       up_down_time,
                       operation_cost,
                       bus,
                       zone,
                       FOR,
                       MTTR)

        return  HydroGenEMIS{BuildPhase}(name,
                                 tech,
                                 decision_year,
                                 construction_year,
                                 retirement_year,
                                 end_life_year,
                                 products,
                                 finance_data)

    elseif type == "BA"
        storage_capacity = parse(Float64, projectdata["Duration Hr"]) * size
        tech = BatteryTech(type,
                        (min = 0.0, max = parse(Float64, projectdata["Input Power Rating pu"]) * size),
                        (min = 0.0, max = parse(Float64, projectdata["Output Power Rating pu"]) * size),
                        ramp_limits,
                        (min = parse(Float64, projectdata["Min Storgae pu"]) * storage_capacity, max = storage_capacity),
                        0.0,
                        (in = parse(Float64, projectdata["Round Trip Efficiency pu"]), out = 1.0),
                        bus,
                        zone,
                        FOR,
                        MTTR)
        return BatteryEMIS{BuildPhase}(name,
                        tech,
                        decision_year,
                        construction_year,
                        retirement_year,
                        end_life_year,
                        products,
                        finance_data)

    end

end

"""
This function creates the Tech struct for ThermalGen projects
based on the CSV and PSY Device data.
"""
function create_tech_type(name::String,
                          device::T,
                          projectdata::DataFrames.DataFrameRow,
                          size::Float64,
                          base_power::Float64,
                          decision_year::Int64,
                          construction_year::Int64,
                          retirement_year::Int64,
                          end_life_year::Int64,
                          products::Vector{Product},
                          finance_data::Finance
                        ) where T <: PSY.ThermalGen

    min_cap = deepcopy(PSY.get_active_power_limits(device)[:min]) * base_power

    fuel_cost = projectdata["Fuel Price \$/MMBTU"]
    bus = deepcopy(PSY.get_bus(device))
    device_ramp_limits = deepcopy(PSY.get_ramp_limits(device))

    if isnothing(device_ramp_limits)
        ramp_limits = nothing
    else
        ramp_limits = (up = device_ramp_limits[:up] * base_power * 60, down = device_ramp_limits[:down] * base_power * 60)
    end

    prime_mover = string(deepcopy(PSY.get_prime_mover(device)))
    fuel = string(deepcopy(PSY.get_fuel(device)))

    if fuel == "NUCLEAR"
        prime_mover = "NU_ST"
    end

    FOR = projectdata["FOR"]
    MTTR = projectdata["MTTR Hr"]

    output_point_fields = String[]
    heat_rate_fields = String[]
    fields = names(projectdata)

    for field in fields
        if occursin("Output_pct_", field)
            push!(output_point_fields, field)
        elseif occursin("HR_", field)
            push!(heat_rate_fields, field)
        end
    end

    @assert length(output_point_fields) > 0
    cost_colnames = zip(heat_rate_fields, output_point_fields)

    heat_rate = [(projectdata[hr], projectdata[mw]) for (hr, mw) in cost_colnames]

    heat_rate = unique([
        (tryparse(Float64, string(c[1])), tryparse(Float64, string(c[2])))
        for c in heat_rate if !in("NA", c)
    ])
    if length(heat_rate) > 1
        heat_rate[2:end] = [
            (
                heat_rate[i][1] *
                (heat_rate[i][2] - heat_rate[i - 1][2]),
                heat_rate[i][2]
            ) for i in 2:length(heat_rate)
        ]
        heat_rate[1] =
            (heat_rate[1][1] * heat_rate[1][2], heat_rate[1][2])

            fixed = max(
                0.0,
                heat_rate[1][1] -
                (heat_rate[2][1] / (heat_rate[2][2] - heat_rate[1][2]) * heat_rate[1][2]),
            )
            heat_rate[1] = (heat_rate[1][1] / 1000.0, heat_rate[1][2] * size)
            for i in 2:length(heat_rate)
                heat_rate[i] = (heat_rate[i - 1][1] + heat_rate[i][1] / 1000.0, heat_rate[i][2] * size)
            end
            pushfirst!(heat_rate, (fixed / 1000.0, 0.0))
    elseif length(heat_rate) == 1
        # if there is only one point, use it to determine the constant $/MW cost
        heat_rate = (heat_rate[1][1] * heat_rate[1][2] / 1000.0, heat_rate[1][2] * size)
        fixed = 0.0
    else
        heat_rate = [(0.0, 0.0)]
        fixed = 0.0
    end


    tech = ThermalTech(prime_mover,
                       fuel,
                       (min = min_cap, max = size),
                       ramp_limits,
                       deepcopy(PSY.get_time_limits(device)),
                       deepcopy(PSY.get_operation_cost(device)),
                       fuel_cost,
                       heat_rate,
                       deepcopy(PSY.get_number(bus)),
                       projectdata["Zone"],
                       FOR,
                       MTTR)

    return  ThermalGenEMIS{BuildPhase}(name,
                                 tech,
                                 decision_year,
                                 construction_year,
                                 retirement_year,
                                 end_life_year,
                                 products,
                                 finance_data)

end

"""
This function creates the Tech struct for RenewableGen projects
based on the CSV and PSY Device data.
"""
function create_tech_type(name::String,
                          device::T,
                          projectdata::DataFrames.DataFrameRow,
                          size::Float64,
                          base_power::Float64,
                          decision_year::Int64,
                          construction_year::Int64,
                          retirement_year::Int64,
                          end_life_year::Int64,
                          products::Vector{Product},
                          finance_data::Finance
                        ) where T <: PSY.RenewableGen

    bus = deepcopy(PSY.get_bus(device))

    type = string(deepcopy(PSY.get_prime_mover(device)))

    FOR = projectdata["FOR"]
    MTTR = projectdata["MTTR Hr"]

    tech = RenewableTech(type,
                       (min = 0.0, max = size),
                       (up = size, down = size),
                       deepcopy(PSY.get_operation_cost(device)),
                       deepcopy(PSY.get_number(bus)),
                       projectdata["Zone"],
                       FOR,
                       MTTR)

    return  RenewableGenEMIS{BuildPhase}(name,
                                 tech,
                                 decision_year,
                                 construction_year,
                                 retirement_year,
                                 end_life_year,
                                 products,
                                 finance_data)

end

"""
This function creates the Tech struct for HydroGen projects
based on the CSV and PSY Device data.
"""
function create_tech_type(name::String,
                          device::P,
                          projectdata::DataFrames.DataFrameRow,
                          size::Float64,
                          base_power::Float64,
                          decision_year::Int64,
                          construction_year::Int64,
                          retirement_year::Int64,
                          end_life_year::Int64,
                          products::Vector{Product},
                          finance_data::Finance
                        ) where P <: PSY.HydroGen

    min_cap = deepcopy(PSY.get_active_power_limits(device)[:min]) * base_power
    bus = deepcopy(PSY.get_bus(device))

    device_ramp_limits = deepcopy(PSY.get_ramp_limits(device))

    if isnothing(device_ramp_limits)
        ramp_limits = nothing
    else
        ramp_limits = (up = device_ramp_limits[:up] * base_power * 60, down = device_ramp_limits[:down] * base_power * 60)
    end

    FOR = projectdata["FOR"]
    MTTR = projectdata["MTTR Hr"]

    tech = HydroTech("HY",
                       (min = min_cap, max = size),
                       ramp_limits,
                       deepcopy(PSY.get_time_limits(device)),
                       deepcopy(PSY.get_operation_cost(device)),
                       deepcopy( PSY.get_number(bus)),
                       projectdata["Zone"],
                       FOR,
                       MTTR)

        return  HydroGenEMIS{BuildPhase}(name,
                                 tech,
                                 decision_year,
                                 construction_year,
                                 retirement_year,
                                 end_life_year,
                                 products,
                                 finance_data)

end

"""
This function creates the Tech struct for Battery projects
based on the CSV and PSY Device data.
"""
function create_tech_type(name::String,
                          device::PSY.GenericBattery,
                          projectdata::DataFrames.DataFrameRow,
                          size::Float64,
                          base_power::Float64,
                          decision_year::Int64,
                          construction_year::Int64,
                          retirement_year::Int64,
                          end_life_year::Int64,
                          products::Vector{Product},
                          finance_data::Finance
                        )

    bus = deepcopy(PSY.get_bus(device))

    type = string(deepcopy(PSY.get_prime_mover(device)))
    input_active_power_limits = deepcopy(PSY.get_input_active_power_limits(device))
    output_active_power_limits = deepcopy(PSY.get_output_active_power_limits(device))
    storage_capacity = deepcopy(PSY.get_state_of_charge_limits(device))

    FOR = projectdata["FOR"]
    MTTR = projectdata["MTTR Hr"]

    tech = BatteryTech(type,
                       (min = input_active_power_limits[:min] * base_power, max = input_active_power_limits[:max] * base_power),
                       (min = output_active_power_limits[:min] * base_power, max = output_active_power_limits[:max] * base_power),
                       (up = size, down = size),
                       (min = storage_capacity[:min] * base_power, max = storage_capacity[:max] * base_power), #originally has *size
                       PSY.get_initial_energy(device) * size * base_power,
                       PSY.get_efficiency(device),
                       PSY.get_number(bus),
                       projectdata["Zone"],
                       FOR,
                       MTTR)

    return  BatteryEMIS{BuildPhase}(name,
                                 tech,
                                 decision_year,
                                 construction_year,
                                 retirement_year,
                                 end_life_year,
                                 products,
                                 finance_data)

end

"""
This function returns a Project struct based only on project's CSV data.
"""
function create_project(projectdata::DataFrames.DataFrameRow,
                          simulation_data::AgentSimulationData,
                          investor_name::String,
                          investor_dir::String,
                          scenario_names::Vector{String},
                          lag_bool::Bool)

    name = String(projectdata["GEN_UID"])

    unit_type = String(projectdata["Unit Type"])
    size_raw = projectdata["Size"]
    if typeof(size_raw) !== Float64
        size_raw = String(size_raw)
    end

    size = Float64(size_in_MW(investor_dir,
                     unit_type,
                     size_raw))

    base_power = size
    carbon_tax = get_carbon_tax(simulation_data)

    products = create_products(simulation_data,
                               projectdata)

    decision_year,
    construction_year,
    retirement_year,
    end_life_year,
    finance_data = create_investment_data(size, projectdata, investor_dir, simulation_data, products, investor_name, scenario_names, lag_bool)

    sys_UC = get_system_UC(simulation_data)

    project = create_tech_type(name, projectdata, size, base_power, decision_year, construction_year, retirement_year, end_life_year, products, finance_data, sys_UC, carbon_tax)

    return project
end



"""
This function returns a Project struct based on CSV data and PSY data.
"""
function create_project(projectdata::DataFrames.DataFrameRow,
                          device::T,
                          simulation_data::AgentSimulationData,
                          investor_name::String,
                          investor_dir::String,
                          scenario_names::Vector{String},
                          lag_bool::Bool) where T <: Union{PSY.Generator, PSY.Storage}

    name = get_name(device)

    base_power = PSY.get_base_power(device)

    size = get_device_size(device) * base_power

    products = create_products(simulation_data,
                               projectdata)

    decision_year,
    construction_year,
    retirement_year,
    end_life_year,
    finance_data = create_investment_data(size, projectdata, investor_dir, simulation_data, products, investor_name, scenario_names, lag_bool)

    project = create_tech_type(name, device, projectdata, size, base_power, decision_year, construction_year, retirement_year, end_life_year, products, finance_data)

    return project
end
