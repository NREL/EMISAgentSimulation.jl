"""
This functions returns an AxisArray of project parameters
included in CEM for price projection and endogeneous Economic Dispatch.
"""
function make_parameter_vector(
    structs::Vector, id::Symbol, parameter::Symbol)
    ids = getproperty.(structs, id)
    vals = getproperty.(structs, parameter)
    return AxisArrays.AxisArray(vals, ids)
end

"""
This functions returns an array of market parameters
included in CEM for price projection and endogeneous Economic Dispatch.
"""
function make_parameter_vector(
    structs::Vector{<:MarketCollection}, market::Symbol, parameter::Symbol)
    markets = getproperty.(structs, market)
    return getproperty.(markets, parameter)
end

"""
This functions returns a matrix of demand parameters
included in CEM for price projection and endogeneous Economic Dispatch.
"""
function make_demand_matrix(structs::Vector{<:MarketCollection}, market::Symbol)
    markets = getproperty.(structs, market)
    matrix = hcat(getproperty.(markets, :demand)...)
    return permutedims(matrix)
end

"""
This functions returns capacity market demand curve parameters
included in CEM for price projection and endogeneous Economic Dispatch.
"""
function make_capacity_demand_vectors(capmarkets::Vector{CapacityMarket})
   break_points =  getproperty.(capmarkets, :break_points)
   price_points =  getproperty.(capmarkets, :price_points)
   num_segments = length.(break_points) .- 1
   segment_size = [zeros(num_segments[i]) for i in 1:length(capmarkets)]
   segment_grad = [zeros(num_segments[i]) for i in 1:length(capmarkets)]
   for p in 1:length(capmarkets)
        for segment in 1:num_segments[p]
           segment_size[p][segment] = break_points[p][segment + 1] - break_points[p][segment]
           if segment_size[p][segment] != 0
                segment_grad[p][segment] = (price_points[p][segment + 1] - price_points[p][segment]) /  segment_size[p][segment]
           end
        end
    end

   return segment_size, segment_grad, price_points, num_segments
end

"""
This functions returns operating reserve demand curve parameters
included in CEM for price projection and endogeneous Economic Dispatch.
"""
function make_ORDC_vectors(ordc_markets::Vector{Dict{String, ReserveORDCMarket{T}}}) where T

    inv_periods = length(ordc_markets)

    products = unique(keys(ordc_markets[p]) for p in inv_periods)[1]

    num_segments = Dict(product => AxisArrays.AxisArray(Array{Int64, 2}(undef, inv_periods, T), 1:inv_periods, 1:T) for product in products)
    segment_size = Dict(product => AxisArrays.AxisArray(Array{Vector{Float64}, 2}(undef, inv_periods, T), 1:inv_periods, 1:T) for product in products)
    segment_grad = Dict(product => AxisArrays.AxisArray(Array{Vector{Float64}, 2}(undef, inv_periods, T), 1:inv_periods, 1:T) for product in products)
    price_points = Dict(product => AxisArrays.AxisArray(Array{Vector{Float64}, 2}(undef, inv_periods, T), 1:inv_periods, 1:T) for product in products)

    for product in products
        break_points =  [getproperty(ordc_markets[p][product], :break_points) for p in 1:inv_periods]
        price_points_raw = [getproperty(ordc_markets[p][product], :price_points) for p in 1:inv_periods]

        stepped = getproperty(ordc_markets[1][product], :stepped)

        for p in 1:inv_periods
            for t in 1:T
                num_segments[product][p, t] = length(break_points[p][t]) - 1
                segment_size[product][p, t] = zeros(num_segments[product][p, t])
                segment_grad[product][p, t] = zeros(num_segments[product][p, t])
                price_points[product][p, t] = price_points_raw[p][t]

                if stepped

                    for segment in 1:num_segments[product][p, t]
                        segment_size[product][p, t][segment] = break_points[p][t][segment + 1] - break_points[p][t][segment]
                        price_points[product][p, t][segment] = price_points[product][p, t][segment + 1]
                    end

                else

                    for segment in 1:num_segments[product][p, t]
                        segment_size[product][p, t][segment] = break_points[p][t][segment + 1] - break_points[p][t][segment]
                        if segment_size[product][p, t][segment] != 0
                            segment_grad[product][p, t][segment] = (price_points[product][p, t][segment + 1] - price_points[product][p, t][segment]) /  segment_size[product][p, t][segment]
                        end
                    end

                end
            end
        end

    end

    return segment_size, segment_grad, price_points, num_segments
 end

 """
 This functions returns capacity market demand curve parameters
 included in the actual clearing of capacity markets.
 """
function make_capacity_demand(capmarket::CapacityMarket)
    break_points =  getproperty(capmarket, :break_points)
    price_points =  getproperty(capmarket, :price_points)
    num_segments = length(break_points) - 1
    segment_size = zeros(num_segments)
    segment_grad = zeros(num_segments)
         for segment in 1:num_segments
            segment_size[segment] = break_points[segment + 1] - break_points[segment]
            segment_grad[segment] = (price_points[segment + 1] - price_points[segment]) /  segment_size[segment]
         end

    return segment_size, segment_grad, price_points
 end

 """
This functions returns capital cost multiplier curve parameters
included in CEM for price projection and endogeneous Economic Dispatch.
"""
function make_capital_cost_curve(options_by_type::Dict{String, Vector{String}},
                              annualized_cap_cost::AxisArrays.AxisArray{Float64, 2},
                              basecostunits::AxisArrays.AxisArray{Int64, 1},
                              maxnewoptions::AxisArrays.AxisArray{Int64, 1},
                              capital_cost_multiplier::Float64)

    invperiods = size(annualized_cap_cost, 2)
    types = collect(keys(options_by_type))

    num_segments = AxisArrays.AxisArray(round.(Int, 2 * ones(length(types), invperiods)), types, 1:invperiods)

    segment_size = AxisArrays.AxisArray(Array{Vector{Float64}, 2}(undef, length(types), invperiods), types, 1:invperiods)
    segment_grad = AxisArrays.AxisArray(Array{Vector{Float64}, 2}(undef, length(types), invperiods), types, 1:invperiods)
    costpoints = AxisArrays.AxisArray(Array{Vector{Float64}, 2}(undef, length(types), invperiods), types, 1:invperiods)

    for t in types
        for p in 1:invperiods
            segment_size[t, p] = zeros(num_segments[t, p])
            segment_grad[t, p] = zeros(num_segments[t, p])
            costpoints[t, p] = zeros(num_segments[t, p])

            if length(options_by_type[t]) >= 1
                total_max_new_options = sum(maxnewoptions[g] for g in options_by_type[t])
                g = options_by_type[t][1]
                segment_size[t, p][1] = basecostunits[g]
                segment_size[t, p][2] = total_max_new_options - basecostunits[g]
                costpoints[t, p] = [annualized_cap_cost[g, p], annualized_cap_cost[g, p]]
                segment_grad[t, p][1] = 0.0
                segment_grad[t, p][2] = costpoints[t, p][2] * capital_cost_multiplier
            end
        end
    end

    return segment_size, segment_grad, costpoints, num_segments
end

"""
This function updates the expected capacity factors of projects which are not in Option phase.
"""
function update_capacity_factors!(project::P,
                                 scenario_name::String,
                                 capacity_factors::Dict{String, Array{Float64, 2}}) where P <: Project{<: BuildPhase}
            for product in get_products(project)
                set_capacity_factors!(product, scenario_name, capacity_factors[get_name(project)])
            end
    return
end

"""
This function updates the expected capacity factors of Option projects.
"""
function update_capacity_factors!(project::P,
                                 scenario_name::String,
                                 capacity_factors::Dict{String, Array{Float64, 2}}) where P <: Project{Option}

            tech = get_tech(project)
            for product in get_products(project)
                set_capacity_factors!(product, scenario_name, capacity_factors["option_$(get_type(tech))_$(get_zone(tech))"])
            end
    return
end

"""
This function updates the expected capacity factors of projects which are not in Option phase.
"""
function update_total_utilization!(project::P,
                                 scenario_name::String,
                                 total_utilization::Dict{String, Array{Float64, 2}}) where P <: Project{<: BuildPhase}

    project_total_utilization = total_utilization[get_name(project)]
    project.finance_data.scenario_total_utilization[scenario_name] = project_total_utilization

    return
end

"""
This function updates the expected capacity factors of Option projects.
"""
function update_total_utilization!(project::P,
                                 scenario_name::String,
                                 total_utilization::Dict{String, Array{Float64, 2}}) where P <: Project{Option}

    tech = get_tech(project)
    project_total_utilization = total_utilization["option_$(get_type(tech))_$(get_zone(tech))"]
    project.finance_data.scenario_total_utilization[scenario_name] = project_total_utilization

    return
end

"""
This function updates the expected accepted capacity percentage of projects which are not in Option phase.
"""
function update_capacity_accepted_perc!(project::P,
                                 scenario_name::String,
                                 capacity_accepted_perc::Dict{String, Array{Float64, 1}}) where P <: Project{<: BuildPhase}
            for product in get_products(project)
                set_accepted_perc!(product, scenario_name, capacity_accepted_perc[get_name(project)])
            end
    return
end

"""
This function updates the expected accepted capacity percentage of Option projects.
"""
function update_capacity_accepted_perc!(project::P,
                                 scenario_name::String,
                                 capacity_accepted_perc::Dict{String, Array{Float64, 1}}) where P <: Project{Option}

            tech = get_tech(project)
            for product in get_products(project)
                set_accepted_perc!(product, scenario_name, capacity_accepted_perc["option_$(get_type(tech))_$(get_zone(tech))"])
            end
    return
end

"""
This function does nothing if the product is not of Capacity type.
"""
function update_initial_capacity_revenues!(project::P,
                                 product::T,
                                 initial_capacity_prices::Vector{Float64},
                                 year::Int64) where {P <: Project{<: BuildPhase}, T <: Product}
    return
end

"""
This function updates the capacity market revenues of project at the start of the simulation.
"""
function update_initial_capacity_revenues!(project::P,
                                 product::Capacity,
                                 initial_capacity_prices::Vector{Float64},
                                 year::Int64) where P <: Project{<: BuildPhase}
    capacity_revenue = initial_capacity_prices[year] * get_maxcap(project) * get_project_derating(project)

    finance_data = get_finance_data(project)
    set_realized_profit!(get_finance_data(project),
                    get_name(product),
                    year,
                    capacity_revenue)

    for scenario_name in keys(get_scenario_profit(finance_data))
                update_forward_profit!(product, finance_data, scenario_name, year, capacity_revenue)
            end

    return
end

"""
This function does nothing if the product is not of Capacity or REC type.
"""
function update_bid!(product::T,
                     capacity_market_bid::Float64,
                     rec_market_bid::Float64,
                     energy_production::Float64,
                     iteration_year::Int64) where T <: Product
    return
end

"""
This function updates the Capacity market bid of projects.
"""
function update_bid!(product::Capacity,
                     capacity_market_bid::Float64,
                     rec_market_bid::Float64,
                     energy_production::Float64,
                     iteration_year::Int64)

    set_capacity_bid!(product, capacity_market_bid)

    return
end

"""
This function updates the REC market bid of projects.
"""
function update_bid!(product::REC,
                     capacity_market_bid::Float64,
                     rec_market_bid::Float64,
                     energy_production::Float64,
                     iteration_year::Int64)

    set_rec_bid!(product, rec_market_bid)
    set_expected_rec_certificates!(product, energy_production * get_rec_correction_factor(product, iteration_year))
    return
end

function calculate_carbon_cost_ratio(product::Product, carbon_cost_ratio::Float64, carbon_tax::Vector{Float64}, year::Int64)
    return carbon_cost_ratio
end

function calculate_carbon_cost_ratio(product::CarbonTax, carbon_cost_ratio::Float64, carbon_tax::Vector{Float64}, year::Int64)
    carbon_cost = get_emission_intensity(product) * carbon_tax[year]
    ratio = carbon_cost / get_fuel_cost(product)
    if !isnan(ratio)
        carbon_cost_ratio = ratio
    end
    return carbon_cost_ratio
end

function update_device_operation_cost!(project::P, sys_UC::PSY.System, var_cost, fixed::Float64) where P <: Project{<:BuildPhase}
    return
end

function update_device_operation_cost!(project::P, sys_UC::PSY.System, var_cost, fixed::Float64) where P <: Project{Existing}
    name = get_name(project)
    device = PSY.get_components_by_name(PSY.Device, sys_UC, name)[1]
    device.operation_cost.variable.cost = var_cost
    device.operation_cost.fixed = fixed
    return
end

function update_operation_cost!(project::P, sys_UC::PSY.System, carbon_tax::Vector{Float64}, year::Int64) where P <: Project{<:BuildPhase}
    return
end

function update_operation_cost!(project::P, sys_UC::Nothing, carbon_tax::Vector{Float64}, year::Int64) where P <: Project{<:BuildPhase}
    return
end

function update_operation_cost!(project::P, sys_UC::PSY.System, carbon_tax::Vector{Float64}, year::Int64) where P <: Project{Existing}
    operation_cost = get_operation_cost(get_tech(project))
    if !(isnothing(operation_cost))
        products = get_products(project)
        carbon_cost_ratio = 0.0

        for product in products
            carbon_cost_ratio = calculate_carbon_cost_ratio(product, carbon_cost_ratio, carbon_tax, year)
        end

        total_cost_scalar = 1 + carbon_cost_ratio

        var_cost = deepcopy(PSY.get_variable(operation_cost).cost)

        fixed = deepcopy(PSY.get_fixed(operation_cost))

        if length(var_cost) > 1
            for i in 1:length(var_cost)
                var_cost[i] = (var_cost[i][1] * total_cost_scalar, var_cost[i][2])
            end
        elseif length(var_cost) == 1
            var_cost = var_cost * total_cost_scalar
        end

        fixed = fixed * total_cost_scalar

        update_device_operation_cost!(project, sys_UC, var_cost, fixed)
    end

    return
end
