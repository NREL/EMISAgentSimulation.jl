
"""
This function returns all the subtypes at the leaf nodes (i.e. which do not have any further subtypes) of a given type.
"""
function leaftypes(t::Any, level=1, leaves = Vector{Any}())
    for s in InteractiveUtils.subtypes(t)
        if length(InteractiveUtils.subtypes(s)) == 0
            push!(leaves, s)
        else
            leaftypes(s, level+1, leaves)
        end
    end
    return leaves
end


"""
This function gets the capacity market forward years for the simulation.
"""
function get_capacity_forward_years(sim::Union{AgentSimulation, AgentSimulationData})
    forward_years = 1
    if get_markets(sim)[:Capacity]
        capacity_mkt_params = DataFrames.DataFrame(CSV.File(joinpath(get_data_dir(get_case(sim)), "markets_data", "capacity_mkt_param.csv")))
        forward_years = capacity_mkt_params.forward_years[1]
    end
    return forward_years
end



"""
This function returns the size of the project (in MW)
based on the project type and size type.
"""
function size_in_MW(investor_dir::String,
                    projecttype::String,
                    sizetype::String)
    filename = joinpath(investor_dir, "sizedict.csv")
    sizedata = read_data(filename)

    size = sizedata[findfirst(x-> x == projecttype, sizedata.type),
                findfirst(x -> x == sizetype, names(sizedata))]
    return size
end

"""
This function returns the size of the project (in MW)
if already specified in MW.
"""
function size_in_MW(investor_dir::String,
                    projecttype::String,
                    size::Float64)
    return size
end

"""
This function returns the size of the PSY Device.
"""
get_device_size(device::PSY.ThermalStandard) = PSY.get_active_power_limits(device)[:max]
get_device_size(device::PSY.RenewableDispatch) = PSY.get_rating(device)
get_device_size(device::PSY.HydroEnergyReservoir) = PSY.get_active_power_limits(device)[:max]
get_device_size(device::PSY.GenericBattery) = PSY.get_rating(device)

"""
This function returns line rating for PSY Lines.
"""
get_line_rating(line::PSY.Line) = PSY.get_rate(line)
get_line_rating(line::PSY.HVDCLine) = line.active_power_limits_from[:max]

"""
This function removes leap days from the data if year is not leap year.
"""
function remove_leap_day!(df:: DataFrames.DataFrame, start_year::Int64)
    if start_year % 4 != 0
        filter!(row -> row["Month"] != 2 || row["Day"] != 29, df)
    end

    return
end

"""
This function finds the start and end years
for reading price data.
"""
function find_price_years(construction_year::Int64,
                          end_life_year::Int64,
                          iteration_year::Int64,
                          yearly_horizon::Int64)

    price_years = (start_year = max((construction_year -
                                    iteration_year + 1),
                                    1),
                   end_year = min((end_life_year -
                                  iteration_year + 1),
                                  (yearly_horizon)))

    return price_years
end

"""
This function finds the start and end years
for updating project profit and npv data.
"""
function find_update_years(construction_year::Int64,
                          end_life_year::Int64,
                          iteration_year::Int64,
                          yearly_horizon::Int64)

    update_years = (start_year = max(construction_year,
                                     iteration_year),
                    end_year = min(end_life_year,
                                   yearly_horizon + iteration_year - 1))

    return update_years
end

"""
This function returns the number of hours modeled in the market clearing problem.
"""
function get_num_hours(system::MarketClearingProblem{Z, T}) where {Z, T}
    num_hours = T

    return num_hours
end

"""
This function returns an O&M cost array filled with the fixed
annual O&M cost values from update year to the end life year.
"""
function create_OMcost_array(finance::Finance,
                            scenario_name::String,
                            project_life_end::Int64,
                            update_start_year::Int64,
                            iteration_year::Int64)

    OM_cost_array = zeros(size(get_scenario_profit(finance)[scenario_name][iteration_year], 2))

    for i in update_start_year:project_life_end
        OM_cost_array[i] = get_fixed_OM_cost(finance)
    end

    return OM_cost_array

end

"""
This function return a queue cost array filled with the queue cost values
corresponding to each year spent in the queue.
"""
function create_queue_cost_array(size_profit::Int64,
         decision_year::Int64,
         queue_time::Int64,
         queue_cost::Vector{Float64})


    queue_cost_array = zeros(size_profit)

    for i in decision_year:(decision_year + queue_time - 1)
        queue_cost_array[i] = queue_cost[i - decision_year + 1]
    end

    return queue_cost_array

end
"""
This function does nothing if the project is not an Option project.
"""
function update_lifecycle!(project::P,
                           iteration_year::Int64,
                           simulation_years::Int64) where P <: Project{<: BuildPhase}
    return
end

"""
This function updates the decision, construction and endlife years of Option projects.
"""
function update_lifecycle!(project::P,
                           iteration_year::Int64,
                           simulation_years::Int64) where P <: Project{Option}


        finance_data = get_finance_data(project)

        queue_time = length(get_queue_cost(finance_data))

        if queue_time + iteration_year <  simulation_years
            set_decision_year!(project, (iteration_year + 1))
            set_construction_year!(project, (get_decision_year(project) +
                                            queue_time +
                                            get_lag_time(finance_data)))

            end_life_year = get_construction_year(project) + get_life_time(finance_data)  - 1
            set_end_life_year!(project, end_life_year)

            set_effective_investment_cost!(project, get_investment_cost(finance_data)[queue_time + iteration_year + 1])

            products = get_products(project)
            finance_data.realized_profit = AxisArrays.AxisArray(zeros(length(products), end_life_year),
                                                    AxisArrays.Axis{:prod}(get_name.(products)),
                                                    AxisArrays.Axis{:year}(1:1:end_life_year))

            finance_data.annual_cashflow = zeros(end_life_year)

            scenario_profit = get_scenario_profit(finance_data)

            for scenario_name in keys(scenario_profit)
                rows = size(scenario_profit[scenario_name][iteration_year], 1)
                cols = get_end_life_year(project)

                set_expected_npv!(finance_data, zeros(cols))
                set_expected_utility!(finance_data, zeros(cols))

                finance_data.scenario_profit[scenario_name] = [AxisArrays.AxisArray(zeros(length(products), cols),
                                                     AxisArrays.Axis{:prod}(get_name.(products)),
                                                     AxisArrays.Axis{:year}(collect(1:1:cols)))
                                                     for y in 1:cols]
                set_scenario_npv!(finance_data, scenario_name, zeros(cols))
                set_scenario_utility!(finance_data, scenario_name, zeros(cols))
            end
        end

    return
end

"""
This function updates the total installed capacity each year.
"""
function update_installed_cap!(installed_capacity::Vector{Float64},
                               announced_projects::Vector{<: Project{<: BuildPhase}},
                               iteration_year::Int64,
                               simulation_years::Int64)

    for year in iteration_year:simulation_years
        installed_capacity[year] = 0.0
        for project in announced_projects
            if (get_construction_year(project) <= year) && (get_end_life_year(project) >= year)
                installed_capacity[year] += get_maxcap(project)
            end
        end
    end
    return installed_capacity
end

"""
This function does nothing if project is not in Existing phase.
"""
function extrapolate_profits!(project::Project{<:BuildPhase}, final_year::Int64)
    return
end

"""
This function extrapolates the profits of Existing projects beyond the simulation horizon.
"""
function extrapolate_profits!(project::Project{Existing}, final_year::Int64)
    end_life_year = get_end_life_year(project)
    finance_data = get_finance_data(project)
    if end_life_year > final_year
        for y in (final_year + 1):end_life_year
            set_annual_cashflow!(finance_data, y, get_annual_cashflow(finance_data)[final_year])
            for product in get_products(project)
                name = get_name(product)
                set_realized_profit!(finance_data, name, y, get_realized_profit(finance_data)[name, final_year])
            end
        end
    end
    return
end

function month_lookup(str::Union{String, SubString{String}})
    if str == "Jan" || str == "January"
        month = 1
    elseif str == "Feb" || str == "February"
        month = 2
    elseif str == "Mar" || str == "March"
        month = 3
    elseif str == "Apr" || str == "April"
        month = 4
    elseif str == "May"
        month = 5
    elseif str == "Jun" || str == "June"
        month = 6
    elseif str == "Jul" || str == "July"
        month = 7
    elseif str == "Aug" || str == "August"
        month = 8
    elseif str == "Sep" || str == "September"
        month = 9
    elseif str == "Oct" || str == "October"
        month = 10
    elseif str == "Nov" || str == "November"
        month = 11
    elseif str == "Dec" || str == "December"
        month = 12
    end
    return  Int(month)
end

function find_rt_periods(hours::Vector{Int64}, num_rt_intervals::Int64)
    rt_periods = [];
    for hour in hours
        append!(rt_periods, collect(((hour - 1) * num_rt_intervals + 1):(hour * num_rt_intervals)))
    end
    return rt_periods
end
