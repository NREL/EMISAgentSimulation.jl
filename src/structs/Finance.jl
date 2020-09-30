"""
This struct contains the financial data of projects.
    investment_cost: Yearly capex.
    effective_investment_cost: Used to fix yearly capex.
    lag_time: Time (in years) for the project to get built.
    life_time: Total life_time of the project.
    capex_years: Number of years for capital cost recovery calculation.
    fixed_OM_cost: Annual investment cost.
    queue_cost: Vector of yearly queue costs.
    scenario_profit: Array of profits for each product through the project's life_time for each scenario.
    realized_profit: AxisArray of realized profits.
    discount_rate: Discount rate (%) for investment in the project.
    scenario_npv: Net Present Value of the project for each scenario.
    expected_npv: Expected value of NPV.
    scenario_utility: Utility of the project for each scenario.
    expected_utility: Expected value of the utility.
    annual_cashflow: Annual realized cashflow.
    ownedby: Name of investor.
"""
mutable struct Finance
    investment_cost::Vector{Float64}
    effective_investment_cost::Float64
    lag_time::Int64
    life_time::Int64
    capex_years::Int64
    fixed_OM_cost::Float64
    queue_cost::Vector{Float64}
    scenario_profit::Dict{String, Vector{AxisArrays.AxisArray{Float64, 2}}}
    realized_profit::AxisArrays.AxisArray{Float64, 2}
    discount_rate::Float64
    scenario_npv::Dict{String, Vector{Float64}}
    expected_npv::Vector{Float64}
    scenario_utility::Dict{String, Vector{Float64}}
    expected_utility::Vector{Float64}
    annual_cashflow::Vector{Float64}
    ownedby::String
end

get_investment_cost(finance_data::Finance) = finance_data.investment_cost
get_effective_investment_cost(finance_data::Finance) = finance_data.effective_investment_cost
get_lag_time(finance_data::Finance) = finance_data.lag_time
get_life_time(finance_data::Finance) = finance_data.life_time
get_capex_years(finance_data::Finance) = finance_data.capex_years
get_fixed_OM_cost(finance_data::Finance) = finance_data.fixed_OM_cost
get_queue_cost(finance_data::Finance) = finance_data.queue_cost
get_scenario_profit(finance_data::Finance) = finance_data.scenario_profit
get_realized_profit(finance_data::Finance) = finance_data.realized_profit
get_discount_rate(finance_data::Finance) = finance_data.discount_rate
get_scenario_npv(finance_data::Finance) = finance_data.scenario_npv
get_expected_npv(finance_data::Finance) = finance_data.expected_npv
get_scenario_utility(finance_data::Finance) = finance_data.scenario_utility
get_expected_utility(finance_data::Finance) = finance_data.expected_utility
get_annual_cashflow(finance_data::Finance) = finance_data.annual_cashflow
get_ownedby(finance_data::Finance) = finance_data.ownedby

function set_scenario_npv!(finance_data::Finance, scenario_name::String, iteration_year::Int64, npv::Float64)
    finance_data.scenario_npv[scenario_name][iteration_year] = npv
    return
end

function set_scenario_npv!(finance_data::Finance, scenario_name::String, npv::Vector{Float64})
    finance_data.scenario_npv[scenario_name] = npv
    return
end

function set_expected_npv!(finance_data::Finance, year::Int64, npv::Float64)
    finance_data.expected_npv[year] = npv
    return
end

function set_expected_npv!(finance_data::Finance, npv_vector::Vector{Float64})
    finance_data.expected_npv = npv_vector
end

function set_scenario_utility!(finance_data::Finance, scenario_name::String, iteration_year::Int64, utility::Float64)
    finance_data.scenario_utility[scenario_name][iteration_year] = utility
    return
end

function set_scenario_utility!(finance_data::Finance, scenario_name::String, utility::Vector{Float64})
    finance_data.scenario_utility[scenario_name] = utility
    return
end

function set_expected_utility!(finance_data::Finance, year::Int64, utility::Float64)
    finance_data.expected_utility[year] = utility
    return
end

function set_expected_utility!(finance_data::Finance, utility_vector::Vector{Float64})
    finance_data.expected_utility = utility_vector
end

function set_scenario_profit!(finance_data::Finance, scenario_name::String, profit::Vector{AxisArrays.Array{Float64, 2}})
    finance_data.scenario_profit[scenario_name] = profit
    return
end

function set_scenario_profit!(finance_data::Finance, scenario_name::String, iteration_year::Int64,  profit::AxisArrays.AxisArray{Float64, 2})
    finance_data.scenario_profit[scenario_name][iteration_year] = profit
    return
end

function set_scenario_profit!(finance_data::Finance,
                    product_name::Symbol,
                    iteration_year::Int64,
                    year_index::Int64,
                    profit::Nothing)
    return
end

function set_scenario_profit!(finance_data::Finance,
                    scenario_name::String,
                    product_name::Symbol,
                    iteration_year::Int64,
                    year_index::Int64,
                    profit::Float64)

    finance_data.scenario_profit[scenario_name][iteration_year][product_name, year_index] = profit
    return
end

function set_scenario_profit!(finance_data::Finance,
                    scenario_name::String,
                    product_name::Symbol,
                    iteration_year::Int64,
                    start_year::Int64,
                    end_year::Int64,
                    profit::Vector{Float64})

    finance_data.scenario_profit[scenario_name][iteration_year][product_name, start_year:end_year] = profit
    return
end

function set_realized_profit!(finance_data::Finance, profit::AxisArrays.AxisArray{Float64, 2})
    finance_data.realized_profit = profit
    return
end

function set_realized_profit!(finance_data::Finance,
                    product_name::Symbol,
                    year_index::Int64,
                    profit::Nothing)
    return
end

function set_realized_profit!(finance_data::Finance,
                    product_name::Symbol,
                    year_index::Int64,
                    profit::Float64)
    finance_data.realized_profit[product_name, year_index] = profit
    return
end

function set_realized_profit!(finance_data::Finance,
                    product_name::Symbol,
                    start_year::Int64,
                    end_year::Int64,
                    profit::Vector{Float64})
    finance_data.realized_profit[product_name, start_year:end_year] = profit
    return
end

function set_annual_cashflow!(finance_data::Finance, year::Int64, cashflow::Float64)
    finance_data.annual_cashflow[year] = cashflow
end
