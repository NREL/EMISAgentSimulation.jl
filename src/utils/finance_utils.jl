"""
This function checks if value is in input dictionary.
"""
function check_input_for_variable(dict, name, default = 0)
  var = default
  if haskey(dict, name)
    var = dict[name]
  end
  return var
end

"""
This function calculates the adjusted capital cost of a project.
"""
function calculate_capital_costs(inputs)
  #capital costs are the cost of getting the plant to completion, minus incentive benefits, minus tax benefits
  plant_cost = calculate_plant_cost(inputs["cost_inputs"])
  incentives = calculate_incentive_benefits(inputs["incentive_inputs"], plant_cost)
  tax_benefits = calculate_tax_benefits(inputs["tax_inputs"], plant_cost, inputs["incentive_inputs"])
  return(plant_cost - incentives - tax_benefits)
end

"""
This function adds the total construction and development costs of a project.
"""
function calculate_plant_cost(cost_inputs)
  return sum(cost_inputs)
end

"""
This function returns the MACRS depreciation schedule of a project.
"""
function get_macrs_schedule(num_depreciation_years, macrs_df)
  return macrs_df[:, Symbol("Y"*string(num_depreciation_years))]
end

"""
This function calculates the incentive benefits which can come from grants and investment tax credits
   grants are in dollars/KW of capacity
   itc is in terms of percent of plant cost
   itc_recovery_factor scales the percent of capital costs which can be covered by the itc. This can be due to transaction frictions.
"""
function calculate_incentive_benefits(incentive_inputs, plant_cost)
  grants = check_input_for_variable(incentive_inputs, "grants", 0)
  itc = check_input_for_variable(incentive_inputs, "itc", 0)
  itc_recovery_factor = check_input_for_variable(incentive_inputs, "itc_recovery_factor", 1)
  return grants + itc*itc_recovery_factor*plant_cost
end
"""

This function calculates the tax benefits which are from depreciating the capital asset. tax inputs include:
    macrs:                            the depreciation schedule in percent for the asset
    depreciable_percent:              fraction of plant costs which can be depreciated
    discount_rate:                    used to transform future tax benefits to present value
    tax_rate:                         tax percent which is deducted.
"""
function calculate_tax_benefits(tax_inputs, plant_cost, incentive_inputs)
  #
  macrs = check_input_for_variable(tax_inputs, "macrs", [0])
  itc = check_input_for_variable(incentive_inputs, "itc", 0)
  depreciable_percent = check_input_for_variable(tax_inputs, "depreciable_percent", 1)
  tax_rate = check_input_for_variable(tax_inputs, "tax_rate", .25)
  discount_rate = check_input_for_variable(tax_inputs, "discount_rate", 0.08)
  #
  return calculate_present_value(plant_cost*depreciable_percent*(1-itc/2)*tax_rate*macrs, discount_rate)
end

"""
This function calculates the present value which takes yearly values and a discount rate and returns the present value.
"""
function calculate_present_value(values, discount_rate)
  return sum([values[i]*(1+discount_rate)^(-i) for i in (1:length(values))])
end

"""
This function calculates calculates the weighted average cost of capital
 Takes a list with: Tax rate, debt interest rate, debt fraction, and either the equity discount rate or the required values to calculate a CAPM model
"""
function calculate_wacc(wacc_inputs)
  tax_rate =  check_input_for_variable(wacc_inputs, "tax_rate", 0.25)
  debt_fraction =  check_input_for_variable(wacc_inputs, "debt_fraction", 0.6)
  debt_interest_rate =  check_input_for_variable(wacc_inputs, "debt_interest_rate", 0.05)
  equity_discount_rate =  check_input_for_variable(wacc_inputs, "equity_discount_rate", 0.12)
  #
  debt_contribution = (1- tax_rate)*debt_interest_rate*debt_fraction
  equity_contribution = equity_discount_rate*(1-debt_fraction)
  return (debt_contribution + equity_contribution)
end
#_________________________________________________________________________________________________________________________________________________
"""
This function calculates a CAPM model equity rate = risk-free rate + beta*(market return - risk_free rate) + project-specific adder
"""
function calculate_equity_discount_rate(inputs)
  project_specific_adder = check_input_for_variable(inputs, "project_specific_adder", 0)
  risk_free_rate = check_input_for_variable(inputs, "risk_free_rate", 0.04)
  beta = check_input_for_variable(inputs, "beta", 1)
  market_return = check_input_for_variable(inputs, "market_return", .12)
  return (risk_free_rate + beta*(market_return - risk_free_rate) + project_specific_adder)
end

"""
This function utilizes all the finance functions to return the WACC and adjusted capital cost of a project.
"""
function calculate_capcost_and_wacc(project_capex::Vector{Float64},
                                    category::String,
                                    investor_dir::String,
                                    online_year::Int64)
  finance_data =  read_data(joinpath(investor_dir, "finance_params.csv"))
  investor_characteristics = read_data(joinpath(investor_dir, "characteristics.csv"))

  project_row = findfirst(x-> x == category, finance_data[:, "Category"])
  macrs_df = read_data(joinpath(investor_dir, "MACRS Schedule.csv"))

  equity_discount_rate = finance_data[project_row, "EQUITY_RATE"]

  if online_year <= 2020
      itc = finance_data[project_row, "ITC 2020"]
  elseif online_year == 2021
      itc = finance_data[project_row, "ITC 2021"]
  else
      itc = finance_data[project_row, "ITC post 2021"]
  end

  incentive_inputs = Dict("itc" => itc,
                          "grants" => finance_data[project_row, "GRANTS"],
                          "itc_recovery_factor" => finance_data[project_row, "ITC_RECOVERY_FACTOR"])

  macrs_schedule = get_macrs_schedule(finance_data[project_row, "MACRS_DEPRECIATION_YEARS"], macrs_df)
  wacc_inputs = Dict("tax_rate" => finance_data[project_row, "TAX_RATE"],
                      "debt_fration" => finance_data[project_row, "DEBT_FRACTION"],
                      "debt_interest_rate" => finance_data[project_row, "DEBT_INTEREST_RATE"],
                      "equity_discount_rate" => finance_data[project_row, "EQUITY_RATE"])

  discount_rate = calculate_wacc(wacc_inputs) + investor_characteristics[1, "discount_rate_adder"]

  tax_inputs = Dict("macrs" => macrs_schedule,
                    "tax_rate" => finance_data[project_row, "TAX_RATE"],
                    "discount_rate" => discount_rate,
                    "depreciable_percent"=> finance_data[project_row, "DEPRECIABLE_PERCENT"])

  cap_cost = zeros(length(project_capex))
  for i in 1:length(project_capex)
      COSTS = [project_capex[i]]
      capital_cost_inputs = Dict("cost_inputs" => [project_capex[i]], "incentive_inputs" => incentive_inputs, "tax_inputs" => tax_inputs)
      cap_cost[i] = calculate_capital_costs(capital_cost_inputs)
  end

  return cap_cost, discount_rate
end
