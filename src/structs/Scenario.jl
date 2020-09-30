"""
This struct contains the data for each scenario.
    name: Scenario name
    probability: Probability of the scenario
    parameter_multipliers: Parameter multipliers for each scenario - to be multiplied to updated parameters from kalman filters
    parameter_values: Dictionary containing parameter values for the scenario
"""
mutable struct Scenario
    name::String
    probability::Float64
    parameter_multipliers::Union{Nothing, Dict{String, Float64}}
    parameter_values::Vector{AxisArrays.AxisArray{Float64, 2}}
end

get_name(scenario::Scenario) = scenario.name
get_probability(scenario::Scenario) = scenario.probability
get_parameter_multipliers(scenario::Scenario) = scenario.parameter_multipliers
get_parameter_values(scenario::Scenario) = scenario.parameter_values

function set_parameter_values!(scenario::Scenario, iteration_year::Int64, parameter_values::AxisArrays.AxisArray{Float64, 2})
    scenario.parameter_values[iteration_year] = parameter_values
end
