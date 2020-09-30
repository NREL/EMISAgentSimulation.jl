"""
This function runs the actual energy and ancillary services market clearing module using endogeneous ED.
Returns the market clearing prices, capacity factors and reserve percentages.
"""
function energy_mkt_clearing(sys_UC::Nothing,
                             sys_ED::Nothing,
                             sys_local_ED::MarketClearingProblem,
                             load_growth::AxisArrays.AxisArray{Float64, 1},
                             zones::Vector{String},
                             num_days::Int64,
                             solver::JuMP.MOI.OptimizerWithAttributes)
    energy_price,
    reserve_up_price,
    reserve_down_price,
    capacity_factors,
    reserve_up_perc,
    reserve_down_perc = economicdispatch(sys_local_ED,
                                       solver,
                                       "C:/Users/manwar2/Documents/GitRepos/emt-tests/data/simulation_data/results/ed_results.txt")

    return energy_price, reserve_up_price, reserve_down_price, capacity_factors, reserve_up_perc, reserve_down_perc;
end

"""
This function runs the actual energy and ancillary services market clearing module using SIIP PSI Simualtion.
Returns the market clearing prices, capacity factors and reserve percentages.
"""
function energy_mkt_clearing(sys_UC::PSY.System,
                             sys_ED::PSY.System,
                             sys_local_ED::Union{Nothing, MarketClearingProblem},
                             load_growth::AxisArrays.AxisArray{Float64, 1},
                             zones::Vector{String},
                             num_days::Int64,
                             solver::JuMP.MOI.OptimizerWithAttributes)

    update_PSY_timeseries!(sys_UC, load_growth, num_days)

    energy_price,
    reserve_up_price,
    reserve_down_price,
    capacity_factors,
    reserve_up_perc,
    reserve_down_perc = create_simulation(sys_UC, sys_ED, zones, num_days, solver)

    return energy_price, reserve_up_price, reserve_down_price, capacity_factors, reserve_up_perc, reserve_down_perc;
end
