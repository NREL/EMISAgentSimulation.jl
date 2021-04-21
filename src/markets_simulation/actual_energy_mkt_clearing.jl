"""
This function runs the actual energy and ancillary services market clearing module using endogeneous ED.
Returns the market clearing prices, capacity factors and reserve percentages.
"""
function energy_mkt_clearing(sys_UC::Nothing,
                             sys_ED::Nothing,
                             sys_local_ED::MarketClearingProblem,
                             load_growth::AxisArrays.AxisArray{Float64, 1},
                             simulation_dir::String,
                             zones::Vector{String},
                             num_days::Int64,
                             solver::JuMP.MOI.OptimizerWithAttributes,
                             iteration_year::Int64)
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
                             simulation_dir::String,
                             load_growth::AxisArrays.AxisArray{Float64, 1},
                             zones::Vector{String},
                             num_days::Int64,
                             solver::JuMP.MOI.OptimizerWithAttributes,
                             iteration_year::Int64)

    #update_PSY_timeseries!(sys_UC, load_growth, simulation_dir)
    #update_PSY_timeseries!(sys_ED, load_growth, simulation_dir)

    energy_price,
    reserve_price,
    capacity_factors,
    reserve_perc = create_simulation(sys_UC, sys_ED, simulation_dir, zones, num_days, solver, iteration_year)

    return energy_price, reserve_price, capacity_factors, reserve_perc;
end
