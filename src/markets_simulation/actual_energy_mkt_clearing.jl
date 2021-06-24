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
                             iteration_year::Int64,
                             da_resolution::Int64,
                             rt_resolution::Int64,
                             case_name::String,
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
                             simulation_dir::String,
                             load_growth::AxisArrays.AxisArray{Float64, 1},
                             zones::Vector{String},
                             num_days::Int64,
                             iteration_year::Int64,
                             da_resolution::Int64,
                             rt_resolution::Int64,
                             case_name::String,
                             solver::JuMP.MOI.OptimizerWithAttributes)

    update_PSY_timeseries!(sys_UC, load_growth, simulation_dir, "UC", iteration_year, da_resolution, rt_resolution)
    update_PSY_timeseries!(sys_ED, load_growth, simulation_dir, "ED", iteration_year, da_resolution, rt_resolution)

    energy_price,
    reserve_price,
    inertia_price,
    capacity_factors,
    reserve_perc,
    inertia_perc,
    start_up_costs,
    shut_down_costs,
    energy_voll,
    reserve_voll,
    inertia_voll = create_simulation(sys_UC, sys_ED, simulation_dir, zones, num_days, da_resolution, rt_resolution, case_name, solver)

    return energy_price, reserve_price, inertia_price, capacity_factors, reserve_perc, inertia_perc, start_up_costs, shut_down_costs, energy_voll, reserve_voll, inertia_voll;
end
