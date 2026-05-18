"""
This function runs the actual energy and ancillary services market clearing module using endogeneous ED.
Returns the market clearing prices, capacity factors and reserve percentages.
"""
function energy_mkt_clearing(sys_UC::Nothing,
                             sys_ED::Nothing,
                             sys_local_ED::MarketClearingProblem,
                             reserve_penalty::String,
                             rec_requirement::Float64,
                             simulation_dir::String,
                             zones::Vector{String},
                             num_days::Int64,
                             pcm_scenario::String,
                             iteration_year::Int64,
                             da_resolution::Int64,
                             rt_resolution::Int64,
                             case_name::String,
                             solver::JuMP.MOI.OptimizerWithAttributes)
    @warn "This functionality has been deprecated. All PCMs are currently run with Sienna."
   
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
function energy_mkt_clearing(sys_MD::PSY.System,
                             sys_UC::PSY.System,
                             sys_ED::PSY.System,
                             sys_local_ED::Union{Nothing, MarketClearingProblem},
                             simulation_dir::String,
                             reserve_penalty::String,
                             rec_requirement::Float64,
                             zones::Vector{String},
                             num_days::Int64,
                             pcm_scenario::String,
                             iteration_year::Int64,
                             da_resolution::Int64,
                             rt_resolution::Int64,
                             case_name::String,
                             solver::JuMP.MOI.OptimizerWithAttributes,
                             base_dir::String,
                             simulation::AgentSimulation,
                             current_siip_sim,
                             md_market_bool::Bool,
                             ordc_curved::Bool,
                             siip_system)

    @warn "Uncomment the update_PSY_timeseries function calls"
    update_PSY_timeseries!(simulation, sys_MD, rec_requirement, simulation_dir, "MD", pcm_scenario, iteration_year, da_resolution, rt_resolution)
    update_PSY_timeseries!(simulation, sys_UC, rec_requirement, simulation_dir, "UC", pcm_scenario, iteration_year, da_resolution, rt_resolution)
    update_PSY_timeseries!(simulation, sys_ED, rec_requirement, simulation_dir, "ED", pcm_scenario, iteration_year, da_resolution, rt_resolution) 
    # TODO: need to update outage timeseries for MD as well
    # update_PSY_outage_timeseries!(sys_UC, sys_ED,get_results_dir(simulation),base_dir,iteration_year)

    energy_price_ed,
    energy_price_uc,
    energy_price_md,
    reserve_price_ed,
    reserve_price_uc,
    reserve_price_md,
    inertia_price,
    capacity_factors_md,
    capacity_factors_uc,
    capacity_factors_ed,
    reserve_perc_md,
    reserve_perc_uc,
    reserve_perc_ed,
    inertia_perc,
    start_up_costs,
    shut_down_costs,
    energy_voll,
    energy_voll_uc,
    energy_voll_md,
    reserve_voll,
    reserve_voll_uc,
    reserve_voll_md,
    inertia_voll = create_simulation(sys_MD, sys_UC, sys_ED, simulation_dir, reserve_penalty, zones, num_days, da_resolution, rt_resolution, case_name, solver, current_siip_sim, md_market_bool, ordc_curved, siip_system)

    return energy_price_ed, energy_price_uc, energy_price_md, reserve_price_ed, reserve_price_uc, reserve_price_md, inertia_price, capacity_factors_md, capacity_factors_uc, capacity_factors_ed, reserve_perc_md, reserve_perc_uc, reserve_perc_ed, inertia_perc, start_up_costs, shut_down_costs, energy_voll, energy_voll_uc, energy_voll_md, reserve_voll, reserve_voll_uc, reserve_voll_md, inertia_voll;
end
