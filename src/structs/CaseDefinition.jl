"""
    This struct contains all the required data, parameters and solvers
    for creating and running AgentSimulation.
        name: Case name
        base_dir: Directory where Simulation input data for markets and investors are stored.
        sys_dir: Test system directory.
        cem_solver: Solvers used for optimization problems. (The solver should be able to solve QP for price prediction and MILP for SIIP production cost model)
        siip_solver: Solvers used for optimization problems. (The solver should be able to solve QP for price prediction and MILP for SIIP production cost model)
        siip_market_clearing: Whether SIIP production cost model is to be used for energy market clearing. If false, the endogenous Economic Dispatch model will be used for market clearing.
        start_year: Start year for the simulation (default is set to 2020)
        total_horizon: Number of years of data available for the simulation.
        rolling_horizon: Number of years to be used for price prediction. If end of rolling horizon exceeds the years of available data, a receding horizon approach is used.
        simulation_years: Number of years to be simulated.
        rep_period_interval: Total number of hours in representative period selection e.g., 24 for representative days, 168 for representative weeks, etc. Default =  24
        num_rep_periods: Number of representative periods to be used for price prediction. Default = 10.
        avg_block_size: Number of hours in CEM aggregated time blocks. Default = 4. Set to 1 if hourly granularity is needed.
        fixed_block_size: Whether the block size in CEM time block aggregation is fixed, i.e., all blocks are of the same length = avg_block_size. If set to FALSE, the model will select variable size time blocks based on chronological clustering.
        rep_chronology_checkpoint: Number of hours in the checkpoint interval for chronological storage constraints in CEM. Default = 0, i.e., no chronological checkpoints
        da_resolution: Resolution of Day Ahead market clearing (minutes)
        rt_resolution: Resolution of Real Time market clearing (minutes)
        rps_target: High, Mid or Low RPS Target
        markets: Dictionary of which markets are simulated
        ordc_curved: Whether to include the curved part of the ORDC
        ordc_unavailability_method: Which method (Sequential Monte Carlo or Convolution) to use for generating unavailability distribution for ORDCs
        derating_scale: Factor for scaling derating factors
        mopr: Whether Minimum Offer Price Rule is applied
        battery_cap_mkt: Whether Batteries can paritcipate in capacity markets
        vre_reserves: Whether VRE can provide reserves
        heterogeneity: Whether investors' heterogeneous financial characteristics and technology preferences are modeled.
        reserve_penalty: High, Mid or Low penalty prices for reserves
        static_capacity_market: Whether the capacity market demand curve is static or RA-informed
        irm_scalar: Scalar for installed reserve margin to be used for creating the capacity market demand curve.
        accreditation_methodology: RA metric used for ELCC and EFC accreditation methodology. Options: LOLE or EUE. Set to "Nothing" if accreditation methodology is TopNetLoad.
        accreditation_metric: Scalar used for modifying the derating factors of VRE and batteries. Range > 0.0
        forecast_type: "Perfect" or "imperfect" forecasts used for price prediction.
        max_carbon_tax_increase: Maximum annual increase in carbon prices due to under-achievement of Clean Energy Targets.
        info_symmetry: Whether investors have symmetric information about forecast parameters.
        belief_update: Whether investors' beliefs are updated each year after actual market clearing.
        uncertainty: Whether multiple probability weighted scenarios are used instead of a deterministic forecast.
        risk_aversion: Whether investors are risk averse.
        parallel_investors: Whether investors' price prediction is to be parallelized.
        parallel_scenarios: Whether each investor's price prediction scenarios are to be parallelized.
"""

struct CaseDefinition
    name::String
    base_dir::String
    sys_dir::String
    solver::JuMP.MOI.OptimizerWithAttributes
    siip_market_clearing::Bool
    start_year::Int64
    total_horizon::Int64
    rolling_horizon::Int64
    simulation_years::Int64
    rep_period_interval::Int64
    num_rep_periods::Int64
    avg_block_size::Int64
    fixed_block_size::Bool
    rep_chronology_checkpoint::Int64
    da_resolution::Int64
    rt_resolution::Int64
    rps_target::String
    markets::Dict{Symbol, Bool}
    ordc_curved::Bool
    ordc_unavailability_method::String
    reserve_penalty::String
    static_capacity_market::Bool
    irm_scalar::Float64
    accreditation_methodology::String
    accreditation_metric::String
    derating_scale::Float64
    mopr::Bool
    battery_cap_mkt::Bool
    vre_reserves::Bool
    heterogeneity::Bool
    forecast_type::String
    max_carbon_tax_increase::Float64
    info_symmetry::Bool
    belief_update::Bool
    uncertainty::Bool
    risk_aversion::Bool
    parallel_investors::Bool
    parallel_scenarios::Bool

    function CaseDefinition(name,
                            base_dir,
                            sys_dir,
                            solver,
                            siip_market_clearing,
                            start_year,
                            total_horizon,
                            rolling_horizon,
                            simulation_years,
                            rep_period_interval,
                            num_rep_periods,
                            avg_block_size,
                            fixed_block_size,
                            rep_chronology_checkpoint,
                            da_resolution,
                            rt_resolution,
                            rps_target,
                            markets,
                            ordc_curved,
                            ordc_unavailability_method,
                            reserve_penalty,
                            static_capacity_market,
                            irm_scalar,
                            accreditation_methodology,
                            accreditation_metric,
                            derating_scale,
                            mopr,
                            battery_cap_mkt,
                            vre_reserves,
                            heterogeneity,
                            forecast_type,
                            max_carbon_tax_increase,
                            info_symmetry,
                            belief_update,
                            uncertainty,
                            risk_aversion,
                            parallel_investors,
                            parallel_scenarios)

        @assert total_horizon >= simulation_years

        forecast_type = lowercase(forecast_type)
        @assert forecast_type == "perfect" || lowercase(forecast_type) == "imperfect"
        @assert lowercase(rps_target) == "high" || lowercase(rps_target) == "mid" || lowercase(rps_target) == "low"
        @assert lowercase(reserve_penalty) == "high" || lowercase(reserve_penalty) == "mid" || lowercase(reserve_penalty) == "low"

        if forecast_type == "perfect"
            @assert info_symmetry == true
            @assert belief_update == false
            @assert uncertainty == false
            @assert risk_aversion == false
        end

        @assert da_resolution >= rt_resolution
        @assert irm_scalar >= 0.0
        #=
        if !(siip_market_clearing)
            @assert da_resolution == rt_resolution
        end
        =#
        return new(name,
                   base_dir,
                   sys_dir,
                   solver,
                   siip_market_clearing,
                   start_year,
                   total_horizon,
                   rolling_horizon,
                   simulation_years,
                   rep_period_interval,
                   num_rep_periods,
                   avg_block_size,
                   fixed_block_size,
                   rep_chronology_checkpoint,
                   da_resolution,
                   rt_resolution,
                   rps_target,
                   markets,
                   ordc_curved,
                   ordc_unavailability_method,
                   reserve_penalty,
                   static_capacity_market,
                   irm_scalar,
                   accreditation_methodology,
                   accreditation_metric,
                   derating_scale,
                   mopr,
                   battery_cap_mkt,
                   vre_reserves,
                   heterogeneity,
                   forecast_type,
                   max_carbon_tax_increase,
                   info_symmetry,
                   belief_update,
                   uncertainty,
                   risk_aversion,
                   parallel_investors,
                   parallel_scenarios)
    end
end

function CaseDefinition(name::String,
                        base_dir::String,
                        sys_dir::String,
                        solver::JuMP.MOI.OptimizerWithAttributes;
                        siip_market_clearing::Bool = true,
                        start_year::Int64 = 2020,
                        total_horizon::Int64 = 15,
                        rolling_horizon::Int64 = 10,
                        simulation_years::Int64 = 15,
                        rep_period_interval::Int64 = 24,
                        num_rep_periods::Int64 = 10,
                        avg_block_size::Int64 = 4,
                        fixed_block_size::Bool = FALSE,
                        rep_chronology_checkpoint::Int64 = 0,
                        da_resolution::Int64 = 60,
                        rt_resolution::Int64 = 5,
                        rps_target::String = "Mid",
                        markets::Dict{Symbol, Bool} = Dict(:Energy => true, :Synchronous => true, :Primary => true, :Reg_Up => true, :Reg_Down => true,	:Flex_Up => true, :Flex_Down => true, :Capacity => true, :REC => true, :CarbonTax => true),
                        ordc_curved::Bool = true,
                        ordc_unavailability_method::String = "Convolution",
                        reserve_penalty::String = "Mid",
                        static_capacity_market::Bool = true,
                        irm_scalar::Float64 = 1.0,
                        accreditation_methodology::String = "TopNetLoad",
                        accreditation_metric::String = "None",
                        derating_scale::Float64 = 1.0,
                        mopr::Bool = false,
                        battery_cap_mkt::Bool = true,
                        vre_reserves::Bool = true,
                        heterogeneity::Bool = false,
                        forecast_type::String = "perfect",
                        max_carbon_tax_increase::Float64 = 0.0,
                        info_symmetry::Bool = true,
                        belief_update::Bool = false,
                        uncertainty::Bool = false,
                        risk_aversion::Bool = false,
                        parallel_investors::Bool = false,
                        parallel_scenarios::Bool = false)

    CaseDefinition(name,
                   base_dir,
                   sys_dir,
                   solver,
                   siip_market_clearing,
                   start_year,
                   total_horizon,
                   rolling_horizon,
                   simulation_years,
                   rep_period_interval,
                   num_rep_periods,
                   avg_block_size,
                   fixed_block_size,
                   rep_chronology_checkpoint,
                   da_resolution,
                   rt_resolution,
                   rps_target,
                   markets,
                   ordc_curved,
                   ordc_unavailability_method,
                   reserve_penalty,
                   static_capacity_market,
                   irm_scalar,
                   accreditation_methodology,
                   accreditation_metric,
                   derating_scale,
                   mopr,
                   battery_cap_mkt,
                   vre_reserves,
                   heterogeneity,
                   forecast_type,
                   max_carbon_tax_increase,
                   info_symmetry,
                   belief_update,
                   uncertainty,
                   risk_aversion,
                   parallel_investors,
                   parallel_scenarios)
end

get_base_dir(case::CaseDefinition) = case.base_dir
get_sys_dir(case::CaseDefinition) = case.sys_dir
get_solver(case::CaseDefinition) = case.solver
get_siip_market_clearing(case::CaseDefinition) = case.siip_market_clearing
get_start_year(case::CaseDefinition) = case.start_year
get_total_horizon(case::CaseDefinition) = case.total_horizon
get_rolling_horizon(case::CaseDefinition) = case.rolling_horizon
get_simulation_years(case::CaseDefinition) = case.simulation_years
get_rep_period_interval(case::CaseDefinition) = case.rep_period_interval
get_num_rep_periods(case::CaseDefinition) = case.num_rep_periods
get_avg_block_size(case::CaseDefinition) = case.avg_block_size
get_fixed_block_size(case::CaseDefinition) = case.fixed_block_size
get_rep_chronology_checkpoint(case::CaseDefinition) = case.rep_chronology_checkpoint
get_da_resolution(case::CaseDefinition) = case.da_resolution
get_rt_resolution(case::CaseDefinition) = case.rt_resolution
get_rps_target(case::CaseDefinition) = case.rps_target
get_markets(case::CaseDefinition) = case.markets
get_ordc_curved(case::CaseDefinition) = case.ordc_curved
get_reserve_penalty(case::CaseDefinition) = case.reserve_penalty
get_static_capacity_market(case::CaseDefinition) = case.static_capacity_market
get_irm_scalar(case::CaseDefinition) = case.irm_scalar
get_ordc_unavailability_method(case::CaseDefinition) = case.ordc_unavailability_method
get_accreditation_methodology(case::CaseDefinition) = case.accreditation_methodology
get_accreditation_metric(case::CaseDefinition) = case.accreditation_metric
get_derating_scale(case::CaseDefinition) = case.derating_scale
get_mopr(case::CaseDefinition) = case.mopr
get_battery_cap_mkt(case::CaseDefinition) = case.battery_cap_mkt
get_vre_reserves(case::CaseDefinition) = case.vre_reserves
get_heterogeneity(case::CaseDefinition) = case.heterogeneity
get_info_symmetry(case::CaseDefinition) = case.info_symmetry
get_belief_update(case::CaseDefinition) = case.belief_update
get_forecast_type(case::CaseDefinition) = case.forecast_type
get_max_carbon_tax_increase(case::CaseDefinition) = case.max_carbon_tax_increase
get_uncertainty(case::CaseDefinition) = case.uncertainty
get_risk_aversion(case::CaseDefinition) = case.risk_aversion
get_parallel_investors(case::CaseDefinition) = case.parallel_investors
get_parallel_scenarios(case::CaseDefinition) = case.parallel_scenarios

function get_name(case::CaseDefinition)
    #=

    if get_info_symmetry(case)
        information = "InfoSym"
    else
        information = "InfoASym"
    end

    if get_belief_update(case)
        update = "UpdateBelief"
    else
        update = "NoUpdate"
    end

    if get_uncertainty(case)
        uncertainty = "Uncertain"
    else
        uncertainty = "Deterministic"
    end

    if get_risk_aversion(case)
        risk = "RiskAverse"
    else
        risk = "RiskNeutral"
    end

    case_name = "$(investors)_$(information)_Forecast-$(get_forecast_type(case))_$(uncertainty)_$(update)_$(risk)_$(get_simulation_years(case))years"


    if get_heterogeneity(case)
        investors = "Het"
    else
        investors = "Hom"
    end

    #New case name

    rps = "$(get_rps_target(case))_RPS"

    if get_markets(case)[:Capacity]
        capacity = "CapMkt"
    else
        capacity = "NoCapMkt"
    end

    if get_ordc_curved(case)
        ordc = "with_ORDC"
    else
        ordc = "without_ORDC"
    end

    penalty = "$(get_reserve_penalty(case))_penalty"

    if get_markets(case)[:CarbonTax]
        carbon = "Carbon_Tax"
    else
        carbon = "No_Carbon_Tax"
    end


    if get_mopr(case)
        mopr = "MOPR_ON"
    else
        mopr = "MOPR_OFF"
    end

    if get_battery_cap_mkt(case)
        bat_cap = "BAT_Cap_ON"
    else
        bat_cap = "BAT_Cap_OFF"
    end

    derating_scale = replace("$(get_derating_scale(case))", "." => "_")

    derating = "derating_scale_$(derating_scale)"

    if get_vre_reserves(case)
        vre_reserves = "VRE_reserves"
    else
        vre_reserves = "No_VRE_and Bat_reserves"
    end

    if get_markets(case)[:Inertia]
        inertia = "Inertia"
    else
        inertia = "No_Inertia"
    end

    case_name = "$(investors)_$(rps)_$(capacity)_$(ordc)_$(penalty)_$(carbon)_$(derating)_$(mopr)_$(bat_cap)_$(vre_reserves)_$(inertia)"

    return case_name

    =#
    return "$(case.name)_$(get_rps_target(case))_RPS"
end

function get_data_dir(case::CaseDefinition)
    base_dir = get_base_dir(case)
    case_dir = joinpath(base_dir, get_name(case))

    return case_dir
end



