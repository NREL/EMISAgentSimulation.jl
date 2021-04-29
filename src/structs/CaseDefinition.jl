"""
    This struct contains all the required data, parameters and solvers
    for creating and running AgentSimulation.
        base_dir: Directory where Simulation input data for markets and investors are stored.
        sys_dir: Test system directory.
        solver: Solvers used for optimization problems. (The solver should be able to solve QP for price prediction and MILP for SIIP production cost model)
        siip_market_clearing: Whether SIIP production cost model is to be used for energy market clearing. If false, the endogenous Economic Dispatch model will be used for market clearing.
        start_year: Start year for the simulation (default is set to 2020)
        total_horizon: Number of years of data available for the simulation.
        rolling_horizon: Number of years to be used for price prediction. If end of rolling horizon exceeds the years of available data, a receding horizon approach is used.
        simulation_years: Number of years to be simulated.
        num_rep_days: Number of representative days used for price prediction.
        da_resolution: Resolution of Day Ahead market clearing (minutes)
        rt_resolution: Resolution of Real Time market clearing (minutes)
        markets: Dictionary of which markets are simulated
        ordc_curved: Whether to include the curved part of the ORDC
        derating_scale: Factor for scaling derating factors
        mopr: Whether Minimum Offer Price Rule is applied
        heterogeneity: Whether investors' heterogeneous financial characteristics and technology preferences are modeled.
        reserve_penalty: High, Mid or Low penalty prices for reserves
        forecast_type: "Perfect" or "imperfect" forecasts used for price prediction.
        info_symmetry: Whether investors have symmetric information about forecast parameters.
        belief_update: Whether investors' beliefs are updated each year after actual market clearing.
        uncertainty: Whether multiple probability weighted scenarios are used instead of a deterministic forecast.
        risk_aversion: Whether investors are risk averse.
        parallel_investors: Whether investors' price prediction is to be parallelized.
        parallel_scenarios: Whether each investor's price prediction scenarios are to be parallelized.
"""

struct CaseDefinition
    base_dir::String
    sys_dir::String
    solver::JuMP.MOI.OptimizerWithAttributes
    siip_market_clearing::Bool
    start_year::Int64
    total_horizon::Int64
    rolling_horizon::Int64
    simulation_years::Int64
    num_rep_days::Int64
    da_resolution::Int64
    rt_resolution::Int64
    markets::Dict{Symbol, Bool}
    ordc_curved::Bool
    reserve_penalty::String
    derating_scale::Float64
    mopr::Bool
    heterogeneity::Bool
    forecast_type::String
    info_symmetry::Bool
    belief_update::Bool
    uncertainty::Bool
    risk_aversion::Bool
    parallel_investors::Bool
    parallel_scenarios::Bool

    function CaseDefinition(base_dir,
                            sys_dir,
                            solver,
                            siip_market_clearing,
                            start_year,
                            total_horizon,
                            rolling_horizon,
                            simulation_years,
                            num_rep_days,
                            da_resolution,
                            rt_resolution,
                            markets,
                            ordc_curved,
                            reserve_penalty,
                            derating_scale,
                            mopr,
                            heterogeneity,
                            forecast_type,
                            info_symmetry,
                            belief_update,
                            uncertainty,
                            risk_aversion,
                            parallel_investors,
                            parallel_scenarios)

        @assert total_horizon >= simulation_years

        forecast_type = lowercase(forecast_type)
        @assert forecast_type == "perfect" || lowercase(forecast_type) == "imperfect"
        @assert reserve_penalty == "High" || reserve_penalty == "Mid" || reserve_penalty == "Low"

        if forecast_type == "perfect"
            @assert info_symmetry == true
            @assert belief_update == false
            @assert uncertainty == false
            @assert risk_aversion == false
        end

        @assert da_resolution >= rt_resolution
        if !(siip_market_clearing)
            @assert da_resolution == rt_resolution
        end

        return new(base_dir,
                   sys_dir,
                   solver,
                   siip_market_clearing,
                   start_year,
                   total_horizon,
                   rolling_horizon,
                   simulation_years,
                   num_rep_days,
                   da_resolution,
                   rt_resolution,
                   markets,
                   ordc_curved,
                   reserve_penalty,
                   derating_scale,
                   mopr,
                   heterogeneity,
                   forecast_type,
                   info_symmetry,
                   belief_update,
                   uncertainty,
                   risk_aversion,
                   parallel_investors,
                   parallel_scenarios)
    end
end

function CaseDefinition(base_dir::String,
                        sys_dir::String,
                        solver::JuMP.MOI.OptimizerWithAttributes;
                        siip_market_clearing::Bool = true,
                        start_year::Int64 = 2020,
                        total_horizon::Int64 = 20,
                        rolling_horizon::Int64 = 10,
                        simulation_years::Int64 = 10,
                        num_rep_days::Int64 = 12,
                        da_resolution::Int64 = 60,
                        rt_resolution::Int64 = 5,
                        markets::Dict{Symbol, Bool} = Dict(:Energy => true, :Synchronous => true, :Primary => true, :Reg_Up => true, :Reg_Down => true,	:Flex_Up => true, :Flex_Down => true, :Capacity => true, :REC => true, :CarbonTax => true),
                        ordc_curved::Bool = true,
                        reserve_penalty::String = "Mid",
                        derating_scale::Float64 = 1.0,
                        mopr::Bool = false,
                        heterogeneity::Bool = false,
                        forecast_type::String = "perfect",
                        info_symmetry::Bool = true,
                        belief_update::Bool = false,
                        uncertainty::Bool = false,
                        risk_aversion::Bool = false,
                        parallel_investors::Bool = false,
                        parallel_scenarios::Bool = false)

    CaseDefinition(base_dir,
                   sys_dir,
                   solver,
                   siip_market_clearing,
                   start_year,
                   total_horizon,
                   rolling_horizon,
                   simulation_years,
                   num_rep_days,
                   da_resolution,
                   rt_resolution,
                   markets,
                   ordc_curved,
                   reserve_penalty,
                   derating_scale,
                   mopr,
                   heterogeneity,
                   forecast_type,
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
get_num_rep_days(case::CaseDefinition) = case.num_rep_days
get_da_resolution(case::CaseDefinition) = case.da_resolution
get_rt_resolution(case::CaseDefinition) = case.rt_resolution
get_markets(case::CaseDefinition) = case.markets
get_ordc_curved(case::CaseDefinition) = case.ordc_curved
get_reserve_penalty(case::CaseDefinition) = case.reserve_penalty
get_derating_scale(case::CaseDefinition) = case.derating_scale
get_mopr(case::CaseDefinition) = case.mopr
get_heterogeneity(case::CaseDefinition) = case.heterogeneity
get_info_symmetry(case::CaseDefinition) = case.info_symmetry
get_belief_update(case::CaseDefinition) = case.belief_update
get_forecast_type(case::CaseDefinition) = case.forecast_type
get_uncertainty(case::CaseDefinition) = case.uncertainty
get_risk_aversion(case::CaseDefinition) = case.risk_aversion
get_parallel_investors(case::CaseDefinition) = case.parallel_investors
get_parallel_scenarios(case::CaseDefinition) = case.parallel_scenarios

function get_name(case::CaseDefinition)
    if get_heterogeneity(case)
        investors = "Heterogeneous"
    else
        investors = "Homogeneous"
    end

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

    return case_name
end

function get_data_dir(case::CaseDefinition)
    base_dir = get_base_dir(case)
    case_dir = joinpath(base_dir, get_name(case))

    return case_dir
end
