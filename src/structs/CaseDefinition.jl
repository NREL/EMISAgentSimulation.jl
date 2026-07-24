"""
    This struct contains all the required data, parameters and solvers
    for creating and running AgentSimulation.
        name: Case name
        base_dir: Directory where Simulation input data for markets and investors are stored.
        sys_dir: Test system directory.
        scratch_dir: Temporary directory to save Sienna runtime data.
        outage_dir: Location of correlated outages data.
        timeseries_data_dir: Directory where timeseries data for the simulation is stored.
        cem_solver: Solvers used for optimization problems. (The solver should be able to solve QP for price prediction and MILP for SIIP production cost model)
        siip_solver: Solvers used for optimization problems. (The solver should be able to solve QP for price prediction and MILP for SIIP production cost model)
        siip_market_clearing: Whether SIIP production cost model is to be used for energy market clearing. If false, the endogenous Economic Dispatch model will be used for market clearing.
        pcm_scenario: Which scenario to use for PCM timeseries runs in the actual market clearing process.Scenarios could be based on electrification, climate models, etc. Default = scenario_1
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
        marginal_cc_switch: Whether marginal CC is used instead of average CC for new VRE and battery resources.
        forecast_type: "Perfect" or "imperfect" forecasts used for price prediction.
        max_carbon_tax_increase: Maximum annual increase in carbon prices due to under-achievement of Clean Energy Targets.
        info_symmetry: Whether investors have symmetric information about forecast parameters.
        belief_update: Whether investors' beliefs are updated each year after actual market clearing.
        uncertainty: Whether multiple probability weighted scenarios are used instead of a deterministic forecast.
        risk_aversion: Whether investors are risk averse.
        parallel_investors: Whether investors' price prediction is to be parallelized.
        parallel_scenarios: Whether each investor's price prediction scenarios are to be parallelized.
        step_size: Step size in hours for the simulation.
"""
# struct CaseDefinition
#     name::String
#     base_dir::String
#     sys_dir::String
#     scratch_dir::String
#     outage_dir::String
#     timeseries_data_dir::String
#     solver::JuMP.MOI.OptimizerWithAttributes
#     siip_market_clearing::Bool
#     pcm_scenario::String
#     start_year::Int64
#     total_horizon::Int64
#     rolling_horizon::Int64
#     simulation_years::Int64
#     rep_period_interval::Int64
#     num_rep_periods::Int64
#     avg_block_size::Int64
#     fixed_block_size::Bool
#     rep_chronology_checkpoint::Int64
#     da_resolution::Int64
#     rt_resolution::Int64
#     rps_target::String
#     markets::Dict{Symbol, Bool}
#     ordc_curved::Bool
#     ordc_unavailability_method::String
#     reserve_penalty::String
#     static_capacity_market::Bool
#     irm_scalar::Float64
#     accreditation_methodology::String
#     accreditation_metric::String
#     marginal_cc_switch::Bool
#     derating_scale::Float64
#     mopr::Bool
#     battery_cap_mkt::Bool
#     vre_reserves::Bool
#     heterogeneity::Bool
#     forecast_type::String
#     max_carbon_tax_increase::Float64
#     info_symmetry::Bool
#     belief_update::Bool
#     uncertainty::Bool
#     risk_aversion::Bool
#     parallel_investors::Bool
#     parallel_scenarios::Bool
#     md_horizon::Int64
#     md_interval::Int64
#     uc_horizon::Int64
#     uc_interval::Int64
#     ed_horizon::Int64
#     ed_interval::Int64
#     md_market::Bool
#     single_stage::Bool
#     step_size::Int64

#     function CaseDefinition(
#         simulation_settings::Dict{String, String},
#         solver::MathOptInterface.OptimizerWithAttributes,
#         dir_dict::Dict{String, String},
#         markets::Dict{Symbol, Bool},
#         )
        
#         name = simulation_settings["name"]
#         siip_market_clearing = EAS.parsebool(simulation_settings["siip_market_clearing"])
#         pcm_scenario = simulation_settings["pcm_scenario"]
#         start_year = EAS.parseint(simulation_settings["start_year"])
#         total_horizon = EAS.parseint(simulation_settings["total_horizon"])
#         rolling_horizon = EAS.parseint(simulation_settings["rolling_horizon"])
#         simulation_years = EAS.parseint(simulation_settings["simulation_years"])
#         rep_period_interval = EAS.parseint(simulation_settings["rep_period_interval"])
#         num_rep_periods = EAS.parseint(simulation_settings["num_rep_periods"])
#         avg_block_size = EAS.parseint(simulation_settings["avg_block_size"])
#         fixed_block_size = EAS.parsebool(simulation_settings["fixed_block_size"])
#         rep_chronology_checkpoint = EAS.parseint(simulation_settings["rep_chronology_checkpoint"])
#         da_resolution = EAS.parseint(simulation_settings["da_resolution"])
#         rt_resolution = EAS.parseint(simulation_settings["rt_resolution"])
#         rps_target = simulation_settings["rps_target"]
#         ordc_curved = EAS.parsebool(simulation_settings["ordc_curved"])
#         ordc_unavailability_method = simulation_settings["ordc_unavailability_method"]
#         reserve_penalty = simulation_settings["reserve_penalty"]
#         static_capacity_market = EAS.parsebool(simulation_settings["static_capacity_market"])
#         irm_scalar = EAS.parsefloat(simulation_settings["irm_scalar"])
#         accreditation_methodology = simulation_settings["accreditation_methodology"]
#         accreditation_metric = simulation_settings["accreditation_metric"]
#         marginal_cc_switch = EAS.parsebool(simulation_settings["marginal_cc_switch"])
#         derating_scale = EAS.parsefloat(simulation_settings["derating_scale"])
#         mopr = EAS.parsebool(simulation_settings["mopr"])
#         battery_cap_mkt = EAS.parsebool(simulation_settings["battery_cap_mkt"])
#         vre_reserves = EAS.parsebool(simulation_settings["vre_reserves"])
#         heterogeneity = EAS.parsebool(simulation_settings["heterogeneity"])
#         forecast_type = simulation_settings["forecast_type"]
#         max_carbon_tax_increase =  EAS.parsefloat(simulation_settings["max_carbon_tax_increase"])
#         info_symmetry = EAS.parsebool(simulation_settings["info_symmetry"])
#         belief_update =EAS.parsebool(simulation_settings["belief_update"])
#         uncertainty = EAS.parsebool(simulation_settings["uncertainty"])
#         risk_aversion = EAS.parsebool(simulation_settings["risk_aversion"])
#         parallel_investors = EAS.parsebool(simulation_settings["parallel_investors"])
#         parallel_scenarios = EAS.parsebool(simulation_settings["parallel_scenarios"])
#         md_horizon = EAS.parseint(simulation_settings["md_horizon"])
#         md_interval = EAS.parseint(simulation_settings["md_interval"])
#         uc_horizon = EAS.parseint(simulation_settings["uc_horizon"])
#         uc_interval = EAS.parseint(simulation_settings["uc_interval"])
#         ed_horizon = EAS.parseint(simulation_settings["ed_horizon"])
#         ed_interval = EAS.parseint(simulation_settings["ed_interval"])
#         md_market = EAS.parsebool(simulation_settings["multi_day_market"])

#         @assert typeof(name) == String
#         @assert typeof(siip_market_clearing) == Bool

#         @assert total_horizon >= simulation_years
#         @assert lowercase(forecast_type) in ["perfect", "imperfect"]
#         @assert lowercase(rps_target) in ["high", "mid", "low"]
#         @assert lowercase(reserve_penalty) in ["high", "mid", "low"]

#         if lowercase(forecast_type) == "perfect"
#             @assert info_symmetry == true
#             @assert belief_update == false
#             @assert uncertainty == false
#             @assert risk_aversion == false
#         end

#         @assert da_resolution >= rt_resolution
#         @assert irm_scalar >= 0.0

#         md_horizon > md_interval || error("MD horizon must be greater than MD interval")
#         uc_horizon > uc_interval || error("UC horizon must be greater than UC interval")
#         ed_horizon > ed_interval || error("ED horizon must be greater than ED interval")

#         base_dir = dir_dict["base_dir"]
#         sys_dir = dir_dict["test_system_dir"]
#         scratch_dir = dir_dict["scratch_dir"]
#         outage_filepath = dir_dict["outage_filepath"]
#         timeseries_data_dir = dir_dict["timeseries_data_dir"]

#         # TODO: add assertion to MD, DA, RT horizon and interval
#         case = new(name,
#             base_dir,
#             sys_dir,
#             scratch_dir,
#             outage_filepath,
#             timeseries_data_dir,
#             solver,
#             siip_market_clearing,
#             pcm_scenario,
#             start_year,
#             total_horizon,
#             rolling_horizon,
#             simulation_years,
#             rep_period_interval,
#             num_rep_periods,
#             avg_block_size,
#             fixed_block_size,
#             rep_chronology_checkpoint,
#             da_resolution,
#             rt_resolution,
#             rps_target,
#             markets,
#             ordc_curved,
#             ordc_unavailability_method,
#             reserve_penalty,
#             static_capacity_market,
#             irm_scalar,
#             accreditation_methodology,
#             accreditation_metric,
#             marginal_cc_switch,
#             derating_scale,
#             mopr,
#             battery_cap_mkt,
#             vre_reserves,
#             heterogeneity,
#             forecast_type,
#             max_carbon_tax_increase,
#             info_symmetry,
#             belief_update,
#             uncertainty,
#             risk_aversion,
#             parallel_investors,
#             parallel_scenarios,
#             md_horizon,
#             md_interval,
#             uc_horizon,
#             uc_interval,
#             ed_horizon,
#             ed_interval,
#             md_market,
#             single_stage,
#             step_size)

#         make_case_data_dir(case)
#         return case
#     end
# end


mutable struct CaseDefinition
    name::String
    base_dir::String
    sys_dir::String
    scratch_dir::String
    outage_dir::String
    timeseries_data_dir::String
    solver::Union{JuMP.MOI.OptimizerWithAttributes, Nothing}
    siip_market_clearing::Bool
    pcm_scenario::String
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
    seasonal_capacity_market::Bool
    irm_scalar::Float64
    accreditation_methodology::String
    accreditation_metric::String
    marginal_cc_switch::Bool
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
    md_horizon::Int64
    md_interval::Int64
    uc_horizon::Int64
    uc_interval::Int64
    ed_horizon::Int64
    ed_interval::Int64
    md_market::Bool
    single_stage::Bool
    step_size::Int64

    function CaseDefinition(name,
        base_dir,
        sys_dir,
        scratch_dir,
        outage_dir,
        timeseries_data_dir,
        solver,
        siip_market_clearing,
        pcm_scenario,
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
        seasonal_capacity_market,
        irm_scalar,
        accreditation_methodology,
        accreditation_metric,
        marginal_cc_switch,
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
        parallel_scenarios,
        md_horizon,
        md_interval,
        uc_horizon,
        uc_interval,
        ed_horizon,
        ed_interval,
        md_market,
        single_stage,
        step_size)
        @assert total_horizon >= simulation_years

        forecast_type = lowercase(forecast_type)
        @assert forecast_type == "perfect" || lowercase(forecast_type) == "imperfect"
        @assert lowercase(rps_target) == "high" || lowercase(rps_target) == "mid" ||
                lowercase(rps_target) == "low"
        @assert lowercase(reserve_penalty) == "high" ||
                lowercase(reserve_penalty) == "mid" || lowercase(reserve_penalty) == "low"

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

        # TODO: add assertion to MD, DA, RT horizon and interval
        case = new(name,
            base_dir,
            sys_dir,
            scratch_dir,
            outage_dir,
            timeseries_data_dir,
            solver,
            siip_market_clearing,
            pcm_scenario,
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
            seasonal_capacity_market,
            irm_scalar,
            accreditation_methodology,
            accreditation_metric,
            marginal_cc_switch,
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
            parallel_scenarios,
            md_horizon,
            md_interval,
            uc_horizon,
            uc_interval,
            ed_horizon,
            ed_interval,
            md_market,
            single_stage,
            step_size)

        make_case_data_dir(case)
        return case
    end
end

function CaseDefinition(name::String,
    base_dir::String,
    sys_dir::String,
    scratch_dir::String,
    outage_dir::String,
    timeseries_data_dir::String,
    markets_included::Dict{Symbol, Bool},
    solver::JuMP.MOI.OptimizerWithAttributes;
    simulation_settings::Dict{String, Any}
    )

    CaseDefinition(name,
        base_dir,
        sys_dir,
        scratch_dir,
        outage_dir,
        timeseries_data_dir,
        solver,
        siip_market_clearing = parsebool(simulation_settings["siip_market_clearing"]),
        pcm_scenario = simulation_settings["pcm_scenario"],
        start_year = parseint(simulation_settings["start_year"]),
        total_horizon = parseint(simulation_settings["total_horizon"]),
        rolling_horizon = parseint(simulation_settings["rolling_horizon"]),
        simulation_years = parseint(simulation_settings["simulation_years"]),
        rep_period_interval = parseint(simulation_settings["rep_period_interval"]),
        num_rep_periods = parseint(simulation_settings["num_rep_periods"]),
        avg_block_size = parseint(simulation_settings["avg_block_size"]),
        fixed_block_size = parsebool(simulation_settings["fixed_block_size"]),
        rep_chronology_checkpoint = parseint(simulation_settings["rep_chronology_checkpoint"]),
        da_resolution = parseint(simulation_settings["da_resolution"]),
        rt_resolution = parseint(simulation_settings["rt_resolution"]),
        rps_target = simulation_settings["rps_target"],
        markets = markets_included,
        ordc_curved = parsebool(simulation_settings["ordc_curved"]),
        ordc_unavailability_method = simulation_settings["ordc_unavailability_method"],
        reserve_penalty = simulation_settings["reserve_penalty"],
        static_capacity_market = parsebool(simulation_settings["static_capacity_market"]),
        seasonal_capacity_market = parsebool(
            get(simulation_settings, "seasonal_capacity_market", "false"),
        ),
        irm_scalar = parsefloat(simulation_settings["irm_scalar"]),
        accreditation_methodology = simulation_settings["accreditation_methodology"],
        accreditation_metric = simulation_settings["accreditation_metric"],
        marginal_cc_switch = parsebool(simulation_settings["marginal_cc_switch"]),
        derating_scale = parsefloat(simulation_settings["derating_scale"]),
        mopr = parsebool(simulation_settings["mopr"]),
        battery_cap_mkt = parsebool(simulation_settings["battery_cap_mkt"]),
        vre_reserves = parsebool(simulation_settings["vre_reserves"]),
        heterogeneity = parsebool(simulation_settings["heterogeneity"]),
        forecast_type = simulation_settings["forecast_type"],
        max_carbon_tax_increase = parsefloat(simulation_settings["max_carbon_tax_increase"]),
        info_symmetry = parsebool(simulation_settings["info_symmetry"]),
        belief_update = parsebool(simulation_settings["belief_update"]),
        uncertainty = parsebool(simulation_settings["uncertainty"]),
        risk_aversion = parsebool(simulation_settings["risk_aversion"]),
        parallel_investors = parsebool(simulation_settings["parallel_investors"]),
        parallel_scenarios = parsebool(simulation_settings["parallel_scenarios"]),
        md_horizon = parseint(simulation_settings["md_horizon"]),
        md_interval = parseint(simulation_settings["md_interval"]),
        uc_horizon = parseint(simulation_settings["uc_horizon"]),
        uc_interval = parseint(simulation_settings["uc_interval"]),
        ed_horizon = parseint(simulation_settings["ed_horizon"]),
        ed_interval = parseint(simulation_settings["ed_interval"]),
        md_market = parsebool(simulation_settings["md_market"]),
    )
end

function CaseDefinition(name::String,
    base_dir::String,
    sys_dir::String,
    scratch_dir::String,
    outage_dir::String,
    timeseries_data_dir::String,
    solver::JuMP.MOI.OptimizerWithAttributes;
    siip_market_clearing::Bool = true,
    pcm_scenario::String = "scenario_1",
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
    markets::Dict{Symbol, Bool} = Dict(
        :Energy => true,
        :Synchronous => true,
        :Primary => true,
        :Reg_Up => true,
        :Reg_Down => true,
        :Flex_Up => true,
        :Flex_Down => true,
        :Capacity => true,
        :REC => true,
        :CarbonTax => true,
    ),
    ordc_curved::Bool = true,
    ordc_unavailability_method::String = "Convolution",
    reserve_penalty::String = "Mid",
    static_capacity_market::Bool = true,
    seasonal_capacity_market::Bool = false,
    irm_scalar::Float64 = 1.0,
    accreditation_methodology::String = "TopNetLoad",
    accreditation_metric::String = "None",
    marginal_cc_switch::Bool = true,
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
    parallel_scenarios::Bool = false,
    md_horizon::Int64 = 168,
    md_interval::Int64 = 168,
    uc_horizon::Int64 = 36,
    uc_interval::Int64 = 24,
    ed_horizon::Int64 = 2,
    ed_interval::Int64 = 1,
    md_market::Bool = false,
    single_stage::Bool = false,
    step_size::Int64 = 1)
    CaseDefinition(name,
        base_dir,
        sys_dir,
        scratch_dir,
        outage_dir,
        timeseries_data_dir,
        solver,
        siip_market_clearing,
        pcm_scenario,
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
        seasonal_capacity_market,
        irm_scalar,
        accreditation_methodology,
        accreditation_metric,
        marginal_cc_switch,
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
        parallel_scenarios,
        md_horizon,
        md_interval,
        uc_horizon,
        uc_interval,
        ed_horizon,
        ed_interval,
        md_market,
        single_stage,
        step_size)
end

get_base_dir(case::CaseDefinition) = case.base_dir
get_sys_dir(case::CaseDefinition) = case.sys_dir
get_scratch_dir(case::CaseDefinition) = case.scratch_dir
get_outage_dir(case::CaseDefinition) = case.outage_dir
get_timeseries_data_dir(case::CaseDefinition) = case.timeseries_data_dir
get_solver(case::CaseDefinition) = case.solver
get_siip_market_clearing(case::CaseDefinition) = case.siip_market_clearing
get_pcm_scenario(case::CaseDefinition) = case.pcm_scenario
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
get_seasonal_capacity_market(case::CaseDefinition) = case.seasonal_capacity_market
get_irm_scalar(case::CaseDefinition) = case.irm_scalar
get_ordc_unavailability_method(case::CaseDefinition) = case.ordc_unavailability_method
get_accreditation_methodology(case::CaseDefinition) = case.accreditation_methodology
get_accreditation_metric(case::CaseDefinition) = case.accreditation_metric
get_marginal_cc_switch(case::CaseDefinition) = case.marginal_cc_switch
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
get_md_horizon(case::CaseDefinition) = case.md_horizon
get_md_interval(case::CaseDefinition) = case.md_interval
get_uc_horizon(case::CaseDefinition) = case.uc_horizon
get_uc_interval(case::CaseDefinition) = case.uc_interval
get_ed_horizon(case::CaseDefinition) = case.ed_horizon
get_ed_interval(case::CaseDefinition) = case.ed_interval
get_md_market(case::CaseDefinition) = case.md_market
get_single_stage(case::CaseDefinition) = case.single_stage
get_step_size(case::CaseDefinition) = case.step_size

set_solver!(case::CaseDefinition, solver) = (case.solver = solver)

function get_name(case::CaseDefinition)
    return "$(case.name)"
    # return "$(case.name)_$(get_rps_target(case))_RPS"
end

function get_data_dir(case::CaseDefinition)
    return joinpath(get_base_dir(case), get_name(case))
end
