module EMISAgentSimulation

#################################################################################
# Exports

# Export Structs
using Base: Tuple, Float64, project_deps_get
export CaseDefinition

export AgentSimulation

export Investor

export KalmanFilter
export InvestorBelief
export Scenario
export Forecast
export Perfect
export Imperfect

export RiskPreference
export RiskNeutral
export RiskAverse

export MarketPrices

export Finance

export Product
export OperatingProduct
export ReserveDirection
export ReserveUpEMIS
export ReserveDownEMIS
export ReserveUpEMIS
export OperatingReserve
export Energy
export Capacity
export REC
export CarbonTax
export Inertia
export ResourceAdequacy

export BuildPhase
export Existing
export Planned
export Queue
export Option
export Retired

export DeviceEMIS
export Project
export GeneratorEMIS
export StorageEMIS

export BatteryTech
export BatteryEMIS

export HydroTech
export HydroGenEMIS

export RenewableTech
export RenewableGenEMIS

export ThermalTech
export ThermalGenEMIS
export ThermalFastStartSIIP

export ZonalLine

export MarketProject
export EnergyMarket
export ReserveUpMarket
export ReserveDownMarket
export ReserveORDCMarket
export CapacityMarket
export RECMarket
export InertiaMarket
export MarketCollection
export MarketClearingProblem

# Export Functions

# Export Core functions
export calculate_annualized_costs
export calculate_capcost_and_wacc
export calculate_derating_data
export calculate_npv
export capacity_market_clearing
export cem
export cluster
export construct_ordc
export construct_gen_unavail_distribution
export construct_net_load_forecast_error_distribution
export create_agent_simulation
export create_capacity_demand_curve
export create_expected_marketdata
export create_ordc_market
export create_PSY_generator
export create_realized_marketdata
export economicdispatch
export energy_mkt_clearing
export finish_construction!
export make_investments!
export rec_market_clearing
export retire_old!
export retire_unprofitable!
export run_agent_simulation
export start_construction!

# Export Utility functions
export calculate_operating_profit
export calculate_required_processes
export create_parallel_workers
export dir_exists
export leaftypes
export make_case_data_dir
export make_results_dir
export read_data
export remove_leap_day!
export size_in_MW
export update_derating_factor!
export update_installed_cap!
export update_PSY_timeseries!
export update_PSY_outage_timeseries!
export update_simulation_derating_data!
export write_data
export parsebool
export parseint
export parsefloat

# Export Getter Functions
export get_accepted_perc
export get_active_power_limit
export get_active_power_limits
export get_activeprojects
export get_allprojects
export get_all_techs
export get_annual_cashflow
export get_annual_growth
export get_base_dir
export get_belief_update
export get_bus
export get_capacity_forward_years
export get_capacity_price
export get_capacity_bid
export get_capacity_factors
export get_cap_cost_multiplier
export get_capex_years
export get_carbon_tax
export get_case
export get_constant
export get_construction_year
export get_data_dir
export get_decidedprojects
export get_decision_year
export get_derating
export get_derating_data
export get_device_size
export get_discount_rate
export get_effective_investment_cost
export get_efficiency
export get_end_life_year
export get_emission
export get_energy_price
export get_error_covariance_estimate
export get_existing
export get_expected_npv
export get_expected_utility
export get_finance_data
export get_fixed_OM_cost
export get_forecast
export get_forecast_type
export get_from_zone
export get_fuel
export get_heat_rate
export get_heterogeneity
export get_hour_weight
export get_info_symmetry
export get_input_active_power_limits
export get_investment_cost
export get_investor_data
export get_investors
export get_iteration_year
export get_kalman_filter
export get_lag_time
export get_life_time
export get_line_rating
export get_lines
export get_marginal_cost
export get_market_prices
export get_markets
export get_max_annual_projects
export get_maxcap
export get_max_limit
export get_maxperc
export get_measurement_covariance
export get_mincap
export get_multiplier
export get_name
export get_num_rep_periods
export get_operation_cost
export get_options
export get_output_active_power_limits
export get_ownedby
export get_parallel_investors
export get_parallel_scenarios
export get_parameter_multipliers
export get_parameter_values
export get_peak_load
export get_planned
export get_prices
export get_probability
export get_process_covariance
export get_products
export get_projects
export get_queue
export get_queue_cost
export get_ramp_limits
export get_realized_profit
export get_expected_rec_certificates
export get_rec_price
export get_rec_bid
export get_rep_hour_weight
export get_reservedown_price
export get_reserveup_price
export get_results_dir
export get_retired
export get_retirement_year
export get_risk_aversion
export get_risk_coefficient
export get_risk_preference
export get_rolling_horizon
export get_scenario_npv
export get_scenario_profit
export get_scenario_utility
export get_scenario_data
export get_siip_market_clearing
export get_simulation_years
export get_soc
export get_solver
export get_start_year
export get_state_estimate
export get_state_measurement_matrix
export get_state_transition_matrix
export get_storage_capacity
export get_system_ED
export get_system_UC
export get_system_services
export get_tech
export get_sys_dir
export get_time_limits
export get_to_zone
export get_total_horizon
export get_type
export get_uncertainty
export get_zone
export get_zones
export update_belief!

# Export Find functions
export find_active_invested_projects
export find_energy_product
export find_operating_products
export find_option_projects

# Export Setter Functions
export set_investors!
export set_iteration_year!
export set_peak_load!

#################################################################################
# Imports
import AxisArrays
import Clustering
import ClusterManagers
import CSV
import DataFrames
import Dates
import Distributed
import Distributions
import FileIO
import InteractiveUtils
import JuMP
export optimizer_with_attributes
import JLD2
import LinearAlgebra
import PooledArrays
import PowerSystems
import PowerSimulations

using EMISExtensions
using PRAS
using HiGHS
using OrderedCollections
import InfrastructureSystems
import ReliablePowerSimulations

const PSY = PowerSystems
const PSI = PowerSimulations
const EMISEx = EMISExtensions
const IS = InfrastructureSystems
const RPSI = ReliablePowerSimulations

import Random
import UUIDs
import StatsBase
import Statistics
import TimeSeries
const TS = TimeSeries

import TimeZones

import Base.convert
import JuMP.value
import MathOptInterface
import MathOptInterface: AbstractOptimizer

const MOI = MathOptInterface

using Revise

import PowerSystems:
    get_value,
    set_value

################################################################################
# Includes

# Include all structs:

include("structs/products/Product.jl")
include("structs/products/Energy.jl")
include("structs/products/OperatingReserve.jl")
include("structs/products/Capacity.jl")
include("structs/products/REC.jl")
include("structs/products/CarbonTax.jl")
include("structs/products/Inertia.jl")
include("structs/ResourceAdequacy.jl")

include("structs/Finance.jl")

include("structs/devices/BuildPhase.jl")
include("structs/devices/Project.jl")
include("structs/devices/BatteryEMIS.jl")
include("structs/devices/HydroGenEMIS.jl")
include("structs/devices/RenewableGenEMIS.jl")
include("structs/devices/ThermalGenEMIS.jl")
include("structs/devices/ThermalFastStartSIIP.jl")

include("structs/CaseDefinition.jl")
include("structs/MarketPrices.jl")
include("structs/Scenario.jl")

include("structs/InvestorBelief.jl")

include("structs/Forecast.jl")
include("structs/RiskPreference.jl")
include("structs/ZonalLine.jl")

include("structs/Investor.jl")

include("structs/AgentSimulation.jl")

include("structs/market_structs/MarketProject.jl")
include("structs/market_structs/EnergyMarket.jl")
include("structs/market_structs/ReserveUpMarket.jl")
include("structs/market_structs/ReserveDownMarket.jl")
include("structs/market_structs/ReserveORDCMarket.jl")
include("structs/market_structs/CapacityMarket.jl")
include("structs/market_structs/RECMarket.jl")
include("structs/market_structs/InertiaMarket.jl")
include("structs/market_structs/MarketCollection.jl")
include("structs/market_structs/MarketClearingProblem.jl")

# Include utilities:
include("utils/general.jl")                 # General utilities.
include("utils/market_utils.jl")           # Utils for actual and expected market clearing modules.
include("utils/product_utils.jl")          # Utils for products provided by the projects
include("utils/siip_psy_utils.jl")         # Utils for SIIP PowerSystems.
include("utils/conversion_utils.jl")       # Define new convert functions for changing buildphase of projects.
include("utils/read_and_write_utils.jl")   # Read and write utils.
include("utils/parallel_utils.jl")         # Utils for parallelizing price prediction runs.
include("utils/finance_utils.jl")          # Functions for calculating adjusted CAPEX and WACC

#Include files containing functions for creating the simulation structs from the given data.
include("struct_creators/simulation_structs/product_creator.jl")
include("struct_creators/simulation_structs/project_creator.jl")
include("struct_creators/simulation_structs/project_vector_creator.jl")
include("struct_creators/simulation_structs/investor_creator.jl")
include("struct_creators/simulation_structs/rts_psy_creator.jl")
include("struct_creators/simulation_structs/siip_psy_device_creator.jl")
include("struct_creators/simulation_structs/agent_simulation_creator.jl")

#Include files containing functions for creating the structs used for market clearing modules.
include("struct_creators/market_structs/market_project_creator.jl")
include("struct_creators/market_structs/cem_creator.jl")
include("struct_creators/market_structs/economic_dispatch_creator.jl")

#Include files for constructing ORDC
include("markets_simulation/ordc_construction/error_distributions.jl")
include("markets_simulation/ordc_construction/ordc_construction.jl")
include("markets_simulation/ordc_construction/ordc_market_creator.jl")

#Include PRAS Resource adequacy functions
include("resource_adequacy/conv.jl")
include("resource_adequacy/PSY2PRAS.jl")
include("resource_adequacy/parsers/power_system_table_data.jl")
include("resource_adequacy/ra_utils.jl")
include("resource_adequacy/generator_unavailability.jl")

#Include Test System Parsers
include("test_system_parsers/test_system_reader.jl")
include("test_system_parsers/rts_reader.jl")

#Include Expected and Actual Market Simultion functions.
include("markets_simulation/cem.jl")
include("markets_simulation/expected_market_simulation.jl")
include("markets_simulation/economicdispatch.jl")
include("markets_simulation/siip_simulation_definition.jl")
include("markets_simulation/actual_energy_mkt_clearing.jl")
include("markets_simulation/actual_capacity_mkt_clearing.jl")
include("markets_simulation/actual_rec_mkt_clearing.jl")
include("markets_simulation/actual_market_simulation.jl")

#Include Investor functions
    include("investor_functions/investor_iteration.jl")            # Runs investors annual iteration.

    #### Predictions #####################
    include("investor_functions/prediction/prediction_methodology.jl")  # Functions for running investors' price prediction methodology.
    include("investor_functions/prediction/operating_profit.jl")        # Functions for calculating expected operating market profits.
    include("investor_functions/prediction/capacity_profit.jl")         # Functions for calculating expected capacity market profits.
    include("investor_functions/prediction/REC_profit.jl")              # Functions for calculating expected REC market profits.
    include("investor_functions/prediction/total_profit.jl")            # Functions for updating the expected profits from all markets.

    #### Decisions ######################
    include("investor_functions/decisions/buildphase_conversion.jl")         # Functions for evaluating when to convert the buildphase of projects.
    include("investor_functions/decisions/investment_decision.jl")           # Functions for making investment decisions
    include("investor_functions/decisions/retirement_decision.jl")           # Functions for making retirement decisions

    #### Decision Metrics ######
    include("investor_functions/decision_metrics/npv_functions.jl")          # Functions for calculating NPV.
    include("investor_functions/decision_metrics/utility_functions.jl")      # Functions for calculating expected utility.

    ### Realized Profits and Updates ######
    include("investor_functions/realized_profits_calculator.jl")   # Functions for calculating realized profits from different markets
    include("investor_functions/annual_updates.jl")                # Functions for updating investor revenues and forecasts each year.


#Include derating factor updating methodology
include("derating_factor_updates/derating_factor_calculator.jl")  # Derating factor calculation for renewables

#Include k-medoids clustering for representative days selection.
include("representative_days/clustering.jl")

#Include main investment simulation function.
include("investment_simulation.jl")

end
