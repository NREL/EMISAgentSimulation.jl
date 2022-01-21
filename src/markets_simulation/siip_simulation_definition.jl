# This file contains functions for calling, running and post-processing
# SIIP PSI Simulation for actual Energy and A/S Market Clearing

function adjust_reserve_voll!(sys::PSY.System,
                             problem::PSI.OperationsProblem,
                             simulation_dir::String,
                             reserve_penalty::String,
                             zones::Vector{String},
                             default_balance_slack_cost::Float64,
                             default_service_slack_cost::Float64,
                             energy_voll_cost::AxisArrays.AxisArray{Float64, 1}
                             )

    base_power = PSY.get_base_power(sys)
    services = get_system_services(sys)

    optimization_container = PSI.get_optimization_container(problem)
    variables = PSI.get_variables(optimization_container)
    cost_expression = optimization_container.cost_function

    for zone in zones
        bus = find_zonal_bus(zone, sys)
        slack_coefficients = ["γ⁺", "γ⁻"]
        for c in slack_coefficients
            variable_name = Symbol("$(c)__P")
            slack_variables = variables[variable_name][PSY.get_number(bus), :]
            delta_cost = energy_voll_cost[zone] * base_power - default_balance_slack_cost
            for v in slack_variables
                JuMP.add_to_expression!(
                    cost_expression,
                    v * delta_cost,
                )
            end
        end
    end

    for s in services
        name = PSY.get_name(s)
        delta_cost = 0.0
        slack_variables = Symbol[]

        if name == "CleanEnergyConstraint"
            delta_cost = default_service_slack_cost * 5
            #delta_cost = -default_service_slack_cost
            slack_variables = variables[Symbol("γ⁺__CleanEnergyConstraint")]
            println(PSY.get_requirement(s))
        else
            reserve_data = read_data(joinpath(simulation_dir, "markets_data", "$(reserve_penalty)_reserve_penalty", "$(name).csv"))
            variable_name = Symbol("γ⁺__$(PSY.get_name(s))")
            if variable_name in keys(variables)
                slack_variables = variables[variable_name]
                penalty_price = reserve_data[1, "price_cap"] * base_power
                delta_cost = penalty_price - default_service_slack_cost
            end
        end


        for v in slack_variables
            JuMP.add_to_expression!(
                cost_expression,
                v * delta_cost,
            )
        end
    end

end

function scale_voll(price::Union{Array{Float64, 2}, AxisArrays.AxisArray{Float64, 2}}, resolution::Int64)
    for t in 1:size(price, 2)
        if price[1, t] >= 990.0
            price[1, t] = price[1, t] * resolution / 60
        end
    end
end

function scale_voll(price::AxisArrays.AxisArray{Float64, 3}, resolution::Int64)
    for z in 1:size(price, 1)
        for t in 1:size(price, 3)
            if price[z, 1, t] >= 990.0
                price[z, 1, t] = price[z, 1, t] * resolution / 60
            end
        end
    end
end

"""
This function returns start-up costs for Thermal generators from PSI Simulation.
"""
function get_start_costs(device::PSY.ThermalStandard,
                              results::Dict{Symbol, DataFrames.DataFrame},
                              data_length_uc::Int64
                                        )
    start_ups = results[:start__ThermalStandard][:, Symbol(get_name(device))]
    start_up_costs = start_ups * PSY.get_start_up(PSY.get_operation_cost(device))
    return start_up_costs
end

function get_start_costs(device::PSYE.ThermalCleanEnergy,
                              results::Dict{Symbol, DataFrames.DataFrame},
                              data_length_uc::Int64
                                        )
    start_ups = results[:start__ThermalCleanEnergy][:, Symbol(get_name(device))]
    start_up_costs = start_ups * PSY.get_start_up(PSY.get_operation_cost(device))
    return start_up_costs
end

function get_start_costs(device::ThermalFastStartSIIP,
                              results::Dict{Symbol, DataFrames.DataFrame},
                              data_length_uc::Int64
                                        )
    start_ups = results[:start__ThermalFastStartSIIP][:, Symbol(get_name(device))]
    start_up_costs = start_ups * PSY.get_start_up(PSY.get_operation_cost(device))
    return start_up_costs
end

"""
This function returns total start-up costs for other generators from PSI Simulation.
"""
function get_start_costs(device::PSY.Device,
                              results::Dict{Symbol, DataFrames.DataFrame},
                              data_length_uc::Int64
                                        )
    return zeros(data_length_uc)
end

"""
This function returns shut-down costs for Thermal generators from PSI Simulation.
"""
function get_shut_costs(device::PSY.ThermalStandard,
                       results::Dict{Symbol, DataFrames.DataFrame},
                       data_length_uc::Int64
                                        )
    shut_downs = results[:stop__ThermalStandard][:, Symbol(get_name(device))]
    shut_down_costs = shut_downs * PSY.get_shut_down(PSY.get_operation_cost(device))
    return shut_down_costs
end

function get_shut_costs(device::PSYE.ThermalCleanEnergy,
                       results::Dict{Symbol, DataFrames.DataFrame},
                       data_length_uc::Int64
                                        )
    shut_downs = results[:stop__ThermalCleanEnergy][:, Symbol(get_name(device))]
    shut_down_costs = shut_downs * PSY.get_shut_down(PSY.get_operation_cost(device))
    return shut_down_costs
end

function get_shut_costs(device::ThermalFastStartSIIP,
                       results::Dict{Symbol, DataFrames.DataFrame},
                       data_length_uc::Int64
                                        )
    shut_downs = results[:stop__ThermalFastStartSIIP][:, Symbol(get_name(device))]
    shut_down_costs = shut_downs * PSY.get_shut_down(PSY.get_operation_cost(device))
    return shut_down_costs
end
"""
This function returns total shut-down costs for other generators from PSI Simulation.
"""
function get_shut_costs(device::PSY.Device,
                         results::Dict{Symbol, DataFrames.DataFrame},
                         data_length_uc::Int64
                                        )
    return zeros(data_length_uc)
end


"""
This function returns realized capacity factors for ThermalStandard generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.ThermalStandard,
                                        results::Dict{Symbol, DataFrames.DataFrame},
                                        results_uc::Dict{Symbol, DataFrames.DataFrame}
                                        )
    energy_production = results[:P__ThermalStandard][:, Symbol(get_name(device))]
    capacity_factors = energy_production / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for ThermalCleanEnergy generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSYE.ThermalCleanEnergy,
                                        results::Dict{Symbol, DataFrames.DataFrame},
                                        results_uc::Dict{Symbol, DataFrames.DataFrame}
                                        )
    energy_production = results[:P__ThermalCleanEnergy][:, Symbol(get_name(device))]
    capacity_factors = energy_production / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for ThermalFastStartSIIP generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::ThermalFastStartSIIP,
                                        results::Dict{Symbol, DataFrames.DataFrame},
                                        results_uc::Dict{Symbol, DataFrames.DataFrame}
                                        )
    energy_production = results[:P__ThermalFastStartSIIP][:, Symbol(get_name(device))]
    capacity_factors = energy_production / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for Renewable generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.RenewableDispatch,
                                        results::Dict{Symbol, DataFrames.DataFrame},
                                        results_uc::Dict{Symbol, DataFrames.DataFrame}
                                        )
    energy_production = results[:P__RenewableDispatch][:, Symbol(get_name(device))]
    capacity_factors = energy_production / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for Hydropower generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.HydroDispatch,
                                        results::Dict{Symbol, DataFrames.DataFrame},
                                        results_uc::Dict{Symbol, DataFrames.DataFrame}
                                        )
    energy_production = results[:P__HydroDispatch][:, Symbol(get_name(device))]
    capacity_factors = energy_production / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for Hydropower generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.HydroEnergyReservoir,
                                        results::Dict{Symbol, DataFrames.DataFrame},
                                        results_uc::Dict{Symbol, DataFrames.DataFrame}
                                        )
    energy_production = results[:P__HydroEnergyReservoir][:, Symbol(get_name(device))]
    capacity_factors = energy_production / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for Hydropower generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.GenericBattery,
                                        results::Dict{Symbol, DataFrames.DataFrame},
                                        results_uc::Dict{Symbol, DataFrames.DataFrame}
                                        )
    energy_production = results[:Pout__GenericBattery][:, Symbol(get_name(device))] - results[:Pin__GenericBattery][:, Symbol(get_name(device))]
    capacity_factors = energy_production / get_device_size(device)
    generation = filter(x -> x > 0, capacity_factors)

    energy_production_uc = results_uc[:Pout__GenericBattery][:, Symbol(get_name(device))] - results_uc[:Pin__GenericBattery][:, Symbol(get_name(device))]
    capacity_factors_uc = energy_production_uc / get_device_size(device)
    generation_uc = filter(x -> x > 0, capacity_factors_uc)
    return capacity_factors
end

"""
This function returns nothing if Service is not of ReserveUp or ReserveDown type.
"""
function update_realized_reserve_perc!(device::PSY.Device,
                                        service::S,
                                        results_ed::Dict{Symbol, DataFrames.DataFrame},
                                        results_uc::Dict{Symbol, DataFrames.DataFrame},
                                        reserve_perc::Dict{String, Dict{String, Array{Float64, 2}}},
                                        inertia_perc::Dict{String, Array{Float64, 2}},
                                        rt_products::Vector{SubString{String}},
                                        only_da_products::Vector{String}) where S <: PSY.Service
    return
end

"""
This function returns realized reserve up provision percentages from PSI Simulation.
"""
function update_realized_reserve_perc!(device::PSY.Device,
                                        service::PSY.VariableReserve{PSY.ReserveUp},
                                        results_ed::Dict{Symbol, DataFrames.DataFrame},
                                        results_uc::Dict{Symbol, DataFrames.DataFrame},
                                        reserve_perc::Dict{String, Dict{String, Array{Float64, 2}}},
                                        inertia_perc::Dict{String, Array{Float64, 2}},
                                        rt_products::Vector{SubString{String}},
                                        only_da_products::Vector{String})

    service_name = PSY.get_name(service)

    if service_name in rt_products
        reserve_provision = results_ed[Symbol("$(service_name)__VariableReserve_PowerSystems.ReserveUp")][:, Symbol(get_name(device))]
        reserve_perc_value = reserve_provision / get_device_size(device)
        reserve_perc[get_name(device)][service_name][1, :] = reserve_perc_value

    elseif service_name in only_da_products
        reserve_provision = results_uc[Symbol("$(service_name)__VariableReserve_PowerSystems.ReserveUp")][:, Symbol(get_name(device))]
        reserve_perc_value = reserve_provision / get_device_size(device)
        reserve_perc[get_name(device)][service_name][1, :] = reserve_perc_value

    end

    return
end

"""
This function returns realized ordc provision percentages from PSI Simulation.
"""
function update_realized_reserve_perc!(device::PSY.Device,
                                        service::PSY.ReserveDemandCurve{PSY.ReserveUp},
                                        results_ed::Dict{Symbol, DataFrames.DataFrame},
                                        results_uc::Dict{Symbol, DataFrames.DataFrame},
                                        reserve_perc::Dict{String, Dict{String, Array{Float64, 2}}},
                                        inertia_perc::Dict{String, Array{Float64, 2}},
                                        rt_products::Vector{SubString{String}},
                                        only_da_products::Vector{String})

    service_name = PSY.get_name(service)

    if service_name in rt_products
        reserve_provision = results_ed[Symbol("$(service_name)__ReserveDemandCurve_PowerSystems.ReserveUp")][:, Symbol(get_name(device))]
        reserve_perc_value = reserve_provision / get_device_size(device)
        reserve_perc[get_name(device)][service_name][1, :] = reserve_perc_value

    elseif service_name in only_da_products
        reserve_provision = results_uc[Symbol("$(service_name)__ReserveDemandCurve_PowerSystems.ReserveUp")][:, Symbol(get_name(device))]
        reserve_perc_value = reserve_provision / get_device_size(device)
        reserve_perc[get_name(device)][service_name][1, :] = reserve_perc_value

    end
    return
end

"""
This function returns realized reserve down provision percentages from PSI Simulation.
"""
function update_realized_reserve_perc!(device::PSY.Device,
                                        service::PSY.VariableReserve{PSY.ReserveDown},
                                        results_ed::Dict{Symbol, DataFrames.DataFrame},
                                        results_uc::Dict{Symbol, DataFrames.DataFrame},
                                        reserve_perc::Dict{String, Dict{String, Array{Float64, 2}}},
                                        inertia_perc::Dict{String, Array{Float64, 2}},
                                        rt_products::Vector{SubString{String}},
                                        only_da_products::Vector{String})

    service_name = PSY.get_name(service)

    if service_name in rt_products
        reserve_provision = results_ed[Symbol("$(service_name)__VariableReserve_PowerSystems.ReserveDown")][:, Symbol(get_name(device))]
        reserve_perc_value = reserve_provision / get_device_size(device)
        reserve_perc[get_name(device)][service_name][1, :] = reserve_perc_value

    elseif service_name in only_da_products
        reserve_provision = results_uc[Symbol("$(service_name)__VariableReserve_PowerSystems.ReserveDown")][:, Symbol(get_name(device))]
        reserve_perc_value = reserve_provision / get_device_size(device)
        reserve_perc[get_name(device)][service_name][1, :] = reserve_perc_value

    end
    return
end

"""
This function returns realized reserve down provision percentages from PSI Simulation.
"""
function update_realized_reserve_perc!(device::PSY.Device,
                                        service::PSY.VariableReserve{PSY.ReserveSymmetric},
                                        results_ed::Dict{Symbol, DataFrames.DataFrame},
                                        results_uc::Dict{Symbol, DataFrames.DataFrame},
                                        reserve_perc::Dict{String, Dict{String, Array{Float64, 2}}},
                                        inertia_perc::Dict{String, Array{Float64, 2}},
                                        rt_products::Vector{SubString{String}},
                                        only_da_products::Vector{String})

    service_name = PSY.get_name(service)
    if service_name == "Inertia"
        inertia_provision = results_ed[Symbol("$(service_name)__VariableReserve_PowerSystems.ReserveSymmetric")][:, Symbol(get_name(device))]
        inertia_perc_value = inertia_provision / get_device_size(device)
        inertia_perc[get_name(device)][1, :] = inertia_perc_value
    end
    return
end

"""
This function creates the Unit Commitment template for PSI Simulation.
"""
#TODO: Update needed
function create_uc_template(inertia_product)

    if !(isempty(inertia_product))

        template = PSI.OperationsProblemTemplate(PSI.PM.NFAPowerModel)
        PSI.set_device_model!(template, PSY.ThermalStandard, PSIE.ThermalInertiaStandardUnitCommitment)
        PSI.set_device_model!(template, PSYE.ThermalCleanEnergy, PSIE.ThermalEmisStandardUnitCommitment)
        PSI.set_device_model!(template, ThermalFastStartSIIP, PSIE.ThermalInertiaStandardUnitCommitment)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSIE.RenewableEmisDispatch)
        PSI.set_device_model!(template, PSY.RenewableFix, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.PowerLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroEnergyReservoir, PSIE.HydroEmisCommitmentRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, PSIE.HydroEmisCommitmentRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(template, PSY.GenericBattery, PSIE.BookKeepingwInertia)
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.HVDCLine, PSI.HVDCLossless)
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveUp}, PSI.RangeReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveDown}, PSI.RangeReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.ReserveDemandCurve{PSY.ReserveUp}, PSIE.QuadraticCostRampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSYE.InertiaReserve{PSY.ReserveSymmetric}, PSIE.VariableInertiaReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSYE.CleanEnergyReserve{PSY.ReserveSymmetric},PSIE.EnergyRequirementReserve))
    else
        template = PSI.OperationsProblemTemplate(PSI.PM.NFAPowerModel)
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalStandardUnitCommitment)
        PSI.set_device_model!(template, PSYE.ThermalCleanEnergy, PSIE.ThermalCleanStandardUnitCommitment)
        PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalStandardUnitCommitment)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSIE.RenewableCleanEnergyDispatch)
        PSI.set_device_model!(template, PSY.RenewableFix, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.PowerLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroEnergyReservoir, PSIE.HydroCleanEnergyRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, PSIE.HydroCleanEnergyRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(template, PSY.GenericBattery, PSI.BookKeepingwReservation)
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.HVDCLine, PSI.HVDCLossless)
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveUp}, PSI.RangeReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveDown}, PSI.RangeReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.ReserveDemandCurve{PSY.ReserveUp}, PSIE.QuadraticCostRampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSYE.InertiaReserve{PSY.ReserveSymmetric}, PSIE.VariableInertiaReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSYE.CleanEnergyReserve{PSY.ReserveSymmetric},PSIE.EnergyRequirementReserve))
    end

    return template
end

"""
This function creates the Economic Dispatch template for PSI Simulation.
"""
#TODO: Update needed
function create_ed_template(inertia_product)

    if !(isempty(inertia_product))
        template = PSI.OperationsProblemTemplate(PSI.PM.NFAPowerModel)
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalDispatch)
        PSI.set_device_model!(template, PSYE.ThermalCleanEnergy, PSI.ThermalDispatch)
        PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalDispatch)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSIE.RenewableFullDispatchInertia)
        PSI.set_device_model!(template, PSY.RenewableFix, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.PowerLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroEnergyReservoir, PSI.HydroDispatchRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, PSI.HydroDispatchRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(template, PSY.GenericBattery, PSI.BookKeeping)
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.HVDCLine, PSI.HVDCLossless)
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveUp}, PSI.RampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveDown}, PSI.RampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.ReserveDemandCurve{PSY.ReserveUp}, PSIE.QuadraticCostRampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSYE.InertiaReserve{PSY.ReserveSymmetric}, PSIE.VariableInertiaReserve))

    else
        template = PSI.OperationsProblemTemplate(PSI.PM.NFAPowerModel)
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalDispatch)
        PSI.set_device_model!(template, PSYE.ThermalCleanEnergy, PSI.ThermalDispatch)
        PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalDispatch)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSI.RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableFix, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.PowerLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroEnergyReservoir, PSI.HydroDispatchRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, PSI.HydroDispatchRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(template, PSY.GenericBattery, PSI.BookKeeping)
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.HVDCLine, PSI.HVDCLossless)
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveUp}, PSI.RampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveDown}, PSI.RampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.ReserveDemandCurve{PSY.ReserveUp}, PSIE.QuadraticCostRampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSYE.InertiaReserve{PSY.ReserveSymmetric}, PSIE.VariableInertiaReserve))
    end

    return template
end

"""
This function creates the Problem for PSI Simulation.
"""

function create_problem(template::PSI.OperationsProblemTemplate, sys::PSY.System, type::String, solver::JuMP.MOI.OptimizerWithAttributes, inertia_product)

    duals = [:nodal_balance_active__Bus]

    services = get_system_services(sys)
    if in(PSY.VariableReserve{PSY.ReserveUp}, typeof.(services))
        push!(duals, Symbol("requirement__VariableReserve_PowerSystems.ReserveUp"))
    end
    if in(PSY.ReserveDemandCurve{PSY.ReserveUp}, typeof.(services))
        push!(duals, Symbol("requirement__ReserveDemandCurve_PowerSystems.ReserveUp"))
    end
    if in(PSY.VariableReserve{PSY.ReserveDown}, typeof.(services))
        push!(duals, Symbol("requirement__VariableReserve_PowerSystems.ReserveDown"))
    end
    if in(PSY.VariableReserve{PSY.ReserveSymmetric}, typeof.(services))
        push!(duals, Symbol("requirement__VariableReserve_PowerSystems.ReserveSymmetric"))
    end

        if type == "UC"
            problem = PSI.OperationsProblem(
                                    PSIE.MILPDualProblem,
                                    template,
                                    sys;
                                    optimizer = solver,
                                    optimizer_log_print = false,
                                    balance_slack_variables = true,
                                    constraint_duals = duals,
                                    warm_start = true,
                                    services_slack_variables = true,
                                    )

        elseif type == "ED"
            if !(isempty(inertia_product))
                problem = PSI.OperationsProblem(
                                        PSI.EconomicDispatchProblem,
                                        template,
                                        sys;
                                        optimizer = solver,
                                        optimizer_log_print = false,
                                        balance_slack_variables = true,
                                        constraint_duals = duals,
                                        warm_start = true,
                                        services_slack_variables = true,
                                        )
            else
                problem = PSI.OperationsProblem(
                                        PSI.EconomicDispatchProblem,
                                        template,
                                        sys;
                                        optimizer = solver,
                                        optimizer_log_print = false,
                                        balance_slack_variables = true,
                                        constraint_duals = duals,
                                        warm_start = true,
                                        services_slack_variables = true,
                                        )
            end
        else
            error("Type should be either UC or ED")
        end

    return problem
end

"""
This function creates the Sequences for PSI Simulation.
"""
function create_sequence(problems::PSI.SimulationProblems, da_resolution::Int64, feedforward_dict)

    sequence = PSI.SimulationSequence(
        problems = problems,
        feedforward_chronologies = Dict(("UC" => "ED") => PSI.Synchronize(periods = 24 * 60 / da_resolution)),
        intervals = Dict(
            "UC" => (Dates.Hour(24 * 60 / da_resolution), PSI.Consecutive()),
            "ED" => (Dates.Hour(1 * 60 / da_resolution), PSI.Consecutive()),
        ),
        feedforward = feedforward_dict,
        ini_cond_chronology = PSI.InterProblemChronology(),
    )
    return sequence
end


"""
This function creates the PSI Simulation and post-processes the results.
"""
#TODO: Update needed
function create_simulation( sys_UC::PSY.System,
                            sys_ED::PSY.System,
                            simulation_dir::String,
                            reserve_penalty::String,
                            zones::Vector{String},
                            days::Int64,
                            da_resolution::Int64,
                            rt_resolution::Int64,
                            case_name::String,
                            solver::JuMP.MOI.OptimizerWithAttributes;
                            kwargs...)

    inertia_product = collect(PSY.get_components_by_name(PSY.Service, sys_ED, "Inertia"))

    template_uc = create_uc_template(inertia_product)
    template_ed = create_ed_template(inertia_product)

    uc_problem = create_problem(template_uc, sys_UC, "UC", solver, inertia_product)
    ed_problem = create_problem(template_ed, sys_ED, "ED", solver, inertia_product)

    if isempty(inertia_product)
        feedforward_dict = Dict(
            ("ED", :devices, Symbol("PowerSystems.ThermalStandard")) => PSI.SemiContinuousFF(
                binary_source_problem = PSI.ON,
                affected_variables = [PSI.ACTIVE_POWER],
            ),
            ("ED", :devices, Symbol("PowerSystemExtensions.ThermalCleanEnergy")) => PSI.SemiContinuousFF(
                binary_source_problem = PSI.ON,
                affected_variables = [PSI.ACTIVE_POWER],
            ),
            ("ED", :devices, Symbol("ThermalFastStartSIIP")) => PSI.SemiContinuousFF(
                binary_source_problem = PSI.ON,
                affected_variables = [PSI.ACTIVE_POWER],
            ),
            ("ED", :devices, Symbol("PowerSystems.HydroDispatch")) => PSIE.EnergyCommitmentFF(
                variable_source_problem = PSI.ACTIVE_POWER,
                affected_variables = [PSI.ACTIVE_POWER],
                affected_time_periods = 12,
            ),
            ("ED", :devices, Symbol("PowerSystems.HydroEnergyReservoir")) => PSIE.EnergyCommitmentFF(
                variable_source_problem = PSI.ACTIVE_POWER,
                affected_variables = [PSI.ACTIVE_POWER],
                affected_time_periods = 12,
            ),
            ("ED", :devices, Symbol("PowerSystems.RenewableDispatch")) => PSIE.EnergyCommitmentFF(
                variable_source_problem = PSI.ACTIVE_POWER,
                affected_variables = [PSI.ACTIVE_POWER],
                affected_time_periods = 12,
            ),
            ("ED", :devices, Symbol("PowerSystems.GenericBattery")) => PSIE.EnergyTargetFF(
                variable_source_problem = PSI.ENERGY,
                affected_variables = [PSI.ENERGY],
                target_period = 12,
                penalty_cost = 1e5,
            ),
        )
    else
        feedforward_dict = Dict(
            ("ED", :devices, Symbol("PowerSystems.ThermalStandard")) => PSIE.InertiaFF(
                binary_source_problem = PSI.ON,
                affected_variables = [PSI.ACTIVE_POWER],
                service = inertia_product[1]
            ),
            ("ED", :devices, Symbol("PowerSystemExtensions.ThermalCleanEnergy")) => PSIE.EmisFF(
                binary_source_problem = PSI.ON,
                variable_source_problem = PSI.ACTIVE_POWER,
                affected_variables = [PSI.ACTIVE_POWER],
                affected_time_periods = 12,
                service = inertia_product[1],
            ),
            ("ED", :devices, Symbol("ThermalFastStartSIIP")) => PSIE.InertiaFF(
                binary_source_problem = PSI.ON,
                affected_variables = [PSI.ACTIVE_POWER],
                service = inertia_product[1]
            ),
            ("ED", :devices, Symbol("PowerSystems.HydroEnergyReservoir")) => PSIE.EmisFF(
                binary_source_problem = PSI.ON,
                variable_source_problem = PSI.ACTIVE_POWER,
                affected_variables = [PSI.ACTIVE_POWER],
                affected_time_periods = 12,
                service = inertia_product[1],
            ),
            ("ED", :devices, Symbol("PowerSystems.HydroDispatch")) => PSIE.EmisFF(
                binary_source_problem = PSI.ON,
                variable_source_problem = PSI.ACTIVE_POWER,
                affected_variables = [PSI.ACTIVE_POWER],
                affected_time_periods = 12,
                service = inertia_product[1],
            ),
            ("ED", :devices, Symbol("PowerSystems.RenewableDispatch")) => PSIE.EnergyCommitmentFF(
                variable_source_problem = PSI.ACTIVE_POWER,
                affected_variables = [PSI.ACTIVE_POWER],
                affected_time_periods = 12,
            ),
            ("ED", :devices, Symbol("PowerSystems.GenericBattery")) => PSIE.EnergyTargetFF(
                variable_source_problem = PSI.ENERGY,
                affected_variables = [PSI.ENERGY],
                target_period = 12,
                penalty_cost = 1e5,
            ),
         )

    end


    rt_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"rt_products"], "; ")
    da_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"da_products"], "; ")

    default_balance_slack_cost = PSI.BALANCE_SLACK_COST
    default_service_slack_cost = PSI.SERVICES_SLACK_COST

    energy_mkt_data = read_data(joinpath(simulation_dir, "markets_data", "Energy.csv"))
    energy_voll_cost = AxisArrays.AxisArray(energy_mkt_data.price_cap * 1.0, zones)

    problems = PSI.SimulationProblems(
                                    UC = uc_problem,
                                    ED = ed_problem
                                      )

    sequence = create_sequence(problems, da_resolution, feedforward_dict)

    sim = PSI.Simulation(
                    name = "emis_$(case_name)",
                    steps = 360,
                    problems = problems,
                    sequence = sequence,
                    simulation_folder = ".",
                    )

    build_out = PSI.build!(sim; serialize = false)

    adjust_reserve_voll!(sys_UC, uc_problem, simulation_dir, reserve_penalty, zones, default_balance_slack_cost, default_service_slack_cost, energy_voll_cost)
    adjust_reserve_voll!(sys_ED, ed_problem, simulation_dir, reserve_penalty, zones, default_balance_slack_cost, default_service_slack_cost, energy_voll_cost)

    execute_out = PSI.execute!(sim; enable_progress_bar = true)

    sim_results = PSI.SimulationResults(sim);

    base_power = PSY.get_base_power(sys_UC)

    res_uc = PSI.get_problem_results(sim_results, "UC")
    res_ed = PSI.get_problem_results(sim_results, "ED")

    dual_values_ed = PSI.read_realized_duals(res_ed)
    dual_values_uc = PSI.read_realized_duals(res_uc)

    result_variables_ed = PSI.read_realized_variables(res_ed)
    result_variables_uc = PSI.read_realized_variables(res_uc)

    data_length_ed = DataFrames.nrow(dual_values_ed[:nodal_balance_active__Bus])
    data_length_uc = DataFrames.nrow(dual_values_uc[:nodal_balance_active__Bus])

    energy_price = AxisArrays.AxisArray(zeros(length(zones), 1, data_length_ed), zones, 1:1, 1:data_length_ed)
    energy_voll = AxisArrays.AxisArray(zeros(length(zones), 1, data_length_ed), zones, 1:1, 1:data_length_ed)

    reserve_price = Dict(s => zeros(1, data_length_ed) for s in String.(rt_products))
    reserve_voll = Dict(s => zeros(1, data_length_ed) for s in String.(rt_products))

    inertia_price = AxisArrays.AxisArray(zeros(1, data_length_ed), 1:1, 1:data_length_ed)
    inertia_voll = AxisArrays.AxisArray(zeros(1, data_length_ed), 1:1, 1:data_length_ed)

    only_da_products = String[]

    for product in da_products
        if !(product in rt_products)
            push!(only_da_products, product)
            reserve_price[String(product)] = zeros(1, data_length_uc)
            reserve_voll[String(product)] = zeros(1, data_length_uc)
        end
    end

    for zone in zones
            bus = find_zonal_bus(zone, sys_UC)
            zone_num = parse(Int64, last(zone, 1))
            if isnothing(bus)
                energy_price[zone, 1, :] = zeros(data_length_ed)
                energy_voll[zone, 1, :] = zeros(data_length_ed)
            else
                energy_price[zone, 1, :] =
                abs.(round.(dual_values_ed[:nodal_balance_active__Bus][:, PSY.get_name(bus)], digits = 5)) / base_power
                energy_voll[zone, 1, :] = abs.(round.(result_variables_ed[:γ⁺__P][:, string(PSY.get_number(bus))], digits = 5)) / base_power
                energy_voll[zone, 1, :] += abs.(round.(result_variables_ed[:γ⁻__P][:, string(PSY.get_number(bus))], digits = 5)) / base_power
            end
    end

    println(any(isnan, energy_price))
    replace!(energy_price, NaN => 0.0)
    scale_voll(energy_price, rt_resolution)

    println(any(isnan, energy_price))
    println(Statistics.mean(energy_price))


    for service in get_system_services(sys_ED)
        name = PSY.get_name(service)
        #println(name)
        if typeof(service) == PSY.VariableReserve{PSY.ReserveUp}
            reserve_price[name][1, :] = abs.(round.(dual_values_ed[Symbol("requirement__VariableReserve_PowerSystems.ReserveUp")][:, Symbol("$(name)")], digits = 5)) / base_power
            replace!(reserve_price[name], NaN => 0.0)
            scale_voll(reserve_price[name], rt_resolution)
            reserve_voll[name][1, :] = abs.(round.(result_variables_ed[Symbol("γ⁺__$(name)")][:, Symbol("γ⁺__$(name)")], digits = 5)) / base_power
            #println(reserve_price[name])
        elseif typeof(service) == PSY.ReserveDemandCurve{PSY.ReserveUp}
            reserve_price[name][1, :] = abs.(round.(dual_values_ed[Symbol("requirement__ReserveDemandCurve_PowerSystems.ReserveUp")][:, Symbol("$(name)")], digits = 5)) / base_power
            replace!(reserve_price[name], NaN => 0.0)
            scale_voll(reserve_price[name], rt_resolution)
            #println(reserve_price[name])
        elseif typeof(service) == PSY.VariableReserve{PSY.ReserveDown}
            reserve_price[name][1, :] = abs.(round.(dual_values_ed[Symbol("requirement__VariableReserve_PowerSystems.ReserveDown")][:, Symbol("$(name)")], digits = 5)) / base_power
            replace!(reserve_price[name], NaN => 0.0)
            scale_voll(reserve_price[name], rt_resolution)
            reserve_voll[name][1, :] = abs.(round.(result_variables_ed[Symbol("γ⁺__$(name)")][:, Symbol("γ⁺__$(name)")], digits = 5)) / base_power
            #println(reserve_price[name])
        elseif typeof(service) == PSY.VariableReserve{PSY.ReserveSymmetric}
            if name == "Inertia"
                inertia_price[1, :] = abs.(round.(dual_values_ed[Symbol("requirement__VariableReserve_PowerSystems.ReserveSymmetric")][:, Symbol("$(name)")], digits = 5)) / base_power
                replace!(inertia_price, NaN => 0.0)
                scale_voll(inertia_price, rt_resolution)
                inertia_voll[1, :] = abs.(round.(result_variables_ed[Symbol("γ⁺__$(name)")][:, Symbol("γ⁺__$(name)")], digits = 5)) / base_power
            end
            #println(Statistics.mean(inertia_price))
            #println(maximum(inertia_price))
        end
    end

    for service in get_system_services(sys_UC)
        name = PSY.get_name(service)
        if name in only_da_products
            if typeof(service) == PSY.VariableReserve{PSY.ReserveUp}
                reserve_price[name][1, :] = abs.(round.(dual_values_uc[Symbol("requirement__VariableReserve_PowerSystems.ReserveUp")][:, Symbol("$(name)")], digits = 5)) / base_power
                replace!(reserve_price[name], NaN => 0.0)
                reserve_voll[name][1, :] = abs.(round.(result_variables_uc[Symbol("γ⁺__$(name)")][:, Symbol("γ⁺__$(name)")], digits = 5)) / base_power
                scale_voll(reserve_price[name], da_resolution)
            elseif typeof(service) == PSY.ReserveDemandCurve{PSY.ReserveUp}
                reserve_price[name][1, :] = abs.(round.(dual_values_uc[Symbol("requirement__ReserveDemandCurve_PowerSystems.ReserveUp")][:, Symbol("$(name)")], digits = 5)) / base_power
                replace!(reserve_price[name], NaN => 0.0)
                scale_voll(reserve_price[name], da_resolution)
            elseif typeof(service) == PSY.VariableReserve{PSY.ReserveDown}
                reserve_price[name][1, :] = abs.(round.(dual_values_uc[Symbol("requirement__VariableReserve_PowerSystems.ReserveDown")][:, Symbol("$(name)")], digits = 5)) / base_power
                replace!(reserve_price[name], NaN => 0.0)
                scale_voll(reserve_price[name], da_resolution)
                reserve_voll[name][1, :] = abs.(round.(result_variables_uc[Symbol("γ⁺__$(name)")][:, Symbol("γ⁺__$(name)")], digits = 5)) / base_power
            elseif typeof(service) == PSY.VariableReserve{PSY.ReserveSymmetric}
                if name == "Inertia"
                    inertia_price[1, :] = abs.(round.(dual_values_uc[Symbol("requirement__VariableReserve_PowerSystems.ReserveSymmetric")][:, Symbol("$(name)")], digits = 5)) / base_power
                    replace!(inertia_price, NaN => 0.0)
                    scale_voll(inertia_price, da_resolution)
                end
                #println(inertia_price)
            end
        end
    end

    sys_techs = get_all_techs(sys_ED)

    tech_names = get_name.(sys_techs)
    capacity_factors = Dict([g => zeros(1, data_length_ed) for g in tech_names])
    start_up_costs = Dict([g => zeros(1, data_length_uc) for g in tech_names])
    shut_down_costs = Dict([g => zeros(1, data_length_uc) for g in tech_names])

    reserve_perc = Dict(g => Dict(s => zeros(1, data_length_ed) for s in String.(rt_products)) for g in tech_names)

    inertia_perc = Dict([g => zeros(1, data_length_ed) for g in tech_names])

    for g in tech_names
        for product in only_da_products
            reserve_perc[g][string(product)] = zeros(1, data_length_uc)
        end
    end

    for tech in sys_techs
        name = get_name(tech)
        capacity_factors[name][1, :] = get_realized_capacity_factors(tech, result_variables_ed, result_variables_uc)
        start_up_costs[name][1, :] = get_start_costs(tech, result_variables_uc, data_length_uc)
        shut_down_costs[name][1, :] = get_shut_costs(tech, result_variables_uc, data_length_uc)

        services_ED = PSY.get_services(tech)
        services_UC = PSY.get_services(PSY.get_components_by_name(PSY.Device, sys_UC, name)[1])

        for service in services_UC
            update_realized_reserve_perc!(tech,
                                         service,
                                         result_variables_ed,
                                         result_variables_uc,
                                         reserve_perc,
                                         inertia_perc,
                                         rt_products,
                                         only_da_products)

        end
    end

    for g in keys(inertia_perc)
        inertia_perc[g] = inertia_perc[g] / PSY.get_base_power(sys_UC)
    end
    return energy_price, reserve_price, inertia_price, capacity_factors, reserve_perc, inertia_perc, start_up_costs, shut_down_costs, energy_voll, reserve_voll, inertia_voll;
end
