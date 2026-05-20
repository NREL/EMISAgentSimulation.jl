# This file contains functions for calling, running and post-processing
# SIIP PSI Simulation for actual Energy and A/S Market Clearing

function adjust_reserve_voll!(sys::PSY.System,
                             problem::PSI.OperationModel,
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

    for zone in zones
        bus = find_zonal_bus(zone, sys)
        slack_coefficients = [PSI.SystemBalanceSlackUp, PSI.SystemBalanceSlackDown]
        for c in slack_coefficients
            slack_key = PSI.VariableKey{c, PSY.ACBus}("")
            index = findall(x -> x == PSY.get_number(bus), variables[slack_key].axes[1])[1]
            slack_variables = variables[slack_key].data[index, :]
            delta_cost = energy_voll_cost[zone] * base_power - default_balance_slack_cost
            for v in slack_variables
                PSI.add_to_objective_variant_expression!(optimization_container, v * delta_cost)
            end
        end
    end

    for s in services
        name = PSY.get_name(s)
        delta_cost = 0.0
        slack_variables = Symbol[]

        # we don't have slack for Primary and Synchronous reserves?
        if name == "Clean_Energy"
            delta_cost = default_service_slack_cost * 5
            #delta_cost = -default_service_slack_cost
            # slack_variables = variables[PSI.VariableKey{PSI.ReserveRequirementSlack, PSY.VariableReserve{PSY.ReserveUp}}("Clean_Energy")]
            # println(PSY.get_requirement(s))
        elseif name == "Reg_Up" || name == "Inertia"
            reserve_data = read_data(joinpath(simulation_dir, "markets_data", "$(reserve_penalty)_reserve_penalty", "$(name).csv"))
            slack_variables = variables[PSI.VariableKey{PSI.ReserveRequirementSlack, PSY.VariableReserve{PSY.ReserveUp}}("$name")]
            penalty_price = reserve_data[1, "price_cap"] * base_power
            delta_cost = penalty_price - default_service_slack_cost
        elseif name == "Reg_Down"
            reserve_data = read_data(joinpath(simulation_dir, "markets_data", "$(reserve_penalty)_reserve_penalty", "$(name).csv"))
            slack_variables = variables[PSI.VariableKey{PSI.ReserveRequirementSlack, PSY.VariableReserve{PSY.ReserveDown}}("$name")]
            penalty_price = reserve_data[1, "price_cap"] * base_power
            delta_cost = penalty_price - default_service_slack_cost
        end


        for v in slack_variables
            PSI.add_to_objective_variant_expression!(optimization_container, v * delta_cost)
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
                              results::Dict{String, DataFrames.DataFrame},
                              data_length_uc::Int64
                                        )
    start_ups = results["StartVariable__ThermalStandard"][:, Symbol(get_name(device))]
    start_up_costs = start_ups * PSY.get_start_up(PSY.get_operation_cost(device))
    return start_up_costs
end


function get_start_costs(device::ThermalFastStartSIIP,
                              results::Dict{String, DataFrames.DataFrame},
                              data_length_uc::Int64
                                        )
    start_ups = results["StartVariable__ThermalFastStartSIIP"][:, Symbol(get_name(device))]
    start_up_costs = start_ups * PSY.get_start_up(PSY.get_operation_cost(device))
    return start_up_costs
end

"""
This function returns total start-up costs for other generators from PSI Simulation.
"""
function get_start_costs(device::PSY.Device,
                              results::Dict{String, DataFrames.DataFrame},
                              data_length_uc::Int64
                                        )
    return zeros(data_length_uc)
end

"""
This function returns shut-down costs for Thermal generators from PSI Simulation.
"""
function get_shut_costs(device::PSY.ThermalStandard,
                       results::Dict{String, DataFrames.DataFrame},
                       data_length_uc::Int64
                                        )
    shut_downs = results["StopVariable__ThermalStandard"][:, Symbol(get_name(device))]
    shut_down_costs = shut_downs * PSY.get_shut_down(PSY.get_operation_cost(device))
    return shut_down_costs
end


function get_shut_costs(device::ThermalFastStartSIIP,
                       results::Dict{String, DataFrames.DataFrame},
                       data_length_uc::Int64
                                        )
    shut_downs = results["StopVariable__ThermalFastStartSIIP"][:, Symbol(get_name(device))]
    shut_down_costs = shut_downs * PSY.get_shut_down(PSY.get_operation_cost(device))
    return shut_down_costs
end
"""
This function returns total shut-down costs for other generators from PSI Simulation.
"""
function get_shut_costs(device::PSY.Device,
                         results::Dict{String, DataFrames.DataFrame},
                         data_length_uc::Int64
                                        )
    return zeros(data_length_uc)
end


"""
This function returns realized capacity factors for ThermalStandard generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.ThermalStandard,
                                        results::Dict{String, DataFrames.DataFrame},
                                        results_uc::Dict{String, DataFrames.DataFrame},
                                        base_power::Float64
                                        )
    energy_production = results["ActivePowerVariable__ThermalStandard"][:, Symbol(get_name(device))]
    capacity_factors = (energy_production / base_power) / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for ThermalCleanEnergy generators from PSI Simulation.
"""

"""
This function returns realized capacity factors for ThermalFastStartSIIP generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::ThermalFastStartSIIP,
                                        results::Dict{String, DataFrames.DataFrame},
                                        results_uc::Dict{String, DataFrames.DataFrame},
                                        base_power::Float64
                                        )
    energy_production = results["ActivePowerVariable__ThermalFastStartSIIP"][:, Symbol(get_name(device))]
    capacity_factors = (energy_production / base_power) / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for Renewable generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.RenewableDispatch,
                                        results::Dict{String, DataFrames.DataFrame},
                                        results_uc::Dict{String, DataFrames.DataFrame},
                                        base_power::Float64
                                        )
    energy_production = results["ActivePowerVariable__RenewableDispatch"][:, Symbol(get_name(device))]
    capacity_factors = (energy_production / base_power) / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for Hydropower generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.HydroDispatch,
                                        results::Dict{String, DataFrames.DataFrame},
                                        results_uc::Dict{String, DataFrames.DataFrame},
                                        base_power::Float64
                                        )
    energy_production = results["ActivePowerVariable__HydroDispatch"][:, Symbol(get_name(device))]
    capacity_factors = (energy_production / base_power) / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for Hydropower generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.HydroEnergyReservoir,
                                        results::Dict{String, DataFrames.DataFrame},
                                        results_uc::Dict{String, DataFrames.DataFrame},
                                        base_power::Float64
                                        )
    energy_production = results["ActivePowerVariable__HydroEnergyReservoir"][:, Symbol(get_name(device))]
    capacity_factors = (energy_production / base_power) / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for Hydropower generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.GenericBattery,
                                        results::Dict{String, DataFrames.DataFrame},
                                        results_uc::Dict{String, DataFrames.DataFrame},
                                        base_power::Float64
                                        )
    energy_production = results["ActivePowerOutVariable__GenericBattery"][:, Symbol(get_name(device))] - results["ActivePowerInVariable__GenericBattery"][:, Symbol(get_name(device))]
    capacity_factors = (energy_production / base_power) / get_device_size(device)
    generation = filter(x -> x > 0, capacity_factors)

    energy_production_uc = results_uc["ActivePowerOutVariable__GenericBattery"][:, Symbol(get_name(device))] - results_uc["ActivePowerInVariable__GenericBattery"][:, Symbol(get_name(device))]
    capacity_factors_uc = (energy_production_uc / base_power) / get_device_size(device)
    generation_uc = filter(x -> x > 0, capacity_factors_uc)
    return capacity_factors
end

"""
This function returns nothing if Service is not of ReserveUp or ReserveDown type.
"""
function update_realized_reserve_perc!(device::PSY.Device,
                                        service::S,
                                        results_ed::Dict{String, DataFrames.DataFrame},
                                        results_uc::Dict{String, DataFrames.DataFrame},
                                        results_md::Union{Nothing, Dict{String, DataFrames.DataFrame}},
                                        reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
                                        reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
                                        reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
                                        inertia_perc::Dict{String, Array{Float64, 2}},
                                        rt_products::Vector{SubString{String}},
                                        da_products::Vector{SubString{String}},
                                        md_products::Vector{SubString{String}},
                                        base_power::Float64,
                                        md_market_bool::Bool,) where S <: PSY.Service
    return
end

"""
This function returns realized reserve up provision percentages from PSI Simulation.
"""
function update_realized_reserve_perc!(device::PSY.Device,
                                        service::PSY.VariableReserve{PSY.ReserveUp},
                                        results_ed::Dict{String, DataFrames.DataFrame},
                                        results_uc::Dict{String, DataFrames.DataFrame},
                                        results_md::Union{Nothing, Dict{String, DataFrames.DataFrame}},
                                        reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
                                        reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
                                        reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
                                        inertia_perc::Dict{String, Array{Float64, 2}},
                                        rt_products::Vector{SubString{String}},
                                        da_products::Vector{SubString{String}},
                                        md_products::Vector{SubString{String}},
                                        base_power::Float64,
                                        md_market_bool::Bool,)

    service_name = PSY.get_name(service)

    if service_name == "Inertia"
        inertia_provision = results_ed["ActivePowerReserveVariable__VariableReserve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
        inertia_perc_value = inertia_provision / get_device_size(device) / base_power
        inertia_perc[get_name(device)][1, :] = inertia_perc_value
    else
        # if service_name in rt_products
        #     reserve_provision = results_ed["ActivePowerReserveVariable__VariableReserve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
        #     reserve_perc_value = reserve_provision / get_device_size(device) / base_power
        #     reserve_perc[get_name(device)][service_name][1, :] = reserve_perc_value

        # elseif service_name in only_da_products
        #     reserve_provision = results_uc["ActivePowerReserveVariable__VariableReserve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
        #     reserve_perc_value = reserve_provision / get_device_size(device) / base_power
        #     reserve_perc[get_name(device)][service_name][1, :] = reserve_perc_value

        # end
        reserve_provision_ed = results_ed["ActivePowerReserveVariable__VariableReserve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
        reserve_perc_value_ed = reserve_provision_ed / get_device_size(device) / base_power
        reserve_perc_ed[get_name(device)][service_name][1, :] = reserve_perc_value_ed

        reserve_provision_uc = results_uc["ActivePowerReserveVariable__VariableReserve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
        reserve_perc_value_uc = reserve_provision_uc / get_device_size(device) / base_power
        reserve_perc_uc[get_name(device)][service_name][1, :] = reserve_perc_value_uc
    end

    if md_market_bool == true
        reserve_provision_md = results_md["ActivePowerReserveVariable__VariableReserve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
        reserve_perc_value_md = reserve_provision_md / get_device_size(device) / base_power
        reserve_perc_md[get_name(device)][service_name][1, :] = reserve_perc_value_md
    end

    return
end

"""
This function returns realized ordc provision percentages from PSI Simulation.
"""
function update_realized_reserve_perc!(device::PSY.Device,
                                        service::PSY.ReserveDemandCurve{PSY.ReserveUp},
                                        results_ed::Dict{String, DataFrames.DataFrame},
                                        results_uc::Dict{String, DataFrames.DataFrame},
                                        results_md::Union{Nothing, Dict{String, DataFrames.DataFrame}},
                                        reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
                                        reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
                                        reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
                                        inertia_perc::Dict{String, Array{Float64, 2}},
                                        rt_products::Vector{SubString{String}},
                                        da_products::Vector{SubString{String}},
                                        md_products::Vector{SubString{String}},
                                        base_power::Float64,
                                        md_market_bool::Bool,)

    service_name = PSY.get_name(service)

    reserve_provision_ed = results_ed["ActivePowerReserveVariable__ReserveDemandCurve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
    reserve_perc_value_ed = reserve_provision_ed / get_device_size(device) / base_power
    reserve_perc_ed[get_name(device)][service_name][1, :] = reserve_perc_value_ed

    reserve_provision_uc = results_uc["ActivePowerReserveVariable__ReserveDemandCurve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
    reserve_perc_value_uc = reserve_provision_uc / get_device_size(device) / base_power
    reserve_perc_uc[get_name(device)][service_name][1, :] = reserve_perc_value_uc

    if md_market_bool == true
        reserve_provision_md = results_md["ActivePowerReserveVariable__ReserveDemandCurve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
        reserve_perc_value_md = reserve_provision_md / get_device_size(device) / base_power
        reserve_perc_md[get_name(device)][service_name][1, :] = reserve_perc_value_md
    end

    # if service_name in rt_products
    #     reserve_provision = results_ed["ActivePowerReserveVariable__ReserveDemandCurve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
    #     reserve_perc_value = reserve_provision / get_device_size(device) / base_power
    #     reserve_perc[get_name(device)][service_name][1, :] = reserve_perc_value

    # elseif service_name in only_da_products
    #     reserve_provision = results_uc["ActivePowerReserveVariable__ReserveDemandCurve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
    #     reserve_perc_value = reserve_provision / get_device_size(device) / base_power
    #     reserve_perc[get_name(device)][service_name][1, :] = reserve_perc_value

    # end
    return
end

"""
This function returns realized static reserve provision percentages from PSI Simulation.
"""
function update_realized_reserve_perc!(device::PSY.Device,
                                        service::PSY.StaticReserve{PSY.ReserveUp},
                                        results_ed::Dict{String, DataFrames.DataFrame},
                                        results_uc::Dict{String, DataFrames.DataFrame},
                                        results_md::Union{Nothing, Dict{String, DataFrames.DataFrame}},
                                        reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
                                        reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
                                        reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
                                        inertia_perc::Dict{String, Array{Float64, 2}},
                                        rt_products::Vector{SubString{String}},
                                        da_products::Vector{SubString{String}},
                                        md_products::Vector{SubString{String}},
                                        base_power::Float64,
                                        md_market_bool::Bool,)

    service_name = PSY.get_name(service)

    reserve_provision_ed = results_ed["ActivePowerReserveVariable__StaticReserve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
    reserve_perc_value_ed = reserve_provision_ed / get_device_size(device) / base_power
    reserve_perc_ed[get_name(device)][service_name][1, :] = reserve_perc_value_ed

    reserve_provision_uc = results_uc["ActivePowerReserveVariable__StaticReserve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
    reserve_perc_value_uc = reserve_provision_uc / get_device_size(device) / base_power
    reserve_perc_uc[get_name(device)][service_name][1, :] = reserve_perc_value_uc

    if md_market_bool == true
        reserve_provision_md = results_md["ActivePowerReserveVariable__StaticReserve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
        reserve_perc_value_md = reserve_provision_md / get_device_size(device) / base_power
        reserve_perc_md[get_name(device)][service_name][1, :] = reserve_perc_value_md
    end

    return
end

"""
This function returns realized reserve down provision percentages from PSI Simulation.
"""
function update_realized_reserve_perc!(device::PSY.Device,
                                        service::PSY.VariableReserve{PSY.ReserveDown},
                                        results_ed::Dict{String, DataFrames.DataFrame},
                                        results_uc::Dict{String, DataFrames.DataFrame},
                                        results_md::Union{Nothing, Dict{String, DataFrames.DataFrame}},
                                        reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
                                        reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
                                        reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
                                        inertia_perc::Dict{String, Array{Float64, 2}},
                                        rt_products::Vector{SubString{String}},
                                        da_products::Vector{SubString{String}},
                                        md_products::Vector{SubString{String}},
                                        base_power::Float64,
                                        md_market_bool::Bool,)

    service_name = PSY.get_name(service)

    reserve_provision_ed = results_ed["ActivePowerReserveVariable__VariableReserve__ReserveDown__$(service_name)"][:, Symbol(get_name(device))]
    reserve_perc_value_ed = reserve_provision_ed / get_device_size(device) / base_power
    reserve_perc_ed[get_name(device)][service_name][1, :] = reserve_perc_value_ed

    reserve_provision_uc = results_uc["ActivePowerReserveVariable__VariableReserve__ReserveDown__$(service_name)"][:, Symbol(get_name(device))]
    reserve_perc_value_uc = reserve_provision_uc / get_device_size(device) / base_power
    reserve_perc_uc[get_name(device)][service_name][1, :] = reserve_perc_value_uc

    if md_market_bool == true
        reserve_provision_md = results_md["ActivePowerReserveVariable__VariableReserve__ReserveDown__$(service_name)"][:, Symbol(get_name(device))]
        reserve_perc_value_md = reserve_provision_md / get_device_size(device) / base_power
        reserve_perc_md[get_name(device)][service_name][1, :] = reserve_perc_value_md
    end

    # if service_name in rt_products
    #     reserve_provision = results_ed["ActivePowerReserveVariable__VariableReserve__ReserveDown__$(service_name)"][:, Symbol(get_name(device))]
    #     reserve_perc_value = reserve_provision / get_device_size(device) / base_power
    #     reserve_perc[get_name(device)][service_name][1, :] = reserve_perc_value

    # elseif service_name in only_da_products
    #     reserve_provision = results_uc["ActivePowerReserveVariable__VariableReserve__ReserveDown__$(service_name)"][:, Symbol(get_name(device))]
    #     reserve_perc_value = reserve_provision / get_device_size(device) / base_power
    #     reserve_perc[get_name(device)][service_name][1, :] = reserve_perc_value

    # end
    return
end


"""
This function returns realized ordc provision percentages from PSI Simulation.
"""
function update_realized_reserve_perc!(device::PSY.GenericBattery,
    service::PSY.ReserveDemandCurve{PSY.ReserveUp},
    results_ed::Dict{String, DataFrames.DataFrame},
    results_uc::Dict{String, DataFrames.DataFrame},
    results_md::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
    reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
    reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
    inertia_perc::Dict{String, Array{Float64, 2}},
    rt_products::Vector{SubString{String}},
    da_products::Vector{SubString{String}},
    md_products::Vector{SubString{String}},
    base_power::Float64,
    md_market_bool::Bool,)

    service_name = PSY.get_name(service)

    reserve_provision_ed = (results_ed["AncillaryServiceVariableDischarge__GenericBattery__ReserveDemandCurve{ReserveUp}_$(service_name)"][:, Symbol(get_name(device))] .+
        results_ed["AncillaryServiceVariableCharge__GenericBattery__ReserveDemandCurve{ReserveUp}_$(service_name)"][:, Symbol(get_name(device))])
    reserve_perc_value_ed = reserve_provision_ed / get_device_size(device) / base_power
    reserve_perc_ed[get_name(device)][service_name][1, :] = reserve_perc_value_ed

    reserve_provision_uc = (results_uc["AncillaryServiceVariableDischarge__GenericBattery__ReserveDemandCurve{ReserveUp}_$(service_name)"][:, Symbol(get_name(device))] .+
        results_uc["AncillaryServiceVariableCharge__GenericBattery__ReserveDemandCurve{ReserveUp}_$(service_name)"][:, Symbol(get_name(device))])
    reserve_perc_value_uc = reserve_provision_uc / get_device_size(device) / base_power
    reserve_perc_uc[get_name(device)][service_name][1, :] = reserve_perc_value_uc

    if md_market_bool == true
        reserve_provision_md = (results_md["AncillaryServiceVariableDischarge__GenericBattery__ReserveDemandCurve{ReserveUp}_$(service_name)"][:, Symbol(get_name(device))] .+
            results_md["AncillaryServiceVariableCharge__GenericBattery__ReserveDemandCurve{ReserveUp}_$(service_name)"][:, Symbol(get_name(device))])
        reserve_perc_value_md = reserve_provision_md / get_device_size(device) / base_power
        reserve_perc_md[get_name(device)][service_name][1, :] = reserve_perc_value_md
    end

    # if service_name in rt_products
    #     reserve_provision = results_ed["ActivePowerReserveVariable__ReserveDemandCurve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
    #     reserve_perc_value = reserve_provision / get_device_size(device) / base_power
    #     reserve_perc[get_name(device)][service_name][1, :] = reserve_perc_value

    # elseif service_name in only_da_products
    #     reserve_provision = results_uc["ActivePowerReserveVariable__ReserveDemandCurve__ReserveUp__$(service_name)"][:, Symbol(get_name(device))]
    #     reserve_perc_value = reserve_provision / get_device_size(device) / base_power
    #     reserve_perc[get_name(device)][service_name][1, :] = reserve_perc_value

    # end
    return
end

"""
This function creates the Unit Commitment template for PSI Simulation.
"""
#TODO: Update needed

function create_md_template(inertia_product, ordc_curved::Bool)

    if !(isempty(inertia_product))

        template = PSI.ProblemTemplate(
            PSI.NetworkModel(
                PSI.PM.NFAPowerModel,
                duals = [PSI.NodalBalanceActiveConstraint],
                use_slacks = true,
            ),
        )
        # PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalStandardUnitCommitment)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalStandardUnitCommitment)
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicUnitCommitment)
        PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalBasicUnitCommitment)
        # PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicDispatch)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalBasicDispatch)
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalStandardUCOutages)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, RPSI.ThermalStandardUCOutages)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSI.RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableFix, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.PowerLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroEnergyReservoir, HSI.HydroCommitmentRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroCommitmentRunOfRiver) # TODO: check which hydro device we have
        # PSI.set_device_model!(template, PSY.HydroEnergyReservoir, HSI.HydroDispatchRunOfRiver)
        # PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroDispatchRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(template, PSY.GenericBattery, SSI.StorageDispatchWithReserves)
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TwoTerminalHVDCLine, PSI.HVDCTwoTerminalLossless)
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Reg_Up",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveDown},
                PSI.RangeReserve,
                "Reg_Down",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.VariableReserve{PSY.ReserveUp},
        #         EMISEx.InertiaReserve,
        #         "Inertia",
        #         use_slacks=true,
        #         duals = [PSI.RequirementConstraint],
        #     )
        # )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.VariableReserve{PSY.ReserveUp},
        #         EMISEx.CleanEnergyReserve,
        #         "Clean_Energy",
        #         use_slacks=true,
        #         duals = [PSI.RequirementConstraint],
        #     )
        # )
    else
        template = PSI.ProblemTemplate(
            PSI.NetworkModel(
                PSI.PM.NFAPowerModel,
                duals = [PSI.NodalBalanceActiveConstraint],
                use_slacks = true,
            ),
        )
        # PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalStandardUnitCommitment)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalStandardUnitCommitment)
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicUnitCommitment)
        PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalBasicUnitCommitment)
        # PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicDispatch)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalBasicDispatch)
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalStandardUCOutages)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, RPSI.ThermalStandardUCOutages)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSI.RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableFix, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.PowerLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroEnergyReservoir, HSI.HydroCommitmentRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroCommitmentRunOfRiver) # TODO: check which hydro device we have
        # PSI.set_device_model!(template, PSY.HydroEnergyReservoir, HSI.HydroDispatchRunOfRiver)
        # PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroDispatchRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(template, PSY.GenericBattery, SSI.StorageDispatchWithReserves)
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TwoTerminalHVDCLine, PSI.HVDCTwoTerminalLossless)
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Reg_Up",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveDown},
                PSI.RangeReserve,
                "Reg_Down",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.VariableReserve{PSY.ReserveUp},
        #         EMISEx.CleanEnergyReserve,
        #         "Clean_Energy",
        #         use_slacks=true,
        #         duals = [PSI.RequirementConstraint],
        #     )
        # )
    end

    if ordc_curved
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.ReserveDemandCurve{PSY.ReserveUp},
                PSI.StepwiseCostReserve,
                "Synchronous",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.ReserveDemandCurve{PSY.ReserveUp},
                PSI.StepwiseCostReserve,
                "Primary",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
    else
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.StaticReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Synchronous",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.StaticReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Primary",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
    end

    return template
end

function create_uc_template(inertia_product, ordc_curved::Bool)

    if !(isempty(inertia_product))

        template = PSI.ProblemTemplate(
            PSI.NetworkModel(
                PSI.PM.NFAPowerModel,
                duals = [PSI.NodalBalanceActiveConstraint],
                use_slacks = true,
            ),
        )
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicUnitCommitment)
        PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalBasicUnitCommitment)
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalStandardUCOutages)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, RPSI.ThermalStandardUCOutages)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSI.RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableFix, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.PowerLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroEnergyReservoir, HSI.HydroCommitmentRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroCommitmentRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(template, PSY.GenericBattery, SSI.StorageDispatchWithReserves)
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TwoTerminalHVDCLine, PSI.HVDCTwoTerminalLossless)
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Reg_Up",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveDown},
                PSI.RangeReserve,
                "Reg_Down",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.VariableReserve{PSY.ReserveUp},
        #         EMISEx.InertiaReserve,
        #         "Inertia",
        #         use_slacks=true,
        #         duals = [PSI.RequirementConstraint],
        #     )
        # )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.VariableReserve{PSY.ReserveUp},
        #         EMISEx.CleanEnergyReserve,
        #         "Clean_Energy",
        #         use_slacks=true,
        #         duals = [PSI.RequirementConstraint],
        #     )
        # )
    else
        template = PSI.ProblemTemplate(
            PSI.NetworkModel(
                PSI.PM.NFAPowerModel,
                duals = [PSI.NodalBalanceActiveConstraint],
                use_slacks = true,
            ),
        )
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicUnitCommitment)
        PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalBasicUnitCommitment)
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalStandardUCOutages)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, RPSI.ThermalStandardUCOutages)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSI.RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableFix, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.PowerLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroEnergyReservoir, HSI.HydroCommitmentRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroCommitmentRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(template, PSY.GenericBattery, SSI.StorageDispatchWithReserves)
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TwoTerminalHVDCLine, PSI.HVDCTwoTerminalLossless)
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Reg_Up",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveDown},
                PSI.RangeReserve,
                "Reg_Down",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.VariableReserve{PSY.ReserveUp},
        #         EMISEx.CleanEnergyReserve,
        #         "Clean_Energy",
        #         use_slacks=true,
        #         duals = [PSI.RequirementConstraint],
        #     )
        # )
    end

    if ordc_curved
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.ReserveDemandCurve{PSY.ReserveUp},
                PSI.StepwiseCostReserve,
                "Synchronous",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.ReserveDemandCurve{PSY.ReserveUp},
                PSI.StepwiseCostReserve,
                "Primary",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
    else
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.StaticReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Synchronous",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.StaticReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Primary",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
    end

    return template
end

"""
This function creates the Economic Dispatch template for PSI Simulation.
"""
#TODO: Update needed
function create_ed_template(inertia_product, ordc_curved::Bool)

    if !(isempty(inertia_product))
        template = PSI.ProblemTemplate(
            PSI.NetworkModel(
                PSI.PM.NFAPowerModel,
                duals = [PSI.NodalBalanceActiveConstraint],
                use_slacks = true,
            ),
        )
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicDispatch)
        PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalBasicUnitCommitment)
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalDispatchOutages)
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalRampLimitedOutages)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, RPSI.ThermalStandardUCOutages)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSI.RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableFix, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.PowerLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroEnergyReservoir, HSI.HydroDispatchRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroDispatchRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(template, PSY.GenericBattery, SSI.StorageDispatchWithReserves)
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TwoTerminalHVDCLine, PSI.HVDCTwoTerminalLossless)
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Reg_Up",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveDown},
                PSI.RangeReserve,
                "Reg_Down",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.VariableReserve{PSY.ReserveUp},
        #         EMISEx.InertiaReserve,
        #         "Inertia",
        #         use_slacks=true,
        #         duals = [PSI.RequirementConstraint],
        #     )
        # )
    else
        template = PSI.ProblemTemplate(
            PSI.NetworkModel(
                PSI.PM.NFAPowerModel,
                duals = [PSI.NodalBalanceActiveConstraint],
                use_slacks = true,
            ),
        )
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicDispatch)
        PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalBasicUnitCommitment)
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalDispatchOutages)
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalRampLimitedOutages)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, RPSI.ThermalStandardUCOutages)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSI.RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableFix, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.PowerLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroEnergyReservoir, HSI.HydroDispatchRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroDispatchRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(template, PSY.GenericBattery, SSI.StorageDispatchWithReserves)
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TwoTerminalHVDCLine, PSI.HVDCTwoTerminalLossless)
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Reg_Up",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveDown},
                PSI.RangeReserve,
                "Reg_Down",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
    end

    if ordc_curved
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.ReserveDemandCurve{PSY.ReserveUp},
                PSI.StepwiseCostReserve,
                "Synchronous",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.ReserveDemandCurve{PSY.ReserveUp},
                PSI.StepwiseCostReserve,
                "Primary",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
    else
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.StaticReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Synchronous",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.StaticReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Primary",
                use_slacks=true,
                duals = [PSI.RequirementConstraint],
            )
        )
    end

    return template
end

"""
This function creates the Problem for PSI Simulation.
"""

function create_problem(template::PSI.ProblemTemplate, sys::PSY.System, type::String, solver::JuMP.MOI.OptimizerWithAttributes, inertia_product)

        if type == "UC"
            problem = PSI.DecisionModel(
                                    template,
                                    sys;
                                    optimizer = solver,
                                    name = "UC",
                                    optimizer_solve_log_print = false,
                                    warm_start = true,
                                    calculate_conflict = true,
                                    store_variable_names = true,
                                    export_pwl_vars=true,
                                    initialize_model = false,
                                    )

        elseif type == "ED"
            if !(isempty(inertia_product))
                problem = PSI.DecisionModel(
                                        template,
                                        sys;
                                        optimizer = solver,
                                        name = "ED",
                                        optimizer_solve_log_print = false,
                                        warm_start = true,
                                        # horizon =2,
                                        calculate_conflict = true,
                                        store_variable_names = true,
                                        export_pwl_vars=true,
                                        initialize_model = false,
                                        )
            else
                problem = PSI.DecisionModel(
                                        template,
                                        sys;
                                        optimizer = solver,
                                        name = "ED",
                                        optimizer_solve_log_print = false,
                                        warm_start = true,
                                        # horizon =2,
                                        calculate_conflict = true,
                                        store_variable_names = true,
                                        export_pwl_vars=true,
                                        initialize_model = false,
                                        )
            end
        elseif type == "MD"
            problem = PSI.DecisionModel(
                                    template,
                                    sys;
                                    optimizer = solver,
                                    name = "MD",
                                    optimizer_solve_log_print = false,
                                    warm_start = true,
                                    calculate_conflict = true,
                                    store_variable_names = true,
                                    export_pwl_vars=true,
                                    initialize_model = false,
                                    )
        else
            error("Type should be either MD, UC or ED")
        end

    return problem
end

"""
This function creates the Sequences for PSI Simulation.
"""
function create_sequence(problems::PSI.SimulationModels, feedforward_dict)

    sequence = PSI.SimulationSequence(
        models = problems,
        feedforwards = feedforward_dict,
        ini_cond_chronology = PSI.InterProblemChronology(),
    )
    return sequence
end


PSY.get_max_output_fraction(value::PSY.ReserveDemandCurve{PSY.ReserveUp}) = 1.0

"""
This function creates the PSI Simulation and post-processes the results.
"""
#TODO: Update needed
function create_simulation( sys_MD::PSY.System,
                            sys_UC::PSY.System,
                            sys_ED::PSY.System,
                            simulation_dir::String,
                            reserve_penalty::String,
                            zones::Vector{String},
                            days::Int64,
                            da_resolution::Int64,
                            rt_resolution::Int64,
                            case_name::String,
                            solver::JuMP.MOI.OptimizerWithAttributes,
                            current_siip_sim,
                            md_market_bool::Bool,
                            ordc_curved::Bool,
                            siip_system;
                            kwargs...)

    # to_json(sys_MD, "/kfs2/projects/gmlcmarkets/Phase2_EMIS_Analysis/Feb2024_ERCOT_2011_MARKET_Test_NGUO_LDES/HPC_Analysis_Runs/storage_ff_debug/modified_test_sys/MD_sys.json", force=true)
    # to_json(sys_UC, "/kfs2/projects/gmlcmarkets/Phase2_EMIS_Analysis/Feb2024_ERCOT_2011_MARKET_Test_NGUO_LDES/HPC_Analysis_Runs/storage_ff_debug/modified_test_sys/UC_sys.json", force=true)
    # to_json(sys_ED, "/kfs2/projects/gmlcmarkets/Phase2_EMIS_Analysis/Feb2024_ERCOT_2011_MARKET_Test_NGUO_LDES/HPC_Analysis_Runs/storage_ff_debug/modified_test_sys/ED_sys.json", force=true)

    # hacky way to incorporate reserve voll
    base_power = PSY.get_base_power(sys_UC)
    reserve_data = read_data(joinpath(simulation_dir, "markets_data", "$(reserve_penalty)_reserve_penalty", "Reg_Up.csv"))
    penalty_price = reserve_data[1, "price_cap"] * base_power
    PSI.SERVICES_SLACK_COST = penalty_price

    inertia_product = collect(PSY.get_components_by_name(PSY.Service, sys_ED, "Inertia"))

    template_uc = create_uc_template(inertia_product, ordc_curved)
    template_ed = create_ed_template(inertia_product, ordc_curved)
    uc_problem = create_problem(template_uc, sys_UC, "UC", solver, inertia_product)
    ed_problem = create_problem(template_ed, sys_ED, "ED", solver, inertia_product)

    if md_market_bool == true

        template_md = create_md_template(inertia_product, ordc_curved)
        md_problem = create_problem(template_md, sys_MD, "MD", solver, inertia_product)

        if isempty(inertia_product)
            feedforward_dict = Dict(
                # "UC" => [
                #     SSI.EnergyLimitFeedforward(
                #         component_type = PSY.GenericBattery,
                #         source = PSI.ActivePowerOutVariable,
                #         affected_values = [PSI.ActivePowerOutVariable],
                #         number_of_periods = 24,
                #     ),
                # ],
                "UC" => [
                    SSI.EnergyTargetFeedforward(
                        component_type = PSY.GenericBattery,
                        source = PSI.EnergyVariable,
                        affected_values = [PSI.EnergyVariable],
                        target_period = 24,
                        penalty_cost = 5000.0,
                    ),
                ],
                "ED" => [
                    SSI.EnergyTargetFeedforward(
                        component_type = PSY.GenericBattery,
                        source = PSI.EnergyVariable,
                        affected_values = [PSI.EnergyVariable],
                        target_period = 2,
                        penalty_cost = 5000.0,
                    ),
                    PSI.SemiContinuousFeedforward(
                        component_type = PSY.ThermalStandard,
                        source = PSI.OnVariable,
                        affected_values = [PSI.ActivePowerVariable],
                    ),
                ]
                # "ED" => [
                #     SSI.EnergyLimitFeedforward(
                #         component_type = PSY.GenericBattery,
                #         source = PSI.ActivePowerOutVariable,
                #         affected_values = [PSI.ActivePowerOutVariable],
                #         number_of_periods = 1,
                #     ),
                # ],
                    # RPSI.SemiContinuousOutageFeedforward(
                    #     component_type = PSY.ThermalStandard,
                    #     source = PSI.OnVariable,
                    #     affected_values = [PSI.ActivePowerVariable],
                    # ),
            )
        else
            feedforward_dict = Dict(
                # "UC" => [
                #     SSI.EnergyLimitFeedforward(
                #         component_type = PSY.GenericBattery,
                #         source = PSI.ActivePowerOutVariable,
                #         affected_values = [PSI.ActivePowerOutVariable],
                #         number_of_periods = 24,
                #     ),
                # ],
                "UC" => [
                    SSI.EnergyTargetFeedforward(
                        component_type = PSY.GenericBattery,
                        source = PSI.EnergyVariable,
                        affected_values = [PSI.EnergyVariable],
                        target_period = 24,
                        penalty_cost = 5000.0,
                    ),
                ],
                "ED" => [
                    SSI.EnergyTargetFeedforward(
                        component_type = PSY.GenericBattery,
                        source = PSI.EnergyVariable,
                        affected_values = [PSI.EnergyVariable],
                        target_period = 2,
                        penalty_cost = 5000.0,
                    ),
                    PSI.SemiContinuousFeedforward(
                        component_type = PSY.ThermalStandard,
                        source = PSI.OnVariable,
                        affected_values = [PSI.ActivePowerVariable],
                    ),
                ]
                # "ED" => [
                #     SSI.EnergyLimitFeedforward(
                #         component_type = PSY.GenericBattery,
                #         source = PSI.ActivePowerOutVariable,
                #         affected_values = [PSI.ActivePowerOutVariable],
                #         number_of_periods = 1,
                #     ),
                # ],
                    # RPSI.SemiContinuousOutageFeedforward(
                    #     component_type = PSY.ThermalStandard,
                    #     source = PSI.OnVariable,
                    #     affected_values = [PSI.ActivePowerVariable],
                    # ),
             )
        end

        models = PSI.SimulationModels(
            decision_models = [
                md_problem,
                uc_problem,
                ed_problem
            ]
        );

    else

        if isempty(inertia_product)
            feedforward_dict = Dict(
                # "UC" => [
                #     SSI.EnergyLimitFeedforward(
                #         component_type = PSY.GenericBattery,
                #         source = PSI.EnergyVariable,
                #         affected_values = [PSI.EnergyVariable],
                #         number_of_periods = 36,
                #     ),
                # ],
                "ED" => [
                    SSI.EnergyTargetFeedforward(
                        component_type = PSY.GenericBattery,
                        source = PSI.EnergyVariable,
                        affected_values = [PSI.EnergyVariable],
                        target_period = 2,
                        penalty_cost = 5000.0,
                    ),
                    PSI.SemiContinuousFeedforward(
                        component_type = PSY.ThermalStandard,
                        source = PSI.OnVariable,
                        affected_values = [PSI.ActivePowerVariable],
                    ),
                    # RPSI.SemiContinuousOutageFeedforward(
                    #     component_type = PSY.ThermalStandard,
                    #     source = PSI.OnVariable,
                    #     affected_values = [PSI.ActivePowerVariable],
                    # ),
                ],
                # "ED" => [
                #     SSI.EnergyLimitFeedforward(
                #         component_type = PSY.GenericBattery,
                #         source = PSI.ActivePowerOutVariable,
                #         affected_values = [PSI.ActivePowerOutVariable],
                #         number_of_periods = 1,
                #     ),
                # ],
            )
        else
            feedforward_dict = Dict(
                # "UC" => [
                #     SSI.EnergyLimitFeedforward(
                #         component_type = PSY.GenericBattery,
                #         source = PSI.EnergyVariable,
                #         affected_values = [PSI.EnergyVariable],
                #         number_of_periods = 36,
                #     ),
                # ],
                "ED" => [
                    SSI.EnergyTargetFeedforward(
                        component_type = PSY.GenericBattery,
                        source = PSI.EnergyVariable,
                        affected_values = [PSI.EnergyVariable],
                        target_period = 2,
                        penalty_cost = 5000.0,
                    ),
                    PSI.SemiContinuousFeedforward(
                        component_type = PSY.ThermalStandard,
                        source = PSI.OnVariable,
                        affected_values = [PSI.ActivePowerVariable],
                    ),
                    # RPSI.SemiContinuousOutageFeedforward(
                    #     component_type = PSY.ThermalStandard,
                    #     source = PSI.OnVariable,
                    #     affected_values = [PSI.ActivePowerVariable],
                    # ),
                ],
             )
        end

        models = PSI.SimulationModels(
            decision_models = [
                uc_problem,
                ed_problem
            ]
        );

    end

    # TODO: need to define reserve products for MD
    md_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"da_products"], "; ")

    rt_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"rt_products"], "; ")
    da_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"da_products"], "; ")
    
    default_balance_slack_cost = PSI.BALANCE_SLACK_COST
    default_service_slack_cost = PSI.SERVICES_SLACK_COST

    energy_mkt_data = read_data(joinpath(simulation_dir, "markets_data", "Energy.csv"))
    energy_voll_cost = AxisArrays.AxisArray(energy_mkt_data.price_cap * 1.0, zones)

    sequence = create_sequence(models, feedforward_dict);

    if md_market_bool == true
        sim = PSI.Simulation(
                        name = "emis_$(case_name)",
                        steps = sys_MD.data.time_series_params.forecast_params.count, # sys_MD.data.time_series_params.forecast_params.count
                        models = models,
                        sequence = sequence,
                        simulation_folder = ".",
                        # initial_time = Dates.DateTime("2018-02-28T00:00:00")
                        );
    else
        sim = PSI.Simulation(
                        name = "emis_$(case_name)",
                        steps = 365, # sys_MD.data.time_series_params.forecast_params.count
                        models = models,
                        sequence = sequence,
                        simulation_folder = ".",
                        # initial_time = Dates.DateTime("2018-02-28T00:00:00")
                        );
    end

    # current_siip_sim[1] = sim

    build_out = PSI.build!(sim; serialize = false)

    #TODO: not sure if adjust_reserve_voll! function is working.
    # in particular, there does not seem to be additional terms added to the objective function.
    # adjust_reserve_voll!(sys_UC, uc_problem, simulation_dir, reserve_penalty, zones, default_balance_slack_cost, default_service_slack_cost, energy_voll_cost)
    # adjust_reserve_voll!(sys_ED, ed_problem, simulation_dir, reserve_penalty, zones, default_balance_slack_cost, default_service_slack_cost, energy_voll_cost)

    # if md_market_bool == true
    #     adjust_reserve_voll!(sys_MD, md_problem, simulation_dir, reserve_penalty, zones, default_balance_slack_cost, default_service_slack_cost, energy_voll_cost)
    # end

    # current_siip_sim[1] = sim
    # push!(siip_system, sys_MD)
    # push!(siip_system, sys_UC)
    # push!(siip_system, sys_ED)

    execute_out = PSI.execute!(sim; enable_progress_bar = true)
    println("finish sienna run")

    sim_results = PSI.SimulationResults(sim)
    println("get sienna results")

    res_uc = PSI.get_decision_problem_results(sim_results, "UC")
    res_ed = PSI.get_decision_problem_results(sim_results, "ED")

    dual_values_ed = PSI.read_realized_duals(res_ed)
    dual_values_uc = PSI.read_realized_duals(res_uc)

    result_variables_ed = PSI.read_realized_variables(res_ed)
    result_variables_uc = PSI.read_realized_variables(res_uc)

    data_length_ed = DataFrames.nrow(dual_values_ed["NodalBalanceActiveConstraint__ACBus"])
    data_length_uc = DataFrames.nrow(dual_values_uc["NodalBalanceActiveConstraint__ACBus"])

    if md_market_bool == true
        res_md = PSI.get_decision_problem_results(sim_results, "MD")
        dual_values_md = PSI.read_realized_duals(res_md)
        result_variables_md = PSI.read_realized_variables(res_md)
        data_length_md = DataFrames.nrow(dual_values_md["NodalBalanceActiveConstraint__ACBus"])
    else
        result_variables_md = nothing
        data_length_md = DataFrames.nrow(dual_values_uc["NodalBalanceActiveConstraint__ACBus"])
    end

    energy_price_ed = AxisArrays.AxisArray(zeros(length(zones), 1, data_length_ed), zones, 1:1, 1:data_length_ed)
    energy_price_uc = AxisArrays.AxisArray(zeros(length(zones), 1, data_length_uc), zones, 1:1, 1:data_length_uc)
    energy_price_md = AxisArrays.AxisArray(zeros(length(zones), 1, data_length_md), zones, 1:1, 1:data_length_md)
    energy_voll = AxisArrays.AxisArray(zeros(length(zones), 1, data_length_ed), zones, 1:1, 1:data_length_ed)
    energy_voll_uc = AxisArrays.AxisArray(zeros(length(zones), 1, data_length_uc), zones, 1:1, 1:data_length_uc)
    energy_voll_md = AxisArrays.AxisArray(zeros(length(zones), 1, data_length_md), zones, 1:1, 1:data_length_md)

    reserve_price_ed = Dict(s => zeros(1, data_length_ed) for s in String.(rt_products))
    reserve_price_uc = Dict(s => zeros(1, data_length_uc) for s in String.(da_products))
    reserve_price_md = Dict(s => zeros(1, data_length_md) for s in String.(md_products))
    reserve_voll = Dict(s => zeros(1, data_length_ed) for s in String.(rt_products))
    reserve_voll_uc = Dict(s => zeros(1, data_length_uc) for s in String.(da_products))
    reserve_voll_md = Dict(s => zeros(1, data_length_md) for s in String.(md_products))

    inertia_price = AxisArrays.AxisArray(zeros(1, data_length_ed), 1:1, 1:data_length_ed)
    inertia_voll = AxisArrays.AxisArray(zeros(1, data_length_ed), 1:1, 1:data_length_ed)

    # TODO: may need to add only_md_products
    only_da_products = String[]

    for product in da_products
        if !(product in rt_products)
            push!(only_da_products, product)
            reserve_price[String(product)] = zeros(1, data_length_uc)
            # reserve_voll[String(product)] = zeros(1, data_length_uc)
        end
    end

    for zone in zones
        bus = find_zonal_bus(zone, sys_UC)
        zone_num = parse(Int64, last(zone, 1))
        if isnothing(bus)
            energy_price_ed[zone, 1, :] = zeros(data_length_ed)
            energy_price_uc[zone, 1, :] = zeros(data_length_uc)
            energy_price_md[zone, 1, :] = zeros(data_length_md)
            energy_voll[zone, 1, :] = zeros(data_length_ed)
            energy_voll_uc[zone, 1, :] = zeros(data_length_uc)
            energy_voll_md[zone, 1, :] = zeros(data_length_md)
        else
            energy_price_ed[zone, 1, :] =
                abs.(round.(dual_values_ed["NodalBalanceActiveConstraint__ACBus"][:, PSY.get_name(bus)], digits = 5)) / base_power
            energy_price_uc[zone, 1, :] =
                abs.(round.(dual_values_uc["NodalBalanceActiveConstraint__ACBus"][:, PSY.get_name(bus)], digits = 5)) / base_power
            if md_market_bool == true
                energy_price_md[zone, 1, :] =
                    abs.(round.(dual_values_md["NodalBalanceActiveConstraint__ACBus"][:, PSY.get_name(bus)], digits = 5)) / base_power
                
                energy_voll_md[zone, 1, :] = abs.(round.(result_variables_md["SystemBalanceSlackUp__ACBus"][:, string(PSY.get_number(bus))], digits = 5)) / base_power
            else
                energy_price_md[zone, 1, :] = zeros(data_length_md)
            end
            energy_voll[zone, 1, :] = abs.(round.(result_variables_ed["SystemBalanceSlackUp__ACBus"][:, string(PSY.get_number(bus))], digits = 5)) / base_power
            energy_voll_uc[zone, 1, :] = abs.(round.(result_variables_uc["SystemBalanceSlackUp__ACBus"][:, string(PSY.get_number(bus))], digits = 5)) / base_power
            # energy_voll[zone, 1, :] += abs.(round.(result_variables_ed["SystemBalanceSlackDown__ACBus"][:, string(PSY.get_number(bus))], digits = 5)) / base_power
        end
    end

    println("recorded sienna energy prices")

    #println(any(isnan, energy_price))
    replace!(energy_price_ed, NaN => 0.0)
    replace!(energy_price_uc, NaN => 0.0)
    replace!(energy_price_md, NaN => 0.0)
    scale_voll(energy_price_ed, rt_resolution)

    #println(any(isnan, energy_price))
    #println(Statistics.mean(energy_price))


    for service in get_system_services(sys_ED)
        name = PSY.get_name(service)
        #println(name)
        if typeof(service) == PSY.VariableReserve{PSY.ReserveUp}
            if name == "Inertia"
                inertia_price[1, :] = abs.(round.(dual_values_ed["RequirementConstraint__VariableReserve__ReserveUp__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
                replace!(inertia_price, NaN => 0.0)
                scale_voll(inertia_price, rt_resolution)
                inertia_voll[1, :] = abs.(round.(result_variables_ed["ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)"][:, Symbol("ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)")], digits = 5)) / base_power
            elseif name == "Clean_Energy"
            else
                reserve_price_ed[name][1, :] = abs.(round.(dual_values_ed["RequirementConstraint__VariableReserve__ReserveUp__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
                replace!(reserve_price_ed[name], NaN => 0.0)
                scale_voll(reserve_price_ed[name], rt_resolution)
                reserve_voll[name][1, :] = abs.(round.(result_variables_ed["ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)"][:, Symbol("ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)")], digits = 5)) / base_power
                #println(reserve_price[name])
            end
        elseif typeof(service) == PSY.ReserveDemandCurve{PSY.ReserveUp}
            reserve_price_ed[name][1, :] = abs.(round.(dual_values_ed["RequirementConstraint__ReserveDemandCurve__ReserveUp__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
            replace!(reserve_price_ed[name], NaN => 0.0)
            scale_voll(reserve_price_ed[name], rt_resolution)
            #println(reserve_price[name])
        elseif typeof(service) == PSY.StaticReserve{PSY.ReserveUp}
            reserve_price_ed[name][1, :] = abs.(round.(dual_values_ed["RequirementConstraint__StaticReserve__ReserveUp__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
            replace!(reserve_price_ed[name], NaN => 0.0)
            scale_voll(reserve_price_ed[name], rt_resolution)
            reserve_voll[name][1, :] = abs.(round.(result_variables_ed["ReserveRequirementSlack__StaticReserve__ReserveUp__$(name)"][:, Symbol("ReserveRequirementSlack__StaticReserve__ReserveUp__$(name)")], digits = 5)) / base_power
            #println(reserve_price[name])
        elseif typeof(service) == PSY.VariableReserve{PSY.ReserveDown}
            reserve_price_ed[name][1, :] = abs.(round.(dual_values_ed["RequirementConstraint__VariableReserve__ReserveDown__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
            replace!(reserve_price_ed[name], NaN => 0.0)
            scale_voll(reserve_price_ed[name], rt_resolution)
            reserve_voll[name][1, :] = abs.(round.(result_variables_ed["ReserveRequirementSlack__VariableReserve__ReserveDown__$(name)"][:, Symbol("ReserveRequirementSlack__VariableReserve__ReserveDown__$(name)")], digits = 5)) / base_power
            #println(reserve_price[name])
        end
    end

    for service in get_system_services(sys_UC)
        name = PSY.get_name(service)
        # if name in only_da_products
        if typeof(service) == PSY.VariableReserve{PSY.ReserveUp}
            if name == "Inertia"
                inertia_price[1, :] = abs.(round.(dual_values_uc["RequirementConstraint__VariableReserve__ReserveUp__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
                replace!(inertia_price, NaN => 0.0)
                scale_voll(inertia_price, da_resolution)
                # inertia_voll[1, :] = abs.(round.(result_variables_uc["ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)"][:, Symbol("ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)")], digits = 5)) / base_power
            elseif name == "Clean_Energy"
            else
                reserve_price_uc[name][1, :] = abs.(round.(dual_values_uc["RequirementConstraint__VariableReserve__ReserveUp__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
                replace!(reserve_price_uc[name], NaN => 0.0)
                scale_voll(reserve_price_uc[name], da_resolution)
                reserve_voll_uc[name][1, :] = abs.(round.(result_variables_uc["ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)"][:, Symbol("ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)")], digits = 5)) / base_power
                #println(reserve_price[name])
            end
        elseif typeof(service) == PSY.ReserveDemandCurve{PSY.ReserveUp}
            reserve_price_uc[name][1, :] = abs.(round.(dual_values_uc["RequirementConstraint__ReserveDemandCurve__ReserveUp__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
            replace!(reserve_price_uc[name], NaN => 0.0)
            scale_voll(reserve_price_uc[name], da_resolution)
            #println(reserve_price[name])
        elseif typeof(service) == PSY.StaticReserve{PSY.ReserveUp}
            reserve_price_uc[name][1, :] = abs.(round.(dual_values_uc["RequirementConstraint__StaticReserve__ReserveUp__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
            replace!(reserve_price_uc[name], NaN => 0.0)
            scale_voll(reserve_price_uc[name], da_resolution)
            reserve_voll_uc[name][1, :] = abs.(round.(result_variables_uc["ReserveRequirementSlack__StaticReserve__ReserveUp__$(name)"][:, Symbol("ReserveRequirementSlack__StaticReserve__ReserveUp__$(name)")], digits = 5)) / base_power
            #println(reserve_price[name])
        elseif typeof(service) == PSY.VariableReserve{PSY.ReserveDown}
            reserve_price_uc[name][1, :] = abs.(round.(dual_values_uc["RequirementConstraint__VariableReserve__ReserveDown__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
            replace!(reserve_price_uc[name], NaN => 0.0)
            scale_voll(reserve_price_uc[name], da_resolution)
            reserve_voll_uc[name][1, :] = abs.(round.(result_variables_uc["ReserveRequirementSlack__VariableReserve__ReserveDown__$(name)"][:, Symbol("ReserveRequirementSlack__VariableReserve__ReserveDown__$(name)")], digits = 5)) / base_power
            #println(reserve_price[name])
        end
        # end
    end

    if md_market_bool == true
        for service in get_system_services(sys_MD)
            name = PSY.get_name(service)
            # TODO: need to replace "only_da_products"
            # if name in only_da_products
            if typeof(service) == PSY.VariableReserve{PSY.ReserveUp}
                if name == "Inertia"
                    inertia_price[1, :] = abs.(round.(dual_values_md["RequirementConstraint__VariableReserve__ReserveUp__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
                    replace!(inertia_price, NaN => 0.0)
                    scale_voll(inertia_price, da_resolution)
                    # inertia_voll[1, :] = abs.(round.(result_variables_uc["ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)"][:, Symbol("ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)")], digits = 5)) / base_power
                elseif name == "Clean_Energy"
                else
                    reserve_price_md[name][1, :] = abs.(round.(dual_values_md["RequirementConstraint__VariableReserve__ReserveUp__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
                    replace!(reserve_price_md[name], NaN => 0.0)
                    scale_voll(reserve_price_md[name], da_resolution)
                    reserve_voll_md[name][1, :] = abs.(round.(result_variables_md["ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)"][:, Symbol("ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)")], digits = 5)) / base_power
                    #println(reserve_price[name])
                end
            elseif typeof(service) == PSY.ReserveDemandCurve{PSY.ReserveUp}
                reserve_price_md[name][1, :] = abs.(round.(dual_values_md["RequirementConstraint__ReserveDemandCurve__ReserveUp__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
                replace!(reserve_price_md[name], NaN => 0.0)
                scale_voll(reserve_price_md[name], da_resolution)
                #println(reserve_price[name])
            elseif typeof(service) == PSY.StaticReserve{PSY.ReserveUp}
                reserve_price_md[name][1, :] = abs.(round.(dual_values_md["RequirementConstraint__StaticReserve__ReserveUp__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
                replace!(reserve_price_md[name], NaN => 0.0)
                scale_voll(reserve_price_md[name], da_resolution)
                reserve_voll_md[name][1, :] = abs.(round.(result_variables_md["ReserveRequirementSlack__StaticReserve__ReserveUp__$(name)"][:, Symbol("ReserveRequirementSlack__StaticReserve__ReserveUp__$(name)")], digits = 5)) / base_power
                #println(reserve_price[name])
            elseif typeof(service) == PSY.VariableReserve{PSY.ReserveDown}
                reserve_price_md[name][1, :] = abs.(round.(dual_values_md["RequirementConstraint__VariableReserve__ReserveDown__$(name)"][:, Symbol("$(name)")], digits = 5)) / base_power
                replace!(reserve_price_md[name], NaN => 0.0)
                scale_voll(reserve_price_md[name], da_resolution)
                reserve_voll_md[name][1, :] = abs.(round.(result_variables_md["ReserveRequirementSlack__VariableReserve__ReserveDown__$(name)"][:, Symbol("ReserveRequirementSlack__VariableReserve__ReserveDown__$(name)")], digits = 5)) / base_power
                #println(reserve_price[name])
            end
            # end
        end
    end

    println("recorded sienna service prices")

    sys_techs = get_all_techs(sys_ED)

    tech_names = get_name.(sys_techs)
    capacity_factors_md = Dict([g => zeros(1, data_length_md) for g in tech_names])
    capacity_factors_uc = Dict([g => zeros(1, data_length_uc) for g in tech_names])
    capacity_factors_ed = Dict([g => zeros(1, data_length_ed) for g in tech_names])
    start_up_costs = Dict([g => zeros(1, data_length_uc) for g in tech_names])
    shut_down_costs = Dict([g => zeros(1, data_length_uc) for g in tech_names])

    reserve_perc_md = Dict(g => Dict(s => zeros(1, data_length_md) for s in String.(union(rt_products, da_products, md_products))) for g in tech_names)
    reserve_perc_uc = Dict(g => Dict(s => zeros(1, data_length_uc) for s in String.(union(rt_products, da_products, md_products))) for g in tech_names)
    reserve_perc_ed = Dict(g => Dict(s => zeros(1, data_length_ed) for s in String.(union(rt_products, da_products, md_products))) for g in tech_names)

    inertia_perc = Dict([g => zeros(1, data_length_ed) for g in tech_names])

    # for g in tech_names
    #     for product in only_da_products
    #         reserve_perc[g][string(product)] = zeros(1, data_length_uc)
    #     end
    # end

    for tech in sys_techs
        name = get_name(tech)
        if md_market_bool == true
            capacity_factors_md[name][1, :] = get_realized_capacity_factors(tech, result_variables_md, result_variables_uc, base_power)
        end
        capacity_factors_uc[name][1, :] = get_realized_capacity_factors(tech, result_variables_uc, result_variables_uc, base_power)
        capacity_factors_ed[name][1, :] = get_realized_capacity_factors(tech, result_variables_ed, result_variables_ed, base_power)
        start_up_costs[name][1, :] = get_start_costs(tech, result_variables_uc, data_length_uc)
        shut_down_costs[name][1, :] = get_shut_costs(tech, result_variables_uc, data_length_uc)

        services_ED = PSY.get_services(tech)
        services_UC = PSY.get_services(PSY.get_components_by_name(PSY.Device, sys_UC, name)[1])

        for service in services_UC
            # TODO: figure out if we need clean energy or not
            if PSY.get_name(service) != "Clean_Energy"
                update_realized_reserve_perc!(tech,
                                            service,
                                            result_variables_ed,
                                            result_variables_uc,
                                            result_variables_md,
                                            reserve_perc_md,
                                            reserve_perc_uc,
                                            reserve_perc_ed,
                                            inertia_perc,
                                            rt_products,
                                            da_products,
                                            md_products,
                                            base_power,
                                            md_market_bool,)
            end

        end
    end

    for g in keys(inertia_perc)
        inertia_perc[g] = inertia_perc[g] / PSY.get_base_power(sys_UC)
    end
    println("recorded all sienna results")
    return energy_price_ed, energy_price_uc, energy_price_md, reserve_price_ed, reserve_price_uc, reserve_price_md, inertia_price, capacity_factors_md, capacity_factors_uc, capacity_factors_ed, reserve_perc_md, reserve_perc_uc, reserve_perc_ed, inertia_perc, start_up_costs, shut_down_costs, energy_voll, energy_voll_uc, energy_voll_md, reserve_voll, reserve_voll_uc, reserve_voll_md, inertia_voll;
end
