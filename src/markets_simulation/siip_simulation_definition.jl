# This file contains functions for calling, running and post-processing
# SIIP PSI Simulation for actual Energy and A/S Market Clearing

# PSI's ServiceModel dispatch for add_constraint_dual! calls get_available_components which
# returns ALL services of the type (not just the named one), then passes the collection to
# assign_dual_variable! which only accepts a single D<:PSY.Service. Both are bugs: this
# override mirrors what construct_service! does — fetch only the named service and register
# its dual individually.
function PSI.add_constraint_dual!(
    container::PSI.OptimizationContainer,
    sys::PSY.System,
    model::PSI.ServiceModel{T, D},
) where {T <: PSY.Service, D <: PSI.AbstractServiceFormulation}
    if !isempty(PSI.get_duals(model))
        name = PSI.get_service_name(model)
        service = PSY.get_component(T, sys, name)
        (service === nothing || !PSY.get_available(service)) && return
        for constraint_type in PSI.get_duals(model)
            PSI.assign_dual_variable!(container, constraint_type, service, D)
        end
    end
    return
end

# PSI's generic AbstractPowerModel dispatch for add_constraint_dual! fetches ACBus components
# and registers duals keyed on ACBus. For AreaBalancePowerModel the actual constraint is
# CopperPlateBalanceConstraint keyed on Area, so we add a more-specific dispatch here.
function PSI.add_constraint_dual!(
    container::PSI.OptimizationContainer,
    sys::PSY.System,
    model::PSI.NetworkModel{PSI.AreaBalancePowerModel},
)
    if !isempty(PSI.get_duals(model))
        areas = PSY.get_name.(PSI.get_available_components(model, PSY.Area, sys))
        time_steps = PSI.get_time_steps(container)
        for constraint_type in PSI.get_duals(model)
            PSI.add_dual_container!(container, constraint_type, PSY.Area, areas, time_steps)
        end
    end
    return
end

function adjust_reserve_voll!(sys::PSY.System,
    problem::PSI.OperationModel,
    simulation_dir::String,
    reserve_penalty::String,
    zones::Vector{String},
    default_balance_slack_cost::Float64,
    default_service_slack_cost::Float64,
    energy_voll_cost::AxisArrays.AxisArray{Float64, 1},
)
    base_power = PSY.get_base_power(sys)
    services = get_system_services(sys)

    optimization_container = PSI.get_optimization_container(problem)
    variables = PSI.get_variables(optimization_container)

    for zone in zones
        area = find_zonal_area(sys, zone)
        slack_coefficients = [PSI.SystemBalanceSlackUp, PSI.SystemBalanceSlackDown]
        for c in slack_coefficients
            slack_key = PSI.VariableKey{c, PSY.ACBus}("")
            index = findall(x -> x == PSY.get_name(area), variables[slack_key].axes[1])[1]
            slack_variables = variables[slack_key].data[index, :]
            delta_cost = energy_voll_cost[zone] * base_power - default_balance_slack_cost
            for v in slack_variables
                PSI.add_to_objective_variant_expression!(
                    optimization_container,
                    v * delta_cost,
                )
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
            reserve_data = read_data(
                joinpath(
                    simulation_dir,
                    "markets_data",
                    "$(reserve_penalty)_reserve_penalty",
                    "$(name).csv",
                ),
            )
            slack_variables = variables[PSI.VariableKey{
                PSI.ReserveRequirementSlack,
                PSY.VariableReserve{PSY.ReserveUp},
            }(
                "$name",
            )]
            penalty_price = reserve_data[1, "price_cap"] * base_power
            delta_cost = penalty_price - default_service_slack_cost
        elseif name == "Reg_Down"
            reserve_data = read_data(
                joinpath(
                    simulation_dir,
                    "markets_data",
                    "$(reserve_penalty)_reserve_penalty",
                    "$(name).csv",
                ),
            )
            slack_variables = variables[PSI.VariableKey{
                PSI.ReserveRequirementSlack,
                PSY.VariableReserve{PSY.ReserveDown},
            }(
                "$name",
            )]
            penalty_price = reserve_data[1, "price_cap"] * base_power
            delta_cost = penalty_price - default_service_slack_cost
        end

        for v in slack_variables
            PSI.add_to_objective_variant_expression!(optimization_container, v * delta_cost)
        end
    end
end

function scale_voll(
    price::Union{Array{Float64, 2}, AxisArrays.AxisArray{Float64, 2}},
    resolution::Int64,
)
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
    data_length_uc::Int64,
)
    # start_ups = results["StartVariable__ThermalStandard"][:, Symbol(get_name(device))]
    start_ups = results["StartVariable__ThermalStandard"][
        results["StartVariable__ThermalStandard"].name .== PSY.get_name(device),
        :value,
    ]
    start_up_costs = start_ups * PSY.get_start_up(PSY.get_operation_cost(device))
    return start_up_costs
end

function get_start_costs(device::ThermalFastStartSIIP,
    results::Dict{String, DataFrames.DataFrame},
    data_length_uc::Int64,
)
    # start_ups = results["StartVariable__ThermalFastStartSIIP"][:, Symbol(get_name(device))]
    start_ups = results["StartVariable__ThermalFastStartSIIP"][
        results["StartVariable__ThermalFastStartSIIP"].name .== PSY.get_name(device),
        :value,
    ]
    start_up_costs = start_ups * PSY.get_start_up(PSY.get_operation_cost(device))
    return start_up_costs
end

"""
This function returns total start-up costs for other generators from PSI Simulation.
"""
function get_start_costs(device::PSY.Device,
    results::Dict{String, DataFrames.DataFrame},
    data_length_uc::Int64,
)
    return zeros(data_length_uc)
end

"""
This function returns shut-down costs for Thermal generators from PSI Simulation.
"""
function get_shut_costs(device::PSY.ThermalStandard,
    results::Dict{String, DataFrames.DataFrame},
    data_length_uc::Int64,
)
    # shut_downs = results["StopVariable__ThermalStandard"][:, Symbol(get_name(device))]
    shut_downs = results["StopVariable__ThermalStandard"][
        results["StopVariable__ThermalStandard"].name .== get_name(device),
        :value,
    ]
    shut_down_costs = shut_downs * PSY.get_shut_down(PSY.get_operation_cost(device))
    return shut_down_costs
end

function get_shut_costs(device::ThermalFastStartSIIP,
    results::Dict{String, DataFrames.DataFrame},
    data_length_uc::Int64,
)
    # shut_downs = results["StopVariable__ThermalFastStartSIIP"][:, Symbol(get_name(device))]
    shut_downs = results["StopVariable__ThermalFastStartSIIP"][
        results["StopVariable__ThermalFastStartSIIP"].name .== get_name(device),
        :value,
    ]
    shut_down_costs = shut_downs * PSY.get_shut_down(PSY.get_operation_cost(device))
    return shut_down_costs
end
"""
This function returns total shut-down costs for other generators from PSI Simulation.
"""
function get_shut_costs(device::PSY.Device,
    results::Dict{String, DataFrames.DataFrame},
    data_length_uc::Int64,
)
    return zeros(data_length_uc)
end

"""
This function returns realized capacity factors for ThermalStandard generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.ThermalStandard,
    results::Dict{String, DataFrames.DataFrame},
    results_uc::Dict{String, DataFrames.DataFrame},
    base_power::Float64,
)
    key = "ActivePowerVariable__ThermalStandard"
    if !haskey(results, key)
        @warn "Key '$key' not found in results for $(get_name(device)) — returning zeros"
        n = haskey(results_uc, key) ? length(unique(results_uc[key].DateTime)) : 0
        return zeros(n)
    end
    energy_production = results[key][results[key].name .== get_name(device), :value]
    capacity_factors = (energy_production / base_power) / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for ThermalMultiStart generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.ThermalMultiStart,
    results::Dict{String, DataFrames.DataFrame},
    results_uc::Dict{String, DataFrames.DataFrame},
    base_power::Float64,
)
    key = "ActivePowerVariable__ThermalMultiStart"
    if !haskey(results, key)
        @warn "Key '$key' not found in results for $(get_name(device)) — returning zeros"
        n = haskey(results_uc, key) ? length(unique(results_uc[key].DateTime)) : 0
        return zeros(n)
    end
    energy_production = results[key][results[key].name .== get_name(device), :value]
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
    base_power::Float64,
)
    key = "ActivePowerVariable__ThermalFastStartSIIP"
    if !haskey(results, key)
        @warn "Key '$key' not found in results for $(get_name(device)) — returning zeros"
        n = haskey(results_uc, key) ? length(unique(results_uc[key].DateTime)) : 0
        return zeros(n)
    end
    energy_production = results[key][results[key].name .== get_name(device), :value]
    capacity_factors = (energy_production / base_power) / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for Renewable generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.RenewableDispatch,
    results::Dict{String, DataFrames.DataFrame},
    results_uc::Dict{String, DataFrames.DataFrame},
    base_power::Float64,
)
    key = "ActivePowerVariable__RenewableDispatch"
    if !haskey(results, key)
        @warn "Key '$key' not found in results for $(get_name(device)) — returning zeros"
        n = haskey(results_uc, key) ? length(unique(results_uc[key].DateTime)) : 0
        return zeros(n)
    end
    energy_production = results[key][results[key].name .== get_name(device), :value]
    capacity_factors = (energy_production / base_power) / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for Hydropower generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.HydroDispatch,
    results::Dict{String, DataFrames.DataFrame},
    results_uc::Dict{String, DataFrames.DataFrame},
    base_power::Float64,
)
    key = "ActivePowerVariable__HydroDispatch"
    if !haskey(results, key)
        @warn "Key '$key' not found in results for $(get_name(device)) — returning zeros"
        n = haskey(results_uc, key) ? length(unique(results_uc[key].DateTime)) : 0
        return zeros(n)
    end
    energy_production = results[key][results[key].name .== get_name(device), :value]
    capacity_factors = (energy_production / base_power) / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for Hydropower generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.HydroTurbine,
    results::Dict{String, DataFrames.DataFrame},
    results_uc::Dict{String, DataFrames.DataFrame},
    base_power::Float64,
)
    key = "ActivePowerVariable__HydroTurbine"
    if !haskey(results, key)
        @warn "Key '$key' not found in results for $(get_name(device)) — returning zeros"
        n = haskey(results_uc, key) ? length(unique(results_uc[key].DateTime)) : 0
        return zeros(n)
    end
    energy_production = results[key][results[key].name .== get_name(device), :value]
    capacity_factors = (energy_production / base_power) / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for Storage generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.EnergyReservoirStorage,
    results::Dict{String, DataFrames.DataFrame},
    results_uc::Dict{String, DataFrames.DataFrame},
    base_power::Float64,
)
    out_key = "ActivePowerOutVariable__EnergyReservoirStorage"
    in_key = "ActivePowerInVariable__EnergyReservoirStorage"
    if !haskey(results, out_key) || !haskey(results, in_key)
        @warn "Storage result keys not found for $(get_name(device)) — returning zeros"
        n = haskey(results_uc, out_key) ? length(unique(results_uc[out_key].DateTime)) : 0
        return zeros(n)
    end
    energy_production =
        results[out_key][results[out_key].name .== get_name(device), :value] -
        results[in_key][results[in_key].name .== get_name(device), :value]
    capacity_factors = (energy_production / base_power) / get_device_size(device)
    return capacity_factors
end

"""
This function returns nothing if Service is not of ReserveUp or ReserveDown type.
"""
function update_realized_reserve_perc!(device::PSY.Device,
    service::S,
    results_ed::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    results_uc::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    results_md::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
    reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
    reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
    inertia_perc::Dict{String, Array{Float64, 2}},
    rt_products::Vector{SubString{String}},
    da_products::Vector{SubString{String}},
    md_products::Vector{SubString{String}},
    base_power::Float64,
    md_market_bool::Bool,
    single_stage_bool::Bool) where {S <: PSY.Service}
    return
end

"""
This function returns realized reserve up provision percentages from PSI Simulation.
"""
function update_realized_reserve_perc!(device::PSY.Device,
    service::PSY.VariableReserve{PSY.ReserveUp},
    results_ed::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    results_uc::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    results_md::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
    reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
    reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
    inertia_perc::Dict{String, Array{Float64, 2}},
    rt_products::Vector{SubString{String}},
    da_products::Vector{SubString{String}},
    md_products::Vector{SubString{String}},
    base_power::Float64,
    md_market_bool::Bool,
    single_stage_bool::Bool)
    service_name = PSY.get_name(service)

    if service_name == "Inertia"
        ed_key = "ActivePowerReserveVariable__VariableReserve__ReserveUp__$(service_name)"
        if haskey(results_ed, ed_key)
            inertia_provision =
                results_ed[ed_key][results_ed[ed_key].name .== PSY.get_name(device), :value]
            inertia_perc_value = inertia_provision / get_device_size(device) / base_power
            inertia_perc[PSY.get_name(device)][1, :] = inertia_perc_value
        else
            @warn "No ED inertia variable found for $(PSY.get_name(device)), skipping"
        end
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
        ed_key = "ActivePowerReserveVariable__VariableReserve__ReserveUp__$(service_name)"
        uc_key = "ActivePowerReserveVariable__VariableReserve__ReserveUp__$(service_name)"
        if single_stage_bool == false
            if haskey(results_ed, ed_key)
                reserve_provision_ed = results_ed[ed_key][
                    results_ed[ed_key].name .== PSY.get_name(device),
                    :value,
                ]
                reserve_perc_value_ed =
                    reserve_provision_ed / get_device_size(device) / base_power
                reserve_perc_ed[PSY.get_name(device)][service_name][1, :] =
                    reserve_perc_value_ed
            else
                @warn "No ED reserve variable found for $(PSY.get_name(device)) - $(service_name), skipping"
            end
        end

        if haskey(results_uc, uc_key)
            reserve_provision_uc =
                results_uc[uc_key][results_uc[uc_key].name .== PSY.get_name(device), :value]
            reserve_perc_value_uc =
                reserve_provision_uc / get_device_size(device) / base_power
            reserve_perc_uc[PSY.get_name(device)][service_name][1, :] =
                reserve_perc_value_uc
        else
            @warn "No UC reserve variable found for $(PSY.get_name(device)) - $(service_name), skipping"
        end
    end

    if md_market_bool == true
        md_key = "ActivePowerReserveVariable__VariableReserve__ReserveUp__$(service_name)"
        if haskey(results_md, md_key)
            reserve_provision_md =
                results_md[md_key][results_md[md_key].name .== PSY.get_name(device), :value]
            reserve_perc_value_md =
                reserve_provision_md / get_device_size(device) / base_power
            reserve_perc_md[PSY.get_name(device)][service_name][1, :] =
                reserve_perc_value_md
        else
            @warn "No MD reserve variable found for $(PSY.get_name(device)) - $(service_name), skipping"
        end
    end

    return
end

"""
This function returns realized ordc provision percentages from PSI Simulation.
"""
function update_realized_reserve_perc!(device::PSY.Device,
    service::PSY.ReserveDemandCurve{PSY.ReserveUp},
    results_ed::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    results_uc::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    results_md::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
    reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
    reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
    inertia_perc::Dict{String, Array{Float64, 2}},
    rt_products::Vector{SubString{String}},
    da_products::Vector{SubString{String}},
    md_products::Vector{SubString{String}},
    base_power::Float64,
    md_market_bool::Bool,
    single_stage_bool::Bool)
    service_name = PSY.get_name(service)

    if single_stage_bool == false
        @info "Updating reserve provision for $(PSY.get_name(device)) - $(service_name) from ED results"
        key_ed = "ActivePowerReserveVariable__ReserveDemandCurve__ReserveUp__$(service_name)"
        if haskey(results_ed, key_ed)
            reserve_provision_ed =
                results_ed[key_ed][results_ed[key_ed].name .== PSY.get_name(device), :value]
            reserve_perc_value_ed =
                reserve_provision_ed / get_device_size(device) / base_power
            reserve_perc_ed[PSY.get_name(device)][service_name][1, :] =
                reserve_perc_value_ed
        else
            @warn "No ED reserve variable found for $(PSY.get_name(device)) - $(service_name), skipping"
        end
    end

    ### NY_change: need to comment this back once ORDC is enabled again
    # reserve_provision_uc = results_uc["ActivePowerReserveVariable__ReserveDemandCurve__ReserveUp__$(service_name)"][:, Symbol(PSY.get_name(device))]
    # reserve_perc_value_uc = reserve_provision_uc / get_device_size(device) / base_power
    # reserve_perc_uc[PSY.get_name(device)][service_name][1, :] = reserve_perc_value_uc

    if md_market_bool == true
        @info "Updating reserve provision for $(PSY.get_name(device)) - $(service_name) from MD results"
        key_md = "ActivePowerReserveVariable__ReserveDemandCurve__ReserveUp__$(service_name)"
        if haskey(results_md, key_md)
            reserve_provision_md =
                results_md[key_md][results_md[key_md].name .== PSY.get_name(device), :value]
            reserve_perc_value_md =
                reserve_provision_md / get_device_size(device) / base_power
            reserve_perc_md[PSY.get_name(device)][service_name][1, :] =
                reserve_perc_value_md
        else
            @warn "No MD reserve variable found for $(PSY.get_name(device)) - $(service_name), skipping"
        end
    end

    # if service_name in rt_products
    #     reserve_provision = results_ed["ActivePowerReserveVariable__ReserveDemandCurve__ReserveUp__$(service_name)"][:, Symbol(PSY.get_name(device))]
    #     reserve_perc_value = reserve_provision / get_device_size(device) / base_power
    #     reserve_perc[PSY.get_name(device)][service_name][1, :] = reserve_perc_value

    # elseif service_name in only_da_products
    #     reserve_provision = results_uc["ActivePowerReserveVariable__ReserveDemandCurve__ReserveUp__$(service_name)"][:, Symbol(PSY.get_name(device))]
    #     reserve_perc_value = reserve_provision / get_device_size(device) / base_power
    #     reserve_perc[PSY.get_name(device)][service_name][1, :] = reserve_perc_value

    # end
    return
end

"""
This function returns realized reserve down provision percentages from PSI Simulation.
"""
function update_realized_reserve_perc!(device::PSY.Device,
    service::PSY.VariableReserve{PSY.ReserveDown},
    results_ed::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    results_uc::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    results_md::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
    reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
    reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
    inertia_perc::Dict{String, Array{Float64, 2}},
    rt_products::Vector{SubString{String}},
    da_products::Vector{SubString{String}},
    md_products::Vector{SubString{String}},
    base_power::Float64,
    md_market_bool::Bool,
    single_stage_bool::Bool)
    service_name = PSY.get_name(service)

    if single_stage_bool == false
        @info "Updating reserve provision for $(PSY.get_name(device)) - $(service_name) from ED results"
        key_ed = "ActivePowerReserveVariable__VariableReserve__ReserveDown__$(service_name)"
        if haskey(results_ed, key_ed)
            reserve_provision_ed =
                results_ed[key_ed][results_ed[key_ed].name .== PSY.get_name(device), :value]
            reserve_perc_value_ed =
                reserve_provision_ed / get_device_size(device) / base_power
            reserve_perc_ed[PSY.get_name(device)][service_name][1, :] =
                reserve_perc_value_ed
        else
            @warn "No ED reserve variable found for $(PSY.get_name(device)) - $(service_name), skipping"
        end
    end

    key_uc = "ActivePowerReserveVariable__VariableReserve__ReserveDown__$(service_name)"
    if haskey(results_uc, key_uc)
        reserve_provision_uc =
            results_uc[key_uc][results_uc[key_uc].name .== PSY.get_name(device), :value]
        reserve_perc_value_uc = reserve_provision_uc / get_device_size(device) / base_power
        reserve_perc_uc[PSY.get_name(device)][service_name][1, :] = reserve_perc_value_uc
    else
        @warn "No UC reserve variable found for $(PSY.get_name(device)) - $(service_name), skipping"
    end

    if md_market_bool == true
        key_md = "ActivePowerReserveVariable__VariableReserve__ReserveDown__$(service_name)"
        if haskey(results_md, key_md)
            reserve_provision_md =
                results_md[key_md][results_md[key_md].name .== PSY.get_name(device), :value]
            reserve_perc_value_md =
                reserve_provision_md / get_device_size(device) / base_power
            reserve_perc_md[PSY.get_name(device)][service_name][1, :] =
                reserve_perc_value_md
        else
            @warn "No MD reserve variable found for $(PSY.get_name(device)) - $(service_name), skipping"
        end
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
function update_realized_reserve_perc!(device::PSY.EnergyReservoirStorage,
    service::PSY.ReserveDemandCurve{PSY.ReserveUp},
    results_ed::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    results_uc::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    results_md::Union{Nothing, Dict{String, DataFrames.DataFrame}},
    reserve_perc_md::Dict{String, Dict{String, Array{Float64, 2}}},
    reserve_perc_uc::Dict{String, Dict{String, Array{Float64, 2}}},
    reserve_perc_ed::Dict{String, Dict{String, Array{Float64, 2}}},
    inertia_perc::Dict{String, Array{Float64, 2}},
    rt_products::Vector{SubString{String}},
    da_products::Vector{SubString{String}},
    md_products::Vector{SubString{String}},
    base_power::Float64,
    md_market_bool::Bool,
    single_stage_bool::Bool)
    service_name = PSY.get_name(service)

    if single_stage_bool == false
        @info "Updating reserve provision for $(PSY.get_name(device)) - $(service_name) from ED results"
        key_ed_discharge = "AncillaryServiceVariableDischarge__EnergyReservoirStorage__ReserveDemandCurve{ReserveUp}_$(service_name)"
        key_ed_charge = "AncillaryServiceVariableCharge__EnergyReservoirStorage__ReserveDemandCurve{ReserveUp}_$(service_name)"
        if haskey(results_ed, key_ed_discharge) && haskey(results_ed, key_ed_charge)
            reserve_provision_ed =
                results_ed[key_ed_discharge][
                    results_ed[key_ed_discharge].name .== PSY.get_name(device),
                    :value,
                ] +
                results_ed[key_ed_charge][
                    results_ed[key_ed_charge].name .== PSY.get_name(device),
                    :value,
                ]
            reserve_perc_value_ed =
                reserve_provision_ed / get_device_size(device) / base_power
            reserve_perc_ed[PSY.get_name(device)][service_name][1, :] =
                reserve_perc_value_ed
        else
            @warn "No ED reserve variable found for $(PSY.get_name(device)) - $(service_name), skipping"
        end
    end
    ### NY_change: need to comment this back once ORDC is enabled again
    # reserve_provision_uc = (results_uc["AncillaryServiceVariableDischarge__EnergyReservoirStorage__ReserveDemandCurve{ReserveUp}_$(service_name)"][:, Symbol(get_name(device))] .+
    #     results_uc["AncillaryServiceVariableCharge__EnergyReservoirStorage__ReserveDemandCurve{ReserveUp}_$(service_name)"][:, Symbol(get_name(device))])
    # reserve_perc_value_uc = reserve_provision_uc / get_device_size(device) / base_power
    # reserve_perc_uc[get_name(device)][service_name][1, :] = reserve_perc_value_uc

    if md_market_bool == true
        @info "Updating reserve provision for $(PSY.get_name(device)) - $(service_name) from MD results"
        key_md_discharge = "AncillaryServiceVariableDischarge__EnergyReservoirStorage__ReserveDemandCurve{ReserveUp}_$(service_name)"
        key_md_charge = "AncillaryServiceVariableCharge__EnergyReservoirStorage__ReserveDemandCurve{ReserveUp}_$(service_name)"
        if haskey(results_md, key_md_discharge) && haskey(results_md, key_md_charge)
            reserve_provision_md =
                results_md[key_md_discharge][
                    results_md[key_md_discharge].name .== PSY.get_name(device),
                    :value,
                ] +
                results_md[key_md_charge][
                    results_md[key_md_charge].name .== PSY.get_name(device),
                    :value,
                ]
            reserve_provision_md = (
                results_md[key_md_discharge][:, Symbol(PSY.get_name(device))] .+
                results_md[key_md_charge][:, Symbol(PSY.get_name(device))]
            )
            reserve_perc_value_md =
                reserve_provision_md / get_device_size(device) / base_power
            reserve_perc_md[PSY.get_name(device)][service_name][1, :] =
                reserve_perc_value_md
        else
            @warn "No MD reserve variable found for $(PSY.get_name(device)) - $(service_name), skipping"
        end
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
function create_md_template(inertia_product)
    if !(isempty(inertia_product))
        template = PSI.ProblemTemplate(
            PSI.NetworkModel(
                PSI.AreaBalancePowerModel;
                duals = [PSI.CopperPlateBalanceConstraint],
                use_slacks = true,
            ),
        )
        # PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalStandardUnitCommitment)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalStandardUnitCommitment)
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicUnitCommitment)
        PSI.set_device_model!(
            template,
            ThermalFastStartSIIP,
            PSI.ThermalBasicUnitCommitment,
        )
        # PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicDispatch)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalBasicDispatch)
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalStandardUCOutages)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, RPSI.ThermalStandardUCOutages)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSI.RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableNonDispatch, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.StandardLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroTurbine, HSI.HydroCommitmentRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroCommitmentRunOfRiver) # TODO: check which hydro device we have
        # PSI.set_device_model!(template, PSY.HydroEnergyReservoir, HSI.HydroDispatchRunOfRiver)
        # PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroDispatchRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(
            template,
            PSY.EnergyReservoirStorage,
            SSI.StorageDispatchWithReserves,
        )
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(
            template,
            PSY.TwoTerminalGenericHVDCLine,
            PSI.HVDCTwoTerminalLossless,
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Reg_Up";
                use_slacks = true,
                duals = [PSI.RequirementConstraint],
            ),
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveDown},
                PSI.RangeReserve,
                "Reg_Down";
                use_slacks = true,
                duals = [PSI.RequirementConstraint],
            ),
        )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.ReserveDemandCurve{PSY.ReserveUp},
        #         PSI.StepwiseCostReserve,
        #         "Synchronous";
        #         use_slacks = true,
        #         duals = [PSI.RequirementConstraint],
        #     ),
        # )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.ReserveDemandCurve{PSY.ReserveUp},
        #         PSI.StepwiseCostReserve,
        #         "Primary";
        #         use_slacks = true,
        #         duals = [PSI.RequirementConstraint],
        #     ),
        # )

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
                PSI.AreaBalancePowerModel;
                duals = [PSI.CopperPlateBalanceConstraint],
                use_slacks = true,
            ),
        )
        # PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalStandardUnitCommitment)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalStandardUnitCommitment)
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicUnitCommitment)
        PSI.set_device_model!(
            template,
            ThermalFastStartSIIP,
            PSI.ThermalBasicUnitCommitment,
        )
        # PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicDispatch)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, PSI.ThermalBasicDispatch)
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalStandardUCOutages)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, RPSI.ThermalStandardUCOutages)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSI.RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableNonDispatch, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.StandardLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroTurbine, HSI.HydroCommitmentRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroCommitmentRunOfRiver) # TODO: check which hydro device we have
        # PSI.set_device_model!(template, PSY.HydroEnergyReservoir, HSI.HydroDispatchRunOfRiver)
        # PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroDispatchRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(
            template,
            PSY.EnergyReservoirStorage,
            SSI.StorageDispatchWithReserves,
        )
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(
            template,
            PSY.TwoTerminalGenericHVDCLine,
            PSI.HVDCTwoTerminalLossless,
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Reg_Up";
                use_slacks = true,
                duals = [PSI.RequirementConstraint],
            ),
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveDown},
                PSI.RangeReserve,
                "Reg_Down";
                use_slacks = true,
                duals = [PSI.RequirementConstraint],
            ),
        )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.ReserveDemandCurve{PSY.ReserveUp},
        #         PSI.StepwiseCostReserve,
        #         "Synchronous";
        #         use_slacks = true,
        #         duals = [PSI.RequirementConstraint],
        #     ),
        # )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.ReserveDemandCurve{PSY.ReserveUp},
        #         PSI.StepwiseCostReserve,
        #         "Primary";
        #         use_slacks = true,
        #         duals = [PSI.RequirementConstraint],
        #     ),
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
    end

    return template
end

function create_uc_template(inertia_product)
    if !(isempty(inertia_product))
        template = PSI.ProblemTemplate(
            PSI.NetworkModel(
                PSI.AreaBalancePowerModel;
                duals = [PSI.CopperPlateBalanceConstraint],
                use_slacks = true,
            ),
        )
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicUnitCommitment)
        PSI.set_device_model!(
            template,
            ThermalFastStartSIIP,
            PSI.ThermalBasicUnitCommitment,
        )
        PSI.set_device_model!(
            template,
            PSY.ThermalMultiStart,
            PSI.ThermalBasicUnitCommitment,
        )
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalStandardUCOutages)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, RPSI.ThermalStandardUCOutages)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSI.RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableNonDispatch, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.StandardLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroTurbine, HSI.HydroCommitmentRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroCommitmentRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(
            template,
            PSY.EnergyReservoirStorage,
            SSI.StorageDispatchWithReserves,
        )
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(
            template,
            PSY.TwoTerminalGenericHVDCLine,
            PSI.HVDCTwoTerminalLossless,
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Reg_Up";
                use_slacks = true,
                duals = [PSI.RequirementConstraint],
            ),
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveDown},
                PSI.RangeReserve,
                "Reg_Down";
                use_slacks = true,
                duals = [PSI.RequirementConstraint],
            ),
        )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.ReserveDemandCurve{PSY.ReserveUp},
        #         PSI.StepwiseCostReserve,
        #         "Synchronous";
        #         use_slacks = true,
        #         duals = [PSI.RequirementConstraint],
        #     ),
        # )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.ReserveDemandCurve{PSY.ReserveUp},
        #         PSI.StepwiseCostReserve,
        #         "Primary";
        #         use_slacks = true,
        #         duals = [PSI.RequirementConstraint],
        #     ),
        # )
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
                PSI.AreaBalancePowerModel;
                duals = [PSI.CopperPlateBalanceConstraint],
                use_slacks = true,
            ),
        )
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicUnitCommitment)
        PSI.set_device_model!(
            template,
            ThermalFastStartSIIP,
            PSI.ThermalBasicUnitCommitment,
        )
        PSI.set_device_model!(
            template,
            PSY.ThermalMultiStart,
            PSI.ThermalBasicUnitCommitment,
        )
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalStandardUCOutages)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, RPSI.ThermalStandardUCOutages)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSI.RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableNonDispatch, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.StandardLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroTurbine, HSI.HydroCommitmentRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroCommitmentRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(
            template,
            PSY.EnergyReservoirStorage,
            SSI.StorageDispatchWithReserves,
        )
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(
            template,
            PSY.TwoTerminalGenericHVDCLine,
            PSI.HVDCTwoTerminalLossless,
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Reg_Up";
                use_slacks = true,
                duals = [PSI.RequirementConstraint],
            ),
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveDown},
                PSI.RangeReserve,
                "Reg_Down";
                use_slacks = true,
                duals = [PSI.RequirementConstraint],
            ),
        )

        # ## NY_change: need to re-enable these
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.ReserveDemandCurve{PSY.ReserveUp},
        #         PSI.StepwiseCostReserve,
        #         "Synchronous";
        #         use_slacks = true,
        #         duals = [PSI.RequirementConstraint],
        #     ),
        # )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.ReserveDemandCurve{PSY.ReserveUp},
        #         PSI.StepwiseCostReserve,
        #         "Primary";
        #         use_slacks = true,
        #         duals = [PSI.RequirementConstraint],
        #     ),
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
    end

    return template
end

"""
This function creates the Economic Dispatch template for PSI Simulation.
"""
#TODO: Update needed
function create_ed_template(inertia_product)
    if !(isempty(inertia_product))
        template = PSI.ProblemTemplate(
            PSI.NetworkModel(
                PSI.AreaBalancePowerModel;
                duals = [PSI.CopperPlateBalanceConstraint],
                use_slacks = true,
            ),
        )
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicDispatch)
        PSI.set_device_model!(
            template,
            ThermalFastStartSIIP,
            PSI.ThermalBasicUnitCommitment,
        )
        PSI.set_device_model!(
            template,
            PSY.ThermalMultiStart,
            PSI.ThermalBasicUnitCommitment,
        )
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalDispatchOutages)
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalRampLimitedOutages)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, RPSI.ThermalStandardUCOutages)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSI.RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableNonDispatch, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.StandardLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroTurbine, HSI.HydroDispatchRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroDispatchRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(
            template,
            PSY.EnergyReservoirStorage,
            SSI.StorageDispatchWithReserves,
        )
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(
            template,
            PSY.TwoTerminalGenericHVDCLine,
            PSI.HVDCTwoTerminalLossless,
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Reg_Up";
                use_slacks = true,
                duals = [PSI.RequirementConstraint],
            ),
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveDown},
                PSI.RangeReserve,
                "Reg_Down";
                use_slacks = true,
                duals = [PSI.RequirementConstraint],
            ),
        )
        # NOTE: StepwiseCostReserve (ORDC) is intentionally disabled in the ED inertia branch.
        # Including ReserveDemandCurve in ED triggers a PSI HDF5 dimension bug:
        # DecrementalPiecewiseLinearSlopeParameter is initialized as a 1D dataset but
        # written as 2D at runtime → "Number of indices does not match dimension of Dataspace".
        # ORDC clears in UC; ED dispatches committed units via feedforward — ORDC not needed in ED.
        # Re-enable only when the upstream PSI bug is fixed.
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.ReserveDemandCurve{PSY.ReserveUp},
        #         PSI.StepwiseCostReserve,
        #         "Synchronous";
        #         use_slacks = true,
        #         duals = [PSI.RequirementConstraint],
        #     ),
        # )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.ReserveDemandCurve{PSY.ReserveUp},
        #         PSI.StepwiseCostReserve,
        #         "Primary";
        #         use_slacks = true,
        #         duals = [PSI.RequirementConstraint],
        #     ),
        # )
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
                PSI.AreaBalancePowerModel;
                duals = [PSI.CopperPlateBalanceConstraint],
                use_slacks = true,
            ),
        )
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicDispatch)
        PSI.set_device_model!(
            template,
            ThermalFastStartSIIP,
            PSI.ThermalBasicUnitCommitment,
        )
        PSI.set_device_model!(
            template,
            PSY.ThermalMultiStart,
            PSI.ThermalBasicUnitCommitment,
        )
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalDispatchOutages)
        # PSI.set_device_model!(template, PSY.ThermalStandard, RPSI.ThermalRampLimitedOutages)
        # PSI.set_device_model!(template, ThermalFastStartSIIP, RPSI.ThermalStandardUCOutages)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSI.RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableNonDispatch, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.StandardLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroTurbine, HSI.HydroDispatchRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, HSI.HydroDispatchRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(
            template,
            PSY.EnergyReservoirStorage,
            SSI.StorageDispatchWithReserves,
        )
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(
            template,
            PSY.TwoTerminalGenericHVDCLine,
            PSI.HVDCTwoTerminalLossless,
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveUp},
                PSI.RangeReserve,
                "Reg_Up";
                use_slacks = true,
                duals = [PSI.RequirementConstraint],
            ),
        )
        PSI.set_service_model!(
            template,
            PSI.ServiceModel(
                PSY.VariableReserve{PSY.ReserveDown},
                PSI.RangeReserve,
                "Reg_Down";
                use_slacks = true,
                duals = [PSI.RequirementConstraint],
            ),
        )
        # NOTE: StepwiseCostReserve (ORDC) intentionally disabled in ED — see note in inertia branch above.
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.ReserveDemandCurve{PSY.ReserveUp},
        #         PSI.StepwiseCostReserve,
        #         "Synchronous";
        #         use_slacks = true,
        #         duals = [PSI.RequirementConstraint],
        #     ),
        # )
        # PSI.set_service_model!(
        #     template,
        #     PSI.ServiceModel(
        #         PSY.ReserveDemandCurve{PSY.ReserveUp},
        #         PSI.StepwiseCostReserve,
        #         "Primary";
        #         use_slacks = true,
        #         duals = [PSI.RequirementConstraint],
        #     ),
        # )
    end

    return template
end

"""
This function creates the Problem for PSI Simulation.
"""

function create_problem(
    template::PSI.ProblemTemplate,
    sys::PSY.System,
    type::String,
    solver::JuMP.MOI.OptimizerWithAttributes,
    inertia_product,
)
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
            export_pwl_vars = true,
            initialize_model = false,
        )

        ##TODO: Need to check this logic, both branches are identical
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
                export_pwl_vars = true,
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
                export_pwl_vars = true,
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
            export_pwl_vars = true,
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
function create_sequence(
    problems::PSI.SimulationModels,
    feedforward_dict,
    single_stage_bool,
)
    if single_stage_bool == false
        sequence = PSI.SimulationSequence(;
            models = problems,
            feedforwards = feedforward_dict,
            ini_cond_chronology = PSI.InterProblemChronology(),
        )
    else
        sequence = PSI.SimulationSequence(;
            models = problems,
            ini_cond_chronology = PSI.InterProblemChronology(),
        )
    end

    return sequence
end

PSY.get_max_output_fraction(value::PSY.ReserveDemandCurve{PSY.ReserveUp}) = 1.0

"""
This function creates the PSI Simulation and post-processes the results.
"""
#TODO: Update needed
function create_simulation(sys_MD::PSY.System,
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
    single_stage_bool::Bool,
    siip_system;
    kwargs...)

    # to_json(sys_MD, "/kfs2/projects/gmlcmarkets/Phase2_EMIS_Analysis/Feb2024_ERCOT_2011_MARKET_Test_NGUO_LDES/HPC_Analysis_Runs/storage_ff_debug/modified_test_sys/MD_sys.json", force=true)
    # to_json(sys_UC, "/kfs2/projects/gmlcmarkets/Phase2_EMIS_Analysis/Feb2024_ERCOT_2011_MARKET_Test_NGUO_LDES/HPC_Analysis_Runs/storage_ff_debug/modified_test_sys/UC_sys.json", force=true)
    # to_json(sys_ED, "/kfs2/projects/gmlcmarkets/Phase2_EMIS_Analysis/Feb2024_ERCOT_2011_MARKET_Test_NGUO_LDES/HPC_Analysis_Runs/storage_ff_debug/modified_test_sys/ED_sys.json", force=true)

    # hacky way to incorporate reserve voll
    base_power = PSY.get_base_power(sys_UC)
    reserve_data = read_data(
        joinpath(
            simulation_dir,
            "markets_data",
            "$(reserve_penalty)_reserve_penalty",
            "Reg_Up.csv",
        ),
    )
    penalty_price = reserve_data[1, "price_cap"] * base_power

    default_service_slack_cost = penalty_price
    default_balance_slack_cost = BALANCE_SLACK_COST

    inertia_product = collect(PSY.get_components_by_name(PSY.Service, sys_ED, "Inertia"))

    ##TODO: Temporary fix to avoid error when inertia product is not included in the system.
    for sys in [sys_ED, sys_UC, sys_MD]
        reg_down = PSY.get_component(PSY.VariableReserve{PSY.ReserveDown}, sys, "Reg_Down")
        reg_up = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, sys, "Reg_Up")
        for service in filter(!isnothing, [reg_down, reg_up])
            for device in Iterators.flatten([
                PSY.get_components(PSY.ThermalStandard, sys),
                PSY.get_components(PSY.EnergyReservoirStorage, sys),
                PSY.get_components(PSY.HydroDispatch, sys),
                PSY.get_components(PSY.HydroTurbine, sys),
            ])
                if !(service in PSY.get_services(device))
                    # @info "Adding $(service) to $(device)"
                    PSY.add_service!(device, service, sys)
                end
            end
        end
    end

    template_uc = create_uc_template(inertia_product)
    uc_problem = create_problem(template_uc, sys_UC, "UC", solver, inertia_product)

    if !single_stage_bool
        template_ed = create_ed_template(inertia_product)
        ed_problem = create_problem(template_ed, sys_ED, "ED", solver, inertia_product)
    end

    if md_market_bool == true && single_stage_bool == false
        template_md = create_md_template(inertia_product)
        md_problem = create_problem(template_md, sys_MD, "MD", solver, inertia_product)

        if isempty(inertia_product)
            feedforward_dict = Dict(
                # "UC" => [
                #     SSI.EnergyLimitFeedforward(
                #         component_type = PSY.EnergyReservoirStorage,
                #         source = PSI.ActivePowerOutVariable,
                #         affected_values = [PSI.ActivePowerOutVariable],
                #         number_of_periods = 24,
                #     ),
                # ],
                "UC" => [
                    SSI.EnergyTargetFeedforward(;
                        component_type = PSY.EnergyReservoirStorage,
                        source = PSI.EnergyVariable,
                        affected_values = [PSI.EnergyVariable],
                        target_period = 24,
                        penalty_cost = PENALTY_COST,
                    ),
                ],
                "ED" => [
                    SSI.EnergyTargetFeedforward(;
                        component_type = PSY.EnergyReservoirStorage,
                        source = PSI.EnergyVariable,
                        affected_values = [PSI.EnergyVariable],
                        target_period = 2,
                        penalty_cost = PENALTY_COST,
                    ),
                    PSI.SemiContinuousFeedforward(;
                        component_type = PSY.ThermalStandard,
                        source = PSI.OnVariable,
                        affected_values = [PSI.ActivePowerVariable],
                    ),
                ],
                # "ED" => [
                #     SSI.EnergyLimitFeedforward(
                #         component_type = PSY.EnergyReservoirStorage,
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
                #         component_type = PSY.EnergyReservoirStorage,
                #         source = PSI.ActivePowerOutVariable,
                #         affected_values = [PSI.ActivePowerOutVariable],
                #         number_of_periods = 24,
                #     ),
                # ],
                "UC" => [
                    SSI.EnergyTargetFeedforward(;
                        component_type = PSY.EnergyReservoirStorage,
                        source = PSI.EnergyVariable,
                        affected_values = [PSI.EnergyVariable],
                        target_period = 24,
                        penalty_cost = PENALTY_COST,
                    ),
                ],
                "ED" => [
                    SSI.EnergyTargetFeedforward(;
                        component_type = PSY.EnergyReservoirStorage,
                        source = PSI.EnergyVariable,
                        affected_values = [PSI.EnergyVariable],
                        target_period = 2,
                        penalty_cost = PENALTY_COST,
                    ),
                    PSI.SemiContinuousFeedforward(;
                        component_type = PSY.ThermalStandard,
                        source = PSI.OnVariable,
                        affected_values = [PSI.ActivePowerVariable],
                    ),
                ],
                # "ED" => [
                #     SSI.EnergyLimitFeedforward(
                #         component_type = PSY.EnergyReservoirStorage,
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

        models = PSI.SimulationModels(;
            decision_models = [
                md_problem,
                uc_problem,
                ed_problem,
            ],
        );

    elseif md_market_bool == false && single_stage_bool == false
        if isempty(inertia_product)
            feedforward_dict = Dict(
                # "UC" => [
                #     SSI.EnergyLimitFeedforward(
                #         component_type = PSY.EnergyReservoirStorage,
                #         source = PSI.EnergyVariable,
                #         affected_values = [PSI.EnergyVariable],
                #         number_of_periods = 36,
                #     ),
                # ],
                "ED" => [
                    SSI.EnergyTargetFeedforward(;
                        component_type = PSY.EnergyReservoirStorage,
                        source = PSI.EnergyVariable,
                        affected_values = [PSI.EnergyVariable],
                        target_period = 2,
                        penalty_cost = PENALTY_COST,
                    ),
                    PSI.SemiContinuousFeedforward(;
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
                #         component_type = PSY.EnergyReservoirStorage,
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
                #         component_type = PSY.EnergyReservoirStorage,
                #         source = PSI.EnergyVariable,
                #         affected_values = [PSI.EnergyVariable],
                #         number_of_periods = 36,
                #     ),
                # ],
                "ED" => [
                    SSI.EnergyTargetFeedforward(;
                        component_type = PSY.EnergyReservoirStorage,
                        source = PSI.EnergyVariable,
                        affected_values = [PSI.EnergyVariable],
                        target_period = 2,
                        penalty_cost = PENALTY_COST,
                    ),
                    PSI.SemiContinuousFeedforward(;
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

        models = PSI.SimulationModels(;
            decision_models = [
                uc_problem,
                ed_problem,
            ],
        );

    elseif md_market_bool == false && single_stage_bool == true
        feedforward_dict = []

        models = PSI.SimulationModels(;
            decision_models = [
                uc_problem,
            ],
        );

    else
        @error "invalid stage combination"
    end

    # TODO: need to define reserve products for MD
    md_products = split(
        read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[
            1,
            "da_products",
        ],
        "; ",
    )
    rt_products = split(
        read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[
            1,
            "rt_products",
        ],
        "; ",
    )
    da_products = split(
        read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[
            1,
            "da_products",
        ],
        "; ",
    )
    energy_mkt_data = read_data(joinpath(simulation_dir, "markets_data", "Energy.csv"))

    energy_voll_cost = AxisArrays.AxisArray(energy_mkt_data.price_cap * 1.0, zones)

    sequence = create_sequence(models, feedforward_dict, single_stage_bool);

    if md_market_bool == true
        sim = PSI.Simulation(;
            name = "emis_$(case_name)",
            steps = PSY.get_forecast_interval(sys_MD), # sys_MD.data.time_series_params.forecast_params.count
            models = models,
            sequence = sequence,
            simulation_folder = ".",
            # initial_time = Dates.DateTime("2018-02-28T00:00:00")
        );
    else
        sim = PSI.Simulation(;
            name = "emis_$(case_name)",
            steps = 365, # sys_MD.data.time_series_params.forecast_params.count
            models = models,
            sequence = sequence,
            simulation_folder = ".",
            # initial_time = Dates.DateTime("2018-02-28T00:00:00")
        );
    end

    current_siip_sim[1] = sim
    push!(siip_system, sys_UC)

    @info "Building Sienna simulation..."
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

    @info "Executing Sienna simulation..."
    execute_out = PSI.execute!(sim; enable_progress_bar = true)
    @info "Simulation execution completed."

    @info "Getting Sienna results..."
    sim_results = PSI.SimulationResults(sim)

    res_uc = PSI.get_decision_problem_results(sim_results, "UC")
    dual_values_uc = PSI.read_realized_duals(res_uc)
    result_variables_uc = PSI.read_realized_variables(res_uc)
    data_length_uc =
        length(unique(dual_values_uc["CopperPlateBalanceConstraint__Area"].DateTime))

    result_variables_ed = nothing
    if !single_stage_bool
        res_ed = PSI.get_decision_problem_results(sim_results, "ED")
        dual_values_ed = PSI.read_realized_duals(res_ed)
        result_variables_ed = PSI.read_realized_variables(res_ed)
        data_length_ed =
            length(unique(dual_values_ed["CopperPlateBalanceConstraint__Area"].DateTime))
    end

    if md_market_bool
        res_md = PSI.get_decision_problem_results(sim_results, "MD")
        dual_values_md = PSI.read_realized_duals(res_md)
        result_variables_md = PSI.read_realized_variables(res_md)
        data_length_md =
            length(unique(dual_values_md["CopperPlateBalanceConstraint__Area"].DateTime))
    else
        result_variables_md = nothing
        data_length_md =
            length(unique(dual_values_uc["CopperPlateBalanceConstraint__Area"].DateTime))
    end

    energy_price_ed = AxisArrays.AxisArray(
        zeros(length(zones), 1, data_length_ed),
        zones,
        1:1,
        1:data_length_ed,
    )
    energy_price_uc = AxisArrays.AxisArray(
        zeros(length(zones), 1, data_length_uc),
        zones,
        1:1,
        1:data_length_uc,
    )
    energy_price_md = AxisArrays.AxisArray(
        zeros(length(zones), 1, data_length_md),
        zones,
        1:1,
        1:data_length_md,
    )
    energy_voll = AxisArrays.AxisArray(
        zeros(length(zones), 1, data_length_ed),
        zones,
        1:1,
        1:data_length_ed,
    )
    energy_voll_uc = AxisArrays.AxisArray(
        zeros(length(zones), 1, data_length_uc),
        zones,
        1:1,
        1:data_length_uc,
    )
    energy_voll_md = AxisArrays.AxisArray(
        zeros(length(zones), 1, data_length_md),
        zones,
        1:1,
        1:data_length_md,
    )

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
        # bus = find_zonal_bus(String(zone), sys_UC)
        area = find_zonal_area(String(zone), sys_UC)
        # zone_num = parse(Int64, last(zone, 1))
        if isnothing(zone)
            energy_price_ed[zone, 1, :] = zeros(data_length_ed)
            energy_price_uc[zone, 1, :] = zeros(data_length_uc)
            energy_price_md[zone, 1, :] = zeros(data_length_md)
            energy_voll[zone, 1, :] = zeros(data_length_ed)
            energy_voll_uc[zone, 1, :] = zeros(data_length_uc)
            energy_voll_md[zone, 1, :] = zeros(data_length_md)
        else
            if single_stage_bool == false
                energy_price_ed[zone, 1, :] =
                    abs.(
                        round.(
                            dual_values_ed["CopperPlateBalanceConstraint__Area"][
                                dual_values_ed["CopperPlateBalanceConstraint__Area"].name .== PSY.get_name(
                                    area,
                                ),
                                :value,
                            ];
                            digits = 5,
                        ),
                    ) / base_power
                energy_price_uc[zone, 1, :] =
                    abs.(
                        round.(
                            dual_values_uc["CopperPlateBalanceConstraint__Area"][
                                dual_values_uc["CopperPlateBalanceConstraint__Area"].name .== PSY.get_name(
                                    area,
                                ),
                                :value,
                            ];
                            digits = 5,
                        ),
                    ) / base_power
                if md_market_bool == true
                    energy_price_md[zone, 1, :] =
                        abs.(
                            round.(
                                dual_values_md["CopperPlateBalanceConstraint__Area"][
                                    dual_values_md["CopperPlateBalanceConstraint__Area"].name .== PSY.get_name(
                                        area,
                                    ),
                                    :value,
                                ];
                                digits = 5,
                            ),
                        ) / base_power
                    energy_voll_md[zone, 1, :] =
                        abs.(
                            round.(
                                result_variables_md["SystemBalanceSlackUp__Area"][
                                    result_variables_md["SystemBalanceSlackUp__Area"].name .== PSY.get_name(
                                        area,
                                    ),
                                    :value,
                                ];
                                digits = 5,
                            ),
                        ) / base_power
                else
                    energy_price_md[zone, 1, :] = zeros(data_length_md)
                end
                energy_voll[zone, 1, :] =
                    abs.(
                        round.(
                            result_variables_ed["SystemBalanceSlackUp__Area"][
                                result_variables_ed["SystemBalanceSlackUp__Area"].name .== PSY.get_name(
                                    area,
                                ),
                                :value,
                            ];
                            digits = 5,
                        ),
                    ) / base_power
                energy_voll_uc[zone, 1, :] =
                    abs.(
                        round.(
                            result_variables_uc["SystemBalanceSlackUp__Area"][
                                result_variables_uc["SystemBalanceSlackUp__Area"].name .== PSY.get_name(
                                    area,
                                ),
                                :value,
                            ];
                            digits = 5,
                        ),
                    ) / base_power
                # energy_voll[zone, 1, :] += abs.(round.(result_variables_ed["SystemBalanceSlackDown__ACBus"][:, string(PSY.get_number(bus))], digits = 5)) / base_power
            else
                energy_price_uc[zone, 1, :] =
                    abs.(
                        round.(
                            dual_values_uc["CopperPlateBalanceConstraint__Area"][
                                dual_values_uc["CopperPlateBalanceConstraint__Area"].name .== PSY.get_name(
                                    area,
                                ),
                                :value,
                            ];
                            digits = 5,
                        ),
                    ) / base_power
                energy_voll_uc[zone, 1, :] =
                    abs.(
                        round.(
                            result_variables_uc["SystemBalanceSlackUp__Area"][
                                result_variables_uc["SystemBalanceSlackUp__Area"].name .== PSY.get_name(
                                    area,
                                ),
                                :value,
                            ];
                            digits = 5,
                        ),
                    ) / base_power
            end
        end
    end

    @info "Recorded sienna energy prices"

    #println(any(isnan, energy_price))
    replace!(energy_price_ed, NaN => 0.0)
    replace!(energy_price_uc, NaN => 0.0)
    replace!(energy_price_md, NaN => 0.0)
    scale_voll(energy_price_ed, rt_resolution)

    #println(any(isnan, energy_price))
    #println(Statistics.mean(energy_price))

    if !single_stage_bool
        for service in get_system_services(sys_ED)
            name = PSY.get_name(service)
            @info "Recording prices for service $(name) - sys_ED"
            if typeof(service) == PSY.VariableReserve{PSY.ReserveUp}
                if name == "Inertia"
                    @info "Recording $(name) product - $(typeof(service)) - sys_ED"
                    dual_key = "RequirementConstraint__VariableReserve__ReserveUp__$(name)"
                    var_key = "ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)"
                    if haskey(dual_values_ed, dual_key) &&
                       haskey(result_variables_ed, var_key)
                        inertia_price[1, :] =
                            abs.(round.(dual_values_ed[dual_key][:, :value]; digits = 5)) /
                            base_power
                        replace!(inertia_price, NaN => 0.0)
                        scale_voll(inertia_price, rt_resolution)
                        inertia_voll[1, :] =
                            abs.(
                                round.(result_variables_ed[var_key][:, :value]; digits = 5),
                            ) / base_power
                    else
                        @warn "Inertia product found in system but: $(dual_key) $(haskey(dual_values_ed, dual_key)), $(var_key) $(haskey(result_variables_ed, var_key)) - sys_ED"
                    end

                    ##TODO: Check why we are not recording reserve price for clean energy product. 
                elseif name == "Clean_Energy"
                    @info "Not processing for Clean_Energy product - sys_ED"
                else
                    @info "Recording $(name) product - $(typeof(service)) - sys_ED"
                    dual_key = "RequirementConstraint__VariableReserve__ReserveUp__$(name)"
                    var_key = "ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)"
                    if haskey(dual_values_ed, dual_key) &&
                       haskey(result_variables_ed, var_key)
                        reserve_price_ed[name][1, :] =
                            abs.(round.(dual_values_ed[dual_key][:, :value]; digits = 5)) /
                            base_power
                        replace!(reserve_price_ed[name], NaN => 0.0)
                        scale_voll(reserve_price_ed[name], rt_resolution)
                        reserve_voll[name][1, :] =
                            abs.(
                                round.(result_variables_ed[var_key][:, :value]; digits = 5),
                            ) / base_power
                    else
                        @warn "Reserve product $(name) found in system but: $(dual_key) $(haskey(dual_values_ed, dual_key)), $(var_key) $(haskey(result_variables_ed, var_key)) - sys_ED"
                    end
                end
            elseif typeof(service) == PSY.ReserveDemandCurve{PSY.ReserveUp}
                @info "Recording $(name) product - $(typeof(service)) - sys_ED"
                dual_key = "RequirementConstraint__ReserveDemandCurve__ReserveUp__$(name)"
                if haskey(dual_values_ed, dual_key)
                    reserve_price_ed[name][1, :] =
                        abs.(round.(dual_values_ed[dual_key][:, :value]; digits = 5)) /
                        base_power
                    replace!(reserve_price_ed[name], NaN => 0.0)
                    scale_voll(reserve_price_ed[name], rt_resolution)
                else
                    @warn "Reserve demand curve product $(name) found in system but: $(dual_key) $(haskey(dual_values_ed, dual_key)) - sys_ED"
                end

            elseif typeof(service) == PSY.VariableReserve{PSY.ReserveDown}
                @info "Recording $(name) product - $(typeof(service)) - sys_ED"
                dual_key = "RequirementConstraint__VariableReserve__ReserveDown__$(name)"
                var_key = "ReserveRequirementSlack__VariableReserve__ReserveDown__$(name)"
                if haskey(dual_values_ed, dual_key) && haskey(result_variables_ed, var_key)
                    reserve_price_ed[name][1, :] =
                        abs.(round.(dual_values_ed[dual_key][:, :value]; digits = 5)) /
                        base_power
                    replace!(reserve_price_ed[name], NaN => 0.0)
                    scale_voll(reserve_price_ed[name], rt_resolution)
                    reserve_voll[name][1, :] =
                        abs.(round.(result_variables_ed[var_key][:, :value]; digits = 5)) /
                        base_power
                else
                    @warn "Reserve product $(name) found in system but: $(dual_key) $(haskey(dual_values_ed, dual_key)), $(var_key) $(haskey(result_variables_ed, var_key)) - sys_ED"
                end
            end
        end

        for service in get_system_services(sys_UC)
            name = PSY.get_name(service)
            # if name in only_da_products
            if typeof(service) == PSY.VariableReserve{PSY.ReserveUp}
                if name == "Inertia"
                    @info "Recording $(name) product - $(typeof(service)) - sys_UC"
                    dual_key = "RequirementConstraint__VariableReserve__ReserveUp__$(name)"
                    var_key = "ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)"
                    if haskey(dual_values_uc, dual_key) &&
                       haskey(result_variables_uc, var_key)
                        inertia_price[1, :] =
                            abs.(round.(dual_values_uc[dual_key][:, :value]; digits = 5)) /
                            base_power
                        replace!(inertia_price, NaN => 0.0)
                        scale_voll(inertia_price, da_resolution)
                        inertia_voll[1, :] =
                            abs.(
                                round.(result_variables_uc[var_key][:, :value]; digits = 5),
                            ) / base_power
                    else
                        @warn "Inertia product found in system but: $(dual_key) $(haskey(dual_values_uc, dual_key)), $(var_key) $(haskey(result_variables_uc, var_key)) - sys_UC"
                    end

                    ##TODO: Check why we are not recording reserve price for clean energy product.
                elseif name == "Clean_Energy"
                    @info "Not processing prices for Clean_Energy product - sys_UC"
                else
                    @info "Recording $(name) product - $(typeof(service)) - sys_UC"
                    dual_key = "RequirementConstraint__VariableReserve__ReserveUp__$(name)"
                    var_key = "ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)"
                    if haskey(dual_values_uc, dual_key) &&
                       haskey(result_variables_uc, var_key)
                        reserve_price_uc[name][1, :] =
                            abs.(round.(dual_values_uc[dual_key][:, :value]; digits = 5)) /
                            base_power
                        replace!(reserve_price_uc[name], NaN => 0.0)
                        scale_voll(reserve_price_uc[name], da_resolution)
                        reserve_voll_uc[name][1, :] =
                            abs.(
                                round.(result_variables_uc[var_key][:, :value]; digits = 5),
                            ) / base_power
                    else
                        @warn "Reserve product $(name) found in system but: $(dual_key) $(haskey(dual_values_uc, dual_key)), $(var_key) $(haskey(result_variables_uc, var_key)) - sys_UC"
                    end
                end
            elseif typeof(service) == PSY.ReserveDemandCurve{PSY.ReserveUp}
                @info "Recording $(name) product - $(typeof(service)) - sys_UC"
                dual_key = "RequirementConstraint__ReserveDemandCurve__ReserveUp__$(name)"
                if haskey(dual_values_uc, dual_key)
                    reserve_price_uc[name][1, :] =
                        abs.(round.(dual_values_uc[dual_key][:, :value]; digits = 5)) /
                        base_power
                    replace!(reserve_price_uc[name], NaN => 0.0)
                    scale_voll(reserve_price_uc[name], da_resolution)
                else
                    @warn "Reserve demand curve product $(name) found in system but: $(dual_key) $(haskey(dual_values_uc, dual_key)) - sys_UC"
                end
            elseif typeof(service) == PSY.VariableReserve{PSY.ReserveDown}
                @info "Recording $(name) product - $(typeof(service)) - sys_UC"
                dual_key = "RequirementConstraint__VariableReserve__ReserveDown__$(name)"
                var_key = "ReserveRequirementSlack__VariableReserve__ReserveDown__$(name)"
                if haskey(dual_values_uc, dual_key) && haskey(result_variables_uc, var_key)
                    reserve_price_uc[name][1, :] =
                        abs.(round.(dual_values_uc[dual_key][:, :value]; digits = 5)) /
                        base_power
                    replace!(reserve_price_uc[name], NaN => 0.0)
                    scale_voll(reserve_price_uc[name], da_resolution)
                    reserve_voll_uc[name][1, :] =
                        abs.(round.(result_variables_uc[var_key][:, :value]; digits = 5)) /
                        base_power
                else
                    @warn "Reserve product $(name) found in system but: $(dual_key) $(haskey(dual_values_uc, dual_key)), $(var_key) $(haskey(result_variables_uc, var_key)) - sys_UC"
                end
            end
            # end
        end

    else
        for service in get_system_services(sys_UC)
            name = PSY.get_name(service)
            # if name in only_da_products
            if typeof(service) == PSY.VariableReserve{PSY.ReserveUp}
                if name == "Inertia"
                    @info "Recording $(name) product - $(typeof(service)) - sys_UC"
                    dual_key = "RequirementConstraint__VariableReserve__ReserveUp__$(name)"
                    var_key = "ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)"
                    if haskey(dual_values_uc, dual_key) &&
                       haskey(result_variables_uc, var_key)
                        inertia_price[1, :] =
                            abs.(round.(dual_values_uc[dual_key][:, :value]; digits = 5)) /
                            base_power
                        replace!(inertia_price, NaN => 0.0)
                        scale_voll(inertia_price, da_resolution)
                        inertia_voll[1, :] =
                            abs.(
                                round.(result_variables_uc[var_key][:, :value]; digits = 5),
                            ) / base_power
                    else
                        @warn "Inertia product found in system but: $(dual_key) $(haskey(dual_values_uc, dual_key)), $(var_key) $(haskey(result_variables_uc, var_key)) - sys_UC"
                    end
                elseif name == "Clean_Energy"
                    @info "Not processing prices for Clean_Energy product - sys_UC"
                else
                    @info "Recording $(name) product - $(typeof(service)) - sys_UC"
                    dual_key = "RequirementConstraint__VariableReserve__ReserveUp__$(name)"
                    var_key = "ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)"
                    if haskey(dual_values_uc, dual_key) &&
                       haskey(result_variables_uc, var_key)
                        reserve_price_uc[name][1, :] =
                            abs.(round.(dual_values_uc[dual_key][:, :value]; digits = 5)) /
                            base_power
                        replace!(reserve_price_uc[name], NaN => 0.0)
                        scale_voll(reserve_price_uc[name], da_resolution)
                        reserve_voll_uc[name][1, :] =
                            abs.(
                                round.(result_variables_uc[var_key][:, :value]; digits = 5),
                            ) / base_power
                    else
                        @warn "Reserve product $(name) found in system but: $(dual_key) $(haskey(dual_values_uc, dual_key)), $(var_key) $(haskey(result_variables_uc, var_key)) - sys_UC"
                    end
                end
            elseif typeof(service) == PSY.ReserveDemandCurve{PSY.ReserveUp}
                @info "Recording $(name) product - $(typeof(service)) - sys_UC"
                dual_key = "RequirementConstraint__ReserveDemandCurve__ReserveUp__$(name)"
                if haskey(dual_values_uc, dual_key)
                    reserve_price_uc[name][1, :] =
                        abs.(round.(dual_values_uc[dual_key][:, :value]; digits = 5)) /
                        base_power
                    replace!(reserve_price_uc[name], NaN => 0.0)
                    scale_voll(reserve_price_uc[name], da_resolution)
                else
                    @warn "Reserve demand curve product $(name) found in system but: $(dual_key) $(haskey(dual_values_uc, dual_key)) - sys_UC"
                end

            elseif typeof(service) == PSY.VariableReserve{PSY.ReserveDown}
                @info "Recording $(name) product - $(typeof(service)) - sys_UC"
                dual_key = "RequirementConstraint__VariableReserve__ReserveDown__$(name)"
                var_key = "ReserveRequirementSlack__VariableReserve__ReserveDown__$(name)"
                if haskey(dual_values_uc, dual_key) && haskey(result_variables_uc, var_key)
                    reserve_price_uc[name][1, :] =
                        abs.(round.(dual_values_uc[dual_key][:, :value]; digits = 5)) /
                        base_power
                    replace!(reserve_price_uc[name], NaN => 0.0)
                    scale_voll(reserve_price_uc[name], da_resolution)
                    reserve_voll_uc[name][1, :] =
                        abs.(round.(result_variables_uc[var_key][:, :value]; digits = 5)) /
                        base_power
                else
                    @warn "Reserve product $(name) found in system but: $(dual_key) $(haskey(dual_values_uc, dual_key)), $(var_key) $(haskey(result_variables_uc, var_key)) - sys_UC"
                end
            end
            # end
        end
    end

    if md_market_bool == true
        for service in get_system_services(sys_MD)
            name = PSY.get_name(service)
            # TODO: need to replace "only_da_products"
            # if name in only_da_products
            if typeof(service) == PSY.VariableReserve{PSY.ReserveUp}
                if name == "Inertia"
                    @info "Recording $(name) product - $(typeof(service)) - sys_MD"
                    dual_key = "RequirementConstraint__VariableReserve__ReserveUp__$(name)"
                    var_key = "ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)"
                    if haskey(dual_values_md, dual_key) &&
                       haskey(result_variables_md, var_key)
                        inertia_price[1, :] =
                            abs.(round.(dual_values_md[dual_key][:, :value]; digits = 5)) /
                            base_power
                        replace!(inertia_price, NaN => 0.0)
                        scale_voll(inertia_price, da_resolution)
                        inertia_voll[1, :] =
                            abs.(
                                round.(result_variables_md[var_key][:, :value]; digits = 5),
                            ) / base_power
                    else
                        @warn "Inertia product found in system but: $(dual_key) $(haskey(dual_values_md, dual_key)), $(var_key) $(haskey(result_variables_md, var_key)) - sys_MD"
                    end

                elseif name == "Clean_Energy"
                    @info "Not processing prices for Clean_Energy product - sys_MD"
                else
                    dual_key = "RequirementConstraint__VariableReserve__ReserveUp__$(name)"
                    var_key = "ReserveRequirementSlack__VariableReserve__ReserveUp__$(name)"
                    if haskey(dual_values_md, dual_key) &&
                       haskey(result_variables_md, var_key)
                        reserve_price_md[name][1, :] =
                            abs.(round.(dual_values_md[dual_key][:, :value]; digits = 5)) /
                            base_power
                        replace!(reserve_price_md[name], NaN => 0.0)
                        scale_voll(reserve_price_md[name], da_resolution)
                        reserve_voll_md[name][1, :] =
                            abs.(
                                round.(result_variables_md[var_key][:, :value]; digits = 5),
                            ) / base_power
                    else
                        @warn "Reserve product found in system but: $(dual_key) $(haskey(dual_values_md, dual_key)), $(var_key) $(haskey(result_variables_md, var_key)) - sys_MD"
                    end
                end
            elseif typeof(service) == PSY.ReserveDemandCurve{PSY.ReserveUp}
                @info "Recording $(name) product - $(typeof(service)) - sys_MD"
                dual_key = "RequirementConstraint__ReserveDemandCurve__ReserveUp__$(name)"
                if haskey(dual_values_md, dual_key)
                    reserve_price_md[name][1, :] =
                        abs.(round.(dual_values_md[dual_key][:, :value]; digits = 5)) /
                        base_power
                    replace!(reserve_price_md[name], NaN => 0.0)
                    scale_voll(reserve_price_md[name], da_resolution)
                else
                    @warn "Reserve demand curve product $(name) found in system but: $(dual_key) $(haskey(dual_values_md, dual_key)) - sys_MD"
                end
            elseif typeof(service) == PSY.VariableReserve{PSY.ReserveDown}
                @info "Recording $(name) product - $(typeof(service)) - sys_MD"
                dual_key = "RequirementConstraint__VariableReserve__ReserveDown__$(name)"
                var_key = "ReserveRequirementSlack__VariableReserve__ReserveDown__$(name)"
                if haskey(dual_values_md, dual_key) && haskey(result_variables_md, var_key)
                    reserve_price_md[name][1, :] =
                        abs.(round.(dual_values_md[dual_key][:, :value]; digits = 5)) /
                        base_power
                    replace!(reserve_price_md[name], NaN => 0.0)
                    scale_voll(reserve_price_md[name], da_resolution)
                    reserve_voll_md[name][1, :] =
                        abs.(round.(result_variables_md[var_key][:, :value]; digits = 5)) /
                        base_power
                else
                    @warn "Reserve product $(name) found in system but: $(dual_key) $(haskey(dual_values_md, dual_key)), $(var_key) $(haskey(result_variables_md, var_key)) - sys_MD"
                end
            elseif typeof(service) == PSY.ReserveDemandCurve{PSY.ReserveUp}
                @info "Recording $(name) product - $(typeof(service)) - sys_MD"
                dual_key = "RequirementConstraint__ReserveDemandCurve__ReserveUp__$(name)"
                var_key = "ReserveRequirementSlack__ReserveDemandCurve__ReserveUp__$(name)"
                if haskey(dual_values_md, dual_key) && haskey(result_variables_md, var_key)
                    reserve_price_md[name][1, :] =
                        abs.(round.(dual_values_md[dual_key][:, :value]; digits = 5)) /
                        base_power
                    replace!(reserve_price_md[name], NaN => 0.0)
                    scale_voll(reserve_price_md[name], da_resolution)
                    reserve_voll_md[name][1, :] =
                        abs.(round.(result_variables_md[var_key][:, :value]; digits = 5)) /
                        base_power
                else
                    @warn "Reserve demand curve product $(name) found in system but: $(dual_key) $(haskey(dual_values_md, dual_key)), $(var_key) $(haskey(result_variables_md, var_key)) - sys_MD"
                end
            end
            # end
        end
    end

    @info "Recorded sienna service prices"

    sys_techs = get_all_techs(sys_ED)

    tech_names = get_name.(sys_techs)
    capacity_factors_md = Dict([g => zeros(1, data_length_md) for g in tech_names])
    capacity_factors_uc = Dict([g => zeros(1, data_length_uc) for g in tech_names])
    capacity_factors_ed = Dict([g => zeros(1, data_length_ed) for g in tech_names])
    start_up_costs = Dict([g => zeros(1, data_length_uc) for g in tech_names])
    shut_down_costs = Dict([g => zeros(1, data_length_uc) for g in tech_names])

    reserve_perc_md = Dict(
        g => Dict(
            s => zeros(1, data_length_md) for
            s in String.(union(rt_products, da_products, md_products))
        ) for g in tech_names
    )
    reserve_perc_uc = Dict(
        g => Dict(
            s => zeros(1, data_length_uc) for
            s in String.(union(rt_products, da_products, md_products))
        ) for g in tech_names
    )
    reserve_perc_ed = Dict(
        g => Dict(
            s => zeros(1, data_length_ed) for
            s in String.(union(rt_products, da_products, md_products))
        ) for g in tech_names
    )

    inertia_perc = Dict([g => zeros(1, data_length_ed) for g in tech_names])

    # for g in tech_names
    #     for product in only_da_products
    #         reserve_perc[g][string(product)] = zeros(1, data_length_uc)
    #     end
    # end

    for tech in sys_techs
        name = get_name(tech)
        @info "Recording results for technology $(name)"
        if md_market_bool
            capacity_factors_md[name][1, :] = get_realized_capacity_factors(
                tech,
                result_variables_md,
                result_variables_uc,
                base_power,
            )
        end
        if !single_stage_bool
            capacity_factors_ed[name][1, :] = get_realized_capacity_factors(
                tech,
                result_variables_ed,
                result_variables_ed,
                base_power,
            )
        end
        capacity_factors_uc[name][1, :] = get_realized_capacity_factors(
            tech,
            result_variables_uc,
            result_variables_uc,
            base_power,
        )
        start_up_costs[name][1, :] =
            get_start_costs(tech, result_variables_uc, data_length_uc)
        shut_down_costs[name][1, :] =
            get_shut_costs(tech, result_variables_uc, data_length_uc)

        services_ED = PSY.get_services(tech)
        services_UC =
            PSY.get_services(PSY.get_components_by_name(PSY.Device, sys_UC, name)[1])

        for service in services_UC
            # TODO: figure out if we need clean energy or not
            if PSY.get_name(service) != "Clean_Energy"
                @info "Updating reserve percentage for technology $(name) - service $(PSY.get_name(service)) - sys_UC"
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
                    md_market_bool,
                    single_stage_bool)
            end
        end
    end

    for g in keys(inertia_perc)
        inertia_perc[g] = inertia_perc[g] / PSY.get_base_power(sys_UC)
    end

    @info "Recorded all Sienna results"

    return energy_price_ed,
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
    inertia_voll;
end
