# This file contains functions for calling, running and post-processing
# SIIP PSI Simulation for actual Energy and A/S Market Clearing

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
This function returns realized capacity factors for Thermal generators from PSI Simulation.
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
    println(maximum(generation))
    println(sum(generation) * get_device_size(device) * 100 / 6)

    energy_production_uc = results_uc[:Pout__GenericBattery][:, Symbol(get_name(device))] - results_uc[:Pin__GenericBattery][:, Symbol(get_name(device))]
    capacity_factors_uc = energy_production_uc / get_device_size(device)
    generation_uc = filter(x -> x > 0, capacity_factors_uc)
    println(maximum(generation_uc))
    println(sum(generation_uc) * get_device_size(device) * 2 * 100)
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
        PSI.set_device_model!(template, PSY.ThermalStandard, PSIE.ThermalInertiaBasicUnitCommitment)
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSIE.RenewableFullDispatchInertia)
        PSI.set_device_model!(template, PSY.RenewableFix, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.PowerLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroEnergyReservoir, PSIE.HydroInertiaCommitmentRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, PSIE.HydroInertiaCommitmentRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(template, PSY.GenericBattery, PSIE.BookKeepingwInertia)
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.HVDCLine, PSI.HVDCLossless)
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveUp}, PSI.RangeReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveDown}, PSI.RangeReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.ReserveDemandCurve{PSY.ReserveUp}, PSIE.QuadraticCostRampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveSymmetric}, PSIE.InertiaReserve))
    else
        template = PSI.OperationsProblemTemplate(PSI.PM.NFAPowerModel)
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalStandardUnitCommitment)
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
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveUp}, PSI.RangeReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveDown}, PSI.RangeReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.ReserveDemandCurve{PSY.ReserveUp}, PSIE.QuadraticCostRampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveSymmetric}, PSIE.InertiaReserve))
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
        PSI.set_device_model!(template, PSY.RenewableDispatch, PSIE.RenewableFullDispatchInertia)
        PSI.set_device_model!(template, PSY.RenewableFix, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.PowerLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroEnergyReservoir, PSI.HydroDispatchRunOfRiver)
        PSI.set_device_model!(template, PSY.HydroDispatch, PSI.HydroDispatchRunOfRiver) # TODO: check which hydro device we have
        PSI.set_device_model!(template, PSY.GenericBattery, PSIE.BookKeepingwInertia)
        PSI.set_device_model!(template, PSY.Line, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.Transformer2W, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.TapTransformer, PSI.StaticBranch)
        PSI.set_device_model!(template, PSY.HVDCLine, PSI.HVDCLossless)
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveUp}, PSIE.RampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveDown}, PSIE.RampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.ReserveDemandCurve{PSY.ReserveUp}, PSIE.QuadraticCostRampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveSymmetric}, PSIE.InertiaReserve))

    else
        template = PSI.OperationsProblemTemplate(PSI.PM.NFAPowerModel)
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalDispatch)
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
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveUp}, PSIE.RampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveDown}, PSIE.RampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.ReserveDemandCurve{PSY.ReserveUp}, PSIE.QuadraticCostRampReserve))
        PSI.set_service_model!(template, PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveSymmetric}, PSIE.InertiaReserve))
    end

    return template
end

"""
This function creates the Problem for PSI Simulation.
"""

function create_problem(template::PSI.OperationsProblemTemplate, sys::PSY.System, type::String, solver::JuMP.MOI.OptimizerWithAttributes)

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
                                    optimizer_log_print = true,
                                    balance_slack_variables = true,
                                    constraint_duals = duals,
                                    warm_start = true,
                                    services_slack_variables = true,
                                    )

        elseif type == "ED"
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

    uc_problem = create_problem(template_uc, sys_UC, "UC", solver)
    ed_problem = create_problem(template_ed, sys_ED, "ED", solver)

    if isempty(inertia_product)
        feedforward_dict = Dict(
            ("ED", :devices, Symbol("PowerSystems.ThermalStandard")) => PSI.SemiContinuousFF(
                binary_source_problem = PSI.ON,
                affected_variables = [PSI.ACTIVE_POWER],
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
            ("ED", :devices, Symbol("PowerSystems.HydroDispatch")) => PSIE.InertiaFF(
                binary_source_problem = PSI.ON,
                affected_variables = [PSI.ACTIVE_POWER],
                service = inertia_product[1]
            ),
            ("ED", :devices, Symbol("PowerSystems.GenericBattery")) => PSIE.EnergyTargetFF(
                variable_source_problem = PSI.ENERGY,
                affected_variables = [PSI.ENERGY],
                target_period = 12,
                penalty_cost = 1e5,
            ),
         )

    end



    problems = PSI.SimulationProblems(
                                    UC = uc_problem,
                                    ED = ed_problem
                                      )

    sequence = create_sequence(problems, da_resolution, feedforward_dict)

    sim = PSI.Simulation(
                    name = "emis_$(case_name)",
                    steps = 30,
                    problems = problems,
                    sequence = sequence,
                    simulation_folder = ".",
                    )

    build_out = PSI.build!(sim; serialize = false)
    execute_out = PSI.execute!(sim)

    base_power = PSY.get_base_power(sys_UC)
    sim_results = PSI.SimulationResults(sim);

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

    rt_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"rt_products"], "; ")
    da_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"da_products"], "; ")

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
