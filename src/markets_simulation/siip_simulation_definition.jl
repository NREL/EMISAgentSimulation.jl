# This file contains functions for calling, running and post-processing
# SIIP PSI Simulation for actual Energy and A/S Market Clearing

"""
This function returns realized capacity factors for Thermal generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.ThermalStandard,
                                        results::Dict{Symbol, DataFrames.DataFrame}
                                        )
    energy_production = results[:P__ThermalStandard][:, Symbol(get_name(device))]
    capacity_factors = energy_production / PSY.get_active_power_limits(device)[:max]
    return capacity_factors
end

"""
This function returns realized capacity factors for Renewable generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.RenewableDispatch,
                                        results::Dict{Symbol, DataFrames.DataFrame}
                                        )
    energy_production = results[:P__RenewableDispatch][:, Symbol(get_name(device))]
    capacity_factors = energy_production / get_device_size(device)
    return capacity_factors
end

"""
This function returns realized capacity factors for Hydropower generators from PSI Simulation.
"""
function get_realized_capacity_factors(device::PSY.HydroEnergyReservoir,
                                        results::Dict{Symbol, DataFrames.DataFrame}
                                        )
    energy_production = results[:P__HydroEnergyReservoir][:, Symbol(get_name(device))]
    capacity_factors = energy_production / get_device_size(device)
    return capacity_factors
end

"""
This function returns nothing if Service is not of ReserveUp or ReserveDown type.
"""
function update_realized_reserve_perc!(device::PSY.Device,
                                        service::S,
                                        results::Dict{Symbol, DataFrames.DataFrame},
                                        reserve_up_perc::Dict{String, Array{Float64, 2}},
                                        reserve_down_perc::Dict{String, Array{Float64, 2}},
                                        data_length::Int64) where S <: PSY.Service
    return
end

"""
This function returns realized reserve up provision percentages from PSI Simulation.
"""
function update_realized_reserve_perc!(device::PSY.Device,
                                        service::PSY.VariableReserve{PSY.ReserveUp},
                                        results::Dict{Symbol, DataFrames.DataFrame},
                                        reserve_up_perc::Dict{String, Array{Float64, 2}},
                                        reserve_down_perc::Dict{String, Array{Float64, 2}})

    reserve_provision = results[Symbol("$(PSY.get_name(service))__VariableReserve_PowerSystems.ReserveUp")][:, Symbol(get_name(device))]
    reserve_perc = reserve_provision / get_device_size(device)
    reserve_up_perc[get_name(device)][1, :] = reserve_perc
    return
end

"""
This function returns realized reserve down provision percentages from PSI Simulation.
"""
function update_realized_reserve_perc!(device::PSY.Device,
                                        service::PSY.VariableReserve{PSY.ReserveDown},
                                        results::Dict{Symbol, DataFrames.DataFrame},
                                        reserve_up_perc::Dict{String, Array{Float64, 2}},
                                        reserve_down_perc::Dict{String, Array{Float64, 2}})

    reserve_provision = results[Symbol("$(PSY.get_name(service))__VariableReserve_PowerSystems.ReserveDown")][:, Symbol(get_name(device))]
    reserve_perc = reserve_provision / get_device_size(device)
    reserve_down_perc[get_name(device)][1, :] = reserve_perc
    return
end

"""
This function creates the Unit Commitment template for PSI Simulation.
"""
function create_uc_template()
    service = Dict(
        :ReserveUp =>
            PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveUp}, PSI.RangeReserve),
        :ReserveDown =>
            PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveDown}, PSI.RangeReserve),
    )
    devices = Dict(
        :Generators => PSI.DeviceModel(PSY.ThermalStandard, PSI.ThermalStandardUnitCommitment),
        :Ren => PSI.DeviceModel(PSY.RenewableDispatch, PSI.RenewableFullDispatch),
        :Loads => PSI.DeviceModel(PSY.PowerLoad, PSI.StaticPowerLoad),
        :HydroROR => PSI.DeviceModel(PSY.HydroEnergyReservoir, PSI.HydroDispatchRunOfRiver),
        :RenFx => PSI.DeviceModel(PSY.RenewableFix, PSI.FixedOutput),
        :Batt => PSI.DeviceModel(PSY.GenericBattery, PSI.BookKeepingwReservation),
    )
    template_uc =
        PSI.template_unit_commitment(network = PSI.PM.NFAPowerModel, devices = devices, services = service)
    return template_uc
end

"""
This function creates the Economic Dispatch template for PSI Simulation.
"""
function create_ed_template()
    service = Dict(
        :ReserveUp =>
            PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveUp}, PSI.RangeReserve),
        :ReserveDown =>
            PSI.ServiceModel(PSY.VariableReserve{PSY.ReserveDown}, PSI.RangeReserve),
    )
    devices = Dict(
        :Generators => PSI.DeviceModel(PSY.ThermalStandard, PSI.ThermalDispatch),
        :Ren => PSI.DeviceModel(PSY.RenewableDispatch, PSI.RenewableFullDispatch),
        :Loads => PSI.DeviceModel(PSY.PowerLoad, PSI.StaticPowerLoad),
        :HydroROR => PSI.DeviceModel(PSY.HydroEnergyReservoir, PSI.HydroDispatchRunOfRiver),
        :RenFx => PSI.DeviceModel(PSY.RenewableFix, PSI.FixedOutput),
        :Batt => PSI.DeviceModel(PSY.GenericBattery, PSI.BookKeepingwReservation),
    )
    template_ed =
        PSI.template_economic_dispatch(network = PSI.PM.NFAPowerModel, devices = devices, services = service)

    return template_ed
end

"""
This function creates the Stages for PSI Simulation.
"""
function create_stage(template::PSI.OperationsProblemTemplate, sys::PSY.System, solver::JuMP.MOI.OptimizerWithAttributes)
    duals = [:nodal_balance_active__Bus]

    services = get_system_services(sys)
    if in(PSY.VariableReserve{PSY.ReserveUp}, typeof.(services))
        push!(duals, Symbol("requirement__VariableReserve_PowerSystems.ReserveUp"))
    end
    if in(PSY.VariableReserve{PSY.ReserveDown}, typeof.(services))
        push!(duals, Symbol("requirement__VariableReserve_PowerSystems.ReserveDown"))
    end

    stage = PSI.Stage(
        PSI.GenericOpProblem,
        template,
        sys,
        solver;
        constraint_duals = duals)

    return stage
end

"""
This function creates the Sequences for PSI Simulation.
"""
function create_sequence(stage_1::AbstractString, stage_2::AbstractString)
    sequence = PSI.SimulationSequence(
        step_resolution = Dates.Hour(24),
        order = Dict(1 => "UC", 2 => "ED"),
        feedforward_chronologies = Dict(("UC" => "ED") => PSI.Synchronize(periods = 24)),
        horizons = Dict("UC" => 24, "ED" => 1),
        intervals = Dict(
            "UC" => (Dates.Hour(24), PSI.Consecutive()),
            "ED" => (Dates.Hour(1), PSI.Consecutive()),
        ),
        feedforward = Dict(
            ("ED", :devices, :Generators) => PSI.SemiContinuousFF(
                binary_source_stage = PSI.ON,
                affected_variables = [PSI.ACTIVE_POWER],
            ),
        ),
        ini_cond_chronology = PSI.InterStageChronology(),
    )
    return sequence
end


"""
This function creates the PSI Simulation and post-processes the results.
"""
function create_simulation( sys_UC::PSY.System,
                            sys_ED::PSY.System,
                            zones::Vector{String},
                            days::Int64,
                            solver::JuMP.MOI.OptimizerWithAttributes;
                            kwargs...)


    template_uc = create_uc_template()
    template_ed = create_ed_template()

    uc_stage = create_stage(template_uc, sys_UC, solver)
    ed_stage = create_stage(template_ed, sys_UC, solver)



    stages_definition = Dict(
        "UC" => uc_stage,
        "ED" => ed_stage,
    )

    sequence = create_sequence("UC", "ED")

    sim = PSI.Simulation(
        name = "aggregation",
        steps = days,
        stages = stages_definition,
        stages_sequence = sequence,
        simulation_folder = "./",
    )
    PSI.build!(sim)
    model = sim.stages["UC"].internal.psi_container.JuMPmodel

    open("SIIP_UC_Model.txt", "w") do f

        println(f, "UC Model:")
        println(f, model)
    end

    sim_results = PSI.execute!(sim; );

    base_power = PSY.get_base_power(sys_UC)

    res_uc = PSI.load_simulation_results(sim_results, "UC")
    res_ed = PSI.load_simulation_results(sim_results, "ED")

    dual_values = res_ed.dual_values

    duals = keys(dual_values)
    data_length = DataFrames.nrow(dual_values[:dual_nodal_balance_active__Bus])


    energy_price = AxisArrays.AxisArray(zeros(length(zones), 1, data_length), zones, 1:1, 1:data_length)
    reserve_up_price = AxisArrays.AxisArray(zeros(length(zones), 1, data_length), zones, 1:1, 1:data_length)
    reserve_down_price = AxisArrays.AxisArray(zeros(length(zones), 1, data_length), zones, 1:1, 1:data_length)


    for zone in zones
            bus = find_zonal_bus(zone, sys_UC)
            zone_num = parse(Int64, last(zone, 1))
            if isnothing(bus)
                energy_price[zone, 1, :] = zeros(data_length)
            else
                energy_price[zone, 1, :] =
                abs.(round.(dual_values[:dual_nodal_balance_active__Bus][:, PSY.get_name(bus)], digits = 5)) / base_power
            end
            if in(Symbol("dual_requirement__VariableReserve_PowerSystems.ReserveUp"), keys(dual_values))
                reserve_up_price[zone, 1, :] =
                abs.(round.(dual_values[Symbol("dual_requirement__VariableReserve_PowerSystems.ReserveUp")][:, Symbol("Reg_Up_R$(zone_num)")], digits = 5)) / base_power
            end
            if in(Symbol("dual_requirement__VariableReserve_PowerSystems.ReserveDown"), keys(dual_values))
                reserve_down_price[zone, 1, :] =
                abs.(round.(dual_values[Symbol("dual_requirement__VariableReserve_PowerSystems.ReserveDown")][:, Symbol("Reg_Down_R$(zone_num)")], digits = 5)) / base_power
            end
    end

    sys_techs = get_all_techs(sys_UC)
    tech_names = get_name.(sys_techs)
    capacity_factors = Dict([g => zeros(1, data_length) for g in tech_names])
    reserve_up_perc = Dict([g => zeros((1, data_length)) for g in tech_names])
    reserve_down_perc = Dict([g => zeros((1, data_length)) for g in tech_names])

    for tech in sys_techs
        name = get_name(tech)
        services = PSY.get_services(tech)
        capacity_factors[name][1, :] = get_realized_capacity_factors(tech, res_ed.variable_values)

        for service in services
            update_realized_reserve_perc!(tech,
                                         service,
                                         res_ed.variable_values,
                                         reserve_up_perc,
                                         reserve_down_perc)
        end
    end

    return energy_price, reserve_up_price, reserve_down_price, capacity_factors, reserve_up_perc, reserve_down_perc;
end
