
function add_inertia_constant!(device::PSY.Device, product::T) where T <: Product
    return
end

function add_inertia_constant!(device::PSY.Device, product::Inertia)
    device.ext["inertia"] = get_h_constant(product)
    return
end

"""
This function creates a PowerSystems ThermalStandard unit.
"""
function create_PSY_generator(gen::ThermalGenEMIS{<: BuildPhase}, sys::PSY.System)

    tech = get_tech(gen)
    base_power = get_maxcap(gen)

    buses = PSY.get_components(PSY.Bus, sys)
    project_bus = PSY.Bus[]

    for b in buses
        if PSY.get_number(b) == get_bus(tech)
            push!(project_bus, b)
        end
    end

    type = deepcopy(get_type(tech))
    if type == "RE_CT"
        type = "CT"
    end

    PSY_gen =  PSY.ThermalStandard(
        get_name(gen), # name
        true,   # available
        true,   # status
        project_bus[1], # bus
        get_maxcap(gen) / base_power,    # active power
        1.0,                # reactive power
        get_maxcap(gen) / base_power,    # rating
        (min = get_mincap(gen) / base_power, max = get_maxcap(gen) / base_power), # active power limits
        nothing,      # reactive power limits
        (up = get_ramp_limits(tech)[:up] / (base_power * 60), down = get_ramp_limits(tech)[:down] / (base_power * 60)), # ramp limits
        get_operation_cost(tech), # operation cost
        base_power, # base power
        get_time_limits(tech), # up and down time limits
        PSY.PrimeMovers(findfirst(x -> Symbol(x) == Symbol(type), collect(instances(PSY.PrimeMovers)))), # primemover
        PSY.ThermalFuels(findfirst(x -> Symbol(x) == Symbol(get_fuel(tech)), collect(instances(PSY.ThermalFuels)))), # fuel type
    )
    for product in get_products(gen)
        add_inertia_constant!(PSY_gen, product)
    end

    add_outage_info!(PSY_gen, tech)

    return PSY_gen
end

"""
This function creates a PowerSystems RenewableDispatch unit.
"""
function create_PSY_generator(gen::RenewableGenEMIS{<: BuildPhase}, sys::PSY.System)
    tech = get_tech(gen)
    base_power = get_maxcap(gen)

    if get_type(tech) == "WT"
        primemover = PSY.PrimeMovers.WT
    elseif get_type(tech) == "PVe"
        primemover = PSY.PrimeMovers.PVe
    end

    buses = PSY.get_components(PSY.Bus, sys)
    project_bus = PSY.Bus[]

    for b in buses
        if PSY.get_number(b) == get_bus(tech)
            push!(project_bus, b)
        end
    end

    PSY_gen =  PSY.RenewableDispatch(
        get_name(gen),  # name
        true,           # available
        project_bus[1], # bus
        get_maxcap(gen) / base_power, # active power
        1.0,             # reactive power
        get_maxcap(gen) / base_power,    # rating
        primemover,     # primemover
        nothing,        # reactivepower limits
        1.0,            # power factor
        get_operation_cost(tech),
        base_power, # base power
    )
    for product in get_products(gen)
        add_inertia_constant!(PSY_gen, product)
    end

    add_outage_info!(PSY_gen, tech)

    return PSY_gen
end

"""
This function creates a PowerSystems GenericBattery unit.
"""
function create_PSY_generator(gen::BatteryEMIS{<: BuildPhase}, sys::PSY.System)
    tech = get_tech(gen)
    base_power = get_maxcap(gen)

    buses = PSY.get_components(PSY.Bus, sys)
    project_bus = PSY.Bus[]

    for b in buses
        if PSY.get_number(b) == get_bus(tech)
            push!(project_bus, b)
        end
    end

    PSY_gen =  PSY.GenericBattery(
        get_name(gen),  # name
        true,           # available
        project_bus[1], # bus
        PSY.PrimeMovers.BA, # primemover
        get_soc(tech) / (base_power), #initial state of charge
        (min = get_storage_capacity(tech)[:min] / base_power, max = get_storage_capacity(tech)[:max] / base_power), # state of charge limits
        get_maxcap(gen) / base_power, # rating
        get_maxcap(gen) / base_power, # active power
        (min = get_input_active_power_limits(tech)[:min] / base_power, max = get_input_active_power_limits(tech)[:max] / base_power), # input active power limits
        (min = get_output_active_power_limits(tech)[:min] / base_power, max = get_output_active_power_limits(tech)[:max] / base_power), # output active power limits
        get_efficiency(tech), # efficiency
        1.0,             # reactive power
        nothing,      # reactive power limits
        base_power, # base power
        nothing,    #operation_cost
    )
    for product in get_products(gen)
        add_inertia_constant!(PSY_gen, product)
    end

    add_outage_info!(PSY_gen, tech)

    return PSY_gen
end
