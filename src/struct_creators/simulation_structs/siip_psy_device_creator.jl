
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
        PSY.PrimeMovers.Symbol(get_type(tech)), # primemover
        PSY.ThermalFuels.Symbol(get_fuel(tech)), # fuel type
    )
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
        PSY.PrimeMovers.Symbol(get_type(tech)), # primemover
        get_soc(tech) / (get_maxcap(gen) * base_power), #initial state of charge
        (min = get_storage_capacity(tech)[:min] / (get_maxcap(gen) * base_power), max = get_storage_capacity(tech)[:max] / (get_maxcap(gen) * base_power)), # state of charge limits
        get_maxcap(gen) / base_power, # rating
        get_maxcap(gen) / base_power, # active power
        (min = get_input_active_power_limits(tech)[:min] / base_power, max = get_input_active_power_limits(tech)[:max] / base_power), # input active power limits
        (min = get_output_active_power_limits(tech)[:min] / base_power, max = get_output_active_power_limits(tech)[:max] / base_power), # output active power limits
        get_efficiency(tech), # efficiency
        1.0,             # reactive power
        nothing,      # reactive power limits
        base_power, # base power
    )
    return PSY_gen
end
