
function add_inertia_constant!(device::PSY.Device, product::T) where {T <: Product}
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
    gen_name = get_name(gen)

    buses = PSY.get_components(PSY.Bus, sys)
    bus = filter(b -> string(PSY.get_number(b)) == get_bus(tech), collect(buses))

    if isempty(bus)
        bus = filter(b -> PSY.get_name(b) == get_bus(tech), collect(buses))
        if isempty(bus)
            error("No matching bus found for generator $(gen_name) with bus name $(get_bus(tech))")
        else
            gen_bus = only(bus)
        end
    else
        gen_bus = only(bus)
    end

    type = deepcopy(get_type(tech))
    if type == "RE_CT"
        type = "CT"
    end

    PSY_gen = PSY.ThermalStandard(
        get_name(gen), # name
        true,   # available
        true,   # status
        gen_bus, # bus
        get_maxcap(gen) / base_power,    # active power
        1.0,                # reactive power
        get_maxcap(gen) / base_power,    # rating
        (min = get_mincap(gen) / base_power, max = get_maxcap(gen) / base_power), # active power limits
        nothing,      # reactive power limits
        (
            up = get_ramp_limits(tech)[:up] / (base_power * 60),
            down = get_ramp_limits(tech)[:down] / (base_power * 60),
        ), # ramp limits
        get_operation_cost(tech), # operation cost
        base_power, # base power
        get_time_limits(tech), # up and down time limits
        false, # must run
        PSY.PrimeMovers(
            findfirst(x -> Symbol(x) == Symbol(type), collect(instances(PSY.PrimeMovers))),
        ), # primemover
        PSY.ThermalFuels(
            findfirst(
                x -> Symbol(x) == Symbol(get_fuel(tech)),
                collect(instances(PSY.ThermalFuels)),
            ),
        ), # fuel type
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
    gen_name = get_name(gen)

    if get_type(tech) == "WT"
        primemover = PSY.PrimeMovers.WT
    elseif get_type(tech) == "PVe"
        primemover = PSY.PrimeMovers.PVe
    end

    buses = PSY.get_components(PSY.Bus, sys)
    bus = filter(b -> string(PSY.get_number(b)) == get_bus(tech), collect(buses))

    if isempty(bus)
        bus = filter(b -> PSY.get_name(b) == get_bus(tech), collect(buses))
        if isempty(bus)
            error("No matching bus found for generator $(gen_name) with bus name $(get_bus(tech))")
        else
            gen_bus = only(bus)
        end
    else
        gen_bus = only(bus)
    end

    PSY_gen = PSY.RenewableDispatch(
        gen_name,  # name
        true,           # available
        gen_bus, # bus
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
This function creates a PowerSystems EnergyReservoirStorage unit.
EnergyReservoirStorage in PSY: https://nrel-sienna.github.io/PowerSystems.jl/stable/model_library/generated_EnergyReservoirStorage/#EnergyReservoirStorage
"""
function create_PSY_generator(gen::BatteryEMIS{<: BuildPhase}, sys::PSY.System)
    # @info "Creating PSY EnergyReservoirStorage for battery $(get_name(gen))"
    tech = get_tech(gen)
    base_power = get_base_power(gen)
    gen_name = get_name(gen)
    rating = get_storage_rating(gen)
    maxcap = get_maxcap(gen)
    storage_capacity = get_storage_capacity(gen)
    storage_level_limits = get_storage_level_limits(gen)
    input_active_power_limits = get_input_active_power_limits(gen)
    output_active_power_limits = get_output_active_power_limits(gen)
    initial_storage_capacity_level = get_initial_storage_capacity_level(gen)
    efficiency = get_efficiency(gen)

    buses = PSY.get_components(PSY.Bus, sys)
    bus = filter(b -> string(PSY.get_number(b)) == get_bus(tech), collect(buses))

    if isempty(bus)
        bus = filter(b -> PSY.get_name(b) == get_bus(tech), collect(buses))
        if isempty(bus)
            error("No matching bus found for generator $(gen_name) with bus name $(get_bus(tech))")
        else
            gen_bus = only(bus)
        end
    else
        gen_bus = only(bus)
    end

    PSY_gen = PSY.EnergyReservoirStorage(
        gen_name,  # name
        true,           # available
        gen_bus,        # bus
        PSY.PrimeMovers.BA, # primemover
        StorageTech.LIB, # storage technology
        storage_capacity[:max], # storage capacity MWh
        storage_level_limits, # state of charge limits
        initial_storage_capacity_level, # initial state of charge
        rating, # rating
        maxcap / base_power, # active power
        (
            min = input_active_power_limits[:min],
            max = input_active_power_limits[:max],
        ), # input active power limits
        (
            min = output_active_power_limits[:min],
            max = output_active_power_limits[:max],
        ), # output active power limits
        efficiency, # in/out efficiency
        1.0,             # reactive power
        nothing,      # reactive power limits
        base_power, # base power
        PSY.StorageCost(),    #operation_cost
        1.0, # conversion factor
        0.0, # storage_target
        1e4, # cycle_limits
        PSY.Device[], # services
        nothing, # dynamic_injector
    )
    for product in get_products(gen)
        add_inertia_constant!(PSY_gen, product)
    end

    add_outage_info!(PSY_gen, tech)

    return PSY_gen
end
