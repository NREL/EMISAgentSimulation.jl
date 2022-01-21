"""
This function finds the hourly availability data for projects
to be passed to expected and actual market clearing modules.
"""
function find_project_availability_data(project::P,
                                     availability_df::DataFrames.DataFrame,
                                     num_invperiods::Int64,
                                     num_hours::Int64) where P <: Project{<:BuildPhase}
    type = get_type(get_tech(project))
    zone = get_zone(get_tech(project))
        if in(get_name(project), names(availability_df))
            availability_raw = availability_df[:, Symbol(get_name(project))]
        elseif in("$(type)_$(zone)", names(availability_df))
            availability_raw = availability_df[:, Symbol("$(type)_$(zone)")]
        else
            error("project availiability data not found!")
        end

        data_years = unique(availability_df.Year)

        data_hours = count(i -> i == data_years[1], availability_df.Year)

        @assert data_hours == num_hours

        availability_input = zeros(num_invperiods, data_hours)

        # Fill raw availability data into yearly and hourly values
        for y in 1:num_invperiods
            availability_input[y, 1:data_hours] = availability_raw[1:data_hours]
        end

        return availability_input
end

"""
This function extracts the technical details for generators
to be passed to expected and actual market clearing modules.
"""
function get_technical_details(project::P) where P <: GeneratorEMIS{<: BuildPhase}

    project_type = "generator"
    min_input = 0.0
    max_input = 0.0
    efficiency_in = 0.0
    efficiency_out = 0.0
    min_storage = 0.0
    max_storage = 0.0
    init_storage = 0.0

    return project_type, min_input, max_input, efficiency_in, efficiency_out, min_storage, max_storage, init_storage
end

"""
This function extracts the technical details for storage units
to be passed to expected and actual market clearing modules.
"""
function get_technical_details(project::P) where P <: StorageEMIS{<: BuildPhase}

    tech = get_tech(project)

    project_type = "storage"

    input_power_limits = get_input_active_power_limits(tech)
    min_input = input_power_limits[:min]
    max_input = input_power_limits[:max]

    efficiency = get_efficiency(tech)
    efficiency_in = efficiency[:in]
    efficiency_out = efficiency[:out]

    storage_capacity = get_storage_capacity(tech)
    min_storage = storage_capacity[:min]
    max_storage = storage_capacity[:max]

    init_storage = get_soc(tech)

    return project_type, min_input, max_input, efficiency_in, efficiency_out, min_storage, max_storage, init_storage
end


function get_marginal_cost_energy(product::T) where T <: Product
    return
end

function get_marginal_cost_energy(product::Energy)
    marginal_cost_energy = get_marginal_cost(product)
    return marginal_cost_energy
end

"""
This function returns the project's marginal cost of energy product to be passed to economic dispatch and CEM.
"""
function get_project_energy_cost(project::P) where P <: Project{<: BuildPhase}
    marginal_cost_energy = 0.0
    for product in find_operating_products(get_products(project))
        marginal_cost_temp = get_marginal_cost_energy(product)
        if !isnothing(marginal_cost_temp)
            marginal_cost_energy = marginal_cost_temp
        end
    end
    return marginal_cost_energy
end

"""
This function returns the project's marginal cost of reeserve product to be passed to economic dispatch and CEM.
"""
function get_project_reserve_cost(product::Product, project_reserve_cost::Dict{String, Float64})
    return project_reserve_cost
end

function get_project_reserve_cost(product::P, project_reserve_cost::Dict{String, Float64}) where P <: OperatingReserve
    project_reserve_cost[String(get_name(product))] = get_marginal_cost(product)
    return project_reserve_cost
end

"""
This function returns the project's maximum participation in reeserve product to be passed to economic dispatch and CEM.
"""
function get_project_reserve_limit(product::Product, max_cap::Float64, project_reserve_limit::Dict{String, Float64})
    return project_reserve_limit
end

function get_project_reserve_limit(product::P, max_cap::Float64, project_reserve_limit::Dict{String, Float64}) where P <: OperatingReserve
    project_reserve_limit[String(get_name(product))] = get_max_limit(product) * max_cap
    return project_reserve_limit
end

"""
This function returns the project's derating factor to be passed to CEM and capacity market clearing module.
Returns 0 if there is no capacity market participation.
"""
function get_project_derating(project::P) where P <: Project{<: BuildPhase}
    derating_factor = 0.
    for product in get_products(project)
        derating_temp = get_derating(product)
        if !isnothing(derating_temp)
            derating_factor = derating_temp
        end
    end
    return derating_factor
end

"""
This function returns the project's capacity market bid to be passed to the capacity market clearing module.
Returns 0 if there is no capacity market participation.
"""
function get_project_capacity_market_bid(project::P) where P <: Project{<: BuildPhase}
    capacity_bid = 0.
    for product in get_products(project)
        capacity_bid_temp = get_capacity_bid(product)
        if !isnothing(capacity_bid_temp)
            capacity_bid = capacity_bid_temp
        end
    end
    return capacity_bid
end


"""
This function returns the project's REC energy output to be passed the REC market market clearing module.
Returns 0 if there is no REC market participation.
"""
function get_project_expected_rec_output(project::P) where P <: Project{<: BuildPhase}
    rec_output = 0.
    for product in get_products(project)
        rec_output_temp = get_expected_rec_certificates(product)
        if !isnothing(rec_output_temp)
            rec_output = rec_output_temp
        end
    end
    return rec_output
end

"""
This function returns the project's REC market bid to be passed the REC market market clearing module.
Returns 0 if there is no REC market participation.
"""
function get_project_rec_market_bid(project::P) where P <: Project{<: BuildPhase}
    rec_bid = 0.
    for product in get_products(project)
        rec_bid_temp = get_rec_bid(product)
        if !isnothing(rec_bid_temp)
            rec_bid = rec_bid_temp
        end
    end
    return rec_bid
end

function calculate_carbon_emissions(project_carbon_emissions::Float64, product::Product)
    return project_carbon_emissions
end

function calculate_carbon_emissions(project_carbon_emissions::Float64, product::CarbonTax)
    project_carbon_emissions = get_emission_intensity(product) * get_avg_heat_rate(product)
    return project_carbon_emissions
end

function get_emission_intensity(project::Project)
    intensity = 0.0
    for product in get_products(project)
        temp = get_emission_intensity(product)
        if !(isnothing(temp))
            intensity += temp
        end
    end
    return intensity
end

function get_inertia_constant(project::Project)
    constant = 0.0
    for product in get_products(project)
        temp = get_h_constant(product)
        if !(isnothing(temp))
            constant += temp
        end
    end
    return constant
end

function get_synchronous_inertia(project::Project)
    synchronous_inertia = true
    for product in get_products(project)
        temp = get_synchronous(product)
        if !(isnothing(temp))
            synchronous_inertia = temp
        end
    end
    return synchronous_inertia
end

function get_rec_correction_factor(project::Project, iteration_year::Int64)
    value = 1.0
    for product in get_products(project)
        temp = get_rec_correction_factor(product, iteration_year)
        if !(isnothing(temp))
            value = temp
        end
    end
    return value
end

"""
This function creates the MarketProject struct to be passed to CEM price projection and endogeneous Economic Dispatch models.
"""
function populate_market_project(project::P,
                                project_type::String,
                                min_input::Float64,
                                max_input::Float64,
                                efficiency_in::Float64,
                                efficiency_out::Float64,
                                min_storage::Float64,
                                max_storage::Float64,
                                init_storage::Float64,
                                availability_input::Array{Float64, 2},
                                existing_units::Int64,
                                units_inqueue::Vector{Float64},
                                remaining_lag_time::Int64,
                                remaining_life_time::Int64,
                                iteration_year::Int64,
                                num_invperiods::Int64) where P <: Project{<: BuildPhase}

    finance_data = get_finance_data(project)

    products = get_products(project)

    project_reserve_cost = Dict{String, Float64}()
    project_reserve_limit = Dict{String, Float64}()
    project_carbon_emissions = 0.0

    max_cap = get_maxcap(project)

    for product in products
        project_reserve_cost = get_project_reserve_cost(product, project_reserve_cost)
        project_reserve_limit = get_project_reserve_limit(product, max_cap, project_reserve_limit)
        project_carbon_emissions = calculate_carbon_emissions(project_carbon_emissions, product)
    end

    inertia_constant = get_inertia_constant(project)
    synchronous_inertia = get_synchronous_inertia(project)

    market_project = MarketProject(
        get_name(project),                                                                       # name
        project_type,                                                                            # is project of storage type
        get_type(get_tech(project)),                                                             # technology type
        get_fixed_OM_cost(finance_data),                                                           # annualfixed O&M costs
        get_queue_cost(finance_data),                                                             # Queue cost
        get_project_energy_cost(project),                                                        # marginal cost of energy
        project_reserve_cost,                                                                    # marginal cost of reserves
        project_carbon_emissions,                                                                # project emission intensity
        get_investment_cost(finance_data)[iteration_year:iteration_year + num_invperiods - 1],    # yearly investment cost
        get_discount_rate(finance_data),                                                          # discount rate
        get_mincap(project),                                                                     # minimum capacity
        max_cap,                                                                                 # maximum capacity
        min_input,                                                                               # minimum input power
        max_input,                                                                               # maximum input power
        efficiency_in,                                                                           # input efficiency
        efficiency_out,                                                                          # output efficiency
        min_storage,                                                                             # minimum storage capacity
        max_storage,                                                                             # maximum storage capacity
        init_storage,                                                                            # initial storage level
        availability_input,                                                                      # Hourly availability
        get_project_derating(project),                                                           # de-rating factor
        get_ramp_limits(get_tech(project)),                                                       # ramp limits
        project_reserve_limit,                                                                    # reserve participation limits
        existing_units,                                                                          # existing units
        units_inqueue,                                                                           # units in queue
        get_lag_time(finance_data),                                                               # construction lead time
        remaining_lag_time,                                                                      # remaining construction time
        1,                                                                                       # maximum units
        1,                                                                                       # base cost units
        get_capex_years(finance_data),                                                            # capital cost recovery years
        get_life_time(finance_data),                                                              # total life_time
        remaining_life_time,                                                                      # remaining life_time
        in(:Capacity, get_name.(get_products(project))),                                         # eligible for capacity markets
        in(:REC, get_name.(get_products(project))),                                              # eligible for rps compliance
        get_rec_correction_factor(project, iteration_year),                                      # rec correction factor
        inertia_constant,                                                                        # inertia H-constant,
        synchronous_inertia,                                                                     # whether inertia is synchronous
        get_zone(get_tech(project)),                                                             # project zone
        [get_ownedby(finance_data)])                                                             # owned by

    return market_project
end

"""
This function processes data of Existing projects to create the MarketProject struct
passed to CEM price projection and endogeneous Economic Dispatch models.
"""
function create_market_project(project::P,
                              pricecap_energy::Float64,
                              max_peak_loads::AxisArrays.AxisArray{Float64, 1},
                              iteration_year::Int64,
                              num_hours::Int64,
                              num_invperiods::Int64,
                              availability_df::DataFrames.DataFrame) where P <: Project{Existing}

    finance_data = get_finance_data(project)

    queue_time = length(get_queue_cost(finance_data))

    # Get technical characteristics based on project type
    project_type,
    min_input,
    max_input,
    efficiency_in,
    efficiency_out,
    min_storage,
    max_storage,
    init_storage = get_technical_details(project)

    availability_input = find_project_availability_data(project, availability_df, num_invperiods, num_hours)

    # Calculate remaining life_time
    remaining_life_time = min(get_life_time(finance_data), get_end_life_year(project) - iteration_year + 1)

    units_inqueue = zeros(queue_time)

    # Existing project have completed their queue time
    units_inqueue[queue_time] = 1.0

    existing_units = 1
    remaining_lag_time = 0


    market_project = populate_market_project(project,
                                             project_type,
                                             min_input,
                                             max_input,
                                             efficiency_in,
                                             efficiency_out,
                                             min_storage,
                                             max_storage,
                                             init_storage,
                                             availability_input,
                                             existing_units,
                                             units_inqueue,
                                             remaining_lag_time,
                                             remaining_life_time,
                                             iteration_year,
                                             num_invperiods)
    return market_project

end

"""
This function processes data of Option projects to create the MarketProject struct
passed to CEM price projection and endogeneous Economic Dispatch models.
"""
function create_market_project(project::P,
                              pricecap_energy::Float64,
                              max_peak_loads::AxisArrays.AxisArray{Float64, 1},
                              iteration_year::Int64,
                              num_hours::Int64,
                              num_invperiods::Int64,
                              availability_df::DataFrames.DataFrame) where P <: Project{Option}

    # Get technical characteristics based on project type
    project_type,
    min_input,
    max_input,
    efficiency_in,
    efficiency_out,
    min_storage,
    max_storage,
    init_storage = get_technical_details(project)

    availability_input = find_project_availability_data(project, availability_df, num_invperiods, num_hours)

    finance_data = get_finance_data(project)

    queue_time = length(get_queue_cost(finance_data))

    # End of life_time for option projects can be up to the end of horizon.
    remaining_life_time = num_invperiods
    units_inqueue = zeros(queue_time)
    existing_units = 0
    remaining_lag_time = get_lag_time(finance_data) + queue_time

    market_project = populate_market_project(project,
                                             project_type,
                                             min_input,
                                             max_input,
                                             efficiency_in,
                                             efficiency_out,
                                             min_storage,
                                             max_storage,
                                             init_storage,
                                             availability_input,
                                             existing_units,
                                             units_inqueue,
                                             remaining_lag_time,
                                             remaining_life_time,
                                             iteration_year,
                                             num_invperiods)

    market_project.name = "option_$(market_project.tech_type)_$(market_project.zone)"

    if market_project.derating_factor == 0.0
        derating = 0.5                              # Dummy derating for maximum capacity if capacity product doesn't exist
    else
        derating = market_project.derating_factor
    end

    modified_max_gen = max_peak_loads[market_project.zone] / derating
    market_project.max_new_options = round(modified_max_gen / market_project.max_gen)

    return market_project

end

"""
This function processes data of Planned projects to create the MarketProject struct
passed to CEM price projection and endogeneous Economic Dispatch models.
"""
function create_market_project(project::P,
                              pricecap_energy::Float64,
                              max_peak_loads::AxisArrays.AxisArray{Float64, 1},
                              iteration_year::Int64,
                              num_hours::Int64,
                              num_invperiods::Int64,
                              availability_df::DataFrames.DataFrame) where P <: Project{Planned}

    # Get technical characteristics based on project type
    project_type,
    min_input,
    max_input,
    efficiency_in,
    efficiency_out,
    min_storage,
    max_storage,
    init_storage = get_technical_details(project)

    availability_input = find_project_availability_data(project, availability_df, num_invperiods, num_hours)

    finance_data = get_finance_data(project)

    queue_time = length(get_queue_cost(finance_data))

    # Number of construction years left
    remaining_lag_time = get_construction_year(project) - iteration_year

    units_inqueue = zeros(queue_time)

    # Planned units have completed their queue time
    units_inqueue[queue_time] = 1.0

    # Calculate remaining life_time
    remaining_life_time = get_life_time(finance_data) + remaining_lag_time

    existing_units = 0

    market_project = populate_market_project(project,
                                             project_type,
                                             min_input,
                                             max_input,
                                             efficiency_in,
                                             efficiency_out,
                                             min_storage,
                                             max_storage,
                                             init_storage,
                                             availability_input,
                                             existing_units,
                                             units_inqueue,
                                             remaining_lag_time,
                                             remaining_life_time,
                                             iteration_year,
                                             num_invperiods)
    return market_project

end

"""
This function processes data of Queue projects to create the MarketProject struct
passed to CEM price projection and endogeneous Economic Dispatch models.
"""
function create_market_project(project::P,
                              pricecap_energy::Float64,
                              max_peak_loads::AxisArrays.AxisArray{Float64, 1},
                              iteration_year::Int64,
                              num_hours::Int64,
                              num_invperiods::Int64,
                              availability_df::DataFrames.DataFrame) where P <: Project{Queue}

    # Get technical characteristics based on project type
    project_type,
    min_input,
    max_input,
    efficiency_in,
    efficiency_out,
    min_storage,
    max_storage,
    init_storage = get_technical_details(project)

    availability_input = find_project_availability_data(project, availability_df, num_invperiods, num_hours)

    finance_data = get_finance_data(project)

    queue_time = length(get_queue_cost(finance_data))

    units_inqueue = zeros(queue_time)

    queue_year = iteration_year - get_decision_year(project)

    # Calculate number of years left in queue
    if queue_year <= queue_time
        units_inqueue[queue_year] = 1.0
    end

    remaining_lag_time = get_lag_time(finance_data) + queue_time - queue_year

    # Calculate remaining life_time
    remaining_life_time = get_life_time(finance_data) + get_lag_time(finance_data) + queue_time - queue_year

    existing_units = 0

     market_project = populate_market_project(project,
                                             project_type,
                                             min_input,
                                             max_input,
                                             efficiency_in,
                                             efficiency_out,
                                             min_storage,
                                             max_storage,
                                             init_storage,
                                             availability_input,
                                             existing_units,
                                             units_inqueue,
                                             remaining_lag_time,
                                             remaining_life_time,
                                             iteration_year,
                                             num_invperiods)
    return market_project

end
