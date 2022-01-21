"""
This function creates the MarketClearingProblem struct for CEM.
"""
function create_cem_mkt_clr_problem(investor_dir::String,
                                    market_names::Vector{Symbol},
                                    carbon_tax::Vector{Float64},
                                    reserve_products::Vector{String},
                                    ordc_products::Vector{String},
                                    rps_target::String,
                                    reserve_penalty::String,
                                    resource_adequacy::ResourceAdequacy,
                                    irm_scalar::Float64,
                                    expected_portfolio::Vector{<: Project{<: BuildPhase}},
                                    zones::Vector{String},
                                    lines::Vector{ZonalLine},
                                    peak_load::Float64,
                                    rep_hour_weight::Vector{Float64},
                                    average_capital_cost_multiplier::Float64,
                                    scenario::Scenario,
                                    iteration_year::Int64,
                                    yearly_horizon::Int64)

    num_invperiods = yearly_horizon

    #Get updated load growth belief
    load_growth = AxisArrays.AxisArray(zeros(length(zones), num_invperiods), zones, collect(1:num_invperiods))
    for (idx, zone) in enumerate(zones)
        load_growth[idx, :] = get_parameter_values(scenario)[iteration_year]["load_$(zone)", iteration_year:iteration_year + num_invperiods - 1]
    end

    average_load_growth = [Statistics.mean(load_growth[:, p]) for p in 1:num_invperiods]

    carbon_tax_vector = carbon_tax[iteration_year:iteration_year + num_invperiods - 1]

    # Gather markets data-------------------------------------------------------------------------------

    ######################################### Create Energy Markets ####################################################
    load_data = read_data(joinpath(investor_dir, "timeseries_data_files", "Load", "load_$(iteration_year - 1).csv"))

    num_hours = DataFrames.nrow(load_data)

    energy_mkt_params = read_data(joinpath(investor_dir, "markets_data", "Energy.csv"))
    price_cap_energy = AxisArrays.AxisArray(energy_mkt_params.price_cap * 1.0, zones)
    zonal_load = AxisArrays.AxisArray(zeros(length(zones), num_hours), zones, (1:num_hours))

    for (zone_num, zone) in enumerate(zones)
        for h in 1:num_hours
            zonal_load[zone, h] = load_data[:, Symbol(zone_num)][h]
        end
    end

    energy_annual_increment = AxisArrays.AxisArray(ones(length(zones), num_invperiods), zones, collect(1:num_invperiods))

    energy_markets = Vector{EnergyMarket}(undef, num_invperiods)

    for p in 1:num_invperiods
        for z in zones
            for i in 1:p
                energy_annual_increment[z, p] = energy_annual_increment[z, p] * (1 + load_growth[z, i])
            end
        end

        energy_markets[p] = EnergyMarket(AxisArrays.AxisArray(zonal_load .* energy_annual_increment[:, p], zones, (1:num_hours)),
                                     price_cap_energy)
    end

    ######################################### Create Reserve Markets ####################################################

    average_annual_increment = ones(num_invperiods)

    reserve_up_markets = Vector{Dict{String, ReserveUpMarket{num_hours}}}(undef, num_invperiods)
    reserve_down_markets = Vector{Dict{String, ReserveDownMarket{num_hours}}}(undef, num_invperiods)
    reserve_ordc_markets = Vector{Dict{String, ReserveORDCMarket{num_hours}}}(undef, num_invperiods)

    reserve_timeseries_data = Dict(r => read_data(joinpath(investor_dir, "timeseries_data_files", "Reserves", "$(r)_$(iteration_year - 1).csv"))[:, r] for r in reserve_products)
    reserve_parameter_data = Dict(r => read_data(joinpath(investor_dir, "markets_data", "$(reserve_penalty)_reserve_penalty", "$(r).csv")) for r in reserve_products)

    reserve_eligible_projects = Dict(product => String[] for product in reserve_products)

    for p in 1:num_invperiods
        reserve_up_market = Dict{String, ReserveUpMarket{num_hours}}()
        reserve_down_market = Dict{String, ReserveDownMarket{num_hours}}()
        reserve_ordc_market = Dict{String, ReserveORDCMarket{num_hours}}()

        for i in 1:p
            average_annual_increment[p] = average_annual_increment[p] * (1 + average_load_growth[i])
        end

        for product in reserve_products

            if Symbol(product) in market_names
                timeseries_data = reserve_timeseries_data[product]
                parameter_data = reserve_parameter_data[product]

                if product in ordc_products
                    market = create_ordc_market(timeseries_data, parameter_data, reserve_eligible_projects[product])
                    reserve_ordc_market[product] = market
                else
                    direction = lowercase(parameter_data[1, "direction"])
                    price_cap = Float64(parameter_data[1, "price_cap"])
                    zones = ["zone_$(n)" for n in split(parameter_data[1, "eligible_zones"], ";")]
                    if direction == "up"
                        market = ReserveUpMarket(timeseries_data * average_annual_increment[p], price_cap, zones, reserve_eligible_projects[product])
                        reserve_up_market[product] = market
                    elseif direction == "down"
                        market = ReserveDownMarket(timeseries_data * average_annual_increment[p], price_cap, zones, reserve_eligible_projects[product])
                        reserve_down_market[product] = market
                    end
                end
            end
        end
        reserve_up_markets[p] = reserve_up_market
        reserve_down_markets[p] = reserve_down_market
        reserve_ordc_markets[p] = reserve_ordc_market
    end

    ######################################### Create Capacity and REC Markets ####################################################

    capacity_market_bool = false
    rec_market_bool = false
    inertia_market_bool = false

   if in(:Capacity, market_names)
        capacity_market_bool = true
    end

    if in(:REC, market_names)
        rec_market_bool = true
    end

    delta_irm = get_delta_irm(resource_adequacy, iteration_year)

    capacity_mkt_param_file = joinpath(investor_dir, "markets_data", "Capacity.csv")

    REC_mkt_params = read_data(joinpath(investor_dir, "markets_data", "REC_$(rps_target)_RPS.csv"))
    price_cap_rec = REC_mkt_params[1, "price_cap"]
    rec_req = REC_mkt_params[1, "rec_req"] * rec_market_bool
    rec_annual_increment = REC_mkt_params[1, "annual_increment"] * rec_market_bool
    rec_non_binding_years = REC_mkt_params[1, "non_binding_years"] * rec_market_bool

    rec_binding_array = Int.(zeros(num_invperiods))
    binding_years_from_now = max((rec_non_binding_years - iteration_year + 1), 0) + 1

    for j in max(4, binding_years_from_now):num_invperiods
        rec_binding_array[j] = 1
    end

    if in(:Inertia, market_names)
        inertia_market_bool = true
    end

    inertia_mkt_params = read_data(joinpath(investor_dir, "markets_data", "$(reserve_penalty)_reserve_penalty", "Inertia.csv"))
    price_cap_inertia = inertia_mkt_params[1, "price_cap"]
    inertia_req_multiplier = inertia_mkt_params[1, "requirement_multiplier"] * inertia_market_bool

    capacity_markets = Vector{CapacityMarket}(undef, num_invperiods)
    rec_markets = Vector{RECMarket}(undef, num_invperiods)
    inertia_markets = Vector{InertiaMarket}(undef, num_invperiods)

    for p in 1:num_invperiods
        system_peak_load = average_annual_increment[p] * peak_load

        capacity_markets[p] = create_capacity_demand_curve(capacity_mkt_param_file, system_peak_load, irm_scalar, delta_irm, capacity_market_bool)
        rec_markets[p] = RECMarket(min(rec_req + rec_annual_increment * (p + iteration_year - 1), 1), price_cap_rec, !(iszero(rec_binding_array[p])))
        inertia_markets[p] = InertiaMarket(system_peak_load * inertia_req_multiplier, price_cap_inertia)
    end

    max_peak_loads = AxisArrays.AxisArray([maximum([maximum(market.demand[z, :]) for market in energy_markets]) for z in zones], zones)

    markets = MarketCollection.(capacity_markets,
                                energy_markets,
                                reserve_up_markets,
                                reserve_down_markets,
                                reserve_ordc_markets,
                                rec_markets,
                                inertia_markets)
    #-----------------------------------------------------------------------------------------------------------------------------
    availability_df = read_data(joinpath(investor_dir, "timeseries_data_files", "Availability", "DAY_AHEAD_availability.csv"))

    invested_portfolio = find_active_invested_projects(expected_portfolio)
    option_portfolio = find_option_projects(expected_portfolio)

    cem_projects = MarketProject[]

    for project in invested_portfolio
        tech = get_tech(project)
        zone = get_zone(tech)
        cem_project = create_market_project(project,
                                          price_cap_energy[zone],
                                          max_peak_loads,
                                          iteration_year,
                                          num_hours,
                                          num_invperiods,
                                          availability_df)

        if !isnothing(cem_project)
            push!(cem_projects, cem_project)
        end

        products =  get_products(project)
        for product in products
            product_name = String(get_name(product))
            if product_name in reserve_products
                push!(reserve_eligible_projects[product_name], cem_project.name)
            end
        end
    end

    aggregated_options = MarketProject[]
    for option in option_portfolio
        tech = get_tech(option)
        zone = get_zone(tech)
        similar_option = filter(p -> (get_type(tech) == p.tech_type && zone == p.zone), aggregated_options)

        if length(similar_option) < 1
            aggregated_option = create_market_project(option,
                                          price_cap_energy[zone],
                                          max_peak_loads,
                                          iteration_year,
                                          num_hours,
                                          num_invperiods,
                                          availability_df)
            push!(aggregated_options, aggregated_option)

            products =  get_products(option)
            for product in products
                product_name = String(get_name(product))
                if product_name in reserve_products
                    push!(reserve_eligible_projects[product_name], aggregated_option.name)
                end
            end
        else
            investment_cost = get_investment_cost(get_finance_data(option))[iteration_year:iteration_year + num_invperiods - 1]
            for (i, cost) in enumerate(similar_option[1].expansion_cost)
                if investment_cost[i] > cost
                    similar_option[1].expansion_cost[i] = investment_cost[i]
                end
            end

            discount_rate = get_discount_rate(get_finance_data(option))
            if discount_rate > similar_option[1].discount_rate
                similar_option[1].discount_rate = discount_rate
            end

            owner = get_ownedby(get_finance_data(option))
            if !in(owner, similar_option[1].ownedby)
               similar_option[1].base_cost_units = similar_option[1].base_cost_units + 1
               push!(similar_option[1].ownedby, owner)
            end
        end

    end

    append!(cem_projects, aggregated_options)

    system = MarketClearingProblem(zones, lines, average_capital_cost_multiplier, markets, carbon_tax_vector, cem_projects, rep_hour_weight)

    return system
end
