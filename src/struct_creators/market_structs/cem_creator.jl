"""
This function creates the MarketClearingProblem struct for CEM.
"""
function create_cem_mkt_clr_problem(investor_dir::String,
                                    market_names::Vector{Symbol},
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

    reserve_up_market_bool = false
    reserve_down_market_bool = false
    capacity_market_bool = false
    rec_market_bool = false

    # Find if the investor participates in the following markets:
    if in(:ReserveUp, market_names)
        reserve_up_market_bool = true
    end

    if in(:ReserveDown, market_names)
        reserve_down_market_bool = true
    end

    if in(:Capacity, market_names)
        capacity_market_bool = true
    end

    if in(:REC, market_names)
        rec_market_bool = true
    end

    # Gather markets data-------------------------------------------------------------------------------
    load_data = read_data(joinpath(investor_dir, "timeseries_data_files", "Load", "load_$(iteration_year - 1).csv"))
    reserve_up_demand_data = read_data(joinpath(investor_dir, "timeseries_data_files", "Reserves", "reserve_up_$(iteration_year - 1).csv"))
    reserve_down_demand_data = read_data(joinpath(investor_dir, "timeseries_data_files", "Reserves", "reserve_down_$(iteration_year - 1).csv"))

    num_hours = DataFrames.nrow(load_data)

    energy_mkt_params = read_data(joinpath(investor_dir, "markets_data", "energy_mkt_param.csv"))
    price_cap_energy = AxisArrays.AxisArray(energy_mkt_params.price_cap * 1.0, zones)

    reserve_up_mkt_params = read_data(joinpath(investor_dir, "markets_data", "reserve_up_mkt_param.csv"))
    reserve_up_num_points = AxisArrays.AxisArray(zeros(Int64,length(zones)), zones)
    reserve_up_break_points = AxisArrays.AxisArray([Float64[] for z in zones], zones)
    reserve_up_price_points = AxisArrays.AxisArray([Float64[] for z in zones], zones)

    reserve_down_mkt_params = read_data(joinpath(investor_dir, "markets_data", "reserve_down_mkt_param.csv"))
    price_cap_reservedown = AxisArrays.AxisArray(reserve_down_mkt_params.price_cap * 1.0, zones)
    perc_req_reserve_down = AxisArrays.AxisArray(reserve_down_mkt_params.perc_req, zones)

    zonal_load = AxisArrays.AxisArray(zeros(length(zones), num_hours), zones, (1:num_hours))
    zonal_reserve_up = AxisArrays.AxisArray(zeros(length(zones), num_hours), zones, (1:num_hours))
    zonal_reserve_down = AxisArrays.AxisArray(zeros(length(zones), num_hours), zones, (1:num_hours))

    for (zone_num, zone) in enumerate(zones)
        idx = findfirst(x -> x == "$(zone)", reserve_up_mkt_params.zones)
        if reserve_up_mkt_params.ORDC[idx]
            reserve_up_num_points[zone] = reserve_up_mkt_params.ORDC_points[idx]
            for point in 1:reserve_up_num_points[zone]
                push!(reserve_up_break_points[zone], reserve_up_mkt_params[idx, Symbol("ORDC_x$(point)")] * reserve_up_market_bool)
                push!(reserve_up_price_points[zone], reserve_up_mkt_params[idx, Symbol("ORDC_y$(point)")] * reserve_up_market_bool)
            end
        else
            reserve_up_break_points[zone] = [0.0, 1.0, 1.0, 10.0] * reserve_up_market_bool
            reserve_up_price_points[zone] = [reserve_up_mkt_params.price_cap[idx], reserve_up_mkt_params.price_cap[idx], 0.0, 0.0] * reserve_up_market_bool
        end

        price_cap_reservedown[zone] = reserve_down_mkt_params.price_cap[idx]

        for h in 1:num_hours
            zonal_load[zone, h] = load_data[:, Symbol(zone_num)][h]
            zonal_reserve_up[zone, h] = reserve_up_demand_data[:, Symbol(zone_num)][h] * reserve_up_market_bool
            zonal_reserve_down[zone, h] = reserve_down_demand_data[:, Symbol(zone_num)][h] * reserve_down_market_bool
        end
    end

    energy_annual_increment = AxisArrays.AxisArray(ones(length(zones), num_invperiods), zones, collect(1:num_invperiods))
    reserveup_annual_increment = AxisArrays.AxisArray(ones(length(zones), num_invperiods), zones, collect(1:num_invperiods))
    reservedown_annual_increment = AxisArrays.AxisArray(ones(length(zones), num_invperiods), zones, collect(1:num_invperiods))

    for z in zones
        for p in 1:num_invperiods
            for i in 1:p
                energy_annual_increment[z, p] = energy_annual_increment[z, p] * (1 + load_growth[z, i])
                reserveup_annual_increment[z, p] = energy_annual_increment[z, p] * (1 + load_growth[z, i])
                reservedown_annual_increment[z, p] = energy_annual_increment[z, p] * (1 + load_growth[z, i])
            end
        end
    end

    capacity_mkt_param_file = joinpath(investor_dir, "markets_data", "capacity_mkt_param.csv")
    capacity_annual_increment = load_growth

    REC_mkt_params = read_data(joinpath(investor_dir, "markets_data", "REC_mkt_param.csv"))
    price_cap_rec = REC_mkt_params.price_cap[1]
    rec_req = REC_mkt_params.rec_req[1] * rec_market_bool
    rec_annual_increment = REC_mkt_params.annual_increment[1] * rec_market_bool

    capacity_markets = Vector{CapacityMarket}(undef, num_invperiods)
    energy_markets = Vector{EnergyMarket}(undef, num_invperiods)
    reserve_up_markets = Vector{ReserveUpMarket}(undef, num_invperiods)
    reserve_down_markets = Vector{ReserveDownMarket}(undef, num_invperiods)

    rec_markets = Vector{RECMarket}(undef, num_invperiods)


    average_capacity_growth = Statistics.mean(capacity_annual_increment)

    # Create market products data for the horizon
    for p in 1:num_invperiods
        system_peak_load = (1 + average_capacity_growth) ^ (p) * peak_load
        capacity_markets[p] = create_capacity_demand_curve(capacity_mkt_param_file, system_peak_load, capacity_market_bool)

        energy_markets[p] = EnergyMarket(AxisArrays.AxisArray(zonal_load .* energy_annual_increment[:, p], zones, (1:num_hours)),
                                     price_cap_energy)
        reserve_up_markets[p] = ReserveUpMarket(reserve_up_break_points,
                                                reserve_up_price_points,
                                                AxisArrays.AxisArray(zonal_reserve_up .* reserveup_annual_increment[:, p], zones, (1:num_hours)))
        reserve_down_markets[p] = ReserveDownMarket(AxisArrays.AxisArray(zonal_reserve_down .* reservedown_annual_increment[:, p], zones, (1:num_hours)),
                                                  price_cap_reservedown)
        rec_markets[p] = RECMarket(min(rec_req + rec_annual_increment * (p + iteration_year - 1), 1),
                             price_cap_rec)
    end

    max_peak_loads = AxisArrays.AxisArray([maximum([maximum(market.demand[z, :]) for market in energy_markets]) for z in zones], zones)

    markets = MarketCollection.(capacity_markets,
                                energy_markets,
                                reserve_up_markets,
                                reserve_down_markets,
                                rec_markets)
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
                                          maximum(reserve_up_break_points[zone]),
                                          price_cap_reservedown[zone],
                                          max_peak_loads,
                                          iteration_year,
                                          num_hours,
                                          num_invperiods,
                                          availability_df)

        if !isnothing(cem_project)
            push!(cem_projects, cem_project)
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
                                          maximum(reserve_up_break_points[zone]),
                                          price_cap_reservedown[zone],
                                          max_peak_loads,
                                          iteration_year,
                                          num_hours,
                                          num_invperiods,
                                          availability_df)
            push!(aggregated_options, aggregated_option)
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

    system = MarketClearingProblem(zones, lines, average_capital_cost_multiplier, markets, cem_projects, rep_hour_weight)

    return system
end
