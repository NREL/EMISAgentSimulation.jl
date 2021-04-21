"""
This function creates the MarketClearingProblem struct for endogenous Economic Dispatch.
"""
function create_economic_dispatch_problem(simulation::AgentSimulation,
                                          sys_UC::PSY.System,
                                          market_names::Vector{Symbol},
                                          num_invperiods::Int64,
                                          load_growth::AxisArrays.AxisArray{Float64, 1},
                                          existing_projects::Vector{<: Project{<: BuildPhase}},
                                          iteration_year::Int64)

    simulation_dir = get_data_dir(get_case(simulation))
    load_data = read_data(joinpath(simulation_dir, "timeseries_data_files", "Load", "load_$(iteration_year - 1).csv"))
    num_hours = DataFrames.nrow(load_data)
    zones = get_zones(simulation)
    zonal_load = AxisArrays.AxisArray(zeros(length(zones), num_hours), zones, (1:num_hours))

    for (zone_num, zone) in enumerate(zones)
        for h in 1:num_hours
            zonal_load[zone, h] = load_data[:, Symbol(zone_num)][h]
        end
    end

    if isnothing(sys_UC)
        lines = get_lines(simulation)

        reserve_up_market_bool = false
        reserve_down_market_bool = false
        capacity_market_bool = false
        rec_market_bool = false

        #Find if the following markets are being simulated:
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
        reserve_up_demand_data = read_data(joinpath(simulation_dir,  "timeseries_data_files", "Reserves", "reserve_up_$(iteration_year - 1).csv"))
        reserve_down_demand_data = read_data(joinpath(simulation_dir,  "timeseries_data_files", "Reserves", "reserve_down_$(iteration_year - 1).csv"))

        energy_mkt_params = read_data(joinpath(simulation_dir, "markets_data", "energy_mkt_param.csv"))
        price_cap_energy = AxisArrays.AxisArray(energy_mkt_params.price_cap * 1.0, zones)

        reserve_up_mkt_params = read_data(joinpath(simulation_dir, "markets_data", "reserve_up_mkt_param.csv"))
        reserve_up_num_points = AxisArrays.AxisArray(zeros(length(zones)), zones)
        reserve_up_break_points = AxisArrays.AxisArray([Float64[] for z in zones], zones)
        reserve_up_price_points = AxisArrays.AxisArray([Float64[] for z in zones], zones)

        reserve_down_mkt_params = read_data(joinpath(simulation_dir, "markets_data", "reserve_down_mkt_param.csv"))
        price_cap_reservedown = AxisArrays.AxisArray(reserve_down_mkt_params.price_cap * 1.0, zones)

        zonal_reserve_up = AxisArrays.AxisArray(zeros(length(zones), num_hours), zones, (1:num_hours))
        zonal_reserve_down = AxisArrays.AxisArray(zeros(length(zones), num_hours), zones, (1:num_hours))

        for (zone_num, zone) in enumerate(zones)
            idx = findfirst(x -> x == "$(zone)", reserve_up_mkt_params.zones)
            if reserve_up_mkt_params.ORDC[idx]
                reserve_up_num_points[zone] = reserve_up_mkt_params.ORDC_points[idx]
                for point in 1:reserve_up_num_points[zone]
                    push!(reserve_up_break_points[zone], reserve_up_mkt_params[idx, Symbol("ORDC_x$(1)")]) * reserve_up_market_bool
                    push!(reserve_up_price_points[zone], reserve_up_mkt_params[idx, Symbol("ORDC_y$(1)")]) * reserve_up_market_bool
                end
            else
                reserve_up_break_points[zone] = [0.0, 1.0, 1.0, 10.0] * reserve_up_market_bool
                reserve_up_price_points[zone] = [reserve_up_mkt_params.price_cap[idx], reserve_up_mkt_params.price_cap[idx], 0.0, 0.0] * reserve_up_market_bool
            end

            price_cap_reservedown[zone] = reserve_down_mkt_params.price_cap[idx]

            for h in 1:num_hours
                zonal_reserve_up[zone, h] = reserve_up_demand_data[:, Symbol(zone_num)][h] * reserve_up_market_bool
                zonal_reserve_down[zone, h] = reserve_down_demand_data[:, Symbol(zone_num)][h] * reserve_down_market_bool
            end
        end

        energy_annual_increment = load_growth
        reserveup_annual_increment = load_growth
        reservedown_annual_increment = load_growth

        capacity_mkt_param_file = joinpath(simulation_dir, "markets_data", "capacity_mkt_param.csv")
        peak_load = get_peak_load(simulation)
        capacity_annual_increment = load_growth

        REC_mkt_params = read_data(joinpath(simulation_dir, "markets_data", "REC_mkt_param.csv"))
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

            energy_markets[p] = EnergyMarket(AxisArrays.AxisArray(zonal_load .* ((1 .+ energy_annual_increment) .^ p), zones, (1:num_hours)),
                                        price_cap_energy)
            reserve_up_markets[p] = ReserveUpMarket(reserve_up_break_points,
                                                    reserve_up_price_points,
                                                    AxisArrays.AxisArray(zonal_reserve_up .* ((1 .+ reserveup_annual_increment) .^ p), zones, (1:num_hours)))
            reserve_down_markets[p] = ReserveDownMarket(AxisArrays.AxisArray(zonal_reserve_down .* ((1 .+ reservedown_annual_increment) .^ p), zones, (1:num_hours)),
                                                    price_cap_reservedown)
            rec_markets[p] = RECMarket(min(rec_req + rec_annual_increment * (p + iteration_year - 1), 1),
                                price_cap_rec)
        end

        max_peak_loads = AxisArrays.AxisArray([maximum([maximum(market.demand[z, :]) for market in energy_markets]) for z in zones], zones)

        hour_weight = get_hour_weight(simulation)

        markets = MarketCollection.(capacity_markets,
                                    energy_markets,
                                    reserve_up_markets,
                                    reserve_down_markets,
                                    rec_markets)
        #--------------------------------------------------------------------------------------------------------------

        availability_df = read_data(joinpath(simulation_dir, "timeseries_data_files", "Availability", "DAY_AHEAD_availability.csv"))

        ed_projects = MarketProject[]

        for project in existing_projects
                tech = get_tech(project)
                ed_project = create_market_project(project,
                                                price_cap_energy[get_zone(tech)],
                                                maximum(reserve_up_break_points[get_zone(tech)]),
                                                price_cap_reservedown[get_zone(tech)],
                                                max_peak_loads,
                                                iteration_year,
                                                num_hours,
                                                num_invperiods,
                                                availability_df)

                push!(ed_projects, ed_project)
        end

        system = MarketClearingProblem(zones, lines, 0.0, markets, ed_projects, hour_weight)

        return zonal_load, system
    else
        return zonal_load, nothing
    end
end
