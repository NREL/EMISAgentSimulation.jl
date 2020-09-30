
"""
This function creates the realized market data based on actual market clearing.
"""
function create_realized_marketdata(simulation::AgentSimulation,
                                    sys_UC::Union{Nothing, PSY.System},
                                    sys_ED::Union{Nothing, PSY.System},
                                    market_names::Vector{Symbol},
                                    existing_projects::Vector{<: Project{<: BuildPhase}},
                                    capacity_market_projects::Vector{<: Project{<: BuildPhase}},
                                    capacity_forward_years::Int64,
                                    iteration_year::Int64,
                                    simulation_years::Int64,
                                    solver::JuMP.MOI.OptimizerWithAttributes,
                                    results_dir::String)

    num_invperiods = 1

    simulation_dir = get_data_dir(get_case(simulation))
    zones = get_zones(simulation)
    lines = get_lines(simulation)
    annual_growth = get_annual_growth(simulation)[:, iteration_year]

    hour_weight = get_hour_weight(simulation)

    # Get actual value of load growth
    load_growth = AxisArrays.AxisArray(zeros(length(zones)), zones)
    for (idx, zone) in enumerate(zones)
        load_growth[idx] = annual_growth["load_$(zone)"]
    end

    energy_annual_increment = load_growth
    reserveup_annual_increment = load_growth
    reservedown_annual_increment = load_growth

    ######## Energy and Ancillary Services market clearing ############################################

    zonal_load, system = create_economic_dispatch_problem(simulation,
                                                          market_names,
                                                          num_invperiods,
                                                          load_growth,
                                                          existing_projects,
                                                          iteration_year)

    num_hours = get_num_hours(system)
    num_days = Int(num_hours/24)

    energy_price,
    reserve_up_price,
    reserve_down_price,
    capacity_factors,
    reserve_up_perc,
    reserve_down_perc = energy_mkt_clearing(sys_UC, sys_ED, system, load_growth, zones, num_days, solver)

    # Replace energy_mkt_clearing(nothing, nothing, system, load_growth, zones, num_days, solver) with
    # energy_mkt_clearing(sys_UC, sys_ED, system, load_growth, zones, num_days, solver) to run SIIP production cost model

    market_prices = Dict(name => zeros(num_invperiods, num_hours) for name in market_names)

    market_prices = MarketPrices()

    # Fill energy and ancillary services market prices data
    set_energy_price!(market_prices, "realized", energy_price)

    if in(:ReserveUp, market_names)
        set_reserveup_price!(market_prices, "realized", reserve_up_price)
    end

    if in(:ReserveDown, market_names)
        set_reservedown_price!(market_prices, "realized", reserve_down_price)
    end

    ######## Capacity market clearing #####################################################################

    capacity_market_bool = false
    if in(:Capacity, market_names)
        capacity_market_bool = true
    end

    capacity_mkt_param_file = joinpath(simulation_dir, "markets_data", "capacity_mkt_param.csv")
    peak_load = get_peak_load(simulation)
    capacity_annual_increment = load_growth

    average_capacity_growth = Statistics.mean(capacity_annual_increment)

    capacity_supply_curve = Vector{Union{String, Float64}}[]

    for project in capacity_market_projects
        for product in get_products(project)
            capacity_supply_curve = update_capacity_supply_curve!(capacity_supply_curve, product, project)
        end
    end

    capacity_price =  AxisArrays.AxisArray(reshape([0.0], 1,), [1])
    capacity_accepted_bids = Dict("no_accepted_bids" => 0.0)

    if in(:Capacity, market_names) && iteration_year + capacity_forward_years - 1 <= simulation_years
        system_peak_load = (1 + average_capacity_growth) ^ (capacity_forward_years) * peak_load
        capacity_demand_curve = create_capacity_demand_curve(capacity_mkt_param_file, system_peak_load, capacity_market_bool)

        sort!(capacity_supply_curve, by = x -> x[3])      # Sort capacity supply curve by capacity bid

        capacity_price, capacity_accepted_bids = capacity_market_clearing(capacity_demand_curve, capacity_supply_curve, solver)
    end

    set_capacity_price!(market_prices, "realized", capacity_price)

    ######### REC market clearing ############################################################################

    rec_market_bool = false
    if in(:REC, market_names)
        rec_market_bool = true
    end

    REC_mkt_params = read_data(joinpath(simulation_dir, "markets_data", "REC_mkt_param.csv"))
    pricecap_rec = REC_mkt_params.price_cap[1]
    rec_req = REC_mkt_params.rec_req[1] * rec_market_bool
    rec_annual_increment = REC_mkt_params.annual_increment[1] * rec_market_bool

    rec_supply_curve = Vector{Union{String, Float64}}[]

    for project in existing_projects
         # Populate REC market supply curves
         for product in get_products(project)
            rec_supply_curve = update_rec_supply_curve!(rec_supply_curve, product, project)
        end
    end

    rec_price = AxisArrays.AxisArray(reshape([pricecap_rec], 1,), [1])
    rec_accepted_bids = Dict{String, Float64}()

    if in(:REC, market_names)
        if length(rec_supply_curve) >= 1
            rec_energy_requirment = 0.0

            for z in zones
                for t in 1:num_hours
                    rec_energy_requirment += zonal_load[z, t] * (1 + energy_annual_increment[z]) * hour_weight[t]
                end
            end

            rec_energy_requirment =  rec_energy_requirment * min(rec_req + (rec_annual_increment * iteration_year), 1)


            sort!(rec_supply_curve, by = x -> x[3])      # Sort REC supply curve by REC bid

            rec_price, rec_accepted_bids = rec_market_clearing(rec_energy_requirment, pricecap_rec, rec_supply_curve, solver)

            set_rec_price!(market_prices, "realized", rec_price)
        end
    end

    ################# Write actual market clearing data ################################################

    output_file = joinpath(results_dir, "realized_market_data", "year_$(iteration_year).jld2")

    FileIO.save(output_file,
                     "capacity_price", capacity_price,
                     "energy_price", energy_price,
                     "reserve_up_price", reserve_up_price,
                     "reserve_down_price", reserve_down_price,
                     "rec_price", rec_price,
                     "capacity_factors", capacity_factors,
                     "reserve_up_perc", reserve_up_perc,
                     "reserve_down_perc", reserve_down_perc,
                     "capacity_accepted_bids", capacity_accepted_bids,
                     "rec_accepted_bids", rec_accepted_bids
        )

    ################ Update realized load growth and peak load #########################################

    load_data = read_data(joinpath(simulation_dir, "timeseries_data_files", "Load", "load_$(iteration_year - 1).csv"))
    reserve_up_demand_data = read_data(joinpath(simulation_dir,  "timeseries_data_files", "Reserves", "reserve_up_$(iteration_year - 1).csv"))
    reserve_down_demand_data = read_data(joinpath(simulation_dir,  "timeseries_data_files", "Reserves", "reserve_down_$(iteration_year - 1).csv"))
    rep_load_data = read_data(joinpath(simulation_dir, "timeseries_data_files", "Load", "rep_load_$(iteration_year - 1).csv"))
    load_n_vg_data = read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))
    rep_reserve_up_data = read_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "rep_reserve_up_$(iteration_year - 1).csv"))
    rep_reserve_down_data = read_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "rep_reserve_down_$(iteration_year - 1).csv"))

    for (idx, z) in enumerate(zones)

        load_data[:, Symbol(idx)] =  load_data[:, Symbol(idx)] * (1 + load_growth[idx])
        load_data[:, "Year"] = fill(load_data[1, "Year"] + 1, DataFrames.nrow(load_data))

        rep_load_data[:, Symbol(idx)] =  rep_load_data[:, Symbol(idx)] * (1 + load_growth[idx])
        rep_load_data[:, "Year"] = fill(rep_load_data[1, "Year"] + 1, DataFrames.nrow(rep_load_data))

        load_n_vg_data[:, Symbol("load_zone_$(idx)")] = load_n_vg_data[:, Symbol("load_zone_$(idx)")] * (1 + load_growth[idx])
        load_n_vg_data[:, "Year"] = fill(load_n_vg_data[1, "Year"] + 1, DataFrames.nrow(load_n_vg_data))

        reserve_up_demand_data[:, Symbol(idx)] =  reserve_up_demand_data[:, Symbol(idx)] * (1 + reserveup_annual_increment[idx])
        reserve_up_demand_data[:, "Year"] = fill(reserve_up_demand_data[1, "Year"] + 1, DataFrames.nrow(reserve_up_demand_data))

        rep_reserve_up_data[:, Symbol(idx)] =  rep_reserve_up_data[:, Symbol(idx)] * (1 + reserveup_annual_increment[idx])
        rep_reserve_up_data[:, "Year"] = fill(rep_reserve_up_data[1, "Year"] + 1, DataFrames.nrow(rep_reserve_up_data))

        reserve_down_demand_data[:, Symbol(idx)] =  reserve_down_demand_data[:, Symbol(idx)] * (1 + reservedown_annual_increment[idx])
        reserve_down_demand_data[:, "Year"] = fill(reserve_down_demand_data[1, "Year"] + 1, DataFrames.nrow(reserve_down_demand_data))

        rep_reserve_down_data[:, Symbol(idx)] =  rep_reserve_down_data[:, Symbol(idx)] * (1 + reservedown_annual_increment[idx])
        rep_reserve_down_data[:, "Year"] = fill(rep_reserve_down_data[1, "Year"] + 1, DataFrames.nrow(rep_reserve_down_data))
    end

    # Write realized load and reserve demand data in a CSV file
    CSV.write(joinpath(simulation_dir, "timeseries_data_files", "Load", "load_$(iteration_year).csv"), load_data)
    CSV.write(joinpath(simulation_dir, "timeseries_data_files", "Load", "rep_load_$(iteration_year).csv"), rep_load_data)

    CSV.write(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"), load_n_vg_data)

    CSV.write(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "reserve_up_$(iteration_year).csv"), reserve_up_demand_data)
    CSV.write(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "rep_reserve_up_$(iteration_year).csv"), rep_reserve_up_data)

    CSV.write(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "reserve_down_$(iteration_year).csv"), reserve_down_demand_data)
    CSV.write(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "rep_reserve_down_$(iteration_year).csv"), rep_reserve_down_data)

    peak_load_new = (1 + average_capacity_growth) * peak_load
    set_peak_load!(simulation, peak_load_new)

    return market_prices, capacity_factors, reserve_up_perc, reserve_down_perc, capacity_accepted_bids, rec_accepted_bids
end

