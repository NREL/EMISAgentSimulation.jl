
"""
This function creates the realized market data based on actual market clearing.
"""
function create_realized_marketdata(simulation::AgentSimulation,
                                    sys_UC::Union{Nothing, PSY.System},
                                    sys_ED::Union{Nothing, PSY.System},
                                    market_names::Vector{Symbol},
                                    rps_target::String,
                                    reserve_penalty::String,
                                    ordc_curved::Bool,
                                    existing_projects::Vector{<: Project{<: BuildPhase}},
                                    capacity_market_projects::Vector{<: Project{<: BuildPhase}},
                                    capacity_forward_years::Int64,
                                    iteration_year::Int64,
                                    simulation_years::Int64,
                                    solver::JuMP.MOI.OptimizerWithAttributes,
                                    results_dir::String,
                                    current_siip_sim)

    num_invperiods = 1

    simulation_dir = get_data_dir(get_case(simulation))
    zones = get_zones(simulation)
    lines = get_lines(simulation)
    case = get_case(simulation)
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
                                                          sys_UC,
                                                          market_names,
                                                          num_invperiods,
                                                          load_growth,
                                                          existing_projects,
                                                          iteration_year)

    num_hours = size(zonal_load)[2]
    num_days = Int(num_hours/24)

    rec_perc_requirement = get_rec_requirement(simulation)[iteration_year]

    energy_price,
    reserve_price,
    inertia_price,
    capacity_factors,
    reserve_perc,
    inertia_perc,
    start_up_costs,
    shut_down_costs,
    energy_voll,
    reserve_voll,
    inertia_voll = energy_mkt_clearing(sys_UC, sys_ED, system, simulation_dir, load_growth, reserve_penalty, rec_perc_requirement, zones, num_days, iteration_year, get_da_resolution(get_case(simulation)), get_rt_resolution(get_case(simulation)), get_name(get_case(simulation)), solver, get_base_dir(get_case(simulation)),simulation,current_siip_sim)

    println("Clean energy requirement for this year is $(get_rec_requirement(simulation)[iteration_year] * 100) percent")
    total_production = 0.0
    total_cec_production = 0.0
    day = 0
    get_rt_resolution(get_case(simulation))
    for time in 1:Int(24*60/get_rt_resolution(get_case(simulation))):(Int(24*60/get_rt_resolution(get_case(simulation))) * 365)
        day += 1
        daily_total_production = 0.0
        daily_cec_production = 0.0
        for gen in get_all_techs(sys_ED)
            name = PSY.get_name(gen)
            if !(occursin("BA", string(PSY.get_prime_mover(gen)))) #!(occursin("BA", name))
                energy_production = sum(capacity_factors[name][time:time + Int(24*60/get_rt_resolution(get_case(simulation)))-1]) * get_device_size(gen)
                total_production += energy_production
                daily_total_production += energy_production
                if occursin("WT", string(PSY.get_prime_mover(gen))) || occursin("PVe", string(PSY.get_prime_mover(gen))) || occursin("HY", string(PSY.get_prime_mover(gen))) #occursin("WT", name) || occursin("WIND", name) || occursin("PV", name) || occursin("HY", name) || occursin("NU", name) || occursin("RE", name)
                    total_cec_production += energy_production
                    daily_cec_production += energy_production
                end
                if occursin("ST", string(PSY.get_prime_mover(gen)))
                    if occursin("NUCLEAR", string(PSY.get_fuel(gen))) 
                        total_cec_production += energy_production
                        daily_cec_production += energy_production
                    end
                end
            end
        end
        #println("Clean energy contribution for day $(day) is $(round(daily_cec_production * 100.0 / daily_total_production, digits = 2)) percent")
    end

    println("Total Annual clean energy contribution is $(round(total_cec_production * 100.0 / total_production, digits = 2)) percent")

    cet_achieved_ratio = round(total_cec_production / total_production, digits = 2) / get_rec_requirement(simulation)[iteration_year]

    # Replace energy_mkt_clearing(nothing, nothing, system, load_growth, zones, num_days, solver) with
    # energy_mkt_clearing(sys_UC, sys_ED, system, load_growth, zones, num_days, solver) to run SIIP production cost model

    # Create empty market prices struct
    market_prices = MarketPrices()

    set_energy_price!(market_prices, "realized", energy_price)

    set_reserve_price!(market_prices, "realized", reserve_price)

    if in(:Inertia, market_names)
        set_inertia_price!(market_prices, "realized", inertia_price)
    end

    ######## Capacity market clearing #####################################################################

    capacity_market_bool = false
    if in(:Capacity, market_names)
        capacity_market_bool = true
    end

    capacity_mkt_param_file = joinpath(simulation_dir, "markets_data", "Capacity.csv")
    peak_load = get_peak_load(simulation)
    capacity_annual_increment = load_growth

    average_load_growth = Statistics.mean(load_growth)

    capacity_supply_curve = Vector{Union{String, Float64}}[]

    delta_irm = get_delta_irm(get_resource_adequacy(simulation), iteration_year)
    irm_scalar = get_irm_scalar(get_case(simulation))

    for project in capacity_market_projects
        for product in get_products(project)
            capacity_supply_curve = update_capacity_supply_curve!(capacity_supply_curve, product, project)
        end
    end

    capacity_price =  AxisArrays.AxisArray(reshape([0.0], 1,), [1])
    capacity_accepted_bids = Dict("no_accepted_bids" => 0.0)

    if in(:Capacity, market_names) && iteration_year + capacity_forward_years - 1 <= simulation_years
        system_peak_load = (1 + average_load_growth) ^ (capacity_forward_years) * peak_load
        capacity_demand_curve = create_capacity_demand_curve(capacity_mkt_param_file, system_peak_load, irm_scalar, delta_irm, capacity_market_bool)

        sort!(capacity_supply_curve, by = x -> x[3])      # Sort capacity supply curve by capacity bid

        capacity_price, capacity_accepted_bids = capacity_market_clearing(capacity_demand_curve, capacity_supply_curve, solver)
    end

    set_capacity_price!(market_prices, "realized", capacity_price)

    ######### REC market clearing ############################################################################

    rec_market_bool = false
    if in(:REC, market_names)
        rec_market_bool = true
    end

    REC_mkt_params = read_data(joinpath(simulation_dir, "markets_data", "REC_$(rps_target)_RPS.csv"))
    pricecap_rec = REC_mkt_params.price_cap[1]
    rec_req = REC_mkt_params.rec_req[1] * rec_market_bool
    rec_annual_increment = REC_mkt_params.annual_increment[1] * rec_market_bool
    rec_non_binding_years = REC_mkt_params.non_binding_years[1] * rec_market_bool

    rec_price = AxisArrays.AxisArray(reshape([pricecap_rec], 1,), [1])
    rec_accepted_bids = Dict{String, Float64}()

    total_demand = 0.0
    for z in zones
        for t in 1:num_hours
            total_demand += zonal_load[z, t] * (1 + energy_annual_increment[z]) * hour_weight[t]
        end
    end

    total_clean_production = 0.0

    rec_supply_curve = Vector{Union{String, Float64}}[]

    total_storage_consumption = 0.0
    for project in existing_projects
         # Populate REC market supply curves
         clean_production = 0.0
         for product in get_products(project)
            rec_supply_curve = update_rec_supply_curve!(rec_supply_curve, product, project)
            clean_production += find_clean_energy_production(product, project)
            total_storage_consumption += find_storage_energy_consumption(product, project)
        end
        total_clean_production += clean_production
    end

    clean_energy_percentage = min(1.0, (total_clean_production / total_demand))
    #println(clean_energy_percentage)

    if in(:REC, market_names)
        if length(rec_supply_curve) >= 1

            rec_energy_requirment =  total_demand * min(rec_req + (rec_annual_increment * iteration_year), 1)

            sort!(rec_supply_curve, by = x -> x[3])      # Sort REC supply curve by REC bid

            rec_energy_requirment = min(total_clean_production, rec_energy_requirment)
            rec_price, rec_accepted_bids = rec_market_clearing_non_binding(rec_energy_requirment, pricecap_rec, rec_supply_curve, solver)

            # if iteration_year <= rec_non_binding_years

            #     rec_energy_requirment = min(total_clean_production, rec_energy_requirment)
            #     #println(rec_energy_requirment)
            #     rec_price, rec_accepted_bids = rec_market_clearing_non_binding(rec_energy_requirment, pricecap_rec, rec_supply_curve, solver)
            # else
            #     total = 0
            #     for i in rec_supply_curve
            #         total += i[2]
            #     end
            #     #println(total)
            #     rec_energy_requirment = min(total_clean_production, rec_energy_requirment)
            #     #println(rec_energy_requirment)
            #     rec_price, rec_accepted_bids = rec_market_clearing_binding(rec_energy_requirment, pricecap_rec, rec_supply_curve, solver)
            # end

            set_rec_price!(market_prices, "realized", rec_price)
        end
    end

    ################# Write actual market clearing data ################################################

    output_file = joinpath(results_dir, "realized_market_data", "year_$(iteration_year).jld2")

    FileIO.save(output_file,
                     "capacity_price", capacity_price,
                     "energy_price", energy_price,
                     "reserve_price", reserve_price,
                     "rec_price", rec_price,
                     "inertia_price", inertia_price,
                     "capacity_factors", capacity_factors,
                     "reserve_perc", reserve_perc,
                     "capacity_accepted_bids", capacity_accepted_bids,
                     "rec_accepted_bids", rec_accepted_bids,
                     "inertia_perc", inertia_perc,
                     "start_up_costs", start_up_costs,
                     "shut_down_costs", shut_down_costs,
                     "energy_voll", energy_voll,
                     "reserve_voll", reserve_voll,
                     "inertia_voll", inertia_voll,
                     "rec_supply_curve", rec_supply_curve,
                     "rec_energy_requirment", rec_energy_requirment,
                     "cet_achieved_ratio", cet_achieved_ratio
        )

    ################ Update realized load growth and peak load #########################################

    load_data = read_data(joinpath(simulation_dir, "timeseries_data_files", "Load", "load_$(iteration_year - 1).csv"))
    rep_load_data = read_data(joinpath(simulation_dir, "timeseries_data_files", "Load", "rep_load_$(iteration_year - 1).csv"))
    load_n_vg_data = read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"))
    load_n_vg_data_rt = read_data(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data_rt.csv"))


    for (idx, z) in enumerate(zones)

        load_data[:, Symbol(idx)] =  load_data[:, Symbol(idx)] * (1 + load_growth[idx])
        #load_data[:, "Year"] = fill(load_data[1, "Year"] + 1, DataFrames.nrow(load_data))

        rep_load_data[:, Symbol(idx)] =  rep_load_data[:, Symbol(idx)] * (1 + load_growth[idx])
        #rep_load_data[:, "Year"] = fill(rep_load_data[1, "Year"] + 1, DataFrames.nrow(rep_load_data))

        load_n_vg_data[:, Symbol("load_zone_$(idx)")] = load_n_vg_data[:, Symbol("load_zone_$(idx)")] * (1 + load_growth[idx])
        load_n_vg_data_rt[:, Symbol("load_zone_$(idx)")] = load_n_vg_data_rt[:, Symbol("load_zone_$(idx)")] * (1 + load_growth[idx])
        #load_n_vg_data[:, "Year"] = fill(load_n_vg_data[1, "Year"] + 1, DataFrames.nrow(load_n_vg_data))
    end

    # Write realized load and reserve demand data in a CSV file
    CSV.write(joinpath(simulation_dir, "timeseries_data_files", "Load", "load_$(iteration_year).csv"), load_data)
    CSV.write(joinpath(simulation_dir, "timeseries_data_files", "Load", "rep_load_$(iteration_year).csv"), rep_load_data)
    CSV.write(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data.csv"), load_n_vg_data)
    CSV.write(joinpath(simulation_dir, "timeseries_data_files", "Net Load Data", "load_n_vg_data_rt.csv"), load_n_vg_data_rt)

    reserve_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"all_products"], "; ")
    ordc_products = split(read_data(joinpath(simulation_dir, "markets_data", "reserve_products.csv"))[1,"ordc_products"], "; ")
    non_ordc_products = filter(p -> !(p in ordc_products), reserve_products)

    reserve_timeseries_data = Dict(r => read_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "$(r)_$(iteration_year - 1).csv")) for r in non_ordc_products)
    rep_reserve_timeseries_data = Dict(r => read_data(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "rep_$(r)_$(iteration_year - 1).csv")) for r in non_ordc_products)

    for product in non_ordc_products
        reserve_timeseries_data[product][:, product] = reserve_timeseries_data[product][:, product] * (1 + average_load_growth)
        #reserve_timeseries_data[product][:, "Year"] = fill(reserve_timeseries_data[product][1, "Year"] + 1, DataFrames.nrow(reserve_timeseries_data[product]))

        rep_reserve_timeseries_data[product][:, product] = rep_reserve_timeseries_data[product][:, product] * (1 + average_load_growth)
        #rep_reserve_timeseries_data[product][:, "Year"] = fill(rep_reserve_timeseries_data[product][1, "Year"] + 1, DataFrames.nrow(rep_reserve_timeseries_data[product]))

        CSV.write(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "$(product)_$(iteration_year).csv"), reserve_timeseries_data[product])
        CSV.write(joinpath(simulation_dir, "timeseries_data_files", "Reserves", "rep_$(product)_$(iteration_year).csv"), rep_reserve_timeseries_data[product])
    end

    ordc_unavailability_method = get_ordc_unavailability_method(get_case(simulation))

    construct_ordc(deepcopy(sys_UC), simulation_dir, get_investors(simulation), iteration_year, get_rep_days(simulation), ordc_curved, ordc_unavailability_method, reserve_penalty)

    peak_load_new = (1 + average_load_growth) * peak_load
    set_peak_load!(simulation, peak_load_new)

    return market_prices, capacity_factors, reserve_perc, inertia_perc, capacity_accepted_bids, rec_accepted_bids, clean_energy_percentage, cet_achieved_ratio
end
