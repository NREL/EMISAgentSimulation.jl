
function cem(system::MarketClearingProblem{Z, T},
    solver::JuMP.MOI.OptimizerWithAttributes,
    resultfile::String="") where {Z, T}

    lines = [line.name for line in system.lines] # Lines
    projects = [project.name for project in system.projects] # Projects
    invperiods = 1:length(system.inv_periods)     # Investment periods
    opperiods = 1:T
    rep_hour_weight = system.rep_hour_weight
    capital_cost_multiplier = system.capital_cost_multiplier

    end_of_day = collect(24:24:T)

    tech_type,                                # technology type of project
    project_type,                             # is the project of storage type or generator type
    marginal_energy_cost,                      # $/MW
    marginal_reserveup_cost,                   # $/MW
    marginal_reservedown_cost,                 # $/MW
    fixed_cost,                                # $/unit
    queue_cost,                                # $/unit/year
    expansion_cost,                            # $/unit
    project_discount_rate,                    # discount rate
    min_gen, max_gen,                           # MW/unit
    min_input, max_input,                       # MW/unit
    efficiency_in, efficiency_out,            # percentage
    min_storage, max_storage,                   # MWh/unit
    init_storage,                              # MWh/unit
    availability,                             # percentage
    derating_factor,                           # percentage
    ramp_limits,                               # MW/unit/hour
    max_reserveup, max_reservedown,             # MW/unit/hour
    existing_units,                            # existing units
    units_in_queue,                             # units in queue
    build_lead_time, remaining_buildtime,       # construction lag times
    max_new_options,                            # maximum units
    base_cost_units,                            # base cost units (for capital cost multiplier)
    capex_years,                               # for investment cost annualization
    life_time,                                 # total possible life_time
    remaining_life_time,                       # remaining life_time
    capacity_eligible,                        # project eligible for capacity market
    rec_eligible,                             # project eligible for REC market
    location =                                # project zone location
    make_parameter_vector.(
    Ref(system.projects), :name,
    [:tech_type, :project_type, :marginal_energy_cost, :marginal_reserveup_cost, :marginal_reservedown_cost,
    :fixed_cost, :queue_cost, :expansion_cost, :discount_rate,
    :min_gen, :max_gen, :min_input, :max_input,
    :efficiency_in, :efficiency_out, :min_storage, :max_storage, :init_storage,
    :availability, :derating_factor,
    :ramp_limits, :max_reserveup, :max_reservedown,
    :existing_units, :units_in_queue, :build_lead_time, :remaining_build_time,
    :max_new_options, :base_cost_units, :capex_years, :life_time, :remaining_life,
    :rec_eligible, :capacity_eligible, :zone])
    life_range = AxisArrays.AxisArray(convert.(Int64, ones(length(projects))), projects)

    tech_types = unique(tech_type)

    option_projects = String[]
    decided_projects = String[]

    option_projects_by_type = Dict(t => Vector{String}() for t in tech_types)

    for g in projects
        if sum(units_in_queue[g]) <= 0
            push!(option_projects, g)

            for t in tech_types
                if tech_type[g] == t
                    push!(option_projects_by_type[t], g)
                end
            end

            if remaining_life_time[g] >= length(invperiods)
                life_range[g] = length(invperiods) - 1
            else
                life_range[g] = remaining_life_time[g] - 1
            end
        else
            push!(decided_projects, g)
        end
    end

    social_discount_rate = maximum(project_discount_rate)   # Social discount rate is kept as the max of all projects' discount rates

    social_npv_array = [(1 / (1 + social_discount_rate)) ^ p for p in invperiods]

    α = AxisArrays.AxisArray(zeros(length(projects)), projects)
    adjusted_investment_cost = AxisArrays.AxisArray(zeros(length(projects), length(invperiods)), projects, invperiods)
    adjusted_queue_cost = AxisArrays.AxisArray(zeros(length(projects)), projects)

    annualized_cap_cost = AxisArrays.AxisArray(zeros(length(projects), length(invperiods)), projects, invperiods)

    for g in projects
        α[g] = project_discount_rate[g] / (1 - (1 + project_discount_rate[g]) ^ (- (capex_years[g])))

        for q in 1:length(queue_cost[g])
            adjusted_queue_cost[g] += queue_cost[g][q] * (1 + project_discount_rate[g]) ^ (build_lead_time[g] + (length(queue_cost[g]) - q))
        end
        for p in invperiods
            adjusted_investment_cost[g, p] = expansion_cost[g][p] * (1 + project_discount_rate[g]) ^ (build_lead_time[g] - 1)
            annualized_cap_cost[g, p] = α[g] * (adjusted_investment_cost[g, p])

        end
    end

    capcost_segmentsize, capcost_segmentgrad, capcost_price_points, capcost_numsegments = make_capital_cost_curve(option_projects_by_type,
                                                                                                              annualized_cap_cost,
                                                                                                              base_cost_units,
                                                                                                              max_new_options,
                                                                                                              capital_cost_multiplier)

    # Populate markets data
    zones = system.zones
    price_cap_e = AxisArrays.AxisArray(zeros(length(zones), length(invperiods)), zones, invperiods)
    demand_e = AxisArrays.AxisArray(zeros(length(zones), length(invperiods), T), zones, invperiods, opperiods)

    demand_ru = AxisArrays.AxisArray(zeros(length(zones), length(invperiods), T), zones, invperiods, opperiods)

    price_cap_rd = AxisArrays.AxisArray(zeros(length(zones), length(invperiods)), zones, invperiods)
    demand_rd = AxisArrays.AxisArray(zeros(length(zones), length(invperiods), T), zones, invperiods, opperiods)

    price_cap_rec = AxisArrays.AxisArray(zeros(length(invperiods)), invperiods)
    rec_requirement = AxisArrays.AxisArray(zeros(length(invperiods)), invperiods)

    for p in invperiods
        price_cap_rec[p] = getproperty(getproperty(system.inv_periods[p], :rec_market), :price_cap)
        rec_requirement[p] = getproperty(getproperty(system.inv_periods[p], :rec_market), :rec_req)
        for z in zones
            price_cap_e[z, p] = getproperty(getproperty(system.inv_periods[p], :energy_market), :price_cap)[z]
            demand_e[z, p, :] = getproperty(getproperty(system.inv_periods[p], :energy_market), :demand)[z, :]

            demand_ru[z, p, :] = getproperty(getproperty(system.inv_periods[p], :reserveup_market), :demand)[z, :]

            price_cap_rd[z, p] = getproperty(getproperty(system.inv_periods[p], :reservedown_market), :price_cap)[z]
            demand_rd[z, p, :] = getproperty(getproperty(system.inv_periods[p], :reservedown_market), :demand)[z, :]
        end
    end

    ru_segmentsize, ru_segmentgrad, ru_price_points, ru_numsegments = make_ORDC_vectors(getproperty.(system.inv_periods, :reserveup_market))

    cap_segmentsize, cap_segmentgrad, cap_price_points, cap_numsegments = make_capacity_demand_vectors(getproperty.(system.inv_periods, :capacity_market))

    capacity_mkt_projects = Vector{String}()
    rps_compliant_projects = Vector{String}()

    generator_projects = Vector{String}()
    storage_projects = Vector{String}()

    zone_projects = Dict(z => Vector{String}() for z in zones)
    zone_storage = Dict(z => Vector{String}() for z in zones)

    ramp_lim_projects = Vector{String}()

    for project in projects

        if ramp_limits[project] !== nothing
            push!(ramp_lim_projects, project)
        end

        for zone in zones
            if location[project] == zone
                push!(zone_projects[zone], project)

                if project_type[project] == "generator"
                    push!(generator_projects, project)
                elseif project_type[project] == "storage"
                    push!(storage_projects, project)
                    push!(zone_storage[zone], project)
                end
            end
        end

        if capacity_eligible[project]
            push!(capacity_mkt_projects, project) # Find Capacity market eligible projects
        end

        if rec_eligible[project]
            push!(rps_compliant_projects, project) # Find REC eligible projects
        end
    end


    # Populate line data
    from_zone,
    to_zone,
    linepowerlimit =
    make_parameter_vector.(
    Ref(system.lines), :name,
    [:from_zone, :to_zone, :active_power_limit]
    )

    lines_from_zone = Dict(z => Vector{String}() for z in zones)
    lines_to_zone = Dict(z => Vector{String}() for z in zones)

    for line in lines
        for zone in zones
            if from_zone[line] == zone
                push!(lines_from_zone[zone], line)
            elseif to_zone[line] == zone
                push!(lines_to_zone[zone], line)
            end
        end
    end

    m = JuMP.Model(solver)
    # INVESTMENT DECISION VARIABLES

    JuMP.@variable(m, n[g in projects, p in invperiods] >= 0)  # New units built
    JuMP.@variable(m, n_by_type[type in tech_types, p in invperiods, s in 1:capcost_numsegments[type, p]] >= 0) # New units by technology type for capital cost multiplier
    JuMP.@variable(m, r[g in projects, p in invperiods] >= 0)  # Units retired
    JuMP.@variable(m, c[g in projects, p in invperiods] >= 0)  # Unit cleared in capacity market
    JuMP.@variable(m, temp == 0) # For avoiding add_to_expression errors

    # ORDC DEMAND VARIABLE
    JuMP.@variable(m, d_ru[z in zones, p in invperiods, t in opperiods, s in 1:ru_numsegments[z, p]] >= 0)

    # CAPACITY MARKET ELASTIC DEMAND
    JuMP.@variable(m, d_cap[p in invperiods, s in 1:cap_numsegments[p]] >= 0)

    # OPERATION DECISION VARIABLES

    JuMP.@variable(m, p_e[g in projects, p in invperiods, t in opperiods] >= 0) # Unit energy production [MW]
    JuMP.@variable(m, p_in[g in storage_projects, p in invperiods, t in opperiods] >= 0) # Storage charging [MW]
    JuMP.@variable(m, p_ru[g in projects, p in invperiods, t in opperiods] >= 0) # reserve up provided [MW]
    JuMP.@variable(m, p_rd[g in projects, p in invperiods, t in opperiods] >= 0) # reserve down provided [MW]

    JuMP.@variable(m, p_in_ru[g in storage_projects, p in invperiods, t in opperiods] >= 0) # reserve up provided by storage charging [MW]
    JuMP.@variable(m, p_in_rd[g in storage_projects, p in invperiods, t in opperiods] >= 0) # reserve down provided by storage charging [MW]

    JuMP.@variable(m, p_out_ru[g in storage_projects, p in invperiods, t in opperiods] >= 0) # reserve up provided by storage discharging [MW]
    JuMP.@variable(m, p_out_rd[g in storage_projects, p in invperiods, t in opperiods] >= 0) # reserve down provided by storage discharging [MW]

    JuMP.@variable(m, flow[l in lines, p in invperiods, t in opperiods])     # Line flow

    JuMP.@variable(m, v_e[z in zones, p in invperiods, t in opperiods] >= 0) # Load shortfall [MW]
    JuMP.@variable(m, v_rd[z in zones, p in invperiods, t in opperiods] >= 0) # reserve down shortfall [MW]
    JuMP.@variable(m, v_rec[p in invperiods] >= 0) # REC shortfall [MWh]


    # Populate units dispatchable expression based on remaining build/construction time
    JuMP.@expression(m, unitsdispatchable[g in projects, p in invperiods], 0
                            + sum(n[g, i] for i in 1:(p - remaining_buildtime[g]))
                           - r[g, p])

    JuMP.@expression(m, unitsmustretire[g in projects, p in invperiods], 0 + temp)

    JuMP.@expression(m, yearly_cap_cost_decided[g in decided_projects, p in invperiods], 0 + temp)

    JuMP.@expression(m, yearly_cap_cost_options[type in tech_types, p in invperiods], 0 + temp ^ 2)

    JuMP.@expression(m, yearly_queue_cost[g in projects, p in invperiods], 0 + temp)

    # Storage level evolution
    JuMP.@expression(m, storage_level[g in storage_projects, p in invperiods, t in opperiods],
                    p_in[g, p, t] * efficiency_in[g] - p_e[g, p, t] / efficiency_out[g])

    for g in projects

        min_life_time = min(life_time[g], remaining_life_time[g])
        remaining_queuetime = max(0, remaining_buildtime[g] - build_lead_time[g])
        total_queuetime = length(queue_cost[g])

        # Units must retire after their life_time
        for p in invperiods[1:end - min_life_time]
            JuMP.add_to_expression!(unitsmustretire[g, p + min_life_time], unitsdispatchable[g, p] + r[g, p])
        end

        # Calculate yearly capital costs allowing different technology costs each year
        for p in invperiods
            start_year = max(1, p - remaining_buildtime[g] - min_life_time + 1)

                for i in start_year:(p - remaining_buildtime[g])

                    JuMP.add_to_expression!(yearly_queue_cost[g, p], α[g] * adjusted_queue_cost[g] * n[g, i])
                    if in(g, decided_projects)
                        JuMP.add_to_expression!(yearly_cap_cost_decided[g, p], n[g, i] * annualized_cap_cost[g, i + remaining_queuetime])
                    end

                end
            if in(g, storage_projects)
                for t in opperiods
                    if t == 1 || in(t, end_of_day + ones(length(end_of_day)))
                        JuMP.add_to_expression!(storage_level[g, p, t], init_storage[g] * max_storage[g] * unitsdispatchable[g, p])
                    else
                        JuMP.add_to_expression!(storage_level[g, p, t], storage_level[g, p, t - 1])
                    end
                end
            end

        end

    end

    for type in tech_types
        if length(option_projects_by_type[type]) >= 1
            g = option_projects_by_type[type][1]
            min_life_time = min(life_time[g], remaining_life_time[g])
            remaining_queuetime = max(0, remaining_buildtime[g] - build_lead_time[g])
            for p in invperiods
                start_year = max(1, p - remaining_buildtime[g] - min_life_time + 1)
                for i in start_year:(p - remaining_buildtime[g])
                    JuMP.add_to_expression!(yearly_cap_cost_options[type, p],
                        sum(n_by_type[type, i, s] * capcost_price_points[type, i + remaining_queuetime][s] for s in 1:capcost_numsegments[type, i])
                        + 1/2 * sum(capcost_segmentgrad[type, i + remaining_queuetime][s] * (n_by_type[type, i, s]^2) for s in 1:capcost_numsegments[type, i]))
                end
            end
        end
    end

    JuMP.@expression(m, annual_project_costs[g in projects, p in invperiods], 0
                          + yearly_queue_cost[g, p]                               # yearly queue_cost
                          + fixed_cost[g] * unitsdispatchable[g, p]                # fixed operating costs
                          + sum((p_e[g, p, t] * marginal_energy_cost[g]               # project operating costs
                              + p_ru[g, p, t] * marginal_reserveup_cost[g]
                              + p_rd[g, p, t] * marginal_reservedown_cost[g]) * rep_hour_weight[t] for t in opperiods))

    # OBJECTIVE FUNCTION

    JuMP.@objective(m, Max,

    # Capacity market welfare
    sum(sum(d_cap[p, s] * cap_price_points[p][s] for s in 1:cap_numsegments[p]) * social_npv_array[p] for p in invperiods)
    + 1/2 * sum(sum(cap_segmentgrad[p][s] * d_cap[p, s]^2 for s in 1:cap_numsegments[p]) * social_npv_array[p] for p in invperiods)

    # ORDC welfare
    + sum(sum(d_ru[z, p, t, s] * demand_ru[z, p, t] * rep_hour_weight[t] * ru_price_points[z, p][s] for z in zones, t in opperiods, s in 1:ru_numsegments[z, p])
                * social_npv_array[p] for p in invperiods)
    + 1/2 * sum(sum(ru_segmentgrad[z, p][s] * (d_ru[z, p, t, s] * demand_ru[z, p, t])^2 * rep_hour_weight[t] for z in zones, t in opperiods, s in 1:ru_numsegments[z, p])
                * social_npv_array[p] for p in invperiods)

    # REC market Alternative Compliance Payment
    - sum(v_rec[p] * price_cap_rec[p] * social_npv_array[p] for p in invperiods)

    - sum(yearly_cap_cost_decided[g, p] * social_npv_array[p] for g in decided_projects, p in invperiods)

    - sum(yearly_cap_cost_options[type, p] * social_npv_array[p] for type in tech_types, p in invperiods)

    # Annual costs incurred by the project
    - sum(annual_project_costs[g, p] * social_npv_array[p]
        for g in projects, p in invperiods)

    +  (# weigh operations according to actual duration
    # operations market shortfall penalties
             - sum((v_e[z, p, t] * price_cap_e[z, p] + v_rd[z, p, t] * price_cap_rd[z, p]) * rep_hour_weight[t] * social_npv_array[p]
                for z in zones, p in invperiods, t in opperiods)

                        )
        )


    # CONSTRAINTS

    JuMP.@constraint(m, maxoptions[g in projects],
        sum(n[g, p] for p in invperiods) - max_new_options[g] <= 0) # Maximum options purchased

    for p in invperiods

        for type in tech_types
            JuMP.@constraint(m, sum(n_by_type[type, p, s] for s in 1:capcost_numsegments[type, p]) == sum(n[g, p] for g in option_projects_by_type[type]))

            for s in 1:capcost_numsegments[type, p]
                JuMP.@constraint(m, n_by_type[type, p, s] <= capcost_segmentsize[type, p][s])
            end

        end
        # Planning Constraints
        for g in projects
            JuMP.@constraint(m, unitsmustretire[g, p] - r[g, p] <= 0) # Must retire units

            if p < length(invperiods)
                JuMP.@constraint(m, r[g, p] - r[g, p + 1] <= 0) # Units should not become dispatchable after being retired
            end

            if p >= remaining_life_time[g] + 1
                JuMP.@constraint(m, unitsdispatchable[g, p] == 0) # Units should not be dispatchable beyond remaining life_time
            end

            if in(g, capacity_mkt_projects)
                JuMP.@constraint(m, c[g, p] - unitsdispatchable[g, p] <= 0) # Maximum capacity cleared in the capacity market
            end

        end

        # Capacity Demand constraints
        for s in 1:cap_numsegments[p]
            JuMP.@constraint(m, d_cap[p, s] <= cap_segmentsize[p][s]) # Capacity market demand curve segment limit
        end

        for t in opperiods
            # Unit operating constraints (dispatch, reserves and ramping)
            for g in projects
                if in(g, generator_projects)
                    # Generator technical constraints
                    JuMP.@constraint(m, p_e[g, p, t] - p_rd[g, p, t] >= unitsdispatchable[g, p] * 0.0) # Minimum dispatch
                    JuMP.@constraint(m, p_e[g, p, t] + p_ru[g, p, t] <= unitsdispatchable[g, p] * max_gen[g] * availability[g][p, t]) # Maximum dispatch
                    JuMP.@constraint(m, p_rd[g, p, t] <= unitsdispatchable[g, p] * max_reservedown[g]) # Maximum reserve down
                    JuMP.@constraint(m, p_ru[g, p, t] <= unitsdispatchable[g, p] * max_reserveup[g]) # Maximum reserve up

                elseif in(g, storage_projects)
                    # Storage technical constraints
                    JuMP.@constraint(m, p_e[g, p, t] - p_out_rd[g, p, t] >= unitsdispatchable[g, p] * 0.0) # Minimum output dispatch
                    JuMP.@constraint(m, p_e[g, p, t] + p_out_ru[g, p, t] <= unitsdispatchable[g, p] * max_gen[g]) # Maximum output dispatch

                    JuMP.@constraint(m, p_in[g, p, t] - p_in_ru[g, p, t] >= unitsdispatchable[g, p] * 0.0) # Minimum input dispatch
                    JuMP.@constraint(m, p_in[g, p, t] + p_in_rd[g, p, t] <= unitsdispatchable[g, p] * max_input[g]) # Maximum input dispatch

                    JuMP.@constraint(m, p_ru[g, p, t] <= p_out_ru[g, p, t] + p_in_ru[g, p, t]) # Total storage reserve up
                    JuMP.@constraint(m, p_rd[g, p, t] <= p_out_rd[g, p, t] + p_in_rd[g, p, t]) # Total storage reserve down

                    JuMP.@constraint(m, storage_level[g, p, t] >= unitsdispatchable[g, p] * min_storage[g]) # Minimum storage level
                    JuMP.@constraint(m, storage_level[g, p, t] <= unitsdispatchable[g, p] * max_storage[g]) # Maximum storage level

                    if in(t, end_of_day)
                        JuMP.@constraint(m, storage_level[g, p, t] == init_storage[g] * max_storage[g] * unitsdispatchable[g, p])
                    end
                end
            end

            # Line flow constraints
            for l in lines
                JuMP.@constraint(m, flow[l, p, t] >= -linepowerlimit[l])
                JuMP.@constraint(m, flow[l, p, t] <= linepowerlimit[l])
            end

            for z in zones
                for s in 1:ru_numsegments[z, p]
                    JuMP.@constraint(m, d_ru[z, p, t, s] <= ru_segmentsize[z, p][s]) #ORDC segment limit
                end
            end
        end
    end

    # Market Clearing Constraints:

    # Energy market
    JuMP.@constraint(m, energy_market[z in zones, p in invperiods, t in opperiods], # Power balance
        sum(p_e[g ,p, t] for g in zone_projects[z])
        + sum(flow[l, p, t] for l in lines_to_zone[z])
        - sum(flow[l, p, t] for l in lines_from_zone[z])
        + v_e[z, p, t] >= demand_e[z, p, t] + sum(p_in[g, p, t] for g in zone_storage[z]))

    # reserve up market
    JuMP.@constraint(m, reserve_up_market[z in zones, p in invperiods, t in opperiods], # reserve up balance
        sum(p_ru[g ,p, t] for g in zone_projects[z]) == sum(d_ru[z, p, t, s] for s in 1:ru_numsegments[z, p]) * demand_ru[z, p, t] )


    # reserve down market
    JuMP.@constraint(m, reserve_down_market[z in zones, p in invperiods, t in opperiods], # reserve down balance
        sum(p_rd[g ,p, t] for g in zone_projects[z]) + v_rd[z, p, t]  == demand_rd[z, p, t])

    # Capacity market
    JuMP.@constraint(m, capacity_market[p in invperiods],
        sum(c[g, p] * max_gen[g] * derating_factor[g] for g in capacity_mkt_projects) == sum(d_cap[p, s] for s in 1:cap_numsegments[p]))


    # REC market
    JuMP.@constraint(m, rps_compliance[p in invperiods],
        sum(p_e[g, p, t] * rep_hour_weight[t] for g in rps_compliant_projects, t in opperiods) + v_rec[p]
        >= sum(demand_e[z, p, t] * rep_hour_weight[t] for z in zones, t in opperiods) * rec_requirement[p])



    println("Price Projection:")
    @time JuMP.optimize!(m)
    println(JuMP.termination_status(m))
    println(JuMP.objective_value(m))

    #Post processing--------------------------------------------------------------------------------------------------------
    capacity_price = JuMP.dual.(capacity_market)
    energy_price = JuMP.dual.(energy_market)
    reserve_up_price = JuMP.dual.(reserve_up_market)
    reserve_down_price = JuMP.dual.(reserve_down_market)
    REC_price = JuMP.dual.(rps_compliance)

    capacity_factor = Dict([g => zeros(length(invperiods), length(opperiods)) for g in projects])
    capacity_accepted_perc = Dict([g => zeros(length(invperiods)) for g in projects])

    # Number of new options invested in the first year of the CEM.
    new_options = Dict([g => 0 for g in option_projects])

    ϵ = 2e-1

    for g in projects
        for p in invperiods
            if value.(unitsdispatchable[g, p]) < 0.0001
                capacity_accepted_perc[g][p] = 0.0
            else
                capacity_accepted_perc[g][p] = value.(c[g, p]) / value.(unitsdispatchable[g, p])
            end
            for t in opperiods
                if value.(unitsdispatchable[g, p]) < 0.0001
                    capacity_factor[g][p, t] = 0.0
                else
                    capacity_factor[g][p, t] = value.(p_e[g, p, t]) /(value.(unitsdispatchable[g, p]) * max_gen[g])
                end
            end

        end

        if in(g, option_projects)
            for i in 0:base_cost_units[g]:max_new_options[g] + base_cost_units[g]
                if value.(n[g, 1]) > i + ϵ
                    new_options[g] += 1
                end
            end
        end
    end

    nominal_capacity_price = AxisArrays.AxisArray(zeros(length(capacity_price)), invperiods)
    nominal_REC_price = AxisArrays.AxisArray(zeros(length(REC_price)), invperiods)

    nominal_energy_price = AxisArrays.AxisArray(zeros(size(energy_price)), zones, invperiods, opperiods)
    nominal_reserve_up_price = AxisArrays.AxisArray(zeros(size(reserve_up_price)), zones, invperiods, opperiods)
    nominal_reserve_down_price = AxisArrays.AxisArray(zeros(size(reserve_down_price)), zones, invperiods, opperiods)


    for p in invperiods
        nominal_capacity_price[p] = capacity_price[p] / social_npv_array[p]
        nominal_REC_price[p] = REC_price[p] / social_npv_array[p]
        for t in opperiods
            nominal_energy_price[:, p, t] = (energy_price[:, p, t] / rep_hour_weight[t]) / social_npv_array[p]
            nominal_reserve_up_price[:, p, t] = (reserve_up_price[:, p, t] / rep_hour_weight[t]) / social_npv_array[p]
            nominal_reserve_down_price[:, p, t] = (reserve_down_price[:, p, t] / rep_hour_weight[t]) / social_npv_array[p]
        end
    end

    for z in zones
        for p in invperiods
            for t in opperiods
                if nominal_energy_price[z, p, t] <= 1e-5
                    nominal_energy_price[z, p, t] = 1e-5 # replace zeros with 1e-5, to avoid no energy production by renewables when price is 0.
                end
            end
        end
    end

    return nominal_capacity_price,
           nominal_energy_price,
           nominal_reserve_up_price,
           nominal_reserve_down_price,
           nominal_REC_price,
           capacity_factor,
           capacity_accepted_perc,
           new_options;

end
