
function economicdispatch(system::MarketClearingProblem{Z, T},
        solver::JuMP.MOI.OptimizerWithAttributes,
        resultfile::String="") where {Z, T}

        lines = [line.name for line in system.lines] # Lines
        projects = [project.name for project in system.projects] # Projects
        invperiods = 1:length(system.inv_periods)     # Investment periods
        opperiods = 1:T
        rep_hour_weight = system.rep_hour_weight

        end_of_day = collect(24:24:T)

        project_type,                             # is the project of storage type
        marginal_energy_cost,                      # $/MW
        marginal_reserveup_cost,                # $/MW
        marginal_reservedown_cost,                # $/MW
        min_gen, max_gen,                           # MW/unit
        min_input, max_input,                       # MW/unit
        efficiency_in, efficiency_out,            # percentage
        min_storage, max_storage,                   # MWh/unit
        init_storage,                              # MWh/unit
        availability,                             # hourly availability factor
        ramp_limits,                               # MW/unit/hour
        max_reserveup, max_reservedown,             # MW/unit/hour
        location =                                # Project zone location
        make_parameter_vector.(
        Ref(system.projects), :name,
        [:project_type, :marginal_energy_cost, :marginal_reserveup_cost, :marginal_reservedown_cost,
        :min_gen, :max_gen, :min_input, :max_input,
        :efficiency_in, :efficiency_out, :min_storage, :max_storage, :init_storage,
        :availability,
        :ramp_limits, :max_reserveup, :max_reservedown, :zone]
        )

        zones = system.zones
        price_cap_e = AxisArrays.AxisArray(zeros(length(zones), length(invperiods)), zones, invperiods)
        demand_e = AxisArrays.AxisArray(zeros(length(zones), length(invperiods), T), zones, invperiods, opperiods)

        demand_ru = AxisArrays.AxisArray(zeros(length(zones), length(invperiods), T), zones, invperiods, opperiods)

        price_cap_rd = AxisArrays.AxisArray(zeros(length(zones), length(invperiods)), zones, invperiods)
        demand_rd = AxisArrays.AxisArray(zeros(length(zones), length(invperiods), T), zones, invperiods, opperiods)

        for p in invperiods
                for z in zones
                        price_cap_e[z, p] = getproperty(getproperty(system.inv_periods[p], :energy_market), :price_cap)[z]
                        demand_e[z, p, :] = getproperty(getproperty(system.inv_periods[p], :energy_market), :demand)[z, :]

                        demand_ru[z, p, :] = getproperty(getproperty(system.inv_periods[p], :reserveup_market), :demand)[z, :]

                        price_cap_rd[z, p] = getproperty(getproperty(system.inv_periods[p], :reservedown_market), :price_cap)[z]
                        demand_rd[z, p, :] = getproperty(getproperty(system.inv_periods[p], :reservedown_market), :demand)[z, :]
                end
        end

        ru_segmentsize, ru_segmentgrad, ru_price_points, ru_numsegments = make_ORDC_vectors(getproperty.(system.inv_periods, :reserveup_market))

        zone_projects = Dict(z => Vector{String}() for z in zones)

        generator_projects = Vector{String}()
        storage_projects = Vector{String}()

        zone_projects = Dict(z => Vector{String}() for z in zones)
        zone_storage = Dict(z => Vector{String}() for z in zones)


        ramp_lim_projects = Vector{String}()

        for project in projects

                if ramp_limits[project] !== nothing
                        push!(ramp_lim_projects, project)
                end

            for z in zones
                if location[project] == z
                    push!(zone_projects[z], project)

                    if project_type[project] == "generator"
                        push!(generator_projects, project)
                    elseif project_type[project] == "storage"
                        push!(storage_projects, project)
                        push!(zone_storage[z], project)
                    end

                end
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

        # ORDC DEMAND VARIABLE
        JuMP.@variable(m, d_ru[z in zones, p in invperiods, t in opperiods, s in 1:ru_numsegments[z, p]] >= 0)

        JuMP.@variable(m, v_e[z in zones, p in invperiods, t in opperiods] >= 0) # Load shortfall [MW]
        JuMP.@variable(m, v_rd[z in zones, p in invperiods, t in opperiods] >= 0) # Reserve Down shortfall [MW]

        JuMP.@expression(m, annual_project_costs[g in projects, p in invperiods], 0
                              + sum((p_e[g, p, t] * marginal_energy_cost[g]
                                  + p_ru[g, p, t] * marginal_reserveup_cost[g]
                                  + p_rd[g, p, t] * marginal_reservedown_cost[g]) * rep_hour_weight[t]  for t in opperiods))    # project operating costs

        # Storage level evolution
        JuMP.@expression(m, storage_level[g in storage_projects, p in invperiods, t in opperiods],
                    p_in[g, p, t] * efficiency_in[g] - p_e[g, p, t] / efficiency_out[g])


        for g in storage_projects
                for p in invperiods
                        for t in opperiods
                                if t == 1 || in(t, end_of_day + ones(length(end_of_day)))
                                        JuMP.add_to_expression!(storage_level[g, p, t], init_storage[g] * max_storage[g])
                                else
                                        JuMP.add_to_expression!(storage_level[g, p, t], storage_level[g, p, t - 1])
                                end
                        end
                end

        end


        # OBJECTIVE FUNCTION

        JuMP.@objective(m, Max,

        # ORDC welfare
        sum(sum(d_ru[z, p, t, s] * demand_ru[z, p, t] * ru_price_points[z, p][s] * rep_hour_weight[t] for z in zones, t in opperiods, s in 1:ru_numsegments[z, p])
                for p in invperiods)
        + 1/2 * sum(sum(ru_segmentgrad[z, p][s] * (d_ru[z, p, t, s] * demand_ru[z, p, t])^2 * rep_hour_weight[t] for z in zones, t in opperiods, s in 1:ru_numsegments[z, p])
                for p in invperiods)

        # Annual costs incurred by the projects
        - sum(annual_project_costs[g, p]
            for g in projects, p in invperiods)

        # weight operations according to actual duration
        # operations market shortfall penalties
        + ( - sum((v_e[z, p, t] * price_cap_e[z, p] + v_rd[z, p, t] * price_cap_rd[z, p]) * rep_hour_weight[t]
                for z in zones, p in invperiods, t in opperiods)
         )

            )

    # CONSTRAINTS

    # Operations:

        # Energy market
        JuMP.@constraint(m, energy_market[z in zones, p in invperiods, t in opperiods], # Power balance
                sum(p_e[g ,p, t] for g in zone_projects[z])
                + sum(flow[l, p, t] for l in lines_to_zone[z])
                - sum(flow[l, p, t] for l in lines_from_zone[z])
                + v_e[z, p, t] == demand_e[z, p, t] + sum(p_in[g, p, t] for g in zone_storage[z]))

        # reserve up market
        JuMP.@constraint(m, reserve_up_market[z in zones, p in invperiods, t in opperiods], # reserve up balance
                sum(p_ru[g ,p, t] for g in zone_projects[z]) == sum(d_ru[z, p, t, s] for s in 1:ru_numsegments[z, p]) * demand_ru[z, p, t] )

        JuMP.@constraint(m, ru_segment_limit[z in zones, p in invperiods, t in opperiods, s in 1:ru_numsegments[z, p]],
                d_ru[z, p, t, s] <= ru_segmentsize[z, p][s])

        # reserve down market
        JuMP.@constraint(m, reserve_down_market[z in zones, p in invperiods, t in opperiods], # reserve down balance
                sum(p_rd[g ,p, t] for g in zone_projects[z]) + v_rd[z, p, t]  == demand_rd[z, p, t])

        # Line flow constraints
        JuMP.@constraint(m, negative_flow[l in lines, p in invperiods, t in opperiods],
                        flow[l, p, t] >= -linepowerlimit[l])

        JuMP.@constraint(m, positive_flow[l in lines, p in invperiods, t in opperiods],
                        flow[l, p, t] <= linepowerlimit[l])

        # Unit operating constraints (dispatch and ramping)

        # Generator technical constraints
        JuMP.@constraint(m, min_generation[g in generator_projects, p in invperiods, t in opperiods],
                   p_e[g, p, t] - p_rd[g, p, t] >= 0.0) # Minimum dispatch

        JuMP.@constraint(m, max_generation[g in generator_projects, p in invperiods, t in opperiods],
                   p_e[g, p, t] + p_ru[g, p, t] <= max_gen[g] * availability[g][p, t]) # Maximum dispatch

        JuMP.@constraint(m, reserve_down_max[g in generator_projects, p in invperiods, t in opperiods],
                   p_rd[g, p, t] <= max_reservedown[g]) # Maximum Reserve Down

        JuMP.@constraint(m, reserve_up_max[g in generator_projects, p in invperiods, t in opperiods],
                   p_ru[g, p, t] <= max_reserveup[g]) # Maximum Reserve Up

        JuMP.@constraint(m, rampdown[g in ramp_lim_projects, p in invperiods, t in opperiods[2:end]], # Down-ramp limit
                      p_e[g, p ,t] - p_rd[g, p ,t] >= p_e[g, p, t-1] - ramp_limits[g][:down])

        JuMP.@constraint(m, rampup[g in ramp_lim_projects, p in invperiods, t in opperiods[2:end]], # Up-ramp limit
                          p_e[g, p, t] + p_ru[g, p, t] <= p_e[g, p, t-1] + ramp_limits[g][:up])

        # Storage technical constraints
        JuMP.@constraint(m, minoutput[g in storage_projects, p in invperiods, t in opperiods],
                   p_e[g, p, t] - p_out_rd[g, p, t] >= 0.0) # Minimum output dispatch

        JuMP.@constraint(m, maxoutput[g in storage_projects, p in invperiods, t in opperiods],
                   p_e[g, p, t] + p_out_ru[g, p, t] <= max_gen[g]) # Maximum output dispatch

        JuMP.@constraint(m, min_input[g in storage_projects, p in invperiods, t in opperiods],
                   p_in[g, p, t] - p_in_ru[g, p, t] >= 0.0) # Minimum input dispatch

        JuMP.@constraint(m, max_input[g in storage_projects, p in invperiods, t in opperiods],
                   p_in[g, p, t] + p_in_rd[g, p, t] <= max_input[g]) # Maximum input dispatch


        JuMP.@constraint(m, storage_reserve_up[g in storage_projects, p in invperiods, t in opperiods],
                   p_ru[g, p, t] <= p_out_ru[g, p, t] + p_in_ru[g, p, t]) # Maximum storage reserve up

        JuMP.@constraint(m, storage_reserve_down[g in storage_projects, p in invperiods, t in opperiods],
                   p_rd[g, p, t] <= p_out_rd[g, p, t] + p_in_rd[g, p, t]) # Maximum storage reserve down

        JuMP.@constraint(m, min_storage_level[g in storage_projects, p in invperiods, t in opperiods],
                   storage_level[g, p, t] >= min_storage[g]) # Minimum storage level

        JuMP.@constraint(m, max_storage_level[g in storage_projects, p in invperiods, t in opperiods],
                   storage_level[g, p, t] <= max_storage[g]) # Maximum storage level

        JuMP.@constraint(m, end_day_storage_level[g in storage_projects, p in invperiods, t in end_of_day],
                   storage_level[g, p, t] == init_storage[g] * max_storage[g]) # Maximum storage level after reserve provision


        println("Actual Energy and A/S Market Clearing:")
        JuMP.optimize!(m)
        println(JuMP.termination_status(m))
        println(JuMP.objective_value(m))

        energy_dual = JuMP.dual.(energy_market)
        reserve_up_dual = JuMP.dual.(reserve_up_market)
        reserve_down_dual = JuMP.dual.(reserve_down_market)

        energy_price = AxisArrays.AxisArray(zeros(size(energy_dual)), zones, invperiods, opperiods)
        reserve_up_price = AxisArrays.AxisArray(zeros(size(reserve_up_dual)), zones, invperiods, opperiods)
        reserve_down_price = AxisArrays.AxisArray(zeros(size(reserve_down_dual)), zones, invperiods, opperiods)


        for p in invperiods
                for t in opperiods
                        energy_price[:, p, t] = round.((energy_dual[:, p, t] / rep_hour_weight[t]), digits = 5)
                        reserve_up_price[:, p, t] = round.((reserve_up_dual[:, p, t] / rep_hour_weight[t]), digits = 5)
                        reserve_down_price[:, p, t] = round.((reserve_down_dual[:, p, t] / rep_hour_weight[t]), digits = 5)
                end
        end

        capacity_factor = Dict([g => zeros(length(invperiods), length(opperiods)) for g in projects])
        reserve_up_perc = Dict([g => zeros(length(invperiods), length(opperiods)) for g in projects])
        reserve_down_perc = Dict([g => zeros(length(invperiods), length(opperiods)) for g in projects])

        for g in projects
                for p in invperiods
                        for t in opperiods
                                capacity_factor[g][p, t] = value.(p_e[g, p, t]) / max_gen[g]
                                reserve_up_perc[g][p, t] = value.(p_ru[g, p, t]) / max_gen[g]
                                reserve_down_perc[g][p, t] = value.(p_rd[g, p, t]) / max_gen[g]
                        end

                end
        end

        return energy_price, reserve_up_price, reserve_down_price, capacity_factor, reserve_up_perc, reserve_down_perc;

    end
