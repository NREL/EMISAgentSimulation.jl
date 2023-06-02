#REC Market Clearing Module

"""
This function does nothing if the product is not REC
"""
function update_rec_supply_curve!(rec_supply_curve::Vector{Vector{Union{String, Float64}}},
    product::T,
    project::P) where {T <: Product, P <: Project{<: BuildPhase}}

    return rec_supply_curve
end

# This function pushes project bids in the REC supply curve
function update_rec_supply_curve!(rec_supply_curve::Vector{Vector{Union{String, Float64}}},
    product::REC,
    project::P) where P <: Project{<: BuildPhase}

    #Each Element of the REC supply Curve:
    #[1] - Project Name (Symbol)
    #[2] - Project REC output (Float64)
    #[3] - Project REC Bid (Float64)


    push!(rec_supply_curve, [get_name(project),
    get_project_expected_rec_output(project),
    get_project_rec_market_bid(project)])

    return rec_supply_curve
end

"""
This function does nothing if the product is not REC
"""
function find_clean_energy_production(product::T,
                                        project::P) where {T <: Product, P <: Project{<: BuildPhase}}

    return 0.0
end

# This function pushes project bids in the REC supply curve
function find_clean_energy_production(product::REC,
                                        project::P) where {P <: Project{<: BuildPhase}}


    return get_project_expected_rec_output(project)
end

"""
This function returns 0 if the project is not of type BatteryEMIS and product is not of type Energy
"""
function find_storage_energy_consumption(product::T,
                                        project::P) where {T <: Product, P <: Project{<: BuildPhase}}

    return 0.0
end

function find_storage_energy_consumption(product::T,
                                        project::BatteryEMIS{Existing}) where {T <: Product}

    return 0.0
end

function find_storage_energy_consumption(product::Energy,
                                        project::P) where {P <: Project{<: BuildPhase}}

    return 0.0
end

# This function pushes project bids in the REC supply curve
function find_storage_energy_consumption(product::Energy,
                                        project::BatteryEMIS{Existing})


    return get_expected_production(product)
end


function get_expected_rec_certificates(project::Project)
    certificates = 0.0
    for product in get_products(project)
        temp = get_expected_rec_certificates(product)
        if !(isnothing(temp))
            certificates = temp
        end
    end
    return certificates
end

"""
This function models the actual REC market clearing.
"""
function rec_market_clearing_non_binding(rec_requirement::Float64,
                                        ACP::Float64,
                                        supply_curve::Vector{Vector{Union{String, Float64}}},
                                        solver::JuMP.MOI.OptimizerWithAttributes)

    #Each Element of the REC supply Curve:
    #[1] - Project Name (Symbol)
    #[2] - Project REC output (Float64)
    #[3] - Project REC Bid (Float64)

    n_supply_seg = length(supply_curve);

    #----------REC Market Clearing Problem---------------------------------------------------------------------------------
    rec_mkt = JuMP.Model(solver);

    #Define the variables
    JuMP.@variables(rec_mkt, begin
    Q_supply[s=1:n_supply_seg] >= 0 # Quantity of cleared REC supply offers
    v_REC >= 0                      # RPS compliance violation
    end)

    #Functions----------------------------------------------------------------------------------------

    #Cost of procuring REC supply for each segment
    supply_cost(Q_seg, s) = supply_curve[s][3] * Q_seg;

    #Expressions---------------------------------------------------------------------------------------

    #Totoal Cost of procuring REC supply
    JuMP.@expression(rec_mkt, total_supply_cost, sum(supply_cost(Q_supply[s], s) for s = 1:n_supply_seg))

    #Constraints--------------------------------------------------------------------------------------

    #Cleared REC supply limit for each segment
    JuMP.@constraint(rec_mkt, [s=1:n_supply_seg], Q_supply[s] <= supply_curve[s][2])

    #Total cleared REC supply should be greater than or equal to the RPS compliance requirement
    JuMP.@constraint(rec_mkt, mkt_clear, sum(Q_supply[s] for s in 1:n_supply_seg) + v_REC >= rec_requirement)

    #Define Objective Function - Social Welfare Maxmization
    JuMP.@objective(rec_mkt, Min, total_supply_cost + v_REC * ACP)

    println("Actual REC Market Clearing:")
    JuMP.optimize!(rec_mkt)

    println(JuMP.termination_status(rec_mkt))
    println(JuMP.objective_value(rec_mkt))

    #REC Market Clearing Price is the shadow variable of the REC balance constraint
    rec_price = AxisArrays.AxisArray(reshape([JuMP.dual(mkt_clear)], 1,), [1])
    rec_accepted_bid = Dict(supply_curve[s][1] => value.(Q_supply[s]) for s in 1:n_supply_seg)
    println(rec_price)
    #------------------------------------------------------------------------------------------------
    return rec_price, rec_accepted_bid

end

"""
This function models the actual REC market clearing.
"""
function rec_market_clearing_binding(rec_requirement::Float64,
                                    ACP::Float64,
                                    supply_curve::Vector{Vector{Union{String, Float64}}},
                                    solver::JuMP.MOI.OptimizerWithAttributes)

    #Each Element of the REC supply Curve:
    #[1] - Project Name (Symbol)
    #[2] - Project REC output (Float64)
    #[3] - Project REC Bid (Float64)

    n_supply_seg = length(supply_curve);

    #----------REC Market Clearing Problem---------------------------------------------------------------------------------
    rec_mkt = JuMP.Model(solver);

    #Define the variables
    JuMP.@variables(rec_mkt, begin
    Q_supply[s=1:n_supply_seg] >= 0 # Quantity of cleared REC supply offers
    end)

    #Functions----------------------------------------------------------------------------------------

    #Cost of procuring REC supply for each segment
    supply_cost(Q_seg, s) = supply_curve[s][3] * Q_seg;

    #Expressions---------------------------------------------------------------------------------------

    #Totoal Cost of procuring REC supply
    JuMP.@expression(rec_mkt, total_supply_cost, sum(supply_cost(Q_supply[s], s) for s = 1:n_supply_seg))

    #Constraints--------------------------------------------------------------------------------------

    #Cleared REC supply limit for each segment
    JuMP.@constraint(rec_mkt, [s=1:n_supply_seg], Q_supply[s] <= supply_curve[s][2])

    #Total cleared REC supply should be greater than or equal to the RPS compliance requirement
    JuMP.@constraint(rec_mkt, mkt_clear, sum(Q_supply[s] for s in 1:n_supply_seg) >= rec_requirement)

    #Define Objective Function - Social Welfare Maxmization
    JuMP.@objective(rec_mkt, Min, total_supply_cost)

    println("Actual REC Market Clearing:")
    JuMP.optimize!(rec_mkt)

    println(JuMP.termination_status(rec_mkt))
    println(JuMP.objective_value(rec_mkt))

    #REC Market Clearing Price is the shadow variable of the REC balance constraint
    rec_price = AxisArrays.AxisArray(reshape([JuMP.dual(mkt_clear)], 1,), [1])
    rec_accepted_bid = Dict(supply_curve[s][1] => value.(Q_supply[s]) for s in 1:n_supply_seg)
    println(rec_price)
    #------------------------------------------------------------------------------------------------
    return rec_price, rec_accepted_bid

end

function update_rec_correction_factors!(active_projects::Vector{Project},
                                       realized_capacity_factors::Dict{String, Array{Float64, 2}},
                                       rt_resolution::Int64,
                                       iteration_year::Int64)

    zones = unique(get_zone.(get_tech.(active_projects)))
    techs = unique(get_type.(get_tech.(active_projects)))

    zone_tech_expected = Dict(["$(zone)_$(tech)" => 0.0 for zone in zones, tech in techs])
    zone_tech_actual = Dict(["$(zone)_$(tech)" => 0.0 for zone in zones, tech in techs])

    for project in active_projects
        project_name = get_name(project)
        zone = get_zone(get_tech(project))
        tech = get_type(get_tech(project))

        if project_name in keys(realized_capacity_factors)
            #println(project_name)
            expected_rec_production = get_expected_rec_certificates(project)
            #println(expected_rec_production)
            actual_rec_production = sum(realized_capacity_factors[project_name]) * get_maxcap(project) * rt_resolution / 60
            #println(actual_rec_production)
            zone_tech_expected["$(zone)_$(tech)"] += expected_rec_production
            zone_tech_actual["$(zone)_$(tech)"] += actual_rec_production
        end

    end

    for project in active_projects
        zone = get_zone(get_tech(project))
        tech = get_type(get_tech(project))
        expected_rec_production = zone_tech_expected["$(zone)_$(tech)"]
        actual_rec_production = zone_tech_actual["$(zone)_$(tech)"]
        previous_correction_factor = get_rec_correction_factor(project, iteration_year)
        if !(iszero(expected_rec_production))
            correction_factor = min((actual_rec_production / expected_rec_production), 1.0)
        else
            correction_factor = 1.0
        end

        for product in get_products(project)
            set_rec_correction_factor!(product, iteration_year + 1, (correction_factor + previous_correction_factor) / 2)
            # set_rec_correction_factor!(product, iteration_year + 1, 1.0)
        end

    end


    return
end
