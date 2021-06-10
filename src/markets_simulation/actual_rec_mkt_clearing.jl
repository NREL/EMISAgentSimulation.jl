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
    get_project_rec_output(project),
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


    return get_project_rec_output(project)
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
