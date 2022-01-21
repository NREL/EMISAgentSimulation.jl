
#Capacity Market Clearing Module

"""
This function does nothing if the product is not Capacity
"""
function update_capacity_supply_curve!(capacity_supply_curve::Vector{Vector{Union{String, Float64}}},
                                       product::T,
                                       project::P) where {T <: Product, P <: Project{<: BuildPhase}}

   return capacity_supply_curve
end

"""
This function pushes project bids in the capacity supply curve
"""
function update_capacity_supply_curve!(capacity_supply_curve::Vector{Vector{Union{String, Float64}}},
                                       product::Capacity,
                                       project::P) where P <: Project{<: BuildPhase}

   # Each Element of the Supply Curve:
   #   [1] - Project Name
   #   [2] - Project Size
   #   [3] - Project Capacity Bid
   #   [4] - Project Derating Factor

   push!(capacity_supply_curve, [get_name(project),
                                 get_maxcap(project),
                                 get_project_capacity_market_bid(project),
                                 get_project_derating(project)])
   return capacity_supply_curve
end

"""
This function creates the capacity market demand curve for actual capacity market clearing.
"""
function create_capacity_demand_curve(input_file::String,
                                      system_peak_load::Float64,
                                      irm_scalar::Float64,
                                      delta_irm::Float64,
                                      capacity_mkt_bool::Bool)

   # Gather parameter data
   capacity_demand_params = read_data(input_file)[1, :]
   eford = capacity_demand_params["EFORd"] # Equivalent demand forced outage rate
   base_irm = capacity_demand_params["IRM"] # Installed Reserve Margin
   adjusted_irm = (base_irm + delta_irm) * irm_scalar
   fpr = (1 + adjusted_irm) * (1 - eford) #Forecast Pool Req
   rel_req = system_peak_load * fpr * capacity_mkt_bool # Reliability Requirement
   irm_perc_points = parse.(Float64, split(capacity_demand_params["IRM perc points"], ";")) #Installed Reserve Margin percentage points
   net_CONE = capacity_demand_params["Net CONE per day"] * 365 * capacity_mkt_bool # Net CONE
   net_CONE_perc_points = parse.(Float64, split(capacity_demand_params["Net CONE perc points"], ";")) #Net CONE percentage points

   @assert length(net_CONE_perc_points) == length(irm_perc_points)

   gross_CONE_value = capacity_demand_params["Gross CONE per day"] * 365 * capacity_mkt_bool  # Gross CONE
   gross_CONE_points = zeros(length(net_CONE_perc_points))
   gross_CONE_points[1] = gross_CONE_value

   max_clear = capacity_demand_params["Max Clear"]

   # Construct demand curve
   demand_curve_break_points = zeros(length(net_CONE_perc_points) + 2)
   demand_curve_price_points = zeros(length(net_CONE_perc_points) + 2)

   demand_curve_break_points[1] = 0.0
   demand_curve_price_points[1] = max(net_CONE * net_CONE_perc_points[1], gross_CONE_points[1])

   demand_curve_break_points[end] = rel_req * max_clear
   demand_curve_price_points[end] = 0.0

   # Demand curve points based on PJM capacity market
   for i in 1:(length(net_CONE_perc_points))
      demand_curve_break_points[i + 1] = rel_req * ( 1 + adjusted_irm + irm_perc_points[i]) / ( 1 + adjusted_irm)
      demand_curve_price_points[i + 1] = max(net_CONE * net_CONE_perc_points[i], gross_CONE_points[i])
   end

   capacity_demand_curve = CapacityMarket(demand_curve_break_points, demand_curve_price_points)

   return capacity_demand_curve
end

"""
This function models actual capacity market clearing.
"""
function capacity_market_clearing(demand_curve::CapacityMarket,
                                  supply_curve::Vector{Vector{Union{String, Float64}}},
                                  solver::JuMP.MOI.OptimizerWithAttributes)


   # Each Element of the Supply Curve:
   #   [1] - Project Name
   #   [2] - Project Size
   #   [3] - Project Capacity Bid
   #   [4] - Project Derating Factor


   demand_segmentsize, demand_segmentgrad, demand_pricepoints = make_capacity_demand(demand_curve)

   #Number of segments of capacity supply and demand curves
   n_demand_seg = length(demand_segmentsize);
   n_supply_seg = length(supply_curve);

   #----------Capacity Market Clearing Problem---------------------------------------------------------------------------------
   cap_mkt = JuMP.Model(solver);

   #Define the variables
   JuMP.@variables(cap_mkt, begin
          Q_supply[s=1:n_supply_seg] >= 0 # Quantity of cleared capacity supply offers
          Q_demand[d=1:n_demand_seg] >= 0 # Quantity of cleared capacity demand from the demand curve
       end)

   #Functions----------------------------------------------------------------------------------------

   #Cost of procuring capacity supply for each segment
   supply_cost(Q_seg, s) = supply_curve[s][3] * Q_seg;

   #Welfare from meeting the capacity resource requirement for each segment
   demand_welfare(Q_seg, d) = demand_pricepoints[d] * Q_seg + 0.5 * demand_segmentgrad[d] * (Q_seg^2);

   #Expressions---------------------------------------------------------------------------------------

   #Totoal Cost of procuring capacity supply
   JuMP.@expression(cap_mkt, total_supply_cost, sum(supply_cost(Q_supply[s], s) for s = 1:n_supply_seg))

   #Total Welfare from meeting the capacity resource requirement
   JuMP.@expression(cap_mkt, total_dem_welfare, sum(demand_welfare(Q_demand[d], d) for d = 1:n_demand_seg))

   #Constraints--------------------------------------------------------------------------------------

   #Cleared capacity supply limit for each segment
   JuMP.@constraint(cap_mkt, [s=1:n_supply_seg], Q_supply[s] <= supply_curve[s][2])

   #Cleared demand limit for each segment
   JuMP.@constraint(cap_mkt, [d=1:n_demand_seg], Q_demand[d] <= demand_segmentsize[d])

   #Total cleared capacity supply should meet the total cleared demand
   JuMP.@constraint(cap_mkt, mkt_clear, sum(Q_supply[s] * supply_curve[s][4] for s in 1:n_supply_seg) == sum(Q_demand))

   #Define Objective Function - Social Welfare Maxmization
   JuMP.@objective(cap_mkt, Max, total_dem_welfare - total_supply_cost)

   println("Actual Capacity Market Clearing:")
   JuMP.optimize!(cap_mkt)
   println(JuMP.termination_status(cap_mkt))
   println(JuMP.objective_value(cap_mkt))

   #Capacity Market Clearing Price is the shadow variable of the capacity balance constraint
   cap_price = AxisArrays.AxisArray(reshape([JuMP.dual(mkt_clear)], 1,), [1])
   cap_accepted_bid = Dict(supply_curve[s][1] => value.(Q_supply[s]) / supply_curve[s][2]  for s in 1:n_supply_seg)
   println(cap_price)
#------------------------------------------------------------------------------------------------
   return cap_price, cap_accepted_bid

end
