"""
This function does nothing if the product is not of Capacity type.
"""
function update_capacity_revenues!(product::T,
                                scenario_name::String,
                                price_years::NamedTuple{(:start_year, :end_year),
                                    Tuple{Int64, Int64}},
                                capacity_revenues::Vector{Float64},
                                capacity_prices::Dict{String, AxisArrays.AxisArray{Float64, 1}},
                                max_cap::Float64) where T <: Product
    return capacity_revenues
end

"""
This function calculates expected capacity revenues for each project.
Returns capacity revenues.
"""
function update_capacity_revenues!(product::Capacity,
                                scenario_name::String,
                                price_years::NamedTuple{(:start_year, :end_year),
                                    Tuple{Int64, Int64}},
                                capacity_revenues::Vector{Float64},
                                capacity_prices::Dict{String, AxisArrays.AxisArray{Float64, 1}},
                                max_cap::Float64)
    capacity_revenues += capacity_prices[scenario_name][price_years[:start_year]:price_years[:end_year]] .*
                         get_accepted_perc(product)[scenario_name][price_years[:start_year]:price_years[:end_year]] *
                         max_cap *
                         get_derating(product)
    return capacity_revenues
end

"""
This function updates the expected profit array of the project
if the product is a capacity product.
Returns nothing.
"""
function update_profit!(project::P,
                        scenario_name::String,
                        product::Capacity,
                        operating_profit::AxisArrays.AxisArray{Float64, 2},
                        capacity_revenues::Vector{Float64},
                        capacity_forward_years::Int64,
                        rec_revenues::Vector{Float64},
                        price_years::NamedTuple{(:start_year, :end_year),
                                     Tuple{Int64, Int64}},
                        update_years::NamedTuple{(:start_year, :end_year),
                                      Tuple{Int64, Int64}},
                        iteration_year::Int64) where P <: Project{<: BuildPhase}

    capacity_revenue_start_year = max(1, capacity_forward_years - price_years[:start_year] + 1)

        for i in capacity_revenue_start_year:length(capacity_revenues)
                set_scenario_profit!(get_finance_data(project),
                            scenario_name,
                            get_name(product),
                            iteration_year,
                            update_years[:start_year] + i - 1,
                            capacity_revenues[i])
        end

    return
end
