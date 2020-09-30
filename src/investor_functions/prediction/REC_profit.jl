"""
This function does nothing if the product is not of Capacity type.
"""
function update_rec_revenues!(product::T,
                                price_years::NamedTuple{(:start_year, :end_year),
                                    Tuple{Int64, Int64}},
                                rec_revenues::Vector{Float64},
                                rec_prices::AxisArrays.AxisArray{Float64, 1},
                                energy_production::AxisArrays.AxisArray{Float64, 1}) where T <: Product
    return rec_revenues
end

"""
This function calculates expected rec revenues for each project.
Returns REC revenues
"""
function update_rec_revenues!(product::REC,
                                price_years::NamedTuple{(:start_year, :end_year),
                                    Tuple{Int64, Int64}},
                                rec_revenues::Vector{Float64},
                                rec_prices::AxisArrays.AxisArray{Float64, 1},
                                energy_production::AxisArrays.AxisArray{Float64, 1})
    rec_revenues += rec_prices[price_years[:start_year]:price_years[:end_year]] .*
    energy_production[1:end]
    return rec_revenues
end

"""
This function updates the expected profit array of the project
if the product is a REC product.
Returns nothing.
"""
function update_profit!(project::P,
                        scenario_name::String,
                        product::REC,
                        operating_profit::AxisArrays.AxisArray{Float64, 2},
                        capacity_revenues::Vector{Float64},
                        capacity_forward_years::Int64,
                        rec_revenues::Vector{Float64},
                        price_years::NamedTuple{(:start_year, :end_year),
                                     Tuple{Int64, Int64}},
                        update_years::NamedTuple{(:start_year, :end_year),
                                      Tuple{Int64, Int64}},
                        iteration_year::Int64) where P <: Project{<: BuildPhase}

    for i in 1:length(rec_revenues)
            set_scenario_profit!(get_finance_data(project),
                        scenario_name,
                        get_name(product),
                        iteration_year,
                        update_years[:start_year] + i - 1,
                        rec_revenues[i])
    end

    return
end
