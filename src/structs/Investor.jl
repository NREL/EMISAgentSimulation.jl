"""
This struct contains all the data for an investor.
    name: Investor's name.
    data_dir: Directoty where investors' data is stored.
    projects: Vector of all projects owned by the investor.
    markets: Vector of markets in which the investor is participating.
    carbon_tax: Vector of annual carbon tax.
    market_prices: Expected market prices.
    rep_hour_weight: Representative hour weight for operating markets.
    forecast: Forecast parameters.
    cap_cost_multiplier: Capital cost multiplier.
    max_annual_projects: Maximum number of projects.
    risk_preference: Whether investor is Risk Neutral or Risk Averse.
    retirement_lookback: Number of past years' data taken into account for making retirement decision
"""
mutable struct Investor
    name::String
    data_dir::String
    projects::Vector{<: Project{<: BuildPhase}}
    markets::Vector{Symbol}
    carbon_tax::Vector{Float64}
    market_prices::MarketPrices
    rep_hour_weight::Vector{Float64}
    forecast::F where F <: Forecast
    cap_cost_multiplier::Float64
    preference_multiplier_range::NamedTuple{(:min, :max),Tuple{Float64,Float64}}
    portfolio_preference_multipliers::Dict{Tuple{String, String}, Vector{Float64}}
    max_annual_projects::Int64
    risk_preference::T where T <: RiskPreference
    retirement_lookback::Int64
end

get_name(investor::Investor) = investor.name
get_data_dir(investor::Investor) = investor.data_dir
get_projects(investor::Investor) = investor.projects
get_markets(investor::Investor) = investor.markets
get_carbon_tax(investor::Investor) = investor.carbon_tax
get_market_prices(investor::Investor) = investor.market_prices
get_rep_hour_weight(investor::Investor) = investor.rep_hour_weight
get_forecast(investor::Investor) = investor.forecast
get_cap_cost_multiplier(investor::Investor) = investor.cap_cost_multiplier
get_max_annual_projects(investor::Investor) = investor.max_annual_projects
get_risk_preference(investor::Investor) = investor.risk_preference
get_retirement_lookback(investor::Investor) = investor.retirement_lookback
get_preference_multiplier_range(investor::Investor) = investor.preference_multiplier_range
get_portfolio_preference_multipliers(investor::Investor) = investor.portfolio_preference_multipliers


# Get projects based on their buildphases
function get_existing(investor::Investor)
    leaftypes_existing = leaftypes(Project{Existing})
    projects = get_projects(investor)
    existing_projects = filter(project -> in(typeof(project), leaftypes_existing), projects)

    return existing_projects
end

function get_planned(investor::Investor)
    leaftypes_planned = leaftypes(Project{Planned})
    projects = get_projects(investor)
    planned_projects = filter(project -> in(typeof(project), leaftypes_planned), projects)

    return planned_projects
end

function get_queue(investor::Investor)
    leaftypes_queue = leaftypes(Project{Queue})
    projects = get_projects(investor)
    queue_projects = filter(project -> in(typeof(project), leaftypes_queue), projects)

    return queue_projects
end

function get_options(investor::Investor)
    leaftypes_option = leaftypes(Project{Option})
    projects = get_projects(investor)
    option_projects = filter(project -> in(typeof(project), leaftypes_option), projects)

    return option_projects
end

function get_retired(investor::Investor)
    leaftypes_retired = leaftypes(Project{Retired})
    projects = get_projects(investor)
    retired_projects = filter(project -> in(typeof(project), leaftypes_retired), projects)

    return retired_projects
end

function find_active_invested_projects(projects::Vector{<: Project})
    leaftypes_option = leaftypes(Project{Option})
    leaftypes_retired= leaftypes(Project{Retired})
    active_invested_projects = filter(project -> (!in(typeof(project), leaftypes_option) && !in(typeof(project), leaftypes_retired) ), projects)

    return active_invested_projects
end

function find_option_projects(projects::Vector{<: Project})
    leaftypes_option = leaftypes(Project{Option})
    option_projects = filter(project -> in(typeof(project), leaftypes_option), projects)

    return option_projects
end

function set_market_prices!(investor::Investor,
                           market_prices::MarketPrices)
    investor.market_prices = market_prices
end

function get_preference_multiplier(investor::Investor, tech::String, zone::String)
    return investor.portfolio_preference_multipliers[(tech, zone)]
end

function set_preference_multiplier!(investor::Investor, tech::String, zone::String, iteration_year::Int64, value::Float64)
    investor.portfolio_preference_multipliers[(tech, zone)][iteration_year] = value
end
