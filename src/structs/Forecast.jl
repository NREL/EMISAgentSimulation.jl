"""
Forecast types include Perfect and Imperfect forecasts
"""
abstract type Forecast end

"""
Perfect forecasts.
"""
struct Perfect <: Forecast
    scenario_data::Vector{Scenario}
end

"""
Imperfect forecasts with Kalman Filter based updates.
"""
struct Imperfect <: Forecast
    kalman_filter::Union{Nothing, KalmanFilter}
    scenario_data::Vector{Scenario}
end

get_scenario_data(forecast::F) where F<: Forecast = forecast.scenario_data
get_kalman_filter(forecast::Imperfect) = forecast.kalman_filter
