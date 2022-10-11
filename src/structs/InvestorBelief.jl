""" This struct defines the necessary parameters for build a Kalman filter to predict
and update investor expectations/beliefs about the future """

mutable struct InvestorBelief
    process_covariance::AxisArrays.AxisArray{Float64, 1}
    measurement_covariance::AxisArrays.AxisArray{Float64, 1}
    state_transition_matrix::Matrix{Float64}
    state_measurement_matrix::Matrix{Float64}

    function InvestorBelief(process_covariance::AxisArrays.AxisArray{Float64, 1},
        measurement_covariance::AxisArrays.AxisArray{Float64, 1},
        state_transition_matrix::Matrix{Float64},
        state_measurement_matrix::Matrix{Float64}
    )
        n = size(process_covariance, 1)
        (m, )  = size(state_measurement_matrix)
        @assert size(state_transition_matrix) == (n, n)
        @assert size(state_measurement_matrix) == (m, n)
        @assert size(measurement_covariance) == (n,)

        new(process_covariance,
            measurement_covariance,
            state_transition_matrix,
            state_measurement_matrix,
        )
    end
end

function InvestorBelief(process_covariance::AxisArrays.AxisArray{Float64, 1},
    measurement_covariance::AxisArrays.AxisArray{Float64, 1}
)


    n = size(process_covariance, 1)
    state_transition_matrix = Matrix(1.0*LinearAlgebra.I, n, n)
    state_measurement_matrix = Matrix(1.0*LinearAlgebra.I, n, n)

    return InvestorBelief(process_covariance,
        measurement_covariance,
        state_transition_matrix,
        state_measurement_matrix,
    )
end


get_process_covariance(v::InvestorBelief) = v.process_covariance
get_measurement_covariance(v::InvestorBelief) = v.measurement_covariance
get_state_transition_matrix(v::InvestorBelief) = v.state_transition_matrix
get_state_measurement_matrix(v::InvestorBelief) = v.state_measurement_matrix
get_n(v::InvestorBelief) = size(v.process_covariance, 1)
get_m(v::InvestorBelief) = size(v.state_measurement_matrix, 1)


mutable struct KalmanFilter
    investor_data::InvestorBelief
    state_estimate::AxisArrays.AxisArray{Float64, 1}
    error_covariance_estimate::AxisArrays.AxisArray{Float64, 2}
end

get_investor_data(v::KalmanFilter) = v.investor_data
get_state_estimate(v::KalmanFilter) = v.state_estimate
get_error_covariance_estimate(v::KalmanFilter) = v.error_covariance_estimate

set_state_estimate(v::KalmanFilter, value) = v.state_estimate = value
set_error_covariance_estimate(v::KalmanFilter, value) = v.error_covariance_estimate = value

function predict(kf::KalmanFilter)

    investor_data =  kf.investor_data
    predicted_state = investor_data.state_transition_matrix * kf.state_estimate
    predicted_error_covariance = ((investor_data.state_transition_matrix * kf.error_covariance_estimate
                                        * investor_data.state_transition_matrix' ) + LinearAlgebra.Diagonal(investor_data.process_covariance))

    return (predicted_state, predicted_error_covariance)
end

function update_belief!(kf::KalmanFilter, measurement::AxisArrays.AxisArray{Float64, 1})
    investor_data =  kf.investor_data
    n = get_n(investor_data)
    state_estimate_priori, error_covariance_priori = predict(kf)
    kalman_gain = ( (error_covariance_priori * investor_data.state_measurement_matrix) /
                    (investor_data.state_measurement_matrix * error_covariance_priori * (investor_data.state_measurement_matrix)' +  LinearAlgebra.Diagonal(investor_data.measurement_covariance)))

    state_estimate_posteriori = state_estimate_priori + kalman_gain * (measurement - investor_data.state_measurement_matrix * state_estimate_priori)

    error_covariance_posteriori = (LinearAlgebra.I(n) - kalman_gain * investor_data.state_measurement_matrix) * error_covariance_priori

    param_names = AxisArrays.axisvalues(investor_data.measurement_covariance)[1]

    state_estimate_posteriori = AxisArrays.AxisArray(state_estimate_posteriori, param_names)
    error_covariance_posteriori = AxisArrays.AxisArray(error_covariance_posteriori, param_names, param_names)

    set_state_estimate(kf, state_estimate_posteriori)
    set_error_covariance_estimate(kf, error_covariance_posteriori)
    return
end

function update_belief!(kf::Nothing, measurement::AxisArrays.AxisArray{Float64, 1})
    return
end
