# Investor's Risk Preference
abstract type RiskPreference end

"""
Risk Neutrality
"""
struct RiskNeutral <: RiskPreference end

"""
This struct contains the parameters for modeling Risk Aversion,
    where a project utility is calculated as follows:
    U = constant - multiplier * e ^ (- risk_coefficient ^ NPV)
"""
struct RiskAverse <: RiskPreference
    constant::Float64
    multiplier::Float64
    risk_coefficient::Float64
end

get_risk_coefficient(preference::RiskAverse) = preference.risk_coefficient
get_constant(preference::RiskAverse) = preference.constant
get_multiplier(preference::RiskAverse) = preference.multiplier
