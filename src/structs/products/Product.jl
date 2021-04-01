abstract type Product end
abstract type OperatingProduct <: Product end
abstract type ReserveDirection end
abstract type ReserveUpEMIS <: ReserveDirection end
abstract type ReserveDownEMIS <: ReserveDirection end
abstract type ReserveEMIS{T<:ReserveDirection} <: OperatingProduct end

get_name(prod::Product) = prod.name
get_marginal_cost(prod::OperatingProduct) = prod.marginal_cost
get_zonal(prod::Product) = prod.zonal
