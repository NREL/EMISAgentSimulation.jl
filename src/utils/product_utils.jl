
"""
This funtion does nothing if the product is not an operating market product.
"""
function get_maxperc(product::T, scenario_name::String, num_years::Int64, num_hours::Int64) where T <: Product
    return
end

"""
This funtion returns the array of capacity factors if product is energy.
"""
function get_maxperc(product::Energy, scenario_name::String, num_years::Int64, num_hours::Int64)
    return get_capacity_factors(product)[scenario_name]
end


"""
This funtion returns the array of maximum percentage limit if product is operating reserve.
"""
function get_maxperc(product::OperatingReserve, scenario_name::String, num_years::Int64, num_hours::Int64)
    return ones(num_years, num_hours) * get_max_limit(product)
end

"""
This function does nothing if the product is not an operating market product.
"""
function push_operating!(array::Array{OperatingProduct}, product::T) where T <: Product
    return
end

"""
This function pushes the product to the array if product is an operating market product.
"""
function push_operating!(array::Array{OperatingProduct}, product::OperatingProduct)
    push!(array, product)
    return
end


function push_energy_product!(array::Array{Energy}, product::OperatingProduct)
    return
end

function push_energy_product!(array::Array{Energy}, product::Energy)
    push!(array, product)
    return
end

function find_energy_product(products::Array{OperatingProduct})
    array = Energy[]

    for product in products
        push_energy_product!(array, product)
    end

    return array
end

"""
This function returns an array of only operating market products from an array of all products.
"""
function find_operating_products(products::Array{Product})
    operating_products = OperatingProduct[]

    for i = 1:length(products)
        push_operating!(operating_products, products[i])
    end

return  operating_products
end
