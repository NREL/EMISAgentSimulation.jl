"""
This struct contains the collection of all markets modeled
for price forcasts using CEM and endogeneous Economic Dispatch.
"""

struct MarketCollection{Z, T}

    capacity::CapacityMarket
    energy::EnergyMarket
    reserveup::Union{Nothing, Dict{String, ReserveUpMarket{T}}}
    reservedown::Union{Nothing, Dict{String, ReserveDownMarket{T}}}
    ordc::Union{Nothing, Dict{String, ReserveORDCMarket{T}}}
    rec::RECMarket
    inertia::InertiaMarket

    function MarketCollection(
                              c::CapacityMarket,
                              e::EnergyMarket{Z, T},
                              ru::Union{Nothing, Dict{String, ReserveUpMarket{T}}},
                              rd::Union{Nothing, Dict{String, ReserveDownMarket{T}}},
                              ordc::Union{Nothing, Dict{String, ReserveORDCMarket{T}}},
                              rec::RECMarket,
                              inertia::InertiaMarket

    ) where {Z, T}
        new{Z, T}(c, e, ru, rd, ordc, rec, inertia)
    end
end
