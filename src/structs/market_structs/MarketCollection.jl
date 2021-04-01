"""
This struct contains the collection of all markets modeled
for price forcasts using CEM and endogeneous Economic Dispatch.
"""

struct MarketCollection{Z, T}

    capacity_market::CapacityMarket
    energy_market::EnergyMarket
    reserveup_market::ReserveUpMarket
    reservedown_market::ReserveDownMarket
    synchronous_reserve_market::Union{Nothing, ReserveUpMarket, ReserveORDCMarket{T}}
    rec_market::RECMarket

    function MarketCollection(
                              c::CapacityMarket,
                              e::EnergyMarket{Z, T},
                              r::ReserveUpMarket,
                              l::ReserveDownMarket{Z, T},
                              s::Union{ReserveUpMarket, ReserveORDCMarket{T}},
                              rec::RECMarket
    ) where {Z, T}
        new{Z, T}(c, e, r, l, s, rec)
    end
end
