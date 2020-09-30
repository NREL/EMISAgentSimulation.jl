"""
This struct contains the collection of all markets modeled
for price forcasts using CEM and endogeneous Economic Dispatch.
"""

struct MarketCollection{Z, T}

    capacity_market::CapacityMarket
    energy_market::EnergyMarket
    reserveup_market::ReserveUpMarket
    reservedown_market::ReserveDownMarket
    rec_market::RECMarket

    function MarketCollection(
                              c::CapacityMarket,
                              e::EnergyMarket{Z, T},
                              r::ReserveUpMarket,
                              l::ReserveDownMarket{Z, T},
                              rec::RECMarket
    ) where {Z, T}
        new{Z, T}(c, e, r, l, rec)
    end
end
