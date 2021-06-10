"""
This struct contains the Inertia market parameters
    for price forcasts using CEM and endogeneous Economic Dispatch.
The Inertia market involves a minimum inertia requirment included in this struct.
"""

struct InertiaMarket

    inertia_req::Float64    # MW-s
    price_cap::Float64        # MW

    function InertiaMarket(req::Float64, price_cap::Float64)
        @assert req >= 0
        @assert price_cap >= 0
        new(req, price_cap)
    end

end
