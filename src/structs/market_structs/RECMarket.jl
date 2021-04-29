"""
This struct contains the REC market parameters
    for price forcasts using CEM and endogeneous Economic Dispatch.
The REC market assumes a percentage RPS requirment of annual energy consumption
which is parameterized within this struct.
"""

struct RECMarket

    rec_req::Float64    # $/MW/investment period
    price_cap::Float64        # MW
    binding::Bool

    function RECMarket(req::Number, price_cap::Number, binding::Bool)
        @assert 0 <= req <= 1
        @assert price_cap >= 0
        new(req, price_cap, binding)
    end

end
