"""
This struct contains data for the REC market product.
    name: REC.
    certificates: Expected output (number of certificates) for REC market bids
    rec_bid: REC market bidding price
"""
mutable struct REC <: Product
    name::Symbol
    certificates::Float64
    rec_bid::Float64
end

# Certificates and REC bids only returned when product is of type REC
get_rec_certificates(prod::Product) = nothing
get_rec_certificates(prod::REC) = prod.certificates
get_rec_bid(prod::Product) = nothing
get_rec_bid(prod::REC) = prod.rec_bid

# REC bids and certificates only set when product is of type REC
function set_rec_certificates!(prod::T, certificates) where T <: Product
    return
end

function set_rec_certificates!(prod::REC, certificates)
    prod.certificates = certificates
    return
end

function set_rec_bid!(prod::T, rec_bid) where T <: Product
    return
end

function set_rec_bid!(prod::REC, rec_bid)
    prod.rec_bid = rec_bid
    return
end
