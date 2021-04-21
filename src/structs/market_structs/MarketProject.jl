"""
This struct contains all the data for creating a market project
    for price forcasts using CEM and endogeneous Economic Dispatch.
"""

mutable struct MarketProject

    name::String                         # name of project
    project_type::String                 # is the project of storage type
    tech_type::String                    # technology type
    fixed_cost::Float64                   # $/unit/investment period
    queue_cost::Vector{Float64}           # queue cost
    marginal_energy_cost::Float64         # $/MW/hour
    marginal_reserve_cost::Dict{String, Float64}
    emission_intensity::Float64
    expansion_cost::Vector{Float64}        # $/unit
    discount_rate::Float64                # discount rate
    min_gen::Float64                      # MW/unit
    max_gen::Float64                      # MW/unit
    min_input::Float64                    # MW/unit
    max_input::Float64                    # MW/unit
    efficiency_in::Float64               # input efficiency
    efficiency_out::Float64              # output efficiency
    min_storage::Float64                  # minimum storage capacity
    max_storage::Float64                  # maximum storage capacity
    init_storage::Float64                 # initial storage
    availability::Array{Float64, 2}      # Hourly availability factor
    derating_factor::Float64              # Derating factor
    ramp_limits::Union{Nothing, NamedTuple{(:up, :down), Tuple{Float64, Float64}}}    # MW/unit/hour
    max_reserve_limits::Dict{String, Float64}
    existing_units::Float64               # existing units
    units_in_queue::Vector{Float64}        # units in queue
    build_lead_time::Int64                 # build lead time
    remaining_build_time::Int64            # remaining build time
    max_new_options::Int                   # maximum units
    base_cost_units::Int64                 # number of units that can be built without capital cost multiplier
    capex_years::Int64                    # capital cost recovery years
    life_time::Int64                      # life_time
    remaining_life::Int64                # Remaining life_time
    capacity_eligible::Bool              # eligible for capacity market participation
    rec_eligible::Bool                   # eligible for rps compliance
    zone::String                         # project zone
    ownedby::Vector{String}              # onned by which investors

    function MarketProject(
        name::String,
        project_type::String,
        tech_type::String,
        fixed_cost::Number,
        queue_cost::Vector{<: Number},
        marginal_energy_cost::Number,
        marginal_reserve_cost::Dict{String, Float64},
        emission_intensity::Float64,
        expansion_cost::Vector{Float64},
        discount_rate::Number,
        min_gen::Number,
        max_gen::Number,
        min_input::Number,
        max_input::Number,
        efficiency_in::Number,
        efficiency_out::Number,
        min_storage::Number,
        max_storage::Number,
        init_storage::Number,
        availability::Array{<: Number, 2},
        derating_factor::Number,
        ramp_limits::Union{Nothing, NamedTuple{(:up, :down), Tuple{Float64, Float64}}},
        max_reserve_limits::Dict{String, Float64},
        existing_units::Number,
        units_in_queue::Vector{<: Number},
        build_lead_time::Number,
        remaining_build_time::Number,
        max_new_options::Number,
        base_cost_units::Number,
        capex_years::Number,
        life_time::Number,
        remaining_life::Int64,
        capacity_eligible::Bool,
        rec_eligible::Bool,
        zone::String,
        ownedby::Vector{String})

        @assert fixed_cost >= 0
        @assert marginal_energy_cost >= 0
        @assert all(values(marginal_reserve_cost) .>= 0)
        @assert emission_intensity >= 0
        @assert all(expansion_cost .>= 0)
        @assert all(queue_cost .>= 0)
        @assert discount_rate >= 0
        @assert min_gen >= 0
        @assert max_gen >= 0
        @assert min_input >= 0
        @assert max_input >= 0
        @assert all(values(max_reserve_limits) .>= 0)
        @assert efficiency_in >= 0
        @assert efficiency_out >= 0
        @assert min_storage >= 0
        @assert max_storage >= 0
        @assert init_storage >= 0
        @assert all(availability .>= 0)
        @assert derating_factor >= 0
        @assert existing_units >= 0
        @assert all(units_in_queue .>= 0)
        @assert build_lead_time >= 0
        @assert remaining_build_time >= 0
        @assert base_cost_units >= 0
        @assert max_new_options >= base_cost_units
        @assert capex_years >- 0
        @assert life_time >= 0
        @assert remaining_life >= 0

        new(name,
            project_type,
            tech_type,
            fixed_cost,
            queue_cost,
            marginal_energy_cost,
            marginal_reserve_cost,
            emission_intensity,
            expansion_cost,
            discount_rate,
            min_gen,
            max_gen,
            min_input,
            max_input,
            efficiency_in,
            efficiency_out,
            min_storage,
            max_storage,
            init_storage,
            availability,
            derating_factor,
            ramp_limits,
            max_reserve_limits,
            existing_units,
            units_in_queue,
            build_lead_time,
            remaining_build_time,
            max_new_options,
            base_cost_units,
            capex_years,
            life_time,
            remaining_life,
            rec_eligible,
            capacity_eligible,
            zone,
            ownedby
            )

    end

end
