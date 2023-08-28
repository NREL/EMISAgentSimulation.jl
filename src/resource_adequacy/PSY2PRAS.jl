#######################################################
# Surya
# NREL
# January 2021
# SIIP --> PRAS Linkage Module
#######################################################
# Loading the required packages
#######################################################
using PowerSystems
using DataFrames
using CSV
using Statistics
using TimeZones
using PRAS
using Dates
const PSY = PowerSystems
#######################################################
# Outage Information CSV
#######################################################
const OUTAGE_INFO_FILE =
    joinpath(@__DIR__, "descriptors", "outage-rates-ERCOT-modified.csv")

df_outage = DataFrames.DataFrame(CSV.File(OUTAGE_INFO_FILE));
#######################################################
# Structs to parse and store the outage information
#######################################################
struct outage_data
    prime_mover::String
    thermal_fuel::String
    capacity::Int64
    FOR::Float64
    MTTR::Int64

    outage_data(prime_mover  = "PrimeMovers.Default", thermal_fuel ="ThermalFuels.Default", capacity = 100, FOR=0.5,MTTR = 50) =new(prime_mover,thermal_fuel,capacity,FOR,MTTR)
end

outage_values =[]
for row in eachrow(df_outage)
    push!(outage_values, outage_data(row.PrimeMovers,row.ThermalFuels,row.NameplateLimit_MW,(row.FOR/100),row.MTTR))
end
##############################################
# Converting FOR and MTTR to λ and μ
##############################################
function outage_to_rate(outage_data::Tuple{Float64, Int64})
    for_gen = outage_data[1]
    mttr = outage_data[2]
    if (mttr != 0)
        μ = 1 / mttr
    else
        μ = 1.0
    end
    λ = (μ * for_gen) / (1 - for_gen)

    return (λ = λ, μ = μ)
end

#######################################################
# Structs to parse and store the outage information
#######################################################
struct outage_data
    prime_mover::String
    thermal_fuel::String
    capacity::Int64
    FOR::Float64
    MTTR::Int64

    outage_data(prime_mover  = "PrimeMovers.Default", thermal_fuel ="ThermalFuels.Default", capacity = 100, FOR=0.5,MTTR = 50) =new(prime_mover,thermal_fuel,capacity,FOR,MTTR)
end 

outage_values =[]
for row in eachrow(df_outage)
    push!(outage_values, outage_data(row.PrimeMovers,row.ThermalFuels,row.NameplateLimit_MW,(row.FOR/100),row.MTTR))
end
#######################################################
# Aux Functions
# Function to get Line Rating
#######################################################
function line_rating(line::Union{PSY.Line,PSY.MonitoredLine})
    rate = PSY.get_rate(line);
    return(forward_capacity = rate , backward_capacity = rate)
end

function line_rating(line::PSY.HVDCLine)
    forward_capacity = getfield(PSY.get_active_power_limits_from(line), :max)
    backward_capacity = getfield(PSY.get_active_power_limits_to(line), :max)
    return(forward_capacity = forward_capacity, backward_capacity = backward_capacity)
end
#######################################################
# Function to get available components in AggregationTopology
#######################################################
function get_available_components_in_aggregation_topology(type::Type{<:PowerSystems.StaticInjection}, sys::PSY.System, region::PSY.AggregationTopology)
    avail_comps =  [comp for comp in PSY.get_components_in_aggregation_topology(type, sys, region) if (PSY.get_available(comp))]

    return avail_comps
end
#######################################################
# Functions to get generator category
#######################################################
function get_generator_category(gen::GEN) where {GEN <: PSY.RenewableGen}
    return string(PSY.get_prime_mover(gen))
end

function get_generator_category(gen::GEN) where {GEN <: PSY.ThermalGen}
    return string(PSY.get_fuel(gen))
end

function get_generator_category(gen::GEN) where {GEN <: PSY.HydroGen}
    return "Hydro"
end

function get_generator_category(stor::GEN) where {GEN <: PSY.Storage}
    if (occursin("Distributed",PSY.get_name(stor)))
        return "Distributed_Storage"
    elseif (occursin("Battery",PSY.get_name(stor)))
        return "Battery_Storage"
    else
        return "Battery"
    end
end

function get_generator_category(stor::GEN) where {GEN <: PSY.HybridSystem}
    return "Hybrid-System"
end
##############################################
# Converting FOR and MTTR to λ and μ
##############################################
function outage_to_rate(outage_data::Tuple{Float64, Int64})
    for_gen = outage_data[1]
    mttr = outage_data[2]

    if (for_gen >1.0)
        for_gen = for_gen/100
    end

    if (mttr != 0)
        μ = 1 / mttr
    else
        μ = 0.0
    end
    λ = (μ * for_gen) / (1 - for_gen)
    #λ = for_gen

    return (λ = λ, μ = μ)
end

function make_pras_system(sys::PSY.System;
                          system_model::Union{Nothing, String} = nothing,aggregation::Union{Nothing, String} = nothing,
                          period_of_interest::Union{Nothing, UnitRange} = nothing,outage_flag=true,lump_pv_wind_gens=false,availability_flag=false, 
                          outage_csv_location::Union{Nothing, String} = nothing) 
    """
    make_pras_system(psy_sys,system_model)
    PSY System and System Model ("Single-Node","Zonal") are taken as arguments 
    and a PRAS SystemModel is returned.
    
    ...
    # Arguments
    - `psy_sys::PSY.System`: PSY System
    - `system_model::String`: "Single-Node" (or) "Zonal"
    - `aggregation::String`: "Area" (or) "LoadZone" {Optional} 
    - `num_time_steps::UnitRange`: Number of timesteps of PRAS SystemModel {Optional} 
    ...
    # Examples
    ```julia-repl
    julia> make_pras_system(psy_sys,system_model,period_of_interest)
    PRAS SystemModel
    ```
    """
    PSY.set_units_base_system!(sys, PSY.IS.UnitSystem.NATURAL_UNITS); # PRAS needs PSY System to be in NATURAL_UNITS
    #######################################################
    # Double counting of HybridSystem subcomponents
    #######################################################
    dup_uuids =[];
    h_s_comps = availability_flag ? PSY.get_components(PSY.HybridSystem, sys, PSY.get_available) : PSY.get_components(PSY.HybridSystem, sys)
    for h_s in h_s_comps
        h_s_subcomps = PSY._get_components(h_s)
        for subcomp in h_s_subcomps
            push!(dup_uuids,PSY.IS.get_uuid(subcomp))
        end
    end
    #######################################################
    # kwargs Handling
    #######################################################
    if (aggregation === nothing)
        aggregation_topology = "Area";
    end

    if (aggregation=="Area")
        aggregation_topology = PSY.Area;
    
    elseif (aggregation=="LoadZone")
        aggregation_topology = PSY.LoadZone;
    else
        error("Unrecognized PSY AggregationTopology")
    end

    all_ts = PSY.get_time_series_multiple(sys)
    first_ts_temp = first(all_ts);
    sys_ts_types = unique(typeof.(PSY.get_time_series_multiple(sys)));
    # Time series information
    sys_for_int_in_hour = round(Dates.Millisecond(PSY.get_forecast_interval(sys)), Dates.Hour)
    sys_res_in_hour = round(Dates.Millisecond(PSY.get_time_series_resolution(sys)), Dates.Hour)
    interval_len = Int(sys_for_int_in_hour.value/sys_res_in_hour.value)
    sys_horizon =  PSY.get_forecast_horizon(sys)
    #######################################################
    # Function to handle PSY timestamps
    #######################################################
    function get_period_of_interest(ts::TS) where {TS <: PSY.StaticTimeSeries}
        return range(1,length = length(ts.data))
    end

    function get_period_of_interest(ts::TS) where {TS <: PSY.Deterministic}
        return range(1,length = length(ts.data)*interval_len)
    end

    function get_len_ts_data(ts::TS) where {TS <: PSY.StaticTimeSeries}
        return length(ts.data)
    end

    function get_len_ts_data(ts::TS) where {TS <: PSY.Deterministic}
        return length(ts.data)*interval_len
    end
    
    # if (period_of_interest === nothing)
    #     period_of_interest = get_period_of_interest(first_ts_temp)
    # else
    #     if (PSY.Deterministic in sys_ts_types)
    #         if !(period_of_interest.start %  interval_len ==1 && period_of_interest.stop %  interval_len == 0)
    #             error("This PSY System has Determinstic time series data with interval length of $(interval_len). The period of interest should therefore be multiples of $(interval_len) to account for forecast windows.")
    #         end
    #     end
    # end

    N = length(period_of_interest);
    # len_ts_data = get_len_ts_data(first_ts_temp)

    # if ((N+(period_of_interest.start-1))> len_ts_data)
    #     error("Cannot make a PRAS System with $(N) timesteps with a PSY System with only $(length(first_ts_temp.data) - (period_of_interest.start-1)) timesteps of time series data")
    # end
    # if (period_of_interest.start >  len_ts_data || period_of_interest.stop >  len_ts_data)
    #     error("Please check the system period of interest selected")
    # end

    # Check if all time series data has a scaling_factor_multiplier
    # if(!all(.!isnothing.(getfield.(all_ts,:scaling_factor_multiplier))))
    #     error("Not all time series associated with components have scaling factor multipliers. This might lead to discrepancies in time series data in the PRAS System.")
    # end
    # if outage_csv_location is passed, perform some data checks
    outage_ts_flag = false
    if (outage_csv_location !== nothing)
        outage_ts_data,outage_ts_flag = try
            @info "Parsing the CSV with outage time series data ..."
            DataFrames.DataFrame(CSV.File(outage_csv_location)), true
        catch ex
            error("Couldn't parse the CSV with outage data at $(outage_csv_location).") 
            throw(ex)
        end
    end
    if (outage_ts_flag)
        if (N>DataFrames.nrow(outage_ts_data))
            @warn "Outage time series data is not available for all System timestamps in the CSV."
        end
    end
    #######################################################
    # PRAS timestamps
    # Need this to select timeseries values of interest
    #######################################################
    start_datetime = PSY.IS.get_initial_timestamp(first_ts_temp);
    start_datetime = start_datetime + Dates.Hour((period_of_interest.start-1)*sys_res_in_hour);
    start_datetime_tz = TimeZones.ZonedDateTime(start_datetime,TimeZones.tz"UTC");
    finish_datetime_tz = start_datetime_tz +  Dates.Hour((N-1)*sys_res_in_hour);
    my_timestamps = StepRange(start_datetime_tz, Dates.Hour(sys_res_in_hour), finish_datetime_tz);
    @info "The first timestamp of PRAS System being built is : $(start_datetime_tz) and last timestamp is : $(finish_datetime_tz) "
    det_ts_period_of_interest = 
    if (PSY.Deterministic in sys_ts_types)
        strt = 
        if (round(Int,period_of_interest.start/interval_len) ==0)
            1
        else
            round(Int,period_of_interest.start/interval_len)
        end
        stp = round(Int,period_of_interest.stop/interval_len)

        range(strt,length = (stp-strt)+1)
        
    end
    #######################################################
    # Common function to handle getting time series values
    #######################################################
    function get_forecast_values(ts::TS) where {TS <: PSY.Deterministic}
        forecast_vals = []
        for it in collect(keys(PSY.get_data(ts)))[det_ts_period_of_interest]
            append!(forecast_vals,collect(values(PSY.get_window(ts, it; len=interval_len))))
        end
        return forecast_vals
    end
   
    function get_forecast_values(ts::TS) where {TS <: PSY.StaticTimeSeries}
        forecast_vals = values(PSY.get_data(ts))[period_of_interest]
        return forecast_vals
    end
    #######################################################
     # PRAS Regions - Areas in SIIP
    #######################################################
    @info "Processing Regions in PSY System... "
    regions = collect(PSY.get_components(aggregation_topology, sys));
    if (length(regions)!=0)
        @info "The PSY System has $(length(regions)) regions based on PSY AggregationTopology : $(aggregation_topology)."
    else
        error("No regions in the PSY System. Cannot proceed with the process of making a PRAS SystemModel.")
    end 

    region_names = PSY.get_name.(regions);
    num_regions = length(region_names);

    region_load = Array{Int64,2}(undef,num_regions,N);
   
    for (idx,region) in enumerate(regions)
        reg_load_comps = availability_flag ? get_available_components_in_aggregation_topology(PSY.PowerLoad, sys, region) :
                                             PSY.get_components_in_aggregation_topology(PSY.PowerLoad, sys, region)
      
        region_load[idx,:]=floor.(Int,sum(get_forecast_values.(first.(PSY.get_time_series_multiple.(reg_load_comps, name = "max_active_power")))
                        .*PSY.get_max_active_power.(reg_load_comps))); # Any issues with using the first of time_series_multiple?
    end

    new_regions = PRAS.Regions{N,PRAS.MW}(region_names, region_load);

    #######################################################
    # kwargs Handling
    #######################################################
    if (system_model === nothing)
        if (num_regions>1)
            system_model = "Zonal"
        else
            system_model = "Single-Node"
        end
    end
    #######################################################
    # Generator Region Indices
    #######################################################
    gens=Array{PSY.Generator}[];
    start_id = Array{Int64}(undef,num_regions); 
    region_gen_idxs = Array{UnitRange{Int64},1}(undef,num_regions); 
    reg_wind_gens_DA = []
    reg_pv_gens_DA = []

    if (lump_pv_wind_gens)
        for (idx,region) in enumerate(regions)
            reg_ren_comps = availability_flag ? get_available_components_in_aggregation_topology(PSY.RenewableGen, sys, region) :
                                                 PSY.get_components_in_aggregation_topology(PSY.RenewableGen, sys, region)
            wind_gs_DA= [g for g in reg_ren_comps if (PSY.get_prime_mover(g) == PSY.PrimeMovers.WT)] 
            pv_gs_DA= [g for g in reg_ren_comps if (PSY.get_prime_mover(g) == PSY.PrimeMovers.PVe)] 
            reg_gen_comps = availability_flag ? get_available_components_in_aggregation_topology(PSY.Generator, sys, region) :
                                                PSY.get_components_in_aggregation_topology(PSY.Generator, sys, region)
            gs= [g for g in reg_gen_comps if (typeof(g) != PSY.HydroEnergyReservoir && PSY.get_max_active_power(g)!=0 && 
                                              PSY.IS.get_uuid(g) ∉ union(dup_uuids,PSY.IS.get_uuid.(wind_gs_DA),PSY.IS.get_uuid.(pv_gs_DA)))] 
            push!(gens,gs)
            push!(reg_wind_gens_DA,wind_gs_DA)
            push!(reg_pv_gens_DA,pv_gs_DA)

            if (idx==1)
                start_id[idx] = 1
            else 
                if (length(reg_wind_gens_DA[idx-1]) > 0 && length(reg_pv_gens_DA[idx-1]) > 0)
                    start_id[idx] =start_id[idx-1]+length(gens[idx-1])+2
                elseif (length(reg_wind_gens_DA[idx-1]) > 0 || length(reg_pv_gens_DA[idx-1]) > 0)
                    start_id[idx] =start_id[idx-1]+length(gens[idx-1])+1
                else
                    start_id[idx] =start_id[idx-1]+length(gens[idx-1])
                end
            end

            if (length(reg_wind_gens_DA[idx]) > 0 && length(reg_pv_gens_DA[idx]) > 0)
                region_gen_idxs[idx] = range(start_id[idx], length=length(gens[idx])+2)
            elseif (length(reg_wind_gens_DA[idx]) > 0 || length(reg_pv_gens_DA[idx]) > 0)
                region_gen_idxs[idx] = range(start_id[idx], length=length(gens[idx])+1)
            else
                region_gen_idxs[idx] = range(start_id[idx], length=length(gens[idx]))
            end
        end
    else
        for (idx,region) in enumerate(regions)
            reg_gen_comps = availability_flag ? get_available_components_in_aggregation_topology(PSY.Generator, sys, region) :
                                                PSY.get_components_in_aggregation_topology(PSY.Generator, sys, region)
            gs= [g for g in reg_gen_comps if (typeof(g) != PSY.HydroEnergyReservoir && PSY.get_max_active_power(g)!=0 && PSY.IS.get_uuid(g) ∉ dup_uuids)]
            push!(gens,gs)
            idx==1 ? start_id[idx] = 1 : start_id[idx] =start_id[idx-1]+length(gens[idx-1])
            region_gen_idxs[idx] = range(start_id[idx], length=length(gens[idx]))
        end
    end
    #######################################################
    # Storages Region Indices
    #######################################################
    stors=[];
    start_id = Array{Int64}(undef,num_regions);
    region_stor_idxs = Array{UnitRange{Int64},1}(undef,num_regions);

    for (idx,region) in enumerate(regions)
        #push!(stors,[s for s in PSY.get_components_in_aggregation_topology(PSY.Storage, sys, region)])
        reg_stor_comps = availability_flag ? get_available_components_in_aggregation_topology(PSY.Storage, sys, region) :
                                             PSY.get_components_in_aggregation_topology(PSY.Storage, sys, region)
        push!(stors,[s for s in reg_stor_comps if (PSY.IS.get_uuid(s) ∉ dup_uuids)])
        idx==1 ? start_id[idx] = 1 : start_id[idx] =start_id[idx-1]+length(stors[idx-1])
        region_stor_idxs[idx] = range(start_id[idx], length=length(stors[idx]))
    end
    #######################################################
    # GeneratorStorages Region Indices
    #######################################################
    gen_stors=[];
    start_id = Array{Int64}(undef,num_regions);
    region_genstor_idxs = Array{UnitRange{Int64},1}(undef,num_regions);

    for (idx,region) in enumerate(regions)
        reg_gen_stor_comps = availability_flag ? get_available_components_in_aggregation_topology(PSY.StaticInjection, sys, region) :
                                                 PSY.get_components_in_aggregation_topology(PSY.StaticInjection, sys, region)
        gs= [g for g in reg_gen_stor_comps if (typeof(g) == PSY.HydroEnergyReservoir || typeof(g)==PSY.HybridSystem)]
        push!(gen_stors,gs)
        idx==1 ? start_id[idx] = 1 : start_id[idx] =start_id[idx-1]+length(gen_stors[idx-1])
        region_genstor_idxs[idx] = range(start_id[idx], length=length(gen_stors[idx]))
    end
    #######################################################
    # PRAS Generators
    #######################################################
    @info "Processing Generators in PSY System... "
    
    # Lumping Wind and PV Generators per Region
    if (lump_pv_wind_gens)
        for i in 1: num_regions
            if (length(reg_wind_gens_DA[i])>0)
                # Wind
                temp_lumped_wind_gen = PSY.RenewableDispatch(nothing)
                PSY.set_name!(temp_lumped_wind_gen,"Lumped_Wind_"*region_names[i])
                PSY.set_prime_mover!(temp_lumped_wind_gen,PSY.PrimeMovers.WT)
                ext = PSY.get_ext(temp_lumped_wind_gen)
                ext["region_gens"] = reg_wind_gens_DA[i]
                ext["outage_probability"] = 0.0
                ext["recovery_probability"] = 1.0
                push!(gens[i],temp_lumped_wind_gen)
            end
            if (length(reg_pv_gens_DA[i])>0)
                # PV
                temp_lumped_pv_gen = PSY.RenewableDispatch(nothing)
                PSY.set_name!(temp_lumped_pv_gen,"Lumped_PV_"*region_names[i])
                PSY.set_prime_mover!(temp_lumped_pv_gen,PSY.PrimeMovers.PVe)
                ext = PSY.get_ext(temp_lumped_pv_gen)
                ext["region_gens"] = reg_pv_gens_DA[i]
                ext["outage_probability"] = 0.0
                ext["recovery_probability"] = 1.0
                push!(gens[i],temp_lumped_pv_gen)
            end
        end
    end

    gen=[];
    for i in 1: num_regions
        if (length(gens[i]) != 0)
            append!(gen,gens[i])
        end
    end
    
    if(length(gen) ==0)
        gen_names = String[];
    else
        gen_names = PSY.get_name.(gen);
    end

    gen_categories = string.(typeof.(gen));
    n_gen = length(gen_names);

    gen_cap_array = Matrix{Int64}(undef, n_gen, N);
    λ_gen = Matrix{Float64}(undef, n_gen, N);
    μ_gen = Matrix{Float64}(undef, n_gen, N);

    for (idx,g) in enumerate(gen)
        # Nominal outage and recovery rate
        (λ,μ) = (0.0,1.0)
        
        if (lump_pv_wind_gens && (PSY.get_prime_mover(g) == PSY.PrimeMovers.WT || PSY.get_prime_mover(g) == PSY.PrimeMovers.PVe))
            reg_gens_DA = PSY.get_ext(g)["region_gens"];
            gen_cap_array[idx,:] = round.(Int,sum(get_forecast_values.(first.(PSY.get_time_series_multiple.(reg_gens_DA, name = "max_active_power")))
                                   .*PSY.get_max_active_power.(reg_gens_DA)));
        else
            if (PSY.has_time_series(g) && ("max_active_power" in PSY.get_name.(PSY.get_time_series_multiple(g))))
                gen_cap_array[idx,:] = floor.(Int,get_forecast_values(first(PSY.get_time_series_multiple(g, name = "max_active_power")))
                                       *PSY.get_max_active_power(g));
            else
                gen_cap_array[idx,:] = fill.(floor.(Int,PSY.get_max_active_power(g)),1,N);
            end
        end

        if (outage_ts_flag)
            try
                # @info "Using FOR time series data for $(PSY.get_name(g)) of type $(gen_categories[idx]). Assuming the mean time to recover (MTTR) is 24 hours to compute the λ and μ time series data ..."
                g_λ_μ_ts_data = outage_to_rate.(zip(outage_ts_data[!,PSY.get_name(g)],fill(24,length(outage_ts_data[!,PSY.get_name(g)]))))
                λ_gen[idx,:] = getfield.(g_λ_μ_ts_data,:λ)[period_of_interest]
                μ_gen[idx,:] = getfield.(g_λ_μ_ts_data,:μ)[period_of_interest] # This assumes a mean time to recover of 24 hours.
            catch ex
                # @warn "FOR time series data for $(PSY.get_name(g)) of type $(gen_categories[idx]) is not available in the CSV. Using nominal outage and recovery probabilities for this generator." 
                λ_gen[idx,:] = fill.(λ,1,N); 
                μ_gen[idx,:] = fill.(μ,1,N);
            end
        else
            if (~outage_flag)
                if (gen_categories[idx] == "ThermalStandard")
                    p_m = string(PSY.get_prime_mover(g))
                    fl = string(PSY.get_fuel(g))

                    p_m_idx = findall(x -> x == p_m, getfield.(outage_values,:prime_mover))
                    fl_idx =  findall(x -> x == fl, getfield.(outage_values[p_m_idx],:thermal_fuel))
                    
                    if (length(fl_idx) ==0)
                        fl_idx =  findall(x -> x == "NA", getfield.(outage_values[p_m_idx],:thermal_fuel))
                    end

                    temp_range = p_m_idx[fl_idx]

                    temp_cap = floor(Int,PSY.get_max_active_power(g))

                    if (length(temp_range)>1)
                        gen_idx = temp_range[1]
                        for (x,y) in zip(temp_range,getfield.(outage_values[temp_range],:capacity))
                            temp=0
                            if (temp<temp_cap<y)
                                gen_idx = x
                                break
                            else
                                temp = y
                            end
                        end
                        f_or = getfield(outage_values[gen_idx],:FOR)
                        mttr_hr = getfield(outage_values[gen_idx],:MTTR)

                        (λ,μ) = outage_to_rate((f_or,mttr_hr))

                    elseif (length(temp_range)==1)
                        gen_idx = temp_range[1]
                        
                        f_or = getfield(outage_values[gen_idx],:FOR)
                        mttr_hr = getfield(outage_values[gen_idx],:MTTR)

                        (λ,μ) = outage_to_rate((f_or,mttr_hr))
                    else
                        # @warn "No outage information is available for $(PSY.get_name(g)) with a $(p_m) prime mover and $(fl) fuel type. Using nominal outage and recovery probabilities for this generator."
                        #λ = 0.0;
                        #μ = 1.0;
                    end

                elseif (gen_categories[idx] == "HydroDispatch")
                    p_m = string(PSY.get_prime_mover(g))
                    p_m_idx = findall(x -> x == p_m, getfield.(outage_values,:prime_mover))

                    temp_cap = floor(Int,PSY.get_max_active_power(g))
                    
                    if (length(p_m_idx)>1)
                        for (x,y) in zip(p_m_idx,getfield.(outage_values[p_m_idx],:capacity))
                            temp=0
                            if (temp<temp_cap<y)
                                gen_idx = x

                                f_or = getfield(outage_values[gen_idx],:FOR)
                                mttr_hr = getfield(outage_values[gen_idx],:MTTR)

                                (λ,μ) = outage_to_rate((f_or,mttr_hr))
                                break
                            else
                                temp = y
                            end
                        end
                    end
                else
                    # @warn "No outage information is available for $(PSY.get_name(g)) of type $(gen_categories[idx]). Using nominal outage and recovery probabilities for this generator."
                    #λ = 0.0;
                    #μ = 1.0;

                end

            else
                ext = PSY.get_ext(g)
                if (!(haskey(ext,"outage_probability") && haskey(ext,"recovery_probability")))
                    # @warn "No outage information is available in ext field of $(PSY.get_name(g)) of type $(gen_categories[idx]). Using nominal outage and recovery probabilities for this generator."
                    #λ = 0.0;
                    #μ = 1.0;
                else
                    λ = ext["outage_probability"];
                    μ = ext["recovery_probability"];
                end
            end
            λ_gen[idx,:] = fill.(λ,1,N); 
            μ_gen[idx,:] = fill.(μ,1,N); 
        end
    end

    new_generators = PRAS.Generators{N,1,PRAS.Hour,PRAS.MW}(gen_names, get_generator_category.(gen), gen_cap_array , λ_gen ,μ_gen);
        
    #######################################################
    # PRAS Storages
    # **TODO Future : time series for storage devices
    #######################################################
    @info "Processing Storages in PSY System... "

    stor=[];
    for i in 1: num_regions
        if (length(stors[i]) != 0)
            append!(stor,stors[i])
        end
    end

    if(length(stor) ==0)
        stor_names=String[];
    else
        stor_names = PSY.get_name.(stor);
    end

    stor_categories = string.(typeof.(stor));

    n_stor = length(stor_names);

    stor_charge_cap_array = Matrix{Int64}(undef, n_stor, N);
    stor_discharge_cap_array = Matrix{Int64}(undef, n_stor, N);
    stor_energy_cap_array = Matrix{Int64}(undef, n_stor, N);
    stor_chrg_eff_array = Matrix{Float64}(undef, n_stor, N);
    stor_dischrg_eff_array  = Matrix{Float64}(undef, n_stor, N);
    λ_stor = Matrix{Float64}(undef, n_stor, N);   
    μ_stor = Matrix{Float64}(undef, n_stor, N);

    for (idx,s) in enumerate(stor)
        stor_charge_cap_array[idx,:] = fill(floor(Int,getfield(PSY.get_input_active_power_limits(s), :max)),1,N);
        stor_discharge_cap_array[idx,:] = fill(floor(Int,getfield(PSY.get_output_active_power_limits(s), :max)),1,N);
        stor_energy_cap_array[idx,:] = fill(floor(Int,getfield(PSY.get_state_of_charge_limits(s),:max)),1,N);
        stor_chrg_eff_array[idx,:] = fill(getfield(PSY.get_efficiency(s), :in),1,N);
        stor_dischrg_eff_array[idx,:]  = fill.(getfield(PSY.get_efficiency(s), :out),1,N);

        if (~outage_flag)
            # @warn "No outage information is available for $(PSY.get_name(s)) of type $(stor_categories[idx]). Using nominal outage and recovery probabilities for this generator."
            λ = 0.0;
            μ = 1.0;
        else
            ext = PSY.get_ext(s)
            if (!(haskey(ext,"outage_probability") && haskey(ext,"recovery_probability")))
                # @warn "No outage information is available in ext field of $(PSY.get_name(s)) of type $(stor_categories[idx]). Using nominal outage and recovery probabilities for this generator."
                λ = 0.0;
                μ = 1.0;
            else
                λ = ext["outage_probability"];
                μ = ext["recovery_probability"];
            end
        end
        
        λ_stor[idx,:] = fill.(λ,1,N); 
        μ_stor[idx,:] = fill.(μ,1,N); 
    end
    
    stor_cryovr_eff = ones(n_stor,N);   # Not currently available/ defined in PowerSystems
    
    new_storage = PRAS.Storages{N,1,PRAS.Hour,PRAS.MW,PRAS.MWh}(stor_names,get_generator_category.(stor),
                                            stor_charge_cap_array,stor_discharge_cap_array,stor_energy_cap_array,
                                            stor_chrg_eff_array,stor_dischrg_eff_array, stor_cryovr_eff,
                                            λ_stor,μ_stor);

    #######################################################
    # PRAS Generator Storages
    # **TODO Consider all combinations of HybridSystem (Currently only works for DER+ESS)
    #######################################################
    @info "Processing GeneratorStorages in PSY System... "

    gen_stor=[];
    for i in 1: num_regions
        if (length(gen_stors[i]) != 0)
            append!(gen_stor,gen_stors[i])
        end
    end
    
    if(length(gen_stor) == 0)
        gen_stor_names=String[];
    else
        gen_stor_names = PSY.get_name.(gen_stor);
    end

    gen_stor_categories = string.(typeof.(gen_stor)); 
    
    n_genstors = length(gen_stor_names);

    gen_stor_charge_cap_array = Matrix{Int64}(undef, n_genstors, N);
    gen_stor_discharge_cap_array = Matrix{Int64}(undef, n_genstors, N);
    gen_stor_enrgy_cap_array = Matrix{Int64}(undef, n_genstors, N);
    gen_stor_inflow_array = Matrix{Int64}(undef, n_genstors, N);
    gen_stor_gridinj_cap_array = Matrix{Int64}(undef, n_genstors, N);

    λ_genstors = Matrix{Float64}(undef, n_genstors, N);   
    μ_genstors = Matrix{Float64}(undef, n_genstors, N);  

    for (idx,g_s) in enumerate(gen_stor)
        if(typeof(g_s) ==PSY.HydroEnergyReservoir)
            if (PSY.has_time_series(g_s))
                if ("inflow" in PSY.get_name.(PSY.get_time_series_multiple(g_s)))
                    gen_stor_charge_cap_array[idx,:] = floor.(Int,get_forecast_values(first(PSY.get_time_series_multiple(g_s, name = "inflow")))
                                                       *PSY.get_inflow(g_s));
                    gen_stor_discharge_cap_array[idx,:] = floor.(Int,get_forecast_values(first(PSY.get_time_series_multiple(g_s, name = "inflow")))
                                                          *PSY.get_inflow(g_s));
                    gen_stor_inflow_array[idx,:] = floor.(Int,get_forecast_values(first(PSY.get_time_series_multiple(g_s, name = "inflow")))
                                                   *PSY.get_inflow(g_s));
                else
                    gen_stor_charge_cap_array[idx,:] = fill.(floor.(Int,PSY.get_inflow(g_s)),1,N);
                    gen_stor_discharge_cap_array[idx,:] = fill.(floor.(Int,PSY.get_inflow(g_s)),1,N);
                    gen_stor_inflow_array[idx,:] = fill.(floor.(Int,PSY.get_inflow(g_s)),1,N);
                end
                if ("storage_capacity" in PSY.get_name.(PSY.get_time_series_multiple(g_s)))
                    gen_stor_enrgy_cap_array[idx,:] = floor.(Int,get_forecast_values(first(PSY.get_time_series_multiple(g_s, name = "storage_capacity")))
                                                      *PSY.get_storage_capacity(g_s));
                else
                    gen_stor_enrgy_cap_array[idx,:] = fill.(floor.(Int,PSY.get_storage_capacity(g_s)),1,N);
                end
                if ("max_active_power" in PSY.get_name.(PSY.get_time_series_multiple(g_s)))
                    gen_stor_gridinj_cap_array[idx,:] = floor.(Int,get_forecast_values(first(PSY.get_time_series_multiple(g_s, name = "max_active_power")))
                                                        *PSY.get_max_active_power(g_s));
                else
                    gen_stor_gridinj_cap_array[idx,:] = fill.(floor.(Int,PSY.get_max_active_power(g_s)),1,N);
                end
            else
                gen_stor_charge_cap_array[idx,:] = fill.(floor.(Int,PSY.get_inflow(g_s)),1,N);
                gen_stor_discharge_cap_array[idx,:] = fill.(floor.(Int,PSY.get_inflow(g_s)),1,N);
                gen_stor_enrgy_cap_array[idx,:] = fill.(floor.(Int,PSY.get_storage_capacity(g_s)),1,N);
                gen_stor_inflow_array[idx,:] = fill.(floor.(Int,PSY.get_inflow(g_s)),1,N);
                gen_stor_gridinj_cap_array[idx,:] = fill.(floor.(Int,PSY.get_max_active_power(g_s)),1,N);
            end  
        else
            gen_stor_charge_cap_array[idx,:] = fill.(floor.(Int,getfield(PSY.get_input_active_power_limits(PSY.get_storage(g_s)), :max)),1,N);
            gen_stor_discharge_cap_array[idx,:] = fill.(floor.(Int,getfield(PSY.get_output_active_power_limits(PSY.get_storage(g_s)), :max)),1,N);
            gen_stor_enrgy_cap_array[idx,:] = fill.(floor.(Int,getfield(PSY.get_state_of_charge_limits(PSY.get_storage(g_s)), :max)),1,N); 
            gen_stor_gridinj_cap_array[idx,:] = fill.(floor.(Int,PSY.getfield(PSY.get_output_active_power_limits(g_s), :max)),1,N);
            
            if (PSY.has_time_series(PSY.get_renewable_unit(g_s)) && ("max_active_power" in PSY.get_name.(PSY.get_time_series_multiple(PSY.get_renewable_unit(g_s)))))
                gen_stor_inflow_array[idx,:] = floor.(Int,get_forecast_values(first(PSY.get_time_series_multiple(PSY.get_renewable_unit(g_s), name = "max_active_power")))
                                               *PSY.get_max_active_power(PSY.get_renewable_unit(g_s))); 
            else
                gen_stor_inflow_array[idx,:] = fill.(floor.(Int,PSY.get_max_active_power(PSY.get_renewable_unit(g_s))),1,N); 
            end
        end
        
        if (~outage_flag)
            if (typeof(g_s) ==PSY.HydroEnergyReservoir)
                p_m = string(PSY.get_prime_mover(g_s))
                p_m_idx = findall(x -> x == p_m, getfield.(outage_values,:prime_mover))

                temp_cap = floor(Int,PSY.get_max_active_power(g_s))
                
                if (length(p_m_idx)>1)
                    for (x,y) in zip(p_m_idx,getfield.(outage_values[p_m_idx],:capacity))
                        temp=0
                        if (temp<temp_cap<y)
                            gen_idx = x

                            f_or = getfield(outage_values[gen_idx],:FOR)
                            mttr_hr = getfield(outage_values[gen_idx],:MTTR)

                            (λ,μ) = outage_to_rate((f_or,mttr_hr))
                            break
                        else
                            temp = y
                        end
                    end
                end
            else
                # @warn "No outage information is available for $(PSY.get_name(g_s)) of type $(gen_stor_categories[idx]). Using nominal outage and recovery probabilities for this generator."
                λ = 0.0;
                μ = 1.0;
            end

        else
            ext = PSY.get_ext(g_s)
            if (!(haskey(ext,"outage_probability") && haskey(ext,"recovery_probability")))
                # @warn "No outage information is available in ext field of $(PSY.get_name(g_s)) of type $(gen_stor_categories[idx]). Using nominal outage and recovery probabilities for this generator."
                λ = 0.0;
                μ = 1.0;
            else
                λ = ext["outage_probability"];
                μ = ext["recovery_probability"];
            end
        end
        
        λ_genstors[idx,:] = fill.(λ,1,N); 
        μ_genstors[idx,:] = fill.(μ,1,N);
    end
    
    gen_stor_gridwdr_cap_array = zeros(Int64,n_genstors, N); # Not currently available/ defined in PowerSystems
    gen_stor_charge_eff = ones(n_genstors,N);                # Not currently available/ defined in PowerSystems
    gen_stor_discharge_eff = ones(n_genstors,N);             # Not currently available/ defined in PowerSystems
    gen_stor_cryovr_eff = ones(n_genstors,N);                # Not currently available/ defined in PowerSystems

    
    new_gen_stors = PRAS.GeneratorStorages{N,1,PRAS.Hour,PRAS.MW,PRAS.MWh}(gen_stor_names,get_generator_category.(gen_stor),
                                                    gen_stor_charge_cap_array, gen_stor_discharge_cap_array, gen_stor_enrgy_cap_array,
                                                    gen_stor_charge_eff, gen_stor_discharge_eff, gen_stor_cryovr_eff,
                                                    gen_stor_inflow_array, gen_stor_gridwdr_cap_array, gen_stor_gridinj_cap_array,
                                                    λ_genstors, μ_genstors);

    #######################################################
    # PRAS SystemModel
    #######################################################
    if (system_model=="Zonal")
    #######################################################
    # PRAS Lines 
    #######################################################
    @info "Collecting all inter regional lines in PSY System..."

    # Dictionary with topology mapping
        line = availability_flag ? 
        collect(PSY.get_components(PSY.Branch, sys, (x -> ~in(typeof(x), [PSY.TapTransformer, PSY.Transformer2W,PSY.PhaseShiftingTransformer]) && PSY.get_available(x)))) :
        collect(PSY.get_components(PSY.Branch, sys, x -> ~in(typeof(x), [PSY.TapTransformer, PSY.Transformer2W,PSY.PhaseShiftingTransformer])));

        mapping_dict = PSY.get_aggregation_topology_mapping(aggregation_topology,sys); # Dict with mapping from Areas to Bus_Names
        new_mapping_dict=Dict{String,Array{Int64,1}}(); 

        for key in keys(mapping_dict)
            push!(new_mapping_dict, key  => PSY.get_number.(mapping_dict[key]))
        end
        
        #######################################################
        # Finding the inter-regional lines and regions_from
        #######################################################
        regional_lines = []; 
        regions_from = [];
        for i in 1:length(line)
            for key in keys(mapping_dict)
                if(PSY.get_number(line[i].arc.from) in new_mapping_dict[key])
                    if(~(PSY.get_number(line[i].arc.to) in new_mapping_dict[key]))
                        push!(regional_lines,line[i])
                        push!(regions_from,key)
                    end
                end
            end
        end
        
        #######################################################
        # Finding the regions_to
        #######################################################
        n_lines = length(regional_lines);
        regions_to = [];
        for i in 1:n_lines
            for key in keys(mapping_dict)
                if(PSY.get_number(regional_lines[i].arc.to) in new_mapping_dict[key])
                    push!(regions_to,key)
                end
            end
        end
        
        #######################################################
        # If there are lines from region1 --> region3 and 
        # region3 --> region1; lines like these need to be grouped
        # to make interface_line_idxs
        #######################################################
        regions_tuple = [];
        for i in 1:length(regions_from)
            region_from_idx = findfirst(x->x==regions_from[i],region_names)
            region_to_idx = findfirst(x->x==regions_to[i],region_names)
            if (region_from_idx < region_to_idx)
                push!(regions_tuple,(regions_from[i],regions_to[i]))
            else
                push!(regions_tuple,(regions_to[i],regions_from[i]))
            end
        end

        temp_regions_tuple = unique(regions_tuple);
        interface_dict = Dict();

        for i in 1: length(temp_regions_tuple)
            temp = findall(x -> x == temp_regions_tuple[i], regions_tuple);
            push!(interface_dict, temp_regions_tuple[i] => (temp,length(temp)))
        end

        num_interfaces = length(temp_regions_tuple);
        sorted_regional_lines = [];
        interface_line_idxs = Array{UnitRange{Int64},1}(undef,num_interfaces);
        start_id = Array{Int64}(undef,num_interfaces); 
        for i in 1: num_interfaces
            for j in interface_dict[temp_regions_tuple[i]][1]
                push!(sorted_regional_lines, regional_lines[j])
            end
            i==1 ? start_id[i] = 1 : start_id[i] =start_id[i-1]+interface_dict[temp_regions_tuple[i-1]][2]
            interface_line_idxs[i] = range(start_id[i], length=interface_dict[temp_regions_tuple[i]][2])
        end
        
        @info "Processing all inter regional lines in PSY System..."
        line_names = PSY.get_name.(sorted_regional_lines);
        line_categories = string.(typeof.(sorted_regional_lines));
        
        line_forward_capacity_array = Matrix{Int64}(undef, n_lines, N);
        line_backward_capacity_array = Matrix{Int64}(undef, n_lines, N);

        λ_lines = Matrix{Float64}(undef, n_lines, N); # Not currently available/ defined in PowerSystems
        μ_lines = Matrix{Float64}(undef, n_lines, N); # Not currently available/ defined in PowerSystems
        for i in 1:n_lines
            line_forward_capacity_array[i,:] = fill.(floor.(Int,getfield(line_rating(sorted_regional_lines[i]),:forward_capacity)),1,N);
            line_backward_capacity_array[i,:] = fill.(floor.(Int,getfield(line_rating(sorted_regional_lines[i]),:backward_capacity)),1,N);

            λ_lines[i,:] .= 0.0; # Not currently available/ defined in PowerSystems # should change when we have this
            μ_lines[i,:] .= 1.0; # Not currently available/ defined in PowerSystems
        end
        
        new_lines = PRAS.Lines{N,1,PRAS.Hour,PRAS.MW}(line_names, line_categories, line_forward_capacity_array, line_backward_capacity_array, λ_lines ,μ_lines);
        #######################################################
        # PRAS Interfaces
        #######################################################
        interface_regions_from = [findfirst(x->x==temp_regions_tuple[i][1],region_names) for i in 1:num_interfaces];
        interface_regions_to = [findfirst(x->x==temp_regions_tuple[i][2],region_names) for i in 1:num_interfaces];
        
        @info "Processing interfaces from inter regional lines in PSY System..."
        interface_forward_capacity_array = Matrix{Int64}(undef, num_interfaces, N);
        interface_backward_capacity_array = Matrix{Int64}(undef, num_interfaces, N);
        for i in 1:num_interfaces
            interface_forward_capacity_array[i,:] =  sum(line_forward_capacity_array[interface_line_idxs[i],:],dims=1)
            interface_backward_capacity_array[i,:] =  sum(line_backward_capacity_array[interface_line_idxs[i],:],dims=1)
        end

        new_interfaces = PRAS.Interfaces{N,PRAS.MW}(interface_regions_from, interface_regions_to, interface_forward_capacity_array, interface_backward_capacity_array);

        pras_system = PRAS.SystemModel(new_regions, new_interfaces, new_generators, region_gen_idxs, new_storage, region_stor_idxs, new_gen_stors,
                          region_genstor_idxs, new_lines,interface_line_idxs,my_timestamps);
    
    elseif (system_model =="Single-Node")
        load_vector = vec(sum(region_load,dims=1));
        pras_system = PRAS.SystemModel(new_generators, new_storage, new_gen_stors, my_timestamps, load_vector);
    else
        error("Unrecognized SystemModel; Please specify correctly if SystemModel is Single-Node or Zonal.")
    end
    @info "Successfully built a PRAS $(system_model) system of type $(typeof(pras_system))."
    return pras_system
end



function generate_csv_outage_profile(pras_system::PRAS.SystemModel;
                                     location::Union{Nothing, String} = nothing,num_runs::Union{Nothing, Int} = nothing,num_scenarios::Union{Nothing, Int} = nothing) 
    """
    generate_outage_profile(pras_system,num_runs,psy_sys,num_scenarios,location)
    Process the assess results to get timeseries of generator status and include 
    this timeseries data to the corresponding component in PSY System and exported
    using to_json method (serializing the PSY System).
    
    ...
    # Arguments
    - `pras_system::PRAS.SystemModel`: PRAS System
    - `num_runs::Int64`: Number of PRAS runs
    - `psy_sys::PSY.System`: PSY System
    - `num_scenarios::Int64`: Number of scenarios of user interest.
    - `location::String`: Location to store outage profile.
    ...
    # Examples
    ```julia-repl
    julia> generate_outage_profile(results,pras_sys,psy_sys,1)
    PSY System exported using to_json method in InfrastructureSystems
    ```
    """
    #kwargs handling
    if (location === nothing)
        location = dirname(dirname(@__DIR__))
        @warn  "Location to save the exported PSY System not specified. Using the data folder of the module."
    end

    if (num_runs === nothing)
        num_runs = 10000
        @warn  "Number of samples of PRAS Assess not specified. Using 10,000 samples."
    end

    if (num_scenarios === nothing)
        num_scenarios = 1
        @warn  "Number of scenarios to be exported not specified. Only exporting one scenario."
    end

    # Run PRAS Analysis
    @info "Running PRAS SequentialMonteCarlo Resource Adequacy analysis for $(num_runs) runs..."
    shortfall_samples,gens_avail,stors_avail,gen_stors_avail,lines_avail = PRAS.assess(pras_system, PRAS.SequentialMonteCarlo(samples=num_runs),  PRAS.ShortfallSamples(),  PRAS.GeneratorAvailability(),PRAS.StorageAvailability(),PRAS.GeneratorStorageAvailability(),PRAS.LineAvailability())
    @info "Successfully completed PRAS Runs. Exporting outage profiles now..."
    # Setup to save the new PSY System with scenario timeseries data
    #working_dir = pwd();
    working_dir = location;
    dt_now = Dates.format(Dates.now(),"dd-u-yy-H-M-S");
    
    dir_name = joinpath(working_dir,"data","Generated-Outage-Profile-JSON",string(UUIDs.uuid4()),dt_now)
   
    #  **TODO: IF generating systems for multiple scenarios, remember to delete the timeseries data from previous scenario.
    asset_dict = Dict([(:generators, (gens_avail,PSY.Generator)), (:storages, (stors_avail,PSY.Storage)),(:generatorstorages, (gen_stors_avail,PSY.StaticInjection)), (:lines, (lines_avail,PSY.Branch))]);
    asset_keys=[];
    
    for k in keys(asset_dict)
        if (length(getfield(pras_system,k)) !=0)
            push!(asset_keys,k)
        end
    end
    @info "Scenarios of interest will be sorted according to sample unserved energy. The availability for individual asset types will be exported to CSV sheets."
    
    sample_idx = sortperm(shortfall_samples[],rev=true);

    for i in 1:num_scenarios
        
        mkpath(joinpath(dir_name,string(i)));

        for k in asset_keys
            asset_status_all = asset_dict[k][1].available[:,:,sample_idx[i]];
            asset_names = getfield(asset_dict[k][1],k)

            df_outage = DataFrames.DataFrame()

            for (j,asset_name) in enumerate(asset_names)
                df_outage[!,asset_name] = asset_status_all[j,:]
            end

            csv_path = joinpath(dir_name,string(i),string(asset_dict[k][2],".csv"))
            CSV.write(csv_path, df_outage,writeheader = true)
        end
        @info "Succesfully exported the outage profile for scenario $(i)."
    end
    @info "Succesfully exported the outage profiles. The CSV files for all asset types for all scenarios are available here : $(dir_name)."
end

#######################################################
# Main Function to make the PRAS System
#######################################################
function add_csv_time_series!(sys_DA::PSY.System,sys_RT::PSY.System, outage_csv_location::String;days_of_interest::Union{Nothing, UnitRange} = nothing,
                              add_scenario::Union{Nothing, Int} = nothing,iteration_year::Int64)
    #######################################################
    # kwarg handling
    #######################################################
    @info "This script currently only adds availability time series to Generators in PSY System. This can be extended to handle Lines, Storage and Static Injections. "

    if (~isdir(outage_csv_location))
        error("Outage CSV location passed isn't a directory.")
    end

    (root, dirs, files) = first(walkdir(outage_csv_location))

    if (length(dirs) == 0)
        error("No scenario outage data available")
    end

    if (add_scenario === nothing)
        add_scenario = 1
        @warn  "The scenario to be used to get the asset availability time series data wasn't passed. Using availability time series data for scenario '1'."
    end

    dir_location = joinpath(root,string(add_scenario))
    csv_locations = readdir(dir_location)

    # if (~("Generator_year$(iteration_year).csv" in csv_locations))
    #     error("Generator Availability time series data not available for scenario $(add_scenario).")
    # end
    #file_location = joinpath(dir_location,csv_locations[1])
    file_location = joinpath(dir_location,"Generator_year$(iteration_year).csv")

    # Handling days of interest to generate time stamps
    first_ts_temp_DA = first(PSY.get_time_series_multiple(sys_DA));
    len_ts_data_DA = length(first_ts_temp_DA.data);

    if (days_of_interest === nothing)
        days_of_interest = range(1,length = len_ts_data_DA)
        @warn  "No days_of_interest passed. Using all the days for which data is available in PSY System. Please make sure these days of interest align with the availability data in the CSV file."
    end

     #Select Full Days to avoid more Processing
    if (~(1<=days_of_interest.start <= len_ts_data_DA) || ~(1<=days_of_interest.stop <= len_ts_data_DA))
        error("Please checked the passed days of interest (e.g 1:24,25:48,1:365 etc.)")
    end

    num_days = length(days_of_interest);
    period_of_interest = range(days_of_interest.start==1 ? 1 : (24*days_of_interest.start)+1,length=num_days*24)
    N = length(period_of_interest);
    #######################################################
    # Timestamps
    #######################################################
    start_datetime_DA = PSY.IS.get_initial_timestamp(first_ts_temp_DA);
    sys_DA_res_in_hour = PSY.get_time_series_resolution(sys_DA)
    start_datetime_DA = start_datetime_DA + Dates.Hour((period_of_interest.start-1)*sys_DA_res_in_hour);
    finish_datetime_DA = start_datetime_DA +  Dates.Hour((N-1)*sys_DA_res_in_hour);
    all_timestamps = StepRange(start_datetime_DA, sys_DA_res_in_hour, finish_datetime_DA);
    ######################################################
    # Reading the CSV Files
    #######################################################
    df_outage_profile = DataFrames.DataFrame(CSV.File(file_location));

    for asset_name in names(df_outage_profile)
        # Creating TimeSeries Data
        DA_data = Dict()
        RT_data = Dict()

        for (idx,timestamp) in enumerate(all_timestamps) #
            if (rem(idx,24) == 1)
                if (df_outage_profile[idx,asset_name] ==0)
                    push!(DA_data,timestamp => zeros(36))
                else
                    push!(DA_data,timestamp => ones(36))
                end
            end

            if (df_outage_profile[idx,asset_name] ==0)
                # If PRAS actual status is actually 0, just use 0 for RT, ReliablePSY can handle this transition
                push!(RT_data,timestamp => zeros(2))
            else
                # If PRAS actual status is actually 1 and DA says 1, use 1 for RT, if DA says 0, use O because ReliablePSY currently can't handle this
                if (rem(idx,24) == 0)
                    offset = 24
                else
                    offset = rem(idx,24)
                end

                if (df_outage_profile[(idx-offset+1),asset_name] ==1) # If PRAS Actual == 1 and DA ==1
                    push!(RT_data,timestamp => ones(2))
                else
                    push!(RT_data,timestamp => ones(2)) # # If PRAS Actual == 1 and DA ==0 - Special case
                end
            end
        end

        DA_availability_forecast = PSY.Deterministic("outage", DA_data,sys_DA_res_in_hour)
        RT_availability_forecast = PSY.Deterministic("outage", RT_data,sys_DA_res_in_hour)

        # Adding TimeSeries Data to PSY System
        try
            # @info "Adding availability time series data for $(asset_name) Generator to DA System..."
            PSY.add_time_series!(sys_DA,PSY.get_component(PSY.Generator,sys_DA, asset_name), DA_availability_forecast)
        catch ex
            # @warn "Couldn't find a Generator with component name $(asset_name) in the DA PSY System. Proceeding without adding availability time series data for this component in the DA System."
        end

        try
            # @info "Adding availability time series data for $(asset_name) Generator to RT System..."
            PSY.add_time_series!(sys_RT,PSY.get_component(PSY.Generator,sys_RT, asset_name), RT_availability_forecast)
        catch ex
            # @warn "Couldn't find a Generator with component name $(asset_name) in the RT PSY System. Proceeding without adding availability time series data for this component in the RT System."
        end
    end

    @info "Succesfully added availability time series data to Generators in PSY DA and RT Systems."
    return sys_DA, sys_RT
end
