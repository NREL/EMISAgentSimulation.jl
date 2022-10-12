function construct_smc_unavailabilities(sys::PSY.System, ordc_unavailability_method::String)
    if ordc_unavailability_method == "SMC"
        system_period_of_interest = range(1, length = 8784);
        pras_system = make_pras_system(sys,
                                       system_model = "Single-Node",
                                       aggregation = "Area",
                                       period_of_interest = system_period_of_interest,
                                       outage_flag = true);

        nsamples = 100
        timeseries = unavailabilities(pras_system, nsamples)
    else
        timeseries = nothing
    end
    return timeseries
end


function unavailabilities(
    sys::SystemModel{N,L,T,P,E}, nsamples::Int
) where {N,L,T,P,E}

    resultspecs = (GeneratorAvailability(),
                StorageAvailability(), GeneratorStorageAvailability(),
                StorageEnergySamples(), GeneratorStorageEnergySamples())

    gen_av, stor_av, genstor_av, stor_soc, genstor_soc =
        PRAS.assess(sys, SequentialMonteCarlo(samples=nsamples), resultspecs...)

    result = Matrix{Float64}(undef, nsamples, N)

    for t in 1:N
        for s in 1:nsamples

            gen_unavailable =
                sum(sys.generators.capacity[.!gen_av.available[:, t, s], t])

            stor_unavailable = 0
            for st in 1:length(sys.storages)

                efficiency = sys.storages.discharge_efficiency[st,t]
                capacitylimit = sys.storages.discharge_capacity[st,t]
                energylimit = if t > 1
                    energytopower(stor_soc.energy[st,t-1,s], E, L, T, P) * efficiency
                else
                    0
                end

                if !stor_av.available[st, t, s]
                    stor_unavailable += capacitylimit
                elseif energylimit < capacitylimit
                    stor_unavailable += capacitylimit - energylimit
                end

            end

            genstor_unavailable = 0
            for gs in 1:length(sys.generatorstorages)

                efficiency = sys.generatorstorages.discharge_efficiency[gs,t]
                capacitylimit = sys.generatorstorages.gridinjection_capacity[gs,t]
                inflow = sys.generatorstorages.inflow[gs,t]
                energylimit = if t > 1
                    energytopower(genstor_soc.energy[gs,t-1,s], E, L, T, P) * efficiency
                else
                    0
                end
                capacityavailable = energylimit + inflow

                if !stor_av.available[gs, t, s]
                    genstor_unavailable += capacitylimit
                elseif capacityavailable < capacitylimit
                    genstor_unavailable += capacitylimit - capacityavailable
                end

            end

            result[s, t] = gen_unavailable + stor_unavailable + genstor_unavailable

        end
    end

    return result
end
