"""
simplequasi(Seq::Array, Numbers::Array, Fittness::Array, Fmutation::Array, Bmutation::Array, Steps)

Model to explore the impact of fitness, forward and backward mutations on the size of sequences with different fitness. Returns a dictionary with the numbers of the sequences for each time step. Note that Fmutation has to be greater than Bmutation, because Fmutation contains Bmutation.

       # Examples
       ```jldoctest
       julia> simplequasi(["A","B"], [99,1], [1.5,1], [0.7,0.6],[0,0], 5)
       ```
	Dict{Any,Any} with 2 entries:
  		"B" => Any[1, 105.35, 298.217, 636.059, 1207.39, 2149.85]
  		"A" => Any[99, 143.55, 208.148, 301.814, 437.63, 634.564]
"""
function simplequasi(Seq::Array, Numbers::Array, Fittness::Array, Fmutation::Array, Bmutation::Array, Steps)


    # Error if Bmutation > Fmutation
    for i in 1:length(Seq)
        if Fmutation[i] < Bmutation[i]
            throw("Fmutation has to be bigger than Bmutation")
        end
    end

    # Dict with Genotype as key, fitness and mutation rates as values
    quasidict = Dict()

    for i in 1:length(Seq)
        quasidict[Seq[i]] = [Numbers[i],Fittness[i],Fmutation[i],Bmutation[i]]
    end

    # Dict with actual numbers

    actualnumbers = Dict()

    for i in 1:length(Seq)
        actualnumbers[Seq[i]] = Any[Numbers[i]]
    end

    # Array for new numbers
    changes = []


    for i in 1:length(Seq)
        append!(changes, quasidict[Seq[i]][1])
    end

    # just a counter
    counter = 1

    # Loop x times
    while counter <= Steps

        # loop for each seq
        for i in 1:length(Seq)

            # case of fittest sequence (no back mutaion)
            if i == 1
                # increase seq 1
                changes[i] = changes[i] + last(actualnumbers[Seq[i]])* quasidict[Seq[i]][2]*(1-quasidict[Seq[i]][3])

                # increase seq 1+1 due to mutation
                changes[i+1] = changes[i+1] + last(actualnumbers[Seq[i]])* quasidict[Seq[i]][2]*quasidict[Seq[i]][3]

            # case of weakest seq (no forward mutation)
            elseif i == length(Seq)
                #increase in last seq
                changes[i] = changes[i] + last(actualnumbers[Seq[i]])* quasidict[Seq[i]][2]*(1-quasidict[Seq[i]][3])

                # increase in fitter seq
                changes[i-1] = changes[i-1] + last(actualnumbers[Seq[i]])* quasidict[Seq[i]][2]*quasidict[Seq[i]][4]
            # 1 < Seq < length(Seq)
            else
                #increase in last seq
                changes[i] = changes[i] + last(actualnumbers[Seq[i]])* quasidict[Seq[i]][2]*(1-quasidict[Seq[i]][3])

                # increase in fitter seq
                changes[i-1] = changes[i-1] + last(actualnumbers[Seq[i]])* quasidict[Seq[i]][2]*quasidict[Seq[i]][4]

                # increase in weaker seq
                changes[i+1] = changes[i+1] + last(actualnumbers[Seq[i]])* quasidict[Seq[i]][2]*(quasidict[Seq[i]][3]-quasidict[Seq[i]][4])
            end
        end


        # add new numbers to dict
        for i in 1:length(Seq)
            push!(actualnumbers[Seq[i]], changes[i])
        end

        # increase counter
        counter = counter + 1
    end

    # return numbers for each seq
    return(actualnumbers)
end
