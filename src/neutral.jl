"""
neutral(Input::Array,Meta::Array, mutationrate, migration::Bool, save::Bool, file)

Neutral-Model with mutation, but without selection. Calculates and returns the population until fixation for each time step based on migration from metapopulation. No changes in the metapopulation. If migration is set to true, then only sequences from the meta community will be used for replacement. If it is set to false, the population and meta community can be used for replacement.  Works with an input of strings or with DNA sequences from BioSequences. The function returns and array with the amout of a sequence for each time step and a dictionary. The data can be stored in a file if save is set to true. Therefore every time step of the first 10000 will be stored and after that only every 10000 step. If save is set to "false" only the last 10000 steps will be returned. File will be stored in same folder as the julia code.

       # Examples
       ```jldoctest
       julia> neutral([dna"GGG",dna"AAA",dna"GAG"],[dna"GGA",dna"GCA",dna"ACA"],0.001, true, true, "yourfilenamehere.jld2")

       ```
       """
function neutral(Input::Array,Meta::Array, mutationrate, migration::Bool, save::Bool, file)

    # apply needed functions
    function remove!(arr, item)
        deleteat!(arr, findall(x->x==item, arr))
    end


    ### First Part

    # give starting pop
    pop = String[]
    if typeof(Input[1]) == String
        pop = String[Input]
    else
        for i in 1:length(Input)
            push!(pop, Input[i])
        end
    end

    #give starting metapop
    metapop = String[]
    if typeof(Meta[1]) == String
        metapop = String[Meta]
    else
        for i in 1:length(Meta)
            push!(metapop, Meta[i])
        end
    end

    # negativpop to clear one individual out
    negativpop = String[]

    # just a counter for loops
    counter = 1

    # pool of combined meta and pop to choose from
    combinedpop = String[]

    # adding base for mutation
    basesA = ["T","G","C"]
    basesG = ["A","T","C"]
    basesT = ["A","G","C"]
    basesC = ["A","T","G"]

    # mutate yes or no with weights
    mutate = ["Yes", "No"]
    mutation = Weights([mutationrate, 1 - mutationrate])
    mutationsteps = 0

    # pop with potential mutant
    mutatepop = String[]

    # splitted mutant sequence for change of bases
    splitmut = []

    # time steps for adding 0s to later mutants in keys
    timesteps = 1

    # counter for numbers in every loop
    actualnumbers = Dict()

    # counter for overall numbers
    plotcount = Dict()

    #creating array for final sequences
    final = []

    # savecounter
    savecounter = 1
    filecounter = 1

    ### Second Part


    # adding values to counter
    for i in 1:length(pop)
    	actualnumbers[pop[i]] = (count(isequal(pop[i]),pop))
    end


    # adding keys of pop with values of input
    for i in 1:length(pop)
        plotcount[pop[i]] = [(count(isequal(pop[i]),pop))]
    end


    ### Third part


    # loop for sampling until pop is one species
    while actualnumbers[pop[1]][1]!= length(pop)

    # choosing one individual out of the pop by negativ selection coefficent
        append!(negativpop, [sample(pop)])

    # marking the individual of negativpop in pop as "x"
        while counter <= length(pop)
            if pop[counter] == negativpop[1]
                pop[counter] = "x"
                break
            else
                counter = counter + 1

            end
        end
        counter = 1

    # remove the marked individual
        remove!(pop,"x")

    # clear negativpop
        negativpop = []

    # combine pop and metapop in combinedpop
        if migration == true
            combinedpop = metapop
        end
        if migration == false
            append!(combinedpop,pop)
            append!(combinedpop,metapop)
        end

    # append mutatepop with a new individual from pop or metapop
        append!(mutatepop, [sample(combinedpop)])
        #println(mutatepop)
    # splitting mutant into bases
        for i in mutatepop
            append!(splitmut, i)
        end

    # check for bases in mutant to mutate
        for i in 1:length(splitmut)
            if splitmut[i] == 'A'
                if sample(mutate, mutation) == "Yes"
                    splitmut[i] = sample(basesA)
                    mutationsteps = mutationsteps + 1

                end
            elseif splitmut[i] == 'G'
                if sample(mutate, mutation) == "Yes"
                    splitmut[i] = sample(basesG)
                    mutationsteps = mutationsteps + 1

                end
            elseif splitmut[i] == 'C'
                if sample(mutate, mutation) == "Yes"
                    splitmut[i] = sample(basesC)
                    mutationsteps = mutationsteps + 1

                end
            elseif splitmut[i] == 'T'
                if sample(mutate, mutation) == "Yes"
                    splitmut[i] = sample(basesT)
                    mutationsteps = mutationsteps + 1

                end
            end
        end



     # putting bases together
        splitmut = [join(splitmut)]

     # adding mutant to pop
        push!(pop, join(splitmut))

    # clear actual numbers
        actualnumbers = Dict()

    # adding mutant to keys if not already there
        for i in splitmut
            if i in keys(plotcount)

            else
                if timesteps == 1
                    plotcount[i] = [0]
                else
                    plotcount[i] = [0]
                    while counter < timesteps
                        push!(plotcount[i], 0)
                        counter = counter + 1
                    end

                end
            end
        end

    # increase timesteps
        timesteps = timesteps + 1

    # reset counter
        counter = 1

    # calculate new actual numbers
        for i in 1:length(pop)
            actualnumbers[pop[i]] = (count(isequal(pop[i]),pop))
        end

    # put values of keys in plotcounter, if key is not in pop than put 0 in
        for i in keys(plotcount)
            if i in keys(actualnumbers)
                push!(plotcount[i],actualnumbers[i][1])
            else
                push!(plotcount[i], 0)

            end

        end


    # clear mutantpop
        mutatepop = []
        splitmut = []

    # clear combinedpop
        combinedpop = []

    # increase savecounter
        savecounter = savecounter +1

    # save and clear Dict every 10000 steps
        if save == true
            if savecounter == 10000 && filecounter == 1
                f = jldopen("$File", "a+")
                write(f, "$filecounter", plotcount)
                close(f)
                savecounter = 1
                filecounter = filecounter + 1
                plotcount = Dict()
                for i in 1:length(pop)
                    plotcount[pop[i]] = [(count(isequal(pop[i]),pop))]
                end
            elseif savecounter == 10000 && filecounter > 1
                plotcount = Dict()
                for i in 1:length(pop)
                    plotcount[pop[i]] = [(count(isequal(pop[i]),pop))]
                end
                f = jldopen("$File", "a+")
                write(f, "$filecounter", plotcount)
                close(f)
                savecounter = 1
                filecounter = filecounter + 1
            end
        else
            if savecounter == 10000
                savecounter = 1
                filecounter = filecounter + 1
                plotcount = Dict()
                for i in 1:length(pop)
                    plotcount[pop[i]] = [(count(isequal(pop[i]),pop))]
                end
            end
        end


    end

    # create a array with data from plotcounter for plotting
    plotdata = []
    for i in keys(plotcount)
        push!(plotdata, plotcount[i] )
    end

    # creating final sequence array as DNA
    for i in 1:length(pop)
    	push!(final, LongDNASeq(pop[i]))
    end

    # Output

    return [plotdata, plotcount, timesteps, pop[1]]
end
