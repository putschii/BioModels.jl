"""
wrightfisher(Input::Array,Fittest, Fitness, Mutationrate)

Wright-Fisher-Model with mutation. Calculates and returns the population until fixation for 	  each time step. Fitness of other sequences is calculated based on differences to fittest  	sequence. Works with an input of strings or with DNA sequences from BioSequences. The function returns and array with the amout of a sequence for each time step and a dictionary.

       # Examples
       ```jldoctest
       julia> wrightfisher([dna"GGG",dna"AAA",dna"GAG"],dna"GGG",0.5,0.01)

       ```
       """
function wrightfisher(Input::Array,Fittest, Fitness, Mutationrate, Steps, dynamics::Bool, File)


### Error for wrong input

    if typeof(Fittest) != BioSequence{DNAAlphabet{4}}
        throw("Fittest has to be a sequence")
    end


    if (Fitness == NaN) || (typeof(Fitness) != Float64) || (typeof(Fitness) != Int64) && (typeof(Fitness) != Float64)
        throw("Fitness has to be a number")
    end


    if Mutationrate == NaN || typeof(Mutationrate) != Float64 || typeof(Mutationrate) != Int64 && typeof(Mutationrate) != Float64
        throw("Mutationrate has to be a number")
    end

    if typeof(dynamics) != Bool
        throw("Dynamics has to be a Bool")
    end

    if typeof(Steps) != Int64 && typeof(Steps) != Int32
        throw("Steps has to be a Int")
    end

    for i in 1:length(Input)
        if typeof(Input[i]) != BioSequence{DNAAlphabet{4}}
            throw("Input must contain an Array with BioSequence{DNAAlphabet{4}}")
        end
    end

    if typeof(File) != String
        throw("File has to be a String")
    end

### First Part

    ## give starting pop
    pop = Input


    ## starting selection coefficent based on last value of input array
    sel = Float64[]


    # fitnesscounter
    fitnesscounter = 0

    ## create empty variable for new pop
    newpop = []

    # just a counter for loops
    counter = 1

    # adding base for mutation
    basesA = [DNA_T,DNA_G,DNA_C]
    basesG = [DNA_A,DNA_T,DNA_C]
    basesT = [DNA_A,DNA_G,DNA_C]
    basesC = [DNA_A,DNA_T,DNA_G]


    # mutate yes or no with weights
    mutate = ["Yes", "No"]
    mutation = Weights([Mutationrate, 1 - Mutationrate])

    # arrays for potential mutants
    splitmut = []
    mutatepop = []

    # number of mutations in simulation
    mutationsteps = 0

    # timestep counter
    timesteps = 1

    # savecounter
    savecounter = 1
    filecounter = 1

    # counter for numbers in every loop
    actualnumbers = Dict()

    # counter for plot
    plotcount = Dict()

    ### Second Part ###


    # calculate selection for elements in array Input
    for o in 1:length(pop)
        if pop[o] == Fittest
           append!(sel, Fitness)
        else
            for j in 1:length(pop[o])
                if pop[o][j] == Fittest[j]
                    fitnesscounter = fitnesscounter +1

                end
            end
            append!(sel, Fitness *(fitnesscounter / length(pop[o])))
            fitnesscounter = 0
        end
    end

    # give weight to selection coefficent
    sel = Weights(sel)

    # count actual numbers in pop
    for o in 1:length(pop)
        actualnumbers[pop[o]] = (count(isequal(pop[o]),pop))
    end

    # adding values to plotcounter
    for o in 1:length(pop)
    	plotcount[pop[o]] = [(count(isequal(pop[o]),pop))]
    end


    ### Third part


    # loop for sampling until pop is one species
    while actualnumbers[pop[1]][1] != length(pop) || Steps == timesteps

    # create new pop (same size) with samples from old pop with mutation
        while counter <= length(pop)
            push!(mutatepop,[sample(pop,sel)])

            for i in mutatepop
                for j in 1:length(i)
                    if i[j] == DNA_A
                        if sample(mutate, mutation) == "Yes"
                            i[j] = sample(basesA)
                            mutationsteps = mutationsteps + 1

                        end
                    elseif i[j] == DNA_G
                        if sample(mutate, mutation) == "Yes"
                            i[j] = sample(basesG)
                            mutationsteps = mutationsteps + 1

                        end
                    elseif i[j] == DNA_C
                        if sample(mutate, mutation) == "Yes"
                            i[j] = sample(basesC)
                            mutationsteps = mutationsteps + 1

                        end
                    elseif i[j] == DNA_T
                        if sample(mutate, mutation) == "Yes"
                            i[j] = sample(basesT)
                            mutationsteps = mutationsteps + 1

                        end
                    end
                end
            end
            for i in mutatepop
                push!(newpop, i)
            end
            mutatepop = []
            counter = counter +1
        end

    # counter reset
        counter = 1

    # clear selection coefficent
        sel = Float64[]

    # turn new pop into old pop

        pop = newpop

    # find new selection coefficent
        for o in 1:length(pop)
            if pop[o] == Fittest
               push!(sel, Fitness)
            else
                for j in 1:length(pop[o])
                    if pop[o][j] == Fittest[j]
                        fitnesscounter = fitnesscounter +1

                    end
                end
                push!(sel, Fitness * (fitnesscounter / length(pop[o])))
                fitnesscounter = 0
            end
        end

    # clear new pop
        newpop = []

    # turn new selection coefficent into old one and give weights
        sel = Weights(sel)


    # clear actual numbers
        actualnumbers = Dict()

    # calculate new actual numbers
        for o in 1:length(pop)
            actualnumbers[pop[o]] = (count(isequal(pop[o]),pop))
        end

     # adding mutant to keys if not already there
        for o in pop
            if o in keys(plotcount)

            else
                if timesteps == 1
                    plotcount[o] = [0]
                else
                    plotcount[o] = [0]
                    while counter < timesteps
                        push!(plotcount[o], 0)
                        counter = counter + 1
                    end

                end
            end
        end

    # increase timesteps
        timesteps = timesteps + 1
        counter = 1




    # put values of keys in plotcounter, if key is not there than put 0 in
        for o in keys(plotcount)
            if o in keys(actualnumbers)
                push!(plotcount[o],actualnumbers[o][1])
            else
                push!(plotcount[o], 0)

            end

        end

    # increase savecounter
        savecounter = savecounter +1

    # Save and clear every 50000 steps or clear Dict every 1000 steps
        if dynamics == true
            if savecounter == 50000
                f = jldopen("$File", "a+")
                write(f, "$filecounter", plotcount)
                close(f)
                savecounter = 1
                filecounter = filecounter + 1
                plotcount = Dict()
                for i in 1:length(pop)
                    plotcount[pop[i]] = [(count(isequal(pop[i]),pop))]
                end
            end
        else
            if savecounter == 1000
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
    for o in keys(plotcount)
        push!(plotdata, plotcount[o] )
    end

    return [plotdata, plotcount,timesteps, pop[1]]
end
