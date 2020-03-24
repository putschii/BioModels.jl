"""
moran(Input::Array,Fittest, Fitness, Mutationrate)

Moran-Model with mutation. Calculates and returns the population until fixation for each  	time step. Fitness of other sequences is calculated based on differences to fittest sequence. Works with an input of strings or with DNA sequences from BioSequences. The function returns and array with the amout of a sequence for each time step and a dictionary.

       # Examples
       ```jldoctest
       julia> moran([dna"GGG",dna"AAA",dna"GAG"],dna"GGG",0.5,0.01)
       julia> moran(["GGG","AAA","GAG"],"GGG",0.5,0.01)
       ```
       """
function moran(Input::Array,Fittest, Fitness, Mutationrate, Steps, dynamics::Bool, File)

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

### apply needed functions
    function remove!(arr, item)
        deleteat!(arr, findall(x->x==item, arr))
    end

### First Part

    ## give starting pop
    pop = Input


    negativpop = []

    ## starting selection coefficent based on last value of input array
    sel = Float64[]
    negativsel = Float64[]

    # fitnesscounter
    fitnesscounter = 0

    # just a counter for loops
    counter = 1

    # adding base for mutation
    basesA = [dna"T",dna"G",dna"C"]
    basesG = [dna"A",dna"T",dna"C"]
    basesT = [dna"A",dna"G",dna"C"]
    basesC = [dna"A",dna"T",dna"G"]

    # mutate yes or no with weights
    mutate = ["Yes", "No"]
    mutation = Weights([Mutationrate, 1 - Mutationrate])

    # arrays for potential mutants
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
            append!(sel, Fitness * (fitnesscounter / length(pop[o])))
            fitnesscounter = 0
        end
    end

    for o in sel
        append!(negativsel, 1 - o)
    end


    # give weight to selection coefficent
    sel = Weights(sel)
    negativsel = Weights(negativsel)


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

    # choosing one individual out of the pop by negativ selection coefficent
        append!(negativpop, [sample(pop,negativsel)])

    # choose one from pop to reproduce with mutation
        append!(mutatepop,[sample(pop,sel)])

        for i in mutatepop
            for j in 1:length(i)
                if i[j] == dna"A"
                    if sample(mutate, mutation) == "Yes"
                        i[j] = sample(basesA)
                        mutationsteps = mutationsteps + 1

                    end
                elseif i[j] == dna"G"
                    if sample(mutate, mutation) == "Yes"
                        i[j] = sample(basesG)
                        mutationsteps = mutationsteps + 1

                    end
                elseif i[j] == dna"C"
                    if sample(mutate, mutation) == "Yes"
                        i[j] = sample(basesC)
                        mutationsteps = mutationsteps + 1

                    end
                elseif i[j] == dna"T"
                    if sample(mutate, mutation) == "Yes"
                        i[j] = sample(basesT)
                        mutationsteps = mutationsteps + 1

                    end
                end
            end
        end

    # marking the individual of negativpop in pop as N
        while counter <= length(pop)
            if pop[counter] == negativpop[1]
                pop[counter] = dna"N"
                break
            else
                counter = counter + 1

            end
        end
        counter = 1

    # remove the marked individual
        remove!(pop,dna"N")

    # clear negativpop
        negativpop = []

    # add Mutant to pop
        push!(pop, mutatepop[1])

    # clear arrays

        mutatepop = []


    # clear selection coefficent
        sel = Float64[]
        negativsel = Float64[]

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

        for o in sel
            append!(negativsel, 1 - o)
        end

    # turn new selection coefficent into old one and give weights
        sel = Weights(sel)
        negativsel = Weights(negativsel)

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

    # save and clear Dict every 50000 steps
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
        # clear every 1000 steps without saving
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

    return [plotdata, plotcount, timesteps, pop[1]]
end
