"""
moran(Input::Array,Fittest, Mutationrate, Steps, Save::Bool, File)

The function works as a Moran-Model with mutation and the mutation rate can be defined by the user. The model will run until it reached fixation, i.e. until the whole population is based on only one DNA sequence for each individual or until it reached the final amount of time steps defined by the user with the variable Steps.
It calculates and returns the population until fixation for each time step. If save is set to true, every time step will be saved in a jld2 File named by the user. It is important that the file is a string and ends with .jld2. The file will be stored in the same folder than the julia file. If Save is set to false, the function returns the population of the last time step and the needed time steps. For further use of the data, it is recommended to set Save to true. The fitness of each sequences is calculated based on differences in bases to fittest sequence defined by the user. The model works with an input of DNA sequences in an array constructed with BioSequences.

       # Examples
       ```jldoctest
       julia> moran([dna"GGG",dna"AAA",dna"GAG"],dna"GGG",0.0001,500,false,"filenamehere.jld2")
       ```
       """
function moran(Input::Array,Fittest, Mutationrate, Steps, Save::Bool, File)


### Error for wrong input




    if Mutationrate == NaN || typeof(Mutationrate) != Float64 || typeof(Mutationrate) != Int64 && typeof(Mutationrate) != Float64
        throw("Mutationrate has to be a number")
    end

    if typeof(Save) != Bool
       throw("Dynamics has to be a Bool")
    end

    if typeof(Steps) != Int64 && typeof(Steps) != Int32
        throw("Steps has to be a Int")
    end



    if typeof(File) != String
        throw("File has to be a String")
    end

### First Part

    ## give starting pop
    pop = Matrix(undef,length(Input),length(Input[1]))

    for i in 1:length(Input)
        for j in 1:length(Input[i])
            pop[i,j] = Input[i][j]
        end
    end

    # Size of Pop
    poplength = []

    for i in 1:size(pop)[1]
        push!(poplength,i)
    end


    # Turn Fittest to Matrix
    fit = Matrix(undef,1,length(Fittest))

    for i in 1:length(Fittest)
        fit[1,i] = Fittest[i]
    end


    # Matrix for selection
    selmatrix = Matrix(undef,length(Input),length(Input[1]))

    sel = Float64[]

    # matrix for fixation
    fixmatrix = Matrix(undef,length(Input),length(Input[1]))


    # fitnesscounter
    fitnesscounter = 0

    # just a counter for loops
    counter = 1

    # adding base for mutation
    basesA = [DNA_T,DNA_G,DNA_C]
    basesG = [DNA_A,DNA_T,DNA_C]
    basesT = [DNA_A,DNA_G,DNA_C]
    basesC = [DNA_A,DNA_T,DNA_G]


    # number of mutations in simulation
    mutationsteps = 0

    # timestep counter
    timesteps = 1

    # savecounter
    savecounter = 1
    filecounter = 1
    # save pop
    if Save == true
        f = jldopen("$File", "a+")
        write(f, "0", pop)
        close(f)
    end
    ### Second Part ###


    # calculate selection for elements in array Input

    for i in 1:size(pop)[1]
        for j in 1:length(pop[i,:])
            if pop[i,j] == fit[1,j]
                selmatrix[i,j] = 1
            else
                selmatrix[i,j] = 0
            end
        end
    end

    for i in 1:size(selmatrix)[1]
        push!(sel,sum(selmatrix[i,:]))
    end

    # give weight to selection coefficent
    sel = Weights(sel)

    # Fixationmatrix with same bases (1 for missmatch, 0 for match)


    for i in 1:size(pop)[2]
        for j in 1:(length(pop[:,i])-1)
            if pop[:,i][j] == pop[:,i][j+1]
                fixmatrix[j,i] = 0
            else
                fixmatrix[j,i] = 1
            end
            if j == (length(pop[:,i])-1)
                if pop[:,i][j] == pop[:,i][j+1]
                    fixmatrix[j+1,i] = 0
                else
                    fixmatrix[j+1,i] = 1
                end
            end
        end

    end

    finish = sum(fixmatrix)

    ### Third part


    # loop for sampling until pop is one species
    while finish != 0 && timesteps <= Steps


    # Populate
        timesteps = timesteps + 1
        sampler = Matrix(undef,1,length(pop[1,:]))
        keeper = sample(poplength,sel)
        negsel = sel *-1
        leaver = sample(poplength,sel)
        for j in 1:length(pop[keeper,:])
            sampler[1,j] = pop[keeper,j]
        end



    # Mutate
        for i in 1:length(sampler)

            if sampler[i] == DNA_A
                if rand() <= Mutationrate
                    sampler[i] = sample(basesA)
                    mutationsteps = mutationsteps + 1

                end
            elseif sampler[i] == DNA_G
                if rand() <= Mutationrate
                   sampler[i] = sample(basesG)
                    mutationsteps = mutationsteps + 1

                end
            elseif sampler[i] == DNA_C
                if rand() <= Mutationrate
                    pop[i] = sample(basesC)
                    mutationsteps = mutationsteps + 1

                end
            elseif sampler[i] == DNA_T
                if rand() <= Mutationrate
                    sampler[i] = sample(basesT)
                    mutationsteps = mutationsteps + 1
                end
            end
        end

    # remove leaver sequence and replace it with sample
        pop[leaver,:] = sampler[1,:]


    # Reset selection

        selmatrix = Matrix(undef,length(Input),length(Input[1]))
        sel = Float64[]
        fixmatrix = Matrix(undef,length(Input),length(Input[1]))

    # calculate selection
        for i in 1:size(pop)[1]
            for j in 1:length(pop[i,:])
                if pop[i,j] == fit[1,j]
                    selmatrix[i,j] = 1
                else
                    selmatrix[i,j] = 0
                end
            end
        end

        for i in 1:size(selmatrix)[1]
            append!(sel,sum(selmatrix[i,:]))
        end

    # give weight to selection coefficent
        sel = Weights(sel)

    # Fixationmatrix with same bases (1 for missmatch, 0 for match)


        for i in 1:size(pop)[2]
            for j in 1:(length(pop[:,i])-1)
                if pop[:,i][j] == pop[:,i][j+1]
                    fixmatrix[j,i] = 0
                else
                    fixmatrix[j,i] = 1
                end
                if j == (length(pop[:,i])-1)
                    if pop[:,i][j] == pop[:,i][j+1]
                        fixmatrix[j+1,i] = 0
                    else
                        fixmatrix[j+1,i] = 1
                    end
                end
            end
        end

    # Save Matrix for every timestep
        if Save == true
                f = jldopen("$File", "a+")
                write(f, "$timesteps", pop)
                close(f)
        end

    # increase timesteps

        finish = sum(fixmatrix)

    end

    if Save == true
        println("Data is saved in File")
        return(timesteps)
    else
        return [pop, timesteps]
    end
end
