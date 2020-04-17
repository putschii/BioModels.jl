"""
neutral(Input::Array,Meta::Array, Mutationrate, Steps, Migration::Bool, dynamics::Bool, File)

The function works as a Neutral-Model with mutation and the mutation rate can be defined by the user. The model will run until it reached fixation, i.e. until the whole population is based on only one DNA sequence for each individual or until it reached the final amount of time steps defined by the user with the variable Steps.
If save is set to true, every time step will be saved in a jld2 File named by the user. It is important that the file is a string and ends with .jld2. The file will be stored in the same folder than the julia file. If Save is set to false, the function returns the population of the last time step and the needed time steps. For further use of the data, it is recommended to set Save to true. The fitness of each sequences is calculated based on differences in bases to fittest sequence defined by the user. The model works with an input of DNA sequences in an array constructed with BioSequences.

       # Examples
       ```jldoctest
       julia> neutral([dna"GGG",dna"AAA",dna"GAG"],[dna"GGA",dna"GCA",dna"ACA"],0.0001, 500,true, true, "yourfilenamehere.jld2")

       ```
       """
function neutral(Input::Array,Meta::Array,Mutationrate,Steps,Migonly::Bool,Save::Bool,File)


### Error for wrong input



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

    # Meta community
    met = Matrix(undef,length(Meta),length(Meta[1]))

    for i in 1:length(Meta)
        for j in 1:length(Meta[i])
            met[i,j] = Meta[i][j]
        end
    end

    # Size of Pop
    metalength = []

    for i in 1:size(met)[1]
        push!(metalength,i)
    end

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
    timesteps = 0

    # save pop
    if Save == true
        f = jldopen("$File", "a+")
        write(f, "0", pop)
        close(f)
    end

    ### Second Part ###




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
    while finish != 0 && timesteps < Steps


    # Populate
        # increase timesteps
        timesteps = timesteps + 1
        sampler = Matrix(undef,1,length(pop[1,:]))
        if Migonly == true
            keeper = sample(metalength)
            leaver = sample(poplength)
            for j in 1:length(met[keeper,:])
                sampler[1,j] = met[keeper,j]
            end


        elseif Migonly == false
            combipop = vcat(pop,met)
            combilength = []
            for i in 1:size(combipop)[1]
                push!(combilength,i)
            end
            keeper = sample(combilength)
            leaver = sample(poplength)
            for j in 1:length(combipop[keeper,:])
                sampler[1,j] = combipop[keeper,j]
            end

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
                    sampler[i] = sample(basesC)
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

        finish = sum(fixmatrix)

    end

    if Save == true
        println("Data is saved in File")
        return(timesteps)
    else
        return [pop, timesteps]
    end
end
