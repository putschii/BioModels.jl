"""
wrightfisher(Input::Array,Fittest, Mutationrate, Steps, Save::Bool, File)


The function works as a Wright-Fisher-Model with mutation and the mutation rate can be defined by the user. The model will run until it reached fixation, i.e. until the whole population is based on only one DNA sequence for each individual or until it reached the final amount of time steps defined by the user with the variable Steps.
It calculates and returns the population until fixation for each time step. If save is set to true, every time step will be saved in a jld2 File named by the user. It is important that the file is a string and ends with .jld2. The file will be stored in the same folder than the julia file. If Save is set to false, the function returns the population of the last time step and the needed time steps. For further use of the data, it is recommended to set Save to true. The fitness of each sequences is calculated based on differences in bases to fittest sequence defined by the user. The model works with an input of DNA sequences in an array constructed with BioSequences.

       # Examples
       ```jldoctest
       julia> wrightfisher([dna"GGG",dna"AAA",dna"GAG"],dna"GGG",0.0001,500,false,"filenamehere.jld2")

       ```
       """
function wrightfisher(Input::Array,Fittest, Mutationrate, Steps, Save::Bool, File)

    
### Error for wrong input
    
    #if typeof(Fittest) != BioSequence{DNAAlphabet{4}}
    #    throw("Fittest has to be a sequence")
   # end

    
        
    if Mutationrate == NaN || typeof(Mutationrate) != Float64 || typeof(Mutationrate) != Int64 && typeof(Mutationrate) != Float64
        throw("Mutationrate has to be a float")
    end
    
    if typeof(Save) != Bool
       throw("Dynamics has to be a Bool")
    end
    
    if typeof(Steps) != Int64 && typeof(Steps) != Int32
        throw("Steps has to be a Int")
    end 
        
  #  for i in 1:length(Input)
  #      if typeof(Input[i]) != BioSequence{DNAAlphabet{4}}
   #         throw("Input must contain an Array with BioSequence{DNAAlphabet{4}}")
  #      end
  #  end
        
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
    timesteps = 0
    
    # savecounter
    savecounter = 1
    filecounter = 1

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
        sizecounter = size(pop)[1]
        while counter <= sizecounter
            sampler = Matrix(undef,1,length(pop[1,:]))
            c = sample(poplength,sel)
            for j in 1:length(pop[c,:])
                sampler[1,j] = pop[c,j]
            end
            pop = vcat(pop,sampler)
            counter = counter + 1
        end
        
    # Counter reset
        counter = 1 
        
    # Remove old Population
        pop = pop[(length(Input)+1):size(pop)[1],:]

    # Mutate
        for i in 1:length(pop)
            
            if pop[i] == DNA_A
                if rand() <= Mutationrate
                    pop[i] = sample(basesA)
                    mutationsteps = mutationsteps + 1
                        
                end
            elseif pop[i] == DNA_G
                if rand() <= Mutationrate
                    pop[i] = sample(basesG)
                    mutationsteps = mutationsteps + 1
                        
                end
            elseif pop[i] == DNA_C
                if rand() <= Mutationrate
                    pop[i] = sample(basesC)
                    mutationsteps = mutationsteps + 1
                        
                end
            elseif pop[i] == DNA_T
                if rand() <= Mutationrate
                    pop[i] = sample(basesT)
                    mutationsteps = mutationsteps + 1
                end
            end
        end

  ##########  # find new selection coefficent
        
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
    

        finish = sum(fixmatrix)
    # Save Matrix for every timestep
        if Save == true
                f = jldopen("$File", "a+")
                write(f, "$timesteps", pop)
                close(f)
        end

    end
    
    if Save == true
        println("Data is saved in file")
    else
        return [pop, timesteps]
    end
end
