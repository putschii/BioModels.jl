"""
neutral(Input::Array,Meta::Array, mutationrate)

Neutral-Model with mutation, but without selection. Calculates and returns the population  	 until fixation for each time step based on migration from metapopulation. No changes in the metapopulation. Works with an input of strings or with DNA sequences from BioSequences. The function returns and array with the amout of a sequence for each time step and a dictionary.     

       # Examples
       ```jldoctest
       julia> neutral([dna"GGG",dna"AAA",dna"GAG"],[dna"GGA",dna"GCA",dna"ACA"],0.01)
       
       ```
       """
function neutral(Input::Array,Meta::Array, mutationrate) 

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
        for o in 1:length(Input)
            push!(pop, Input[o])
        end
    end

    #give starting metapop
    metapop = String[]
    if typeof(Meta[1]) == String
        metapop = String[Meta]
    else        
        for o in 1:length(Meta)
            push!(metapop, Meta[o])
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

    ### Second Part
    

    # adding values to counter
    for o in 1:length(pop)
    	actualnumbers[pop[o]] = (count(isequal(pop[o]),pop))
    end

    
    # adding keys of pop with values of input
    for o in 1:length(pop)
        plotcount[pop[o]] = [(count(isequal(pop[o]),pop))]
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
        append!(combinedpop,pop)
        append!(combinedpop,metapop)
        
    # append mutatepop with a new individual from pop or metapop
        append!(mutatepop, [sample(combinedpop)])
    
    # splitting mutant into bases
        for o in mutatepop
            append!(splitmut, o)
        end
        
    # check for bases in mutant to mutate
        for o in 1:length(splitmut)
            if splitmut[o] == 'A'
                if sample(mutate, mutation) == "Yes"
                    splitmut[o] = sample(basesA)
                    mutationsteps = mutationsteps + 1
                    
                end
            elseif splitmut[o] == 'G'
                if sample(mutate, mutation) == "Yes"
                    splitmut[o] = sample(basesG)
                    mutationsteps = mutationsteps + 1
                    
                end
            elseif splitmut[o] == 'C'
                if sample(mutate, mutation) == "Yes"
                    splitmut[o] = sample(basesC)
                    mutationsteps = mutationsteps + 1
                    
                end
            elseif splitmut[o] == 'T'
                if sample(mutate, mutation) == "Yes"
                    splitmut[o] = sample(basesT)
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
        for o in splitmut
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
        
    # reset counter
        counter = 1
        
    # calculate new actual numbers   
        for o in 1:length(pop)
            actualnumbers[pop[o]] = (count(isequal(pop[o]),pop))
        end
        
    # put values of keys in plotcounter, if key is not in pop than put 0 in
        for o in keys(plotcount)
            if o in keys(actualnumbers)              
                push!(plotcount[o],actualnumbers[o][1])      
            else
                push!(plotcount[o], 0)

            end

        end

    
    # clear mutantpop
        mutatepop = []
        splitmut = []
        
    # clear combinedpop
        combinedpop = []
        
    end
    
    # create a array with data from plotcounter for plotting
    plotdata = []
    for o in keys(plotcount)
        push!(plotdata, plotcount[o] )
    end

    # creating final sequence array as DNA
    for o in 1:length(pop)
    	push!(final, LongDNASeq(pop[o]))
    end
    
    # Output
    println("Fixation after ", timesteps, " timesteps" )
    println("Mutationsteps: ", mutationsteps)
    return [plotdata, plotcount]
end
