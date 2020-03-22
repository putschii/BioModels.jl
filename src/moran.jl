"""
moran(Input::Array,Fittest, Fitness, Mutationrate)

Moran-Model with mutation. Calculates and returns the population until fixation for each  	time step. Fitness of other sequences is calculated based on differences to fittest sequence. Works with an input of strings or with DNA sequences from BioSequences. The function returns and array with the amout of a sequence for each time step and a dictionary. 

       # Examples
       ```jldoctest
       julia> moran([dna"GGG",dna"AAA",dna"GAG"],dna"GGG",0.5,0.01)
       julia> moran(["GGG","AAA","GAG"],"GGG",0.5,0.01)
       ```
       """
function moran(Input::Array,Fittest, Fitness, Mutationrate)

    # apply needed functions
    function remove!(arr, item)
        deleteat!(arr, findall(x->x==item, arr))
    end          
### First Part
    
    ## give starting pop
    pop = String[]
    if typeof(Input[1]) == String
        pop = String[Input]
    else        
        for o in 1:length(Input)
            push!(pop, Input[o])
        end
    end
    
    if typeof(Fittest) == LongSequence{DNAAlphabet{4}}
        Fittest = convert(String, Fittest)
    end
    
    negativpop = String[]
    
    ## starting selection coefficent based on last value of input array
    sel = Float64[]
    negativsel = Float64[]
    
    # fitnesscounter    
    fitnesscounter = 0

    # just a counter for loops
    counter = 1
    
    # adding base for mutation
    basesA = ["T","G","C"]
    basesG = ["A","T","C"]
    basesT = ["A","G","C"]
    basesC = ["A","T","G"]

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

    
    # counter for numbers in every loop
    actualnumbers = Dict()
    
    # counter for plot
    plotcount = Dict()

    #creating array for final sequences
    final = []


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
    while actualnumbers[pop[1]][1] != length(pop)
    
    # choosing one individual out of the pop by negativ selection coefficent
        append!(negativpop, [sample(pop,negativsel)])    
    
    # choose one from pop to reproduce with mutation
        append!(mutatepop,[sample(pop,sel)])
        for o in mutatepop
            append!(splitmut, o)
        end
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
        splitmut = [join(splitmut)]    
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
            
    # add Mutant to pop
        push!(pop, join(splitmut))

    # clear arrays
        splitmut = []
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
    println("Fixation after ", timesteps, " timesteps" )
    println("Mutationsteps: ", mutationsteps)
    return [plotdata, plotcount]
end
