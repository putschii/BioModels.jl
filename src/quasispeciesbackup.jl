"""
quasispecies(Input::Array,Fittest, Fitness, Mutationrate, Deathrate, Steps, Rep)

Quasispecies-Model with mutation and selection. Calculates and returns the population until defined steps are reached. The amount out reproduction hast to be set at the beginning. The function returns and array with the final sequences at the last time step, and array with the amout of a sequence for each time step and a dictionary.     

       # Examples
       ```jldoctest
       julia> quasispecies([dna"GGG",dna"AAA",dna"GAG"],dna"GGG",0.0001,0.3,20,3)
       
       ```
       """
function quasispecies(Input::Array, Fittest, Fitness, Mutationrate, Deathrate, Steps, Rep)

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
   
    ## convert Fittest into String
    if typeof(Fittest) == LongSequence{DNAAlphabet{4}}
        Fittest = convert(String, Fittest)
    end    
    
    # newpop for offsprings of pop
    newpop = []
    
    # selection coefficient
    sel = Float64[]
    
    # arrays for reproduction
    reppop = []
    repsel = Float64[]
    samplepop = []

    
    # start for n
    n = 1

    # just a counter for loops
    counter = 1
    
    # adding base for mutation
    basesA = ["T","G","C"]
    basesG = ["A","T","C"]
    basesT = ["A","G","C"]
    basesC = ["A","T","G"]
    
    # mutate yes or no with weights
    mutate = ["Yes", "No"]
    mutation = Weights([Mutationrate, 1-Mutationrate])
    
    # fitnesscounter    
    fitnesscounter = 0    
    
    # arrays for potential mutants
    splitmut = []
    mutatepop = []
    
    # number of mutations in simulation
    mutationsteps = 0
    
    # timestep counter
    timesteps = 1
    
    #deathrate
    death = Weights([Deathrate, 1-Deathrate])
    
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
    while n <= Steps
    
    # create new pop (same size) with samples from old pop with mutation
        


        for i in 1:length(sel)
            if sel[i] >= mean(sel)
                push!(reppop, pop[i])
                push!(repsel, sel[i])
            end
        end
        repsel = Weights(repsel)
        while counter <= Rep
            push!(samplepop, [sample(reppop,repsel)])
            counter = counter + 1
        end
        counter = 1
        
        for t in samplepop
            append!(splitmut,t)
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
            push!(newpop, join(splitmut))
            splitmut = []
            mutatepop = []     
            
        end


    # clear actual numbers   
        actualnumbers = Dict()   
    
     # adding mutant to keys if not already there
        for o in newpop
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
                counter = 1
            end
        end
        
    # increase timesteps
        timesteps = timesteps + 1
        
    # death in pop
        for i in 1:length(pop)
            if sample(mutate,death) == "Yes"
                    pop[i] = "x"
            end
        end
        remove!(pop, "x")
  
        
        
    # put newpop in pop
        append!(pop, newpop)
           
    # calculate selection for elements in array Input
        sel = Float64[]
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

      
    
    # calculate new actual numbers   
        for o in 1:length(pop)
            actualnumbers[pop[o]] = (count(isequal(pop[o]),pop))
        end       
    
    #clear newpop
        newpop= []        
        repop = []
        samplepop = []
        repsel = Float64[]
    # put values of keys in plotcounter, if key is not there than put 0 in
        for o in keys(plotcount)
            if o in keys(actualnumbers)              
                push!(plotcount[o],actualnumbers[o][1])      
            else
                push!(plotcount[o], 0)

            end

        end

    # increase n by 1
        n = n+1
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
    
    # output
    println("Mutationsteps: ", mutationsteps)
    return [final, plotdata, plotcount]
end
