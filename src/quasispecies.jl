"""
quasispecies(Input::Array, Fittest, Fitness, Mutationrate, Steps)

Quasispecies-Model with mutation and selection. Calculates and returns the population until defined steps are reached. The function returns and array with the final sequences at the last time step, and array with the amout of a sequence for each time step and a dictionary.     

       # Examples
       ```jldoctest
       julia> quasispecies([dna"GGG",dna"AAA",dna"GAG"],dna"GGG",0.5,0.001,20)
       
       ```
       """
function quasispecies(Input::Array, Fittest, Fitness, Mutationrate, Steps)

    function remove!(arr, item)
        deleteat!(arr, findall(x->x==item, arr))
    end        
        
### First Part
    
    ## give starting pop
    pop = String[]
    if typeof(Input[1]) == String
        pop = String[Input]
    else        
        for i in 1:length(Input)
            push!(pop, Input[i])
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
    println(mutation)
    # fitnesscounter    
    fitnesscounter = 0    
    
    # arrays for potential mutants
    splitmut = []
    
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
    for i in 1:length(pop)
        if pop[i] == Fittest
           append!(sel, Fitness) 
        else
            for j in 1:length(pop[i])
                if pop[i][j] == Fittest[j]
                    fitnesscounter = fitnesscounter +1
                    
                end
            end
            append!(sel, Fitness * (fitnesscounter / length(pop[i])))
            fitnesscounter = 0
        end
    end
    
    # give weight to sel
    sel = Weights(sel)
    

    # count actual numbers in pop
    for i in 1:length(pop)
        actualnumbers[pop[i]] = (count(isequal(pop[i]),pop))
    end
        
    # adding values to plotcounter
    for i in 1:length(pop)
        plotcount[pop[i]] = [(count(isequal(pop[i]),pop))]
    end


    ### Third part
    
    
    # loop for sampling until pop is one species
    while n <= Steps
    
    # create new pop (same size) with samples from old pop with mutation
        
        while counter <= length(pop)
            push!(samplepop,sample(pop,sel))
            counter = counter + 1
        end

        counter = 1

        for t in samplepop
            append!(splitmut,t)
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
            splitmut = [join(splitmut)]
            push!(newpop, join(splitmut))
            splitmut = []   
            
        end


    # clear actual numbers   
        actualnumbers = Dict()   
    
     # adding mutant to keys if not already there
        for i in newpop
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
                counter = 1
            end
        end
        
    # increase timesteps
        timesteps = timesteps + 1
    
    # put newpop in pop
        append!(pop, newpop)    
    

        
        
        
           
    # calculate selection for elements in array Input
        sel = Float64[]
        for i in 1:length(pop)
            if pop[i] == Fittest
               append!(sel, Fitness) 
            else
                for j in 1:length(pop[i])
                    if pop[i][j] == Fittest[j]
                    fitnesscounter = fitnesscounter +1
                    
                    end
                end
                append!(sel, Fitness * (fitnesscounter / length(pop[i])))
                fitnesscounter = 0
            end
        end 
    
    # give weight to sel
        sel = Weights(sel)        
    # create dataframe with Seq and Sel
        popdata = DataFrame(Seq = [], Sel = []) 
        
    # put data into frame
        for i in 1:length(pop)
            push!(popdata, [pop[i] sel[i]])
        end
    
    # sort by selection
        sort!(popdata, [:Sel])        
    
    # remove half of seq with lowest selection
        die = length(pop) / 2
        while counter <= die
            deleterows!(popdata, 1)
            counter = counter  +1 
        end
        counter = 1

    # create pop with remaining seq
        pop = String[]
        popseq = convert(Matrix, popdata[:,1:1])
        for i in 1:length(popseq)
            push!(pop, popseq[i])
        end
    
    # create sel with with remaining sel from dataframe
        sel = Float64[]
        popsel = convert(Matrix, popdata[:,2:2])
        for i in 1:length(popsel)
            push!(sel, popsel[i])
        end

    
    # give weight to sel
        sel = Weights(sel)
    
    # calculate new actual numbers   
        for i in 1:length(pop)
            actualnumbers[pop[i]] = (count(isequal(pop[i]),pop))
        end       
    
    #clear newpop and samplepop
        newpop= []        
        samplepop = []
    
    # put values of keys in plotcounter, if key is not there than put 0 in
        for i in keys(plotcount)
            if i in keys(actualnumbers)              
                push!(plotcount[i],actualnumbers[i][1])      
            else
                push!(plotcount[i], 0)

            end

        end

    # increase n by 1
        n = n+1
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
    
    # output
    return [final, plotdata, plotcount]
end
