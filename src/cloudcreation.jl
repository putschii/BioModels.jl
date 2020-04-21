"""
cloudcreation(Size, Length, Distance)

Function to create an array with a random sequences of type dna with defined length. The first created sequence will represent the master sequence,
all other sequences are copies with defined distance from the master sequence.

       # Examples
       ```jldoctest
       julia> cloudcreation(3,4,2)

       ```
	["ATAG","AGGG","TCAG"]
"""
function cloudcreation(Size, Length, Distance)
    
    
    ### Errors
    
    if typeof(Size) != Int64 && typeof(Size) != Int32
        throw("Size has to be a Int")
    end
    if typeof(Length) != Int64 && typeof(Length) != Int32
        throw("Length has to be a Int")
    end 
    
    # adding base for mutation
    basesA = [DNA_T,DNA_G,DNA_C]
    basesG = [DNA_A,DNA_T,DNA_C]
    basesT = [DNA_A,DNA_G,DNA_C]
    basesC = [DNA_A,DNA_T,DNA_G]
    # Start at second squence
    start = 2
    # Setup array for sequences
    seqstorage = []
    # Mastersequence
    seq = randdnaseq(Length)
    # size of non mutations
    no_mutation = zeros(Length-Distance)
    # size of mutations
    mutation = ones(Distance)
    # Array for mutation
    mutate = []
    # Outputarray
    finaloutput = []
    # Combine mutation and non mutation
    append!(mutate,mutation)
    append!(mutate,no_mutation)
    # Put Mastersequence into storage
    push!(seqstorage,seq)
    # Loop while
    while start <= Size
        # Shuffle array with 0 and 1 for mutation
        mutate = shuffle(mutate)
        seq = seqstorage[1][1:end]
        # Mutate
        for i in 1:length(mutate)
            if mutate[i] == 1.0
                if seq[i] == DNA_A
                    seq[i] = sample(basesA)
                elseif seq[i] == DNA_G
                    seq[i] = sample(basesG)
                elseif seq[i] == DNA_C
                    seq[i] = sample(basesC)
                elseif seq[i] == DNA_T
                    seq[i] = sample(basesT)
                end
            end
        end
        # Add new sequence to storage
        push!(seqstorage, seq)
        # Increase Start
        start = start + 1


	end

    # Output
    return seqstorage
end
