"""
cloudcreation(size, length, distance)

Function to create an array with a random sequences of type dna with defined length and copies with defined distance.

       # Examples
       ```jldoctest
       julia> cloudcreation(3,4,2)

       ```
	["ATAG","AGGG","TCAG"]
"""
function cloudcreation(size, length, distance)
    # adding base for mutation
    basesA = ["T","G","C"]
    basesG = ["A","T","C"]
    basesT = ["A","G","C"]
    basesC = ["A","T","G"]
    # Start at second squence
    start = 2
    # Setup array for sequences
    seqstorage = []
    # Mastersequence
    seq = randdnaseq(length)
    seq = convert(String,seq)
    # Array for splitting bases
    split = []
    # size of non mutations
    no_mutation = zeros(length-distance)
    # size of mutations
    mutation = ones(distance)
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
    while start <= size
        # Shuffle array with 0 and 1 for mutation
        mutate = shuffle(mutate)
	# Split bases
        for i in seq
            append!(split,i)
        end
        # Mutate
        for i in 1:length
            if mutate[i] == 1.0
                if split[i] == 'A'
                    split[i] = sample(basesA)
                elseif split[i] == 'G'
                    split[i] = sample(basesG)
                elseif split[i] == 'C'
                    split[i] = sample(basesC)
                elseif split[i] == 'T'
                    split[i] = sample(basesT)
                end
            end
        end
        # Add new sequence to storage
        push!(seqstorage, join(split))
        # Clear array
        split = []
        # Increase Start
        start = start + 1


	end
        # Create output with sequences of type DNA
	for i in 1:size
    	push!(finaloutput,seqstorage[i])
    end
    # Output
    return finaloutput
end
