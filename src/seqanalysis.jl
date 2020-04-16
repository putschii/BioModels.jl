"""
seqanalysis(File,Start,Steps,Stop)

The function seqanalysis works as a tool to analyse the created files with the wrightfisher-, moran-, neutral- and quasispecies-function. The function will loop through each time step of the file defined by the user. The user can define the start and stop time step and how many steps will be tested between these two values.


       # Examples
       ```jldoctest
       julia> seqanalysis("filenamehere.jld2",1,1,2)

       ```
       """
function seqanalysis(File,Start,Steps,Stop)

    # turn File into variable
    file = jldopen(File, "a+")

    # counter for looping through long sequence
    counter1 = 1
    counter2 = size(file["1"])[2]


    # Matrix for all sequences
    allseq = Matrix(undef,size(file["1"])[1],1)

    # create array to load specific steps
    filenumberarray = [Start:Steps:Stop;]

    # turn array with int to string (string is needed to open data)
    stringarray = []
    for i in filenumberarray
        push!(stringarray, string(i))
    end



    # for-loop for n files
    for element in stringarray


        # Create a empty sequence for all bases
        collect = dna""


        # load data from file
        basematrix  = file[element]

        # create sequence to collect all bases of the file
        sequences = Matrix(undef,size(basematrix)[1],1)
        # create for each row a key
        for j in 1:size(basematrix)[1]
            # loop through rows in first file
            for k in basematrix[j,:]

                # push bases into one long sequence
                push!(collect,k)
            end
        end

        # Loop through big sequence and cut it into single sequences
        if counter1 < length(collect)
            for i in 1:size(file["1"])[1]
                sequences[i,1] = collect[counter1:counter2]
                counter1 = counter1 + size(file["1"])[2]
                counter2 = counter2 + size(file["1"])[2]
            end
        end

        # Reset counter
        counter1 = 1
        counter2 = size(file["1"])[2]

        # combine results
        allseq = hcat(allseq,sequences)

        # empty stored sequence
        clearseq(collect)


    end

    # remove undef from array
    allseq = allseq[:,2:end]

    # create a dict with all possible sequences as key and start value of 0
    seqkeys = Dict()

    for i in allseq
        seqkeys[i] = [0]
    end

    # loop through each column
    for i in 1:size(allseq)[2]

        # if sequences is there than increase number by 1
        for sequences in allseq[:,i]
            seqkeys[sequences][i] =  seqkeys[sequences][i] + 1
        end

        # add new starting value for next column
        for element in keys(seqkeys)
            push!(seqkeys[element], 0)
        end
    end

    # remove unnecessary pushed 0
    for element in keys(seqkeys)
        pop!(seqkeys[element])
    end




    # Output
    return([seqkeys,allseq])

end
