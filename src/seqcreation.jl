"""
seqcreation(size, length)

Function to create an array with random sequences of type DNA with a defined length.

       # Examples
       ```jldoctest
       julia> seqcreation(3,4)

       ```
	["ATAG","AGGG","TCTG"]
"""
function seqcreation(size, length)


### Errors

    if typeof(size) != Int64 && typeof(size) != Int32
        throw("Size has to be a Int")
    end
    if typeof(length) != Int64 && typeof(length) != Int32
        throw("Length has to be a Int")
    end

# Startin at value of 1
	i = 1
	# Storage for sequences
	seqstorage = [randdnaseq(length)]
    i = i +1
	# Adding random sequences of defined length to the storage
	while i <= size
		seq = randdnaseq(length)
		push!(seqstorage,seq)
		i = i +1
	end
	i = 1
	return seqstorage
end
