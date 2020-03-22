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
	# Startin at value of 1
	i = 1
	# Storage for sequences
	seqstorage = []
	# Adding random sequences of defined length to the storage
	while i <= size
		seq = randdnaseq(length)
		push!(seqstorage,seq)
		i = i +1
	end
	i = 1
	return seqstorage
end
