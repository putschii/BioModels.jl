"""

Package with different functions to calculate biological models like Wright-Fisher-, Moran-, Neutral- and Quasispecies-Model. The package works with DNA sequences from the package BioSequences.
The functions cloudcreation and seqcreation can be used to create test setups for the models.

# Functions contained in the package:
```
wrightfisher(Input::Array,Fittest, Mutationrate, Steps, Save::Bool, File)
moran(Input::Array,Fittest, Mutationrate, Steps, Save::Bool, File)
neutral(Input::Array,Meta::Array, Mutationrate, Steps, Migration::Bool, dynamics::Bool, File)
quasispecies(Input::Array,Fittest, Mutationrate, Steps, Save::Bool, File)
simplequasi(Seq::Array, Numbers::Array, Fittness::Array, Fmutation::Array, Bmutation::Array, Steps)
seqanalysis(File,Start,Steps,Stop)
seqcreation(size, length)
cloudcreation(size, length, distance)
```
# For details of the functions call:
```
?functionname
```
       """
module BioModels

using StatsBase, BioSequences, Random, DataFrames, FileIO, JLD2, LinearAlgebra
export neutral, wrightfisher, moran, quasispecies, simplequasi, seqcreation, cloudcreation, seqanalysis

include("neutral.jl")
include("wrightfisher.jl")
include("moran.jl")
include("seqcreation.jl")
include("cloudcreation.jl")
include("quasispecies.jl")
include("simplequasi.jl")
include("seqanalysis.jl")


end # module
