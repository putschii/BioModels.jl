"""

Package with different functions to calculate biological models like Wright-Fisher-, Moran-, Neutral- and Quasispecies-Model. The package works with DNA sequences from the package BioSequences.
The functions cloudcreation and seqcreation can be used to create test setups for the models.

# Functions contained in the package:
```
wrightfisher(Input::Array,Fittest, Fitness, Mutationrate)
moran(Input::Array,Fittest, Fitness, Mutationrate)
neutral(Input::Array,Meta::Array, mutationrate, migration::Bool, save::Bool, file)
quasispecies(Input::Array, Fittest, Fitness, Mutationrate, Deathrate, Steps, Rep)
simplequasi(Seq::Array, Numbers::Array, Fittness::Array, Fmutation::Array, Bmutation::Array, Steps)
seqcreation(size, length)
cloudcreation(size, length, distance)
```
# For details of the functions call:
```
?functionname
```
       """
module BioModels

using StatsBase, BioSequences, Random, DataFrames, FileIO, JLD2
export neutral, wrightfisher, moran, quasispecies, simplequasi, seqcreation, cloudcreation

include("neutral.jl")
include("wrightfisher.jl")
include("moran.jl")
include("seqcreation.jl")
include("cloudcreation.jl")
include("quasispecies.jl")
include("simplequasi.jl")


end # module
