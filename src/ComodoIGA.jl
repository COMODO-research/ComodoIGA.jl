module ComodoIGA

## use the packages
using Comodo
using Comodo.GeometryBasics
using Ferrite
using FerriteIGA
using GLMakie
using FerriteGismo

## export the packages
export Comodo
export GeometryBasics
export Ferrite
export FerriteIGA
export GLMakie
export FerriteGismo

include("functions.jl")

end