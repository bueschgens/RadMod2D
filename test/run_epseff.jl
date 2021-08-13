using Pkg
Pkg.activate(".")
Pkg.instantiate()


using RadMod2D
using DelimitedFiles

using GLMakie
# using CairoMakie

include("./epseff2D.jl")

# epsilon_effective_keilspalt()
epsilon_effective_rundrille()

# include("./test/run_epseff.jl")