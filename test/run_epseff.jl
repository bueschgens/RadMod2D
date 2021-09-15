using Pkg
Pkg.activate(".")
Pkg.instantiate()

using RadMod2D

include("./cases_epseff.jl")

# epsilon_effective_round_groove()
# epsilon_effective_v_groove()
# epsilon_effective_square_groove()
# epsilon_effective_elbow_piece()
# epsilon_effective_trapez()
epsilon_effective_cosinus()


# include("./test/run_epseff.jl")