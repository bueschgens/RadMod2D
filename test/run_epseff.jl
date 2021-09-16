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

# plot_epseff_round_groove()
# plot_epseff_v_groove()
# plot_epseff_square_groove()
# plot_epseff_elbow_piece()
# plot_epseff_trapez()
plot_epseff_cosinus()


# include("./test/run_epseff.jl")