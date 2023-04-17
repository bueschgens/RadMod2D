using Revise

using Pkg
Pkg.activate(".")
Pkg.instantiate()


using RadMod2D

using GLMakie
# using CairoMakie

include("./tests.jl")


# test_mesh2D()
# test_raycast2D()
# test_view2D_blocking_2elem() # not working
# test_view2D_blocking_bf_vs_tiles() 
# test_view2D_vf()
# test_view2D_tiles(n = 30)
# test_therm2D()

# include("./test/run_test.jl")