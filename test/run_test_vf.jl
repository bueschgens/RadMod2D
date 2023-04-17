using Pkg
Pkg.activate(".")
Pkg.instantiate()

using RadMod2D

include("./models2D.jl")

# test parameter
d1 = 0.8
d2 = 1.6
elemsize = 0.005
n = 15
criteria = 1E-5

# analytical
r1 = d1 / 2
r2 = d2 / 2
vfmata = zeros(Float64, 2, 2)
vfmata[1,1] = 0
vfmata[1,2] = 1
vfmata[2,1] = r1 / r2
vfmata[2,2] = 1- vfmata[2,1]

# RadMod2D
m = model_circle_in_circle_centered(0.8, 1.6, elemsize)
vfmat = zeros(Float64, m.no_elements, m.no_elements)
@time existing_vf!(m, vfmat)
dx, dy = get_tile_dimensions(m, n)
@time t_occ = check_tile_occupation(m, dx, dy, n)
@time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
@time calculating_vf!(m, vfmat, normit = false)
vfmatp = compact_vfmat_to_parts(m, vfmat, normit = false)

# compare results to check if algorithm is working properly
println("analytical results:")
for i = 1:size(vfmata,1)
    println(vfmata[i,:])
end
println("RadMod2D results:")
for i = 1:size(vfmatp,1)
    println(vfmatp[i,:])
end
diff = vfmata - vfmatp
# println("difference:")
# for i = 1:size(diff,1)
#     println(diff[i,:])
# end
for i = 1:4
    if diff[i] > criteria
        error("Radmod2D not working properly")
    end
end

# include("./test/run_test_vf.jl")