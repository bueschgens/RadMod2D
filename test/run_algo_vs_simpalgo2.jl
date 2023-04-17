using Pkg
Pkg.activate(".")
Pkg.instantiate()

using RadMod2D

include("./models2D.jl")

#### calculation
elemsize = 0.01
n = 50

m = model_line_to_line_with_obstacles(1.0, 1.0, 0.5, 6, elemsize)

# tiles
dx, dy = get_tile_dimensions(m, n)
t_occ = check_tile_occupation(m, dx, dy, n, blockparts = 3:m.no_parts)
tile_occ_analysis(t_occ, printit = true)

println("elem size / tile size: ", round(elemsize/((dx+dy)/2), digits = 2), " %")

# classic algo
vfmat = zeros(Float64, m.no_elements, m.no_elements)
existing_vf!(m, vfmat)
@time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
calculating_vf!(m, vfmat, normit = false)
vfmatp = compact_vfmat_to_parts(m, vfmat, normit = false)

# simplified algo
vfmat_simp = zeros(Float64, m.no_elements, m.no_elements)
existing_vf!(m, vfmat_simp)
@time blocking_vf_with_tiles_simplified!(m, vfmat_simp, dx, dy, n, t_occ, elem_in_t = 1, skipped_t = 1)
calculating_vf!(m, vfmat_simp, normit = false)
vfmatp_simp = compact_vfmat_to_parts(m, vfmat_simp, normit = false)


#### 2D plot
fig = Figure(resolution = (900, 900))
ax = fig[1, 1] = Axis(fig)
ax.xlabel = "X in m"
ax.ylabel = "Y in m"
ax.aspect = DataAspect()
plot_empty_tiles(fig, ax, dx, dy, n)
plot_occupied_tiles(fig, ax, dx, dy, n, t_occ)
plot_model(fig, ax, m, shownvec = true, shownodes = false, showcom = false, showleg = true)
display(fig)


#### compare results of classic algo with simplified algo

vf12_classic = vfmatp[1,2]
vf12_simp = vfmatp_simp[1,2]

diff = vf12_classic - vf12_simp

println("classic algo: ", vf12_classic)
println("simplified algo: ", vf12_simp)
println("diff: ", diff)

println("error rate: ", round((1-vf12_simp/vf12_classic)*100, digits = 2), " %")


# include("./test/run_algo_vs_simpalgo2.jl")