using Pkg
Pkg.activate(".")
Pkg.instantiate()

using RadMod2D

include("./models2D.jl")

#### calculation
elemsize = 0.01
n = 60

# m = model_circles_in_circle_rand_full(0.1, 1.8, elemsize)
# m = model_circles_in_circle_rand_half(0.1, 1.8, elemsize)
m = model_circles_in_circle_rand_quarter(0.1, 1.8, elemsize)
# m = model_circles_in_circle_cross(0.1, 1.8, 9, elemsize)

# tiles
dx, dy = get_tile_dimensions(m, n)
@time t_occ = check_tile_occupation(m, dx, dy, n, blockparts = 1:(m.no_parts-1))
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
@time blocking_vf_with_tiles_simplified!(m, vfmat_simp, dx, dy, n, t_occ, elem_in_t = 3, skipped_t = 3)
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
plot_model(fig, ax, m, shownvec = true, shownodes = false, showcom = false, showleg = false)
display(fig)


#### compare results of classic algo with simplified algo

control_classic = sum(vfmatp, dims=2)
control = zeros(size(control_classic,1),3)
control[:,1] = control_classic

control_simp = sum(vfmatp_simp, dims=2)
control[:,2] = control_simp

diff = control_classic - control_simp
control[:,3] = diff

println("classic algo / simplified algo / diff")
for i = 1:size(control,1)
    println(control[i,:])
end

sum_classic = sum(control_classic)
sum_simp = sum(control_simp)
println("error rate: ", round((1-sum_simp/sum_classic)*100, digits = 2), " %")


# include("./test/run_algo_vs_simpalgo.jl")