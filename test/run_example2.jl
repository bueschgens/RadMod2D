using Pkg
Pkg.activate(".")

using RadMod2D

include("./models2D.jl")
include("./plot2D.jl")

# backend
# using GLMakie
using CairoMakie
CairoMakie.activate!(type = "svg")
filetype = ".svg"

# parameter
r1x = 0.8
r1y = 0.8
r2x = 1.4
r2y = 1.4
elemsize = 0.05
n = 15

# model
m = model_rect_in_rect(r1x, r1y, r2x, r2y, elemsize)

# 2D plot of model
fig1 = Figure(resolution = (350, 350), font = "Arial", fontsize = 12)
ax1 = fig1[1,1] = Axis(fig1)
ax1.xlabel = "X in m"
ax1.ylabel = "Y in m"
ax1.aspect = DataAspect()
plot_model(fig1, ax1, m, shownvec = true, shownodes = false, showcom = false, showleg = true, legpos = "right", linewidth = 1.5)

# view factors
vfmat = zeros(Float64, m.nelem, m.nelem)
existing_vf!(m, vfmat)
dx, dy = create_tiles(m, n)
t_occ = check_tile_occupation(m, dx, dy, n)
blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
calculating_vf!(m, vfmat, normit = true)

# net radiation method
epsilon = zeros(m.nelem,1)
set_bc_part!(m, epsilon, 1, 0.3)
set_bc_part!(m, epsilon, 2, 0.6)
temp = zeros(m.nelem,1)
set_bc_part!(m, temp, 1:2, 300)
temp[79:106,1] .= 600 # bc for specific elements
@time Q, G = tempsolver(m, vfmat, temp, epsilon)
area = [m.elem[i].area for i = m.elem2par[1].first:m.elem2par[end].last]
q = Q[:] ./ area[:]

# 2D plot of heat flux density
fig2 = Figure(resolution = (350, 350), font = "Arial", fontsize = 12)
ax2 = fig2[1,1] = Axis(fig2)
ax2.xlabel = "X in m"
ax2.ylabel = "Y in m"
ax2.aspect = DataAspect()
cbar = plot_model_with_value(fig2, ax2, m, q, "Heat flux density in W/m", cmap = :jet, showcbar = true, cbarpos = "right", linewidth = 1.5)
cbar.height = Relative(2/3)

# display(fig1)
# display(fig2)
save("fig_rect_in_rect_mesh2D"*filetype, fig1)
save("fig_rect_in_rect_q2D"*filetype, fig2)

# include("./test/run_example2.jl")