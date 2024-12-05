using Pkg
Pkg.activate(".")
Pkg.instantiate()

using RadMod2D

include("./models2D.jl")
include("./plot2D.jl")

using CairoMakie
# CairoMakie.activate!(type = "svg")
# filetype = ".svg"
CairoMakie.activate!(type = "png")
filetype = ".png"

d = 0.0445 # d rohr
l = 3.6 # kanallaenge
w = 2.55 # kanalbreite
nrl = 20 # rohranzahl laengs in l-dir
nrq = 31 # rohranzahl quer in w-dir
lt = 0.069 # laengsteilung
qt = 0.08 # querteilung
elemsize_kanal = 0.05 # diskretisierung
elemsize_rohr = 0.01 # diskretisierung

m = model_reku(d, l, w, nrl, nrq, lt, qt, elemsize_kanal, elemsize_rohr)

fontsize = 24
linewidth = 3.0



function fig_mesh(m)    
    # 2D plot
    fig1 = Figure(resolution = (1000, 1000), font = "Arial", fontsize = fontsize)
    ax1 = fig1[1,1] = Axis(fig1)
    ax1.xlabel = "X in m"
    ax1.ylabel = "Y in m"
    ax1.aspect = DataAspect()
    colors = Vector{Any}(undef,m.npar)
    colors[1:4] .= :red
    colors[5:end] .= :blue
    plot_model(fig1, ax1, m, shownvec = false, shownodes = false, showcom = false, showleg = false, legpos = "right", linewidth = linewidth, colors = colors)
    # display(fig1)
    save("fig_reku_mesh2D"*filetype, fig1)
end


function calc_vf(m)
    # view factors
    vfmat = zeros(Float64, m.nelem, m.nelem)
    existing_vf!(m, vfmat)
    n = 20
    dx, dy = create_tiles(m, n)
    t_occ = check_tile_occupation(m, dx, dy, n)
    blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
    calculating_vf!(m, vfmat, normit = false)
    vfmatp = vfmat_to_parts(m, vfmat, normit = true)
    println("vf mat for parts:")
    for i = 1:size(vfmatp,1)
        println(vfmatp[i,:])
    end
end


function calc_q(m)
    # view factors
    vfmat = zeros(Float64, m.nelem, m.nelem)
    existing_vf!(m, vfmat)
    n = 20
    dx, dy = create_tiles(m, n)
    t_occ = check_tile_occupation(m, dx, dy, n)
    blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
    calculating_vf!(m, vfmat, normit = true)
    # net radiation method
    epsilon = zeros(m.nelem)
    set_bc_part!(m, epsilon, [1,3], 1.0)
    set_bc_part!(m, epsilon, [2,4], 0.6)
    set_bc_part!(m, epsilon, 5:m.npar, 0.4)
    temp = zeros(m.nelem)
    set_bc_part!(m, temp, [1,2,4], 800+273.15)
    set_bc_part!(m, temp, 3, 400+273.15)
    set_bc_part!(m, temp, 5:m.npar, 500+273.15)
    Q, G = tempsolver(m, vfmat, temp, epsilon)
    area = [m.elem[i].area for i = m.elem2par[1].first:m.elem2par[end].last]
    q = Q[:] ./ area[:]
    # 2D plot of temperature
    fig2 = Figure(resolution = (1000, 1000), font = "Arial", fontsize = fontsize)
    ax2 = fig2[1,1] = Axis(fig2)
    ax2.xlabel = "X in m"
    ax2.ylabel = "Y in m"
    ax2.aspect = DataAspect()
    cbar2 = plot_model_with_value(fig2, ax2, m, temp .- 273.15, "Temperatur in °C", cmap = :jet, showcbar = true, cbarpos = "right", linewidth = linewidth)
    cbar2.height = Relative(2/3)
    save("fig_reku_temp2D"*filetype, fig2)
    # 2D plot of heat flux
    fig3 = Figure(resolution = (1000, 1000), font = "Arial", fontsize = fontsize)
    ax3 = fig3[1,1] = Axis(fig3)
    ax3.xlabel = "X in m"
    ax3.ylabel = "Y in m"
    ax3.aspect = DataAspect()
    cbar3 = plot_model_with_value(fig3, ax3, m, q .* (-1), "Wärmestromdichte in W/m", cmap = :jet, showcbar = true, cbarpos = "right", linewidth = linewidth)
    cbar3.height = Relative(2/3)
    save("fig_reku_q2D"*filetype, fig3)
end

# fig_mesh(m)
# calc_vf(m)
calc_q(m)

# include("./test/run_reku.jl")