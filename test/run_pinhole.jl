using Pkg
Pkg.activate(".")
Pkg.instantiate()

using RadMod2D

include("./models2D.jl")


function fig_pinhole_therm2D_eps(eps_lab)
    # two rectangles connected with pinhole
    m = model_two_rectangles_with_holes(0.02)
    # blocking with tilewalk
    vfmat = zeros(Float64, m.no_elements, m.no_elements)
    existing_vf!(m, vfmat)
    n = 30
    dx, dy = get_tile_dimensions(m, n)
    @time t_occ = check_tile_occupation(m, dx, dy, n)
    @time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
    @time calculating_vf!(m, vfmat, normit = true)
    # solve Qp
    temp = zeros(m.no_elements,1)
    set_bc_part!(m, temp, 1:2, 400)
    set_bc_part!(m, temp, 3, 600)
    set_bc_part!(m, temp, 4:5, 400)
    # set_bc_part!(m, temp, 1:5, 600)
    set_bc_part!(m, temp, 6:10, 400)
    epsilon = zeros(m.no_elements,1)
    set_bc_part!(m, epsilon, 1:10, 1.0)
    set_bc_part!(m, epsilon, 8, eps_lab)
    @time Qp, G = tempsolver(m, vfmat, temp, epsilon)
    Qp_parts = [sum(Qp[m.elem2par[i].first:m.elem2par[i].last,1]) for i = 1:m.no_parts]
    # rect2_min = minimum(Qp[m.elem2par[6].first:m.elem2par[10].last,1])
    # rect2_max = maximum(Qp[m.elem2par[6].first:m.elem2par[10].last,1])
    # println("Qp of rect2 between :", rect2_min, " and ", rect2_max)
    area = [m.elem[i].area for i = m.elem2par[1].first:m.elem2par[end].last]
    qp_area = Qp[:] ./ area[:]
    qp_area_kw = qp_area ./ 1000
    rect2_min = minimum(qp_area_kw[m.elem2par[6].first:m.elem2par[10].last,1])
    rect2_max = maximum(qp_area_kw[m.elem2par[6].first:m.elem2par[10].last,1])
    println("qp_area_kw of rect2 between :", rect2_min, " and ", rect2_max)
    # modifying qw for displaying only incoming heat flux on rect2
    qp_area_kw_mod = qp_area_kw
    qp_area_kw_mod[qp_area_kw.>=rect2_max] .= 0.0
    qp_area_kw_mod[qp_area_kw.<=rect2_min] .= rect2_min
    qp_area_kw_mod .*= (-1)
    # 2D plot
    fig = Figure(resolution = (270, 400), font = "Arial", fontsize = 10)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis!(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    cbar = plot_model_with_value(fig, ax, m, qp_area_kw_mod, "Incoming heat flux density in kW/m", cmap = :viridis, linewidth = linewidth*1.5, showcbar = true)
    # cropping here at 0 and only displaying incoming heat flux - also making all positive
    # cbar.ticks = -50:20:50
    # cbar.limits = (-10,0.1)
    cbar.limits = (0, 0.2)
    cbar.ticks = 0.0:0.05:0.2
    cbar.vertical = false
    cbar.width = 120
    cbar.height = 12
    cbar.labelpadding = 1.0 # default 5.0
    cbar.ticklabelpad = 1.0 # default 3.0
    cbar.ticksize = 2.0 # default 6.0
    cbar.tellwidth = false # erlaubt es das colorbar und plot andere groeßen zueinander haben
    # display(fig)
    eps_save = replace(string(eps_lab),"." => "n")
    save("fig_c2_therm2D_eps"*eps_save*"_cb"*filetype, fig)
    # for i = 1:m.nnodes
    #     println(m.nodes[i].x, ";", m.nodes[i].y)
    # end
    # for i = 1:m.no_elements
    #     println(m.elem[i].node1, ";", m.elem[i].node2, ";", m.elem[i].com, ";", m.elem[i].nvec, ";", m.elem[i].area)
    # end
    # for i = 1:m.no_elements
    #     println(temp[i,:], ";", epsilon[i,:], ";", qp_area_kw_mod[i,:])
    # end
end

# fig_pinhole_therm2D_eps(0.5)
fig_pinhole_therm2D_eps(1.0)

# include("./test/run_pinhole.jl")