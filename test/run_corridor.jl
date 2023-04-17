using Pkg
Pkg.activate(".")
Pkg.instantiate()

using RadMod2D

include("./models2D.jl")


filetype = ".svg"
# CairoMakie.activate!(type = "png")
# filetype = ".png"

function fig_corridor_therm2D_eps(eps_lab)
    # test blocking with shadow plot
    m = model_labyrinth([0.4, -0.4, 0.4, 0.2, 0.5, -0.3, -0.2], 0.1, 0.02)
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
    set_bc_part!(m, temp, 1, 600)
    set_bc_part!(m, temp, 2:m.no_parts, 400)
    # eps = 1.0
    epsilon = zeros(m.no_elements,1)
    set_bc_part!(m, epsilon, 1, 1.0)
    set_bc_part!(m, epsilon, 2:(m.no_parts-1), eps_lab)
    set_bc_part!(m, epsilon, m.no_parts, 1.0)
    @time Qp, G = tempsolver(m, vfmat, temp, epsilon)
    Qp_parts = [sum(Qp[m.elem2par[i].first:m.elem2par[i].last,1]) for i = 1:m.no_parts]
    println("heat flux start of labyrith: ", Qp_parts[1])
    println("heat flux end of labyrith: ", Qp_parts[end])
    # lab_min = minimum(Qp[m.elem2par[2].first:m.elem2par[end].last,1])
    # lab_max = maximum(Qp[m.elem2par[2].first:m.elem2par[end].last,1])
    # println("Qp of labyrith between: ", lab_min, " and ", lab_max)
    # heat flux density
    area = [m.elem[i].area for i = m.elem2par[1].first:m.elem2par[end].last]
    qp_area = Qp[:] ./ area[:]
    qp_area_kw = qp_area ./ 1000
    lab_min = minimum(qp_area_kw[m.elem2par[2].first:m.elem2par[end].last,1])
    lab_max = maximum(qp_area_kw[m.elem2par[2].first:m.elem2par[end].last,1])
    println("qp_area_kw of labyrith between: ", lab_min, " and ", lab_max)
    # modifying qw for displaying only incoming heat flux
    qp_area_kw_mod = qp_area_kw
    qp_area_kw_mod[qp_area_kw.>=0.0] .= 0.0
    qp_area_kw_mod .*= (-1)
    # 2D plot
    fig = Figure(resolution = (270, 400), font = "Arial", fontsize = 10)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis!(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    # changing legend in therm2D plotting to appear underneath
    cbar = plot_model_with_value(fig, ax, m, qp_area_kw_mod, "Incoming heat flux density in kW/m", cmap = :viridis, linewidth = linewidth*1.5, showcbar = true)
    # cbar.ticks = -50:20:50
    # cbar.limits = (lab_min,lab_max)
    cbar.limits = (0,0.4)
    cbar.ticks = 0.0:0.1:0.4
    cbar.vertical = false
    cbar.width = 120
    cbar.height = 12
    cbar.labelpadding = 1.0 # default 5.0
    cbar.ticklabelpad = 1.0 # default 3.0
    cbar.ticksize = 2.0 # default 6.0
    cbar.tellwidth = false # erlaubt es das colorbar und plot andere groeßen zueinander haben
    # plot_model_elements(fig, ax, m, m.elem2par[1].first:m.elem2par[1].last, color = :red, linewidth = linewidth*1.5)
    # display(fig)
    eps_save = replace(string(eps_lab),"." => "n")
    save("fig_c3_therm2D_eps"*eps_save*"_cb"*filetype, fig)
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

# fig_corridor_therm2D_eps(0.0001)
# fig_corridor_therm2D_eps(0.05)
# fig_corridor_therm2D_eps(0.75)

# include("./test/run_corridor.jl")