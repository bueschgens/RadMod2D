using Pkg
Pkg.activate(".")
Pkg.instantiate()

using RadMod2D

include("./models2D.jl")
include("./plot2D.jl")

using CairoMakie
CairoMakie.activate!(type = "svg")
filetype = ".svg"
# CairoMakie.activate!(type = "png")
# filetype = ".png"

function fig_c1_mesh2D()
    m = model_circle_in_circle(0.8, 1.6, 0.05)
    # 2D plot
    fig = Figure(resolution = (270, 270), font = "Arial", fontsize = 10)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    plot_model(fig, ax, m, shownvec = true, shownodes = false, showcom = false, showleg = false, colors = [:blue, :red], linewidth = linewidth)
    # display(fig)
    save("fig_c1_mesh2D"*filetype, fig)
end

function fig_c1_mesh2D_legend()
    # settings for leg
    # leg = fig[end+1, 1] = Legend(fig, plt_elem, leg_string, orientation = :horizontal, tellwidth = false, tellheight = true)
    # leg.padding = (5.0f0, 5.0f0, 4.0f0, 4.0f0)
    # leg.linewidth = 1.0
    # leg.patchsize = (5, 3)
    # leg.rowgap = 0.1
    # leg.colgap = 8
    # for copy in mesh2D()
    m = model_circle_in_circle(0.8, 1.6, 0.05)
    # 2D plot
    fig = Figure(resolution = (270, 270), font = "Arial", fontsize = 10)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    # legend horizontal at the bottom
    plot_model(fig, ax, m, shownvec = true, shownodes = false, showcom = false, showleg = true, colors = [:blue, :red], linewidth = linewidth)
    # display(fig)
    save("fig_c1_mesh2D_legend"*filetype, fig)
end

function fig_c1_view2D_existing_2elem()
    m = model_circle_in_circle(0.8, 1.6, 0.1)
    # 2D plot
    fig = Figure(resolution = (270, 270), font = "Arial", fontsize = 10)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    # plot all
    plot_model(fig, ax, m, shownvec = false, shownodes = false, showcom = false, showleg = false, colors = [:black, :black], linewidth = linewidth)
    # elem pair check -> existing
    i1 = 47
    i2 = 33
    plot_model_elements(fig, ax, m, [i1,i2], shownvec = true, shownodes = false, showcom = false, color = :blue, linewidth = linewidth)
    linesegments!(ax, [Point2f0(m.elem[i1].com.x, m.elem[i1].com.y) => Point2f0(m.elem[i2].com.x, m.elem[i2].com.y)], color = :red, linewidth = linewidth)
    # elem pair check -> not existing
    i2 = 20
    plot_model_elements(fig, ax, m, [i1,i2], shownvec = true, shownodes = false, showcom = false, color = :blue, linewidth = linewidth)
    linesegments!(ax, [Point2f0(m.elem[i1].com.x, m.elem[i1].com.y) => Point2f0(m.elem[i2].com.x, m.elem[i2].com.y)], color = :red, linewidth = linewidth, linestyle = :dash)
    # display(fig)
    save("fig_c1_view2D_existing_2elem"*filetype, fig)
end

function fig_c1_view2D_existing1()
    m = model_circle_in_circle(0.8, 1.6, 0.1)
    #### vfmat calculation
    vfmat = zeros(Float64, m.nelem, m.nelem)
    @time existing_vf!(m, vfmat)
    #### 2D plot
    fig = Figure(resolution = (270, 270), font = "Arial", fontsize = 10)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    plot_model(fig, ax, m, shownvec = true, shownodes = false, showcom = false, showleg = false, colors = [:black, :black], linewidth = linewidth)
    plot_existing(fig, ax, m, 75, vfmat, color = :red, linewidth = 0.5*linewidth)
    # display(fig)
    save("fig_c1_view2D_existing1"*filetype, fig)
end

function fig_c1_view2D_existing2()
    m = model_circle_in_circle(0.8, 1.6, 0.1)
    #### vfmat calculation
    vfmat = zeros(Float64, m.nelem, m.nelem)
    @time existing_vf!(m, vfmat)
    n = 10
    dx, dy = create_tiles(m, n)
    @time t_occ = check_tile_occupation(m, dx, dy, n)
    @time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
    #### 2D plot
    fig = Figure(resolution = (270, 270), font = "Arial", fontsize = 10)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    plot_model(fig, ax, m, shownvec = true, shownodes = false, showcom = false, showleg = false, colors = [:black, :black], linewidth = linewidth)
    plot_existing(fig, ax, m, 75, vfmat, color = :red, linewidth = 0.5*linewidth)
    # display(fig)
    save("fig_c1_view2D_existing2"*filetype, fig)
end

function fig_c1_view2D_blocking_occ()
    m = model_circle_in_circle(0.8, 1.6, 0.05)
    # 2D plot
    fig = Figure(resolution = (270, 270), font = "Arial", fontsize = 10)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    # create tiles
    n = 15
    dx, dy = create_tiles(m, n)
    t_occ = check_tile_occupation(m, dx, dy, n)
    plot_occupied_tiles(fig, ax, dx, dy, n, t_occ)
    plot_empty_tiles(fig, ax, dx, dy, n, linewidth = linewidth)
    plot_model(fig, ax, m, shownvec = true, shownodes = false, showcom = false, showleg = false, colors = [:black, :black], linewidth = linewidth)
    # display(fig)
    save("fig_c1_view2D_bl_occ_nv"*filetype, fig)
end

function fig_c1_view2D_blocking_2elem_tiles()
    m = model_circle_in_circle(0.8, 1.6, 0.05)
    # 2D plot
    fig = Figure(resolution = (270, 270), font = "Arial", fontsize = 10)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    # create tiles
    n = 15
    dx, dy = create_tiles(m, n)
    t_occ = check_tile_occupation(m, dx, dy, n)
    # plot_occupied_tiles(fig, ax, dx, dy, n, t_occ)
    # elem pair check
    i1 = 90
    i2 = 150
    isexisting = existing_vf_2elem(m, i1, i2)
    if isexisting
        blocking_vf_with_tiles_2elem_tiles(fig, ax, i1, i2, m, dx, dy, n, t_occ)
    end
    plot_model(fig, ax, m, shownvec = false, shownodes = false, showcom = false, showleg = false, colors = [:black, :black], linewidth = linewidth)
    plot_empty_tiles(fig, ax, dx, dy, n, linewidth = linewidth)
    plot_model_elements(fig, ax, m, [i1,i2], shownvec = true, shownodes = false, showcom = false, color = :blue, linewidth = linewidth)
    if isexisting
        linesegments!(ax, [Point2f0(m.elem[i1].com.x, m.elem[i1].com.y) => Point2f0(m.elem[i2].com.x, m.elem[i2].com.y)], color = :red, linewidth = linewidth)
    end
    # display(fig)
    save("fig_c1_view2D_bl_2elem_tiles"*filetype, fig)
end

function fig_c1_view2D_blocking_2elem()
    m = model_circle_in_circle(0.8, 1.6, 0.05)
    # 2D plot
    fig = Figure(resolution = (270, 270), font = "Arial", fontsize = 10)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    # create tiles
    n = 15
    dx, dy = create_tiles(m, n)
    t_occ = check_tile_occupation(m, dx, dy, n)
    plot_occupied_tiles(fig, ax, dx, dy, n, t_occ)
    # elem pair check
    i1 = 90
    i2 = 150
    isexisting = existing_vf_2elem(m, i1, i2)
    if isexisting
        blocking_vf_with_tiles_2elem(fig, ax, i1, i2, m, dx, dy, n, t_occ)
    end
    plot_model(fig, ax, m, shownvec = false, shownodes = false, showcom = false, showleg = false, colors = [:black, :black], linewidth = linewidth)
    plot_empty_tiles(fig, ax, dx, dy, n, linewidth = linewidth)
    plot_model_elements(fig, ax, m, [i1,i2], shownvec = true, shownodes = false, showcom = false, color = :blue, linewidth = linewidth)
	# elements in bucket
	plot_model_elements(fig, ax, m, t_occ[6,11], shownvec = true, shownodes = false, showcom = false, color = :orange, linewidth = linewidth)
	# hitten element
	plot_model_elements(fig, ax, m, [17], shownvec = true, shownodes = false, showcom = false, color = :green, linewidth = linewidth)
    if isexisting
        linesegments!(ax, [Point2f0(m.elem[i1].com.x, m.elem[i1].com.y) => Point2f0(m.elem[i2].com.x, m.elem[i2].com.y)], color = :red, linewidth = linewidth)
    end
    # display(fig)
    save("fig_c1_view2D_bl_2elem"*filetype, fig)
end

function fig_c1_therm2D_temp()
    m = model_circle_in_circle(0.8, 1.6, 0.02)
    temp = zeros(m.nelem,1)
    set_bc_part!(m, temp, 1, 300)
    set_bc_part!(m, temp, 2, 400)
    # bc for some elements
    temp[189:313,1] .= 600
    # 2D plot
    fig = Figure(resolution = (270, 270), font = "Arial", fontsize = 10)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    cbar = plot_model_with_value(fig, ax, m, temp[:,1], "Temperature in K", cmap = :viridis, linewidth = linewidth*1.5, showcbar = false)
    # cbar.ticks = 300:100:600
    # display(fig)
    save("fig_c1_therm2D_temp"*filetype, fig)
    # for i = 1:m.nnodes
    #     println(m.nodes[i].x, ";", m.nodes[i].y)
    # end
    # for i = 1:m.nelem
    #     println(m.elem[i].node1, ";", m.elem[i].node2, ";", m.elem[i].com, ";", m.elem[i].nvec, ";", m.elem[i].area)
    # end
    # for i = 1:m.nelem
    #     println(temp[i,:])
    # end
end

function fig_c1_therm2D_temp_legend()
    m = model_circle_in_circle(0.8, 1.6, 0.02)
    temp = zeros(m.nelem,1)
    set_bc_part!(m, temp, 1, 300)
    set_bc_part!(m, temp, 2, 400)
    # bc for some elements
    temp[189:313,1] .= 600
    # 2D plot
    fig = Figure(resolution = (500, 270), font = "Arial", fontsize = 10)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    cbar = plot_model_with_value(fig, ax, m, temp[:,1], "Temperature in K", cmap = :viridis, linewidth = linewidth*1.5, showcbar = true)
    cbar.ticks = 300:100:600
    # cbar.flipaxis = true
    cbar.vertical = false
    cbar.width = 120
    cbar.height = 12
    cbar.labelpadding = 1.0 # default 5.0
    cbar.ticklabelpad = 1.0 # default 3.0
    cbar.ticksize = 2.0 # default 6.0
    # display(fig)
    save("fig_c1_therm2D_temp_legend"*filetype, fig)
end

function fig_c1_therm2D_qrad()
    m = model_circle_in_circle(0.8, 1.6, 0.02)
    # vfmat calculation
    vfmat = zeros(Float64, m.nelem, m.nelem)
    existing_vf!(m, vfmat)
    n = 15
    dx, dy = create_tiles(m, n)
    t_occ = check_tile_occupation(m, dx, dy, n)
    blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
    calculating_vf!(m, vfmat, normit = true)
    # solve Qp
    epsilon = zeros(m.nelem,1)
    set_bc_part!(m, epsilon, 1, 0.9)
    set_bc_part!(m, epsilon, 2, 0.5)
    temp = zeros(m.nelem,1)
    set_bc_part!(m, temp, 1, 300)
    set_bc_part!(m, temp, 2, 400)
    # bc for some elements
    temp[189:313,1] .= 600
    @time Qp, G = tempsolver(m, vfmat, temp, epsilon)
    Qp_parts = [sum(Qp[m.elem2par[i].first:m.elem2par[i].last,1]) for i = 1:m.npar]
    area = [m.elem[i].area for i = m.elem2par[1].first:m.elem2par[end].last]
    qp_area = Qp[:] ./ area[:]
    # 2D plot
    fig = Figure(resolution = (270, 270), font = "Arial", fontsize = 10)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    cbar = plot_model_with_value(fig, ax, m, qp_area, "Heat flux density in W/m", cmap = :viridis, linewidth = linewidth*1.5, showcbar = false)
    # cbar.ticks = -80:20:60
    # display(fig)
    save("fig_c1_therm2D_qrad"*filetype, fig)
    # for i = 1:m.nelem
    #     println(qp_area[i,:])
    # end
end

function fig_c1_therm2D_qrad_legend()
    m = model_circle_in_circle(0.8, 1.6, 0.02)
    # vfmat calculation
    vfmat = zeros(Float64, m.nelem, m.nelem)
    existing_vf!(m, vfmat)
    n = 15
    dx, dy = create_tiles(m, n)
    t_occ = check_tile_occupation(m, dx, dy, n)
    blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
    calculating_vf!(m, vfmat, normit = true)
    # solve Qp
    epsilon = zeros(m.nelem,1)
    set_bc_part!(m, epsilon, 1, 0.9)
    set_bc_part!(m, epsilon, 2, 0.5)
    temp = zeros(m.nelem,1)
    set_bc_part!(m, temp, 1, 300)
    set_bc_part!(m, temp, 2, 400)
    # bc for some elements
    temp[189:313,1] .= 600
    @time Qp, G = tempsolver(m, vfmat, temp, epsilon)
    Qp_parts = [sum(Qp[m.elem2par[i].first:m.elem2par[i].last,1]) for i = 1:m.npar]
    area = [m.elem[i].area for i = m.elem2par[1].first:m.elem2par[end].last]
    qp_area = Qp[:] ./ area[:]
    qp_area_kw = qp_area ./ 1000
    # 2D plot
    fig = Figure(resolution = (270, 270), font = "Arial", fontsize = 10)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    cbar = plot_model_with_value(fig, ax, m, qp_area_kw, "Heat flux density in kW/m", cmap = :viridis, linewidth = linewidth*1.5, showcbar = true)
    cbar.limits = (-4, 3)
    cbar.ticks = -4:1:3
    cbar.vertical = false
    cbar.width = 120
    cbar.height = 12
    cbar.labelpadding = 1.0 # default 5.0
    cbar.ticklabelpad = 1.0 # default 3.0
    cbar.ticksize = 2.0 # default 6.0
    # display(fig)
    save("fig_c1_therm2D_qrad_legend"*filetype, fig)
end

function calc_c1_view2D(;elemsize = 0.05)
    m = model_circle_in_circle(0.8, 1.6, elemsize)
    #### vfmat calculation
    vfmat = zeros(Float64, m.nelem, m.nelem)
    @time existing_vf!(m, vfmat)
    n = 15
    dx, dy = create_tiles(m, n)
    @time t_occ = check_tile_occupation(m, dx, dy, n)
    @time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
    @time calculating_vf!(m, vfmat, normit = false)
    vfmatp = vfmat_to_parts(m, vfmat, normit = false)
    vfmatp
end

function calc_c1_analytical()
    d1 = 0.8
    d2 = 1.6
    r1 = d1 / 2
    r2 = d2 / 2
    vfmat = zeros(Float64, 2, 2)
    vfmat[1,1] = 0
    vfmat[1,2] = 1
    vfmat[2,1] = r1 / r2
    vfmat[2,2] = 1- vfmat[2,1]
    vfmat
end

##### cyl in cyl as example for algorithm
# fig_c1_mesh2D()
# fig_c1_mesh2D_legend()
# fig_c1_view2D_existing_2elem()
# fig_c1_view2D_existing1()
# fig_c1_view2D_existing2()
# fig_c1_view2D_blocking_occ()
# fig_c1_view2D_blocking_2elem_tiles()
# fig_c1_view2D_blocking_2elem()
# fig_c1_therm2D_temp()
# fig_c1_therm2D_temp_legend()
# fig_c1_therm2D_qrad()
# fig_c1_therm2D_qrad_legend()

##### vf calculation vs analytical solution
# calc_c1_analytical()
# calc_c1_view2D(elemsize = 0.2)
# calc_c1_view2D(elemsize = 0.1)
# calc_c1_view2D(elemsize = 0.05)
# calc_c1_view2D(elemsize = 0.01)
# calc_c1_view2D(elemsize = 0.005)

# include("./test/run_cyl_in_cyl.jl")