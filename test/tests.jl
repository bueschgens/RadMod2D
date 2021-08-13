#= ######################################################
includes
###################################################### =# 

include("./models2D.jl")

#= ######################################################
mesh2D - testing
###################################################### =# 

function test_mesh2D()
    # which model
    # m = model_two_intesecting_circles(3.0, 120, 24)
    # m = model_rectangle(2.0, 4.0, 0.1)
    m = model_circle_in_circle(0.8, 1.6, 0.05)
    # m = model_circle_with_opening_line(3.0, 120, 15)
    # m = model_isosceles_triangle(2.0, 90, 0.2)
    # m = model_two_rectangles_with_holes(0.02)
    # m = model_lines()
    # m = model_cosinus(0.0, 1.0, 100)
    # m = model_trapez(2.0, 3.0, 10.0, 0.1)
    # m = model_right_triangle(3.0, 2.0, 0.05)
    # m = model_furnace1()
    # m = model_furnace2(0.195,0.5,3,0.7)
    # m = model_labyrinth([0.4, -0.4, 0.4, 0.2, 0.5, -0.3, -0.2], 0.1, 0.02)
    # m = model_square_in_square(0.8, 1.6, 0.01)
    # m = model_rect_in_rect(1.0, 3.0, 2.0, 4.0, 0.03)
    # m = model_circles_in_circle_cross(0.1, 1.8, 9, 0.01)
    # m = model_circles_in_circle_rand(0.1, 1.8, 17, 0.01)
    # m = model_circles_in_circle_rand_full(0.1, 1.8, 0.01)
    # m = model_circles_in_circle_rand_half(0.1, 1.8, 0.01)
    # m = model_circles_in_circle_rand_quarter(0.1, 1.8, 0.01)
    # element analysis
    element_analysis(m, printit = true)
    # 2D plot
    fig = Figure(resolution = (900, 900))
    ax = fig[1, 1] = Axis(fig)
    # ax.xlabelpadding = 0.0
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    plot_model(fig, ax, m, shownvec = true, shownodes = false, showcom = false)
    display(fig)
end

#= ######################################################
raycast2D - testing
###################################################### =# 

function test_raycast2D()
    # 2D plot
    fig = Figure(resolution = (900, 900))
    ax = fig[1, 1] = Axis(fig)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    # create tiles
    n = 10
    dx = 0.5
    dy = 0.5
    plot_empty_tiles(fig, ax, dx, dy, n)
    # create two points by random
    p1 = Point2D(rand()*n*dx,rand()*n*dy)
    p2 = Point2D(rand()*n*dx,rand()*n*dy)
    # p1 = Point2D(1*n*dx,rand()*n*dy)
    # p2 = Point2D(1*n*dx,rand()*n*dy)
    # p1 = Point2D(0.2,0.1)
    # p2 = Point2D(4.8,4.9)
    plot_points_and_connection(fig, ax, p1, p2)
    # tile walk - used for blocking calc
    max_steps = get_max_steps(n)
    tiles = Vector{Index2D{Int64}}(undef,max_steps)
    @time ntiles = tilewalk_with_return!(tiles, p1, p2, dx, dy, n)
    for i = 1:ntiles
        t = tiles[i]
        # print tile in color
        poly!(ax, Rect((t.x-1) * dx, (t.y-1) * dy, dx, dy), color = (:blue, 0.3))
    end
    display(fig)
end

function test_raycast2D_occ(;option = 3)
    # not correctly working because of missing Makie pkg in RadMod2D
    # 2D plot
    fig = Figure(resolution = (900, 900))
    ax = fig[1, 1] = Axis(fig)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    # create tiles
    n = 10
    dx = 0.5
    dy = 0.5
    plot_empty_tiles(fig, ax, dx, dy, n)
    t_occ = create_randomly_occupied_tiles(n, density = 0.2)
    plot_occupied_tiles(fig, ax, dx, dy, n, t_occ)
    # create two points by random
    p1 = Point2D(rand()*n*dx,rand()*n*dy)
    p2 = Point2D(rand()*n*dx,rand()*n*dy)
    # p1 = Point2D(1*n*dx,rand()*n*dy)
    # p2 = Point2D(1*n*dx,rand()*n*dy)
    # p1 = Point2D(0.2,0.1)
    # p2 = Point2D(4.8,4.9)
    plot_points_and_connection(fig, ax, p1, p2)
    # tile walk options
    if option == 1
        @time tilewalk_with_check_for_occ_tiles(fig, ax, p1, p2, dx, dy, n, t_occ)
    elseif option == 2
        @time tiles = tilewalk_with_return(fig, ax, p1, p2, dx, dy, n)
        for t in tiles
            println("tile number: ", t.x, ", ", t.y)
        end
    elseif option == 3
        # used for blocking calc
        max_steps = get_max_steps(n)
        tiles = Vector{Index2D{Int64}}(undef,max_steps)
        @time ntiles = tilewalk_with_return!(tiles, p1, p2, dx, dy, n)
        for i = 1:ntiles
            t = tiles[i]
            # print tile in color
            poly!(ax, Rect((t.x-1) * dx, (t.y-1) * dy, dx, dy), color = (:blue, 0.3))
        end
    else
        println("Option not defined...")
    end
    display(fig)
end

#= ######################################################
view2D - testing
###################################################### =# 

function test_view2D_blocking_2elem()
    # test double open circle case
    m = model_two_intesecting_circles(3.0, 120, 3)
    # m = model_circle_with_opening_line(2.0, 360, 143)
    # 2D plot
    fig = Figure(resolution = (900, 900))
    ax = fig[1, 1] = Axis(fig)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    # create tiles
    n = 30
    dx, dy = create_tiles(m, n)
    t_occ = check_tile_occupation(m, dx, dy, n)
    plot_empty_tiles(fig, ax, dx, dy, n)
    plot_occupied_tiles(fig, ax, dx, dy, n, t_occ)
    plot_model(fig, ax, m, shownvec = false, shownodes = false, showcom = false)
    # elem pair check
    i1 = 50
    i2 = 160
    plot_model_elements(fig, ax, m, [i1,i2], shownvec = true, shownodes = false, showcom = false)
    isexisting = existing_vf_2elem(m, i1, i2)
    if isexisting
        linesegments!(ax, [Point2f0(m.elem[i1].com.x, m.elem[i1].com.y) => Point2f0(m.elem[i2].com.x, m.elem[i2].com.y)], color = :blue, linewidth = 3)
        blocking_vf_with_tiles_2elem(fig, ax, i1, i2, m, dx, dy, n, t_occ)
    end
    display(fig)
end

function test_view2D_blocking_bf_vs_tiles()
    # test blocking: brute force vs tiles
    # m = model_two_intesecting_circles(3.0, 360, 36)
    m = model_two_rectangles_with_holes(0.02)
    # blocking with brute force
    vfmat = zeros(Float64, m.nelem, m.nelem)
    existing_vf!(m, vfmat)
    @time blocking_vf!(m, vfmat)
    # blocking with tilewalk
    vfmat2 = zeros(Float64, m.nelem, m.nelem)
    existing_vf!(m, vfmat2)
    n = 10
    dx, dy = create_tiles(m, n)
    @time t_occ = check_tile_occupation(m, dx, dy, n)
    @time blocking_vf_with_tiles!(m, vfmat2, dx, dy, n, t_occ)
    # check for difference
    diff = vfmat2 .- vfmat
    println(sum(diff))
    # findfirst(diff.==1.0)
    # 2D plot
    # fig = Figure(resolution = (1400, 900))
    # ax = fig[1, 1] = Axis(fig)
    # ax.xlabel = "X in m"
    # ax.ylabel = "Y in m"
    # ax.aspect = DataAspect()
    # plot_empty_tiles(fig, ax, dx, dy, n)
    # plot_occupied_tiles(fig, ax, dx, dy, n, t_occ)
    # plot_model(fig, ax, m, shownvec = false, shownodes = false, showcom = false)
    # plot_existing(fig, ax, m, 40, vfmat2)
    # display(fig)
end

function test_view2D_vf()
    #### which model
    # m = model_two_intesecting_circles(3.0, 360, 45)
    # m = model_rectangle(2.0, 5.0, 0.1)
    # m = model_circle_in_circle(3.3, 5.2, 0.1)
    # m = model_circle_with_opening_line(2.0, 360, 143)
    # m = model_labyrinth(0.02)
    # m = model_square_in_square(0.8, 1.6, 0.01)
    # m = model_rect_in_rect(0.8, 0.8, 1.6, 1.6, 0.05)
    # m = model_circles_in_circle_rand_full(0.1, 1.8, 0.005)
    # m = model_circles_in_circle_rand_half(0.1, 1.8, 0.005)
    # m = model_circles_in_circle_rand_quarter(0.1, 1.8, 0.005)
    m = model_circles_in_circle_cross(0.1, 1.8, 9, 0.01)
    #### vfmat calculation
    vfmat = zeros(Float64, m.nelem, m.nelem)
    @time existing_vf!(m, vfmat)
    n = 20
    dx, dy = create_tiles(m, n)
    @time t_occ = check_tile_occupation(m, dx, dy, n)
    @time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
    @time calculating_vf!(m, vfmat, normit = false)
    vfmatp = vfmat_to_parts(m, vfmat, normit = false)
    #### 2D plot
    fig = Figure(resolution = (1400, 900))
    ax = fig[1, 1] = Axis(fig)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    # plot_empty_tiles(fig, ax, dx, dy, n)
    # plot_occupied_tiles(fig, ax, dx, dy, n, t_occ)
    plot_model(fig, ax, m, shownvec = false, shownodes = false, showcom = false)
    plot_existing(fig, ax, m, 170, vfmat)
    display(fig)
    #### post processing
    # control = sum(vfmatp, dims=2)
    # save_vfmat_part_to_csv(vfmatp, "figures", "cyl_in_cyl_rand_cross_0n005", control = true)
    # vfmatp
end

function test_view2D_tiles(; n = 15)
    #### which model
    # m = model_two_intesecting_circles(3.0, 360, 45)
    # m = model_rectangle(2.0, 5.0, 0.1)
    # m = model_circle_in_circle(0.8, 1.6, 0.001)
    # m = model_circle_with_opening_line(2.0, 360, 143)
    # m = model_labyrinth(0.02)
    m = model_square_in_square(0.8, 1.6, 0.01)
    #### vfmat calculation
    vfmat = zeros(Float64, m.nelem, m.nelem)
    @time existing_vf!(m, vfmat)
    dx, dy = create_tiles(m, n)
    @time t_occ = check_tile_occupation(m, dx, dy, n)
    #### 2D plot
    fig = Figure(resolution = (1400, 900))
    ax = fig[1, 1] = Axis(fig)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    plot_empty_tiles(fig, ax, dx, dy, n)
    plot_occupied_tiles(fig, ax, dx, dy, n, t_occ)
    plot_model(fig, ax, m, shownvec = false, shownodes = false, showcom = false)
    display(fig)
end

function test_view2D_shadow()
    # test blocking with shadow plot
    m = model_two_rectangles_with_holes(0.02)
    # blocking with tilewalk
    vfmat = zeros(Float64, m.nelem, m.nelem)
    existing_vf!(m, vfmat)
    n = 30
    dx, dy = create_tiles(m, n)
    @time t_occ = check_tile_occupation(m, dx, dy, n)
    @time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
    # 2D plot
    fig = Figure(resolution = (1400, 900))
    ax = fig[1, 1] = Axis(fig)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    # plot_empty_tiles(fig, ax, dx, dy, n)
    # plot_occupied_tiles(fig, ax, dx, dy, n, t_occ)
    # plot_model(fig, ax, m, shownvec = true, shownodes = false, showcom = false)
    # plot_existing(fig, ax, m, 150, vfmat)
    # plot_model_shadow_1to1(fig, ax, m, vfmat, 3, 8)
    plot_model_shadow_1toAll(fig, ax, m, vfmat, 3)
    display(fig)
end

function test_tile_occ_analysis(;n = 10)
    m = model_circle_in_circle(0.8, 1.6, 0.05)
    dx, dy = create_tiles(m, n)
    @time t_occ = check_tile_occupation(m, dx, dy, n)
    # tile_occ_analysis(t_occ)
    t_max, t_min, t_mean, ratio = tile_occ_analysis(t_occ, printit = false)
end

#= ######################################################
therm2D - testing
###################################################### =# 

function test_therm2D()
    # which model
    # m = model_two_intesecting_circles(3.0, 120, 24)
    m = model_rectangle(2.0, 2.5, 0.1)
    # m = model_circle_in_circle(3.3, 5.2, 0.1)
    # vfmat calculation
    vfmat = zeros(Float64, m.nelem, m.nelem)
    existing_vf!(m, vfmat)
    n = 30
    dx, dy = create_tiles(m, n)
    t_occ = check_tile_occupation(m, dx, dy, n)
    blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
    calculating_vf!(m, vfmat, normit = true)
    # solve Qp
    epsilon = zeros(m.nelem,1)
    set_bc_part!(m, epsilon, 1, 0.9)
    set_bc_part!(m, epsilon, 2:4, 0.3)
    temp = zeros(m.nelem,1)
    set_bc_part!(m, temp, 1, 602)
    set_bc_part!(m, temp, 2, 600)
    set_bc_part!(m, temp, 3:4, 590)
    @time Qp, G = tempsolver(m, vfmat, temp, epsilon)
    Qp_parts = [sum(Qp[m.elem2par[i].first:m.elem2par[i].last,1]) for i = 1:m.npar]
    # 2D plot
    fig = Figure(resolution = (1400, 900))
    ax = fig[1, 1] = Axis(fig)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    # plot_model(fig, ax, m, shownvec = true, shownodes = false, showcom = false)
    plot_model_with_value(fig, ax, m, Qp, "Wärmestrom in W", showcbar = false)
    display(fig)
end