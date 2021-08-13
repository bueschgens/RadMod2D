#= ######################################################
includes
###################################################### =# 

include("./models2D.jl")

#= ######################################################
inside functions
###################################################### =# 

function get_epsilon_effective(m, mat, p_hole; epsilon_ink = 0.1)
    # calculate epsilon effective
    # hole is part with epsilon effective / encl is enclosure
    p_encl = collect(1:m.npar)
    filter!(x->x != p_hole, p_encl)
    temp_encl = 0
    temp_hole = 1000 # Kelvin
    temp = zeros(m.nelem,1)
    set_bc_part!(m, temp, p_hole, temp_hole)
    set_bc_part!(m, temp, p_encl, temp_encl)
    epsilon_hole = 1.0
    n_epsilon_encl = round(Integer,1.0/epsilon_ink)
    epsilon_encl = Matrix{Float64}(undef,n_epsilon_encl,2)
    epsilon_encl[:,1] = [epsilon_ink*i for i =1:n_epsilon_encl]
    epsilon = zeros(m.nelem,1)
    set_bc_part!(m, epsilon, p_hole, epsilon_hole)
    for i = 1:n_epsilon_encl
        set_bc_part!(m, epsilon, p_encl, epsilon_encl[i,1])
        # calculation of Qp
        Qp, Ga = tempsolver(m, mat, temp, epsilon)
        area_hole = get_area_of_part(m, p_hole)
        G_hole = sum(Ga[m.elem2par[p_hole].first:m.elem2par[p_hole].last])
        Qp_hole = sum(Qp[m.elem2par[p_hole].first:m.elem2par[p_hole].last])
        # calculation of epsilon effective
        epsilon_eff = 1.0 - (G_hole / (area_hole * sigma * (temp_hole^4 - temp_encl^4)))
        epsilon_encl[i,2] = epsilon_eff
    end
    return epsilon_encl
end

#= ######################################################
final functions
###################################################### =# 

function epsilon_effective_rundrille()
    # unendlich lange Rundrille
    case = [0.01 3.56;
            0.05 17.2;
            0.1 33.15;
            0.2 62.55;
            0.5 142.79;
            0.9999 357.19] 
    # cant calculate last position with phi = 1.0 therefore using 0.999
    # colors = [:red, :green, :blue, :orange, :cyan, :purple1]
    colors = to_colormap(:rainbow, size(case,1))
    plt_lin = Vector(undef,(size(case,1)))
    epsilon_ink = 0.01
    eps_eff = Matrix{Float64}(undef,round(Integer,1.0/epsilon_ink),2*size(case,1))
    for i = 1:(size(case,1))
        m = model_circle_with_opening_line(2.0, 360, round(Integer,case[i,2]))
        # vfmat calculation
        vfmat = zeros(Float64, m.nelem, m.nelem)
        existing_vf!(m, vfmat)
        n = 30
        dx, dy = create_tiles(m, n)
        @time t_occ = check_tile_occupation(m, dx, dy, n)
        @time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
        calculating_vf!(m, vfmat, normit = true)
        # epsilon effective
        eps_eff_t = get_epsilon_effective(m, vfmat, 2, epsilon_ink = epsilon_ink)
        eps_eff[:,(2*i-1)] = eps_eff_t[:,1]
        eps_eff[:,(2*i)] = eps_eff_t[:,2]
    end
    # 2D plot
    fig = Figure(resolution = (900, 700))
    ax = fig[1, 1] = Axis(fig)
    ax.xlabel = "epsilon"
    ax.ylabel = "epsilon effective"
    ax.aspect = DataAspect()
    for i = 1:(size(case,1))
        plt_lin[i] = lines!(ax, eps_eff[:,(2*i-1)], eps_eff[:,(2*i)], linewidth = 3, color = colors[i])
    end
    xlims!(ax, (0, 1))
    ylims!(ax, (0, 1))
    ax.xticks = 0:0.2:1.0
    ax.yticks = 0:0.2:1.0
    display(fig)
    # legend
    leg_string = ["phi = $(case[i,1])" for i = 1:size(case,1)]
    leg = fig[1, end+1] = Legend(fig, plt_lin, leg_string)
    # save results 
    header1 = Matrix{String}(undef,1,2*size(case,1))
    header2 = Matrix{String}(undef,1,2*size(case,1))
    for i = 1:size(case,1)
        # replace(string(case[i]),"." => "n")
        header1[1,(2*i-1)] = "phi = "*string(case[i])
        header1[1,(2*i)] = "phi = "*string(case[i])
        header2[1,(2*i-1)] = "eps"
        header2[1,(2*i)] = "eps_eff"
    end
    open("data_eps_eff_rundrille.txt", "w") do io
        writedlm(io, header1)
        writedlm(io, header2)
        writedlm(io, eps_eff)
    end
end

function epsilon_effective_keilspalt()
    # unendlich langer Keilspalt
    case = [1E-3 1;
            0.224 25.9;
            0.316 36.8;
            0.447 53.1;
            0.592 72.6;
            0.707 90
            1.0 0.0]
    # cant calc with 0 first position -> therefore 1 deg
    # cant calculate last position with phi = 1.0
    # adding it manually at the end
    colors = to_colormap(:Dark2_8, size(case,1))
    plt_lin = Vector(undef,(size(case,1)))
    # 2D plot
    fig = Figure(resolution = (900, 700))
    ax = fig[1, 1] = Axis(fig)
    ax.xlabel = "epsilon"
    ax.ylabel = "epsilon effective"
    ax.aspect = DataAspect()
    for i = 1:(size(case,1)-1)
        m = model_isosceles_triangle(3.0, case[i,2], 0.01)
        # vfmat calculation
        vfmat = zeros(Float64, m.nelem, m.nelem)
        existing_vf!(m, vfmat)
        n = 30
        dx, dy = create_tiles(m, n)
        @time t_occ = check_tile_occupation(m, dx, dy, n)
        @time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
        calculating_vf!(m, vfmat, normit = true)
        # epsilon effective
        eps_eff = get_epsilon_effective(m, vfmat, 3, epsilon_ink = 0.01)
        plt_lin[i] = lines!(ax, eps_eff[:,1], eps_eff[:,2], linewidth = 3, color = colors[i])
    end
    # add line for phi = 1.0 manually
    plt_lin[end] = lines!(ax, [0,1], [0,1], linewidth = 3, color = colors[end])
    xlims!(ax, (0, 1))
    ylims!(ax, (0, 1))
    ax.xticks = 0:0.2:1.0
    ax.yticks = 0:0.2:1.0
    display(fig)
    # legend
    leg_string = ["phi = $(case[i,1])" for i = 1:size(case,1)]
    leg = fig[1, end+1] = Legend(fig, plt_lin, leg_string)
end

function epsilon_effective_rechteckigenut()
    # rechteckige Nut
    h = 2.0
    case = Matrix{Float64}(undef,10,2) # L/h, L
    for i = 1:10
        case[i,1] = i
        case[i,2] = i * h
    end
    colors = to_colormap(:Dark2_8, size(case,1))
    plt_lin = Vector(undef,(size(case,1)))
    # 2D plot
    fig = Figure(resolution = (900, 700))
    ax = fig[1, 1] = Axis(fig)
    ax.xlabel = "epsilon"
    ax.ylabel = "epsilon effective"
    ax.aspect = DataAspect()
    for i = 1:(size(case,1))
        m = model_rectangle(h, case[i,2], 0.05)
        # vfmat calculation
        vfmat = zeros(Float64, m.nelem, m.nelem)
        existing_vf!(m, vfmat)
        n = 30
        dx, dy = create_tiles(m, n)
        @time t_occ = check_tile_occupation(m, dx, dy, n)
        @time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
        calculating_vf!(m, vfmat, normit = true)
        # epsilon effective
        eps_eff = get_epsilon_effective(m, vfmat, 3, epsilon_ink = 0.01)
        plt_lin[i] = lines!(ax, eps_eff[:,1], eps_eff[:,2], linewidth = 3, color = colors[i])
    end
    xlims!(ax, (0, 1))
    ylims!(ax, (0, 1))
    ax.xticks = 0:0.2:1.0
    ax.yticks = 0:0.2:1.0
    display(fig)
    # legend
    leg_string = ["L/h = $(case[i,1])" for i = 1:size(case,1)]
    leg = fig[1, end+1] = Legend(fig, plt_lin, leg_string)
end

function epsilon_effective_winkelstueck()
    # Winkelstueck
    case = [0.707 1.0 1.0;
            0.850 1.0 0.2]
    colors = to_colormap(:rainbow, size(case,1)+1)
    plt_lin = Vector(undef,(size(case,1)))
    # 2D plot
    fig = Figure(resolution = (900, 700))
    ax = fig[1, 1] = Axis(fig)
    ax.xlabel = "epsilon"
    ax.ylabel = "epsilon effective"
    ax.aspect = DataAspect()
    for i = 1:(size(case,1))
        m = model_right_triangle(case[i,2], case[i,3], 0.01)
        # vfmat calculation
        vfmat = zeros(Float64, m.nelem, m.nelem)
        existing_vf!(m, vfmat)
        n = 30
        dx, dy = create_tiles(m, n)
        @time t_occ = check_tile_occupation(m, dx, dy, n)
        @time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
        calculating_vf!(m, vfmat, normit = true)
        # epsilon effective
        eps_eff = get_epsilon_effective(m, vfmat, 3, epsilon_ink = 0.01)
        plt_lin[i] = lines!(ax, eps_eff[:,1], eps_eff[:,2], linewidth = 3, color = colors[i])
    end
    # add line for 1.0 manually - not in legend
    lines!(ax, [0,1], [0,1], linewidth = 3, color = colors[end])
    xlims!(ax, (0, 1))
    ylims!(ax, (0, 1))
    ax.xticks = 0:0.2:1.0
    ax.yticks = 0:0.2:1.0
    display(fig)
    # legend
    leg_string = ["phi = $(case[i,1])" for i = 1:size(case,1)]
    leg = fig[1, end+1] = Legend(fig, plt_lin, leg_string)
end

function epsilon_effective_cosinus()
    # Cosinus Wave
    case = [0.0 1.0;
            0.5 1.0;
            1.0 1.0;
            1.5 1.0;
            2.0 1.0;
            2.5 1.0;
            3.0 1.0]
    colors = to_colormap(:rainbow, size(case,1)+1)
    plt_lin = Vector(undef,(size(case,1)))
    # 2D plot
    fig = Figure(resolution = (900, 700))
    ax = fig[1, 1] = Axis(fig)
    ax.xlabel = "epsilon"
    ax.ylabel = "epsilon effective"
    ax.aspect = DataAspect()
    # add line for a = 0.0 manually
    plt_lin[1] = lines!(ax, [0,1], [0,1], linewidth = 3, color = colors[1])
    for i = 2:(size(case,1))
        m = model_cosinus(case[i,1], case[i,2], 100)
        # vfmat calculation
        vfmat = zeros(Float64, m.nelem, m.nelem)
        existing_vf!(m, vfmat)
        n = 30
        dx, dy = create_tiles(m, n)
        @time t_occ = check_tile_occupation(m, dx, dy, n)
        @time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
        calculating_vf!(m, vfmat, normit = true)
        # epsilon effective
        eps_eff = get_epsilon_effective(m, vfmat, 2, epsilon_ink = 0.01)
        plt_lin[i] = lines!(ax, eps_eff[:,1], eps_eff[:,2], linewidth = 3, color = colors[i])
    end
    xlims!(ax, (0, 1))
    ylims!(ax, (0, 1))
    ax.xticks = 0:0.2:1.0
    ax.yticks = 0:0.2:1.0
    display(fig)
    # legend
    leg_string = ["Amp = $(case[i,1])" for i = 1:size(case,1)]
    leg = fig[1, end+1] = Legend(fig, plt_lin, leg_string)
end

function epsilon_effective_trapez()
    # Trapez
    a = 2.0
    b = 3.0
    case = collect(1:10)
    colors = to_colormap(:rainbow, size(case,1))
    plt_lin = Vector(undef,(size(case,1)))
    # 2D plot
    fig = Figure(resolution = (900, 700))
    ax = fig[1, 1] = Axis(fig)
    ax.xlabel = "epsilon"
    ax.ylabel = "epsilon effective"
    ax.aspect = DataAspect()
    for i = 1:(size(case,1))
        m = model_trapez(a, b, case[i,1], 0.05)
        # vfmat calculation
        vfmat = zeros(Float64, m.nelem, m.nelem)
        existing_vf!(m, vfmat)
        n = 30
        dx, dy = create_tiles(m, n)
        @time t_occ = check_tile_occupation(m, dx, dy, n)
        @time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
        calculating_vf!(m, vfmat, normit = true)
        # epsilon effective
        eps_eff = get_epsilon_effective(m, vfmat, 4, epsilon_ink = 0.01)
        plt_lin[i] = lines!(ax, eps_eff[:,1], eps_eff[:,2], linewidth = 3, color = colors[i])
    end
    xlims!(ax, (0, 1))
    ylims!(ax, (0, 1))
    ax.xticks = 0:0.2:1.0
    ax.yticks = 0:0.2:1.0
    display(fig)
    # legend
    leg_string = ["h = $(case[i,1])" for i = 1:size(case,1)]
    leg = fig[1, end+1] = Legend(fig, plt_lin, leg_string)
end