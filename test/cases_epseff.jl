#= ######################################################
includes
###################################################### =# 

using DelimitedFiles


include("./epseff2D.jl")
include("./models2D.jl")


#= ######################################################
inside functions
###################################################### =# 

function setup_legend(leg)
    leg.rowgap = -4 # default: 3
    leg.groupgap = 8 # default: 16
    leg.patchlabelgap = 5 # default: 5
    leg.patchsize = (20.0f0, 20.0f0) # default: (20.0f0, 20.0f0)
    leg.padding = (3.0f0, 3.0f0, 1.0f0, 1.0f0) # default: (10.0f0, 10.0f0, 8.0f0, 8.0f0)
end

function setup_axis!_epseff(ax)
    linewidth = 2.0
    setup_axis!(ax, linewidth)
    ax.xlabel = "epsilon"
    ax.ylabel = "epsilon effective"
    ax.aspect = DataAspect()
    xlims!(ax, (0, 1))
    ylims!(ax, (0, 1))
    ax.xticks = 0:0.2:1.0
    ax.yticks = 0:0.2:1.0
    return linewidth
end

function plot_epseff(case, val, filename)
    # 2D plot
    fig = Figure(resolution = (300, 300), font = "Arial", fontsize = 12)
    ax = fig[1, 1] = Axis(fig)
    linewidth = setup_axis!_epseff(ax)
    # colors = [:red, :green, :blue, :orange, :cyan, :purple1]
    # colors = to_colormap(:rainbow, size(case,1))
    colors = to_colormap(:Dark2_8, size(case,1))
    plt_lin = Vector(undef,(size(case,1)))
    for i = 1:size(case,1)
        plt_lin[i] = lines!(ax, val[:,(2*i-1)], val[:,(2*i)], linewidth = linewidth, color = colors[i])
    end
    # legend
    leg_string = ["phi = $(case[i,1])" for i = 1:size(case,1)]
    leg = axislegend(ax, plt_lin, leg_string, position = :rb, orientation = :vertical)
    setup_legend(leg)
    save(filename*filetype, fig, px_per_unit = 8)
end

function calc_vfmat(m; n = 30)
    vfmat = zeros(Float64, m.no_elements, m.no_elements)
    existing_vf!(m, vfmat)
    dx, dy = get_tile_dimensions(m, n)
    @time t_occ = check_tile_occupation(m, dx, dy, n)
    @time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
    calculating_vf!(m, vfmat, normit = true)
    return vfmat
end

function plot_cases(case, m, arr, filename; lim = (nothing, nothing), refcase = 1)
    fig = Figure(resolution = (700, 700), font = "Arial", fontsize = 12)
    for i = 1:size(case,1)
        ax = fig[arr[i,1], arr[i,2]] = Axis(fig)
        linewidth = 2.0
        setup_axis!(ax, linewidth)
        ax.title = "phi = $(case[i,1])"
        ax.xlabel = "X in m"
        ax.ylabel = "Y in m"
        ax.limits = lim
        ax.aspect = DataAspect()
        colors = Vector{Any}(undef,18)
        colors[:] .= :black
        plot_model(fig, ax, m[refcase], showleg = false, colors = colors, linewidth = 0)
        plot_model(fig, ax, m[i], showleg = false, colors = colors, linewidth = linewidth)
    end
    save(filename*filetype, fig, px_per_unit = 8)
end

function set_arrangement(n_cases; col = 2)
    # arrangement in figure
    arr = Matrix{Int64}(undef, n_cases, 2)
    c_col = 0
    c_row = 1
    for i = 1:n_cases
        if c_col >= col
            c_col = 1
            c_row += 1
        else
            c_col += 1
        end
        arr[i,1] = c_row
        arr[i,2] = c_col
    end
    return arr
end

#= ######################################################
cases calculation
###################################################### =# 

function epsilon_effective_round_groove()
    # infinetely long round groove
    case = [0.01 3.56;
            0.05 17.2;
            0.1 33.15;
            0.2 62.55;
            0.5 142.79;
            0.9999 357.19] # cant calculate phi = 1.0
    epsilon_ink = 0.01
    eps_eff = Matrix{Float64}(undef,round(Integer,1.0/epsilon_ink),2*size(case,1))
    for i = 1:size(case,1)
        m = model_circle_with_opening_line_centered(2.0, 360, round(Integer,case[i,2]))
        vfmat = calc_vfmat(m)
        eps_eff_t = get_epsilon_effective(m, vfmat, 2, epsilon_ink = epsilon_ink)
        eps_eff[:,(2*i-1)] = eps_eff_t[:,1]
        eps_eff[:,(2*i)] = eps_eff_t[:,2]
    end
    plot_epseff(case, eps_eff, "fig_epseff_round_groove")
end

function epsilon_effective_v_groove()
    # infinetely long V groove
    case = [1E-3 0.1; # cant calc with phi = 0
            0.224 25.9;
            0.316 36.8;
            0.447 53.1;
            0.592 72.6;
            0.707 90;
            0.9999 178.4] # cant calculate phi = 1.0
    epsilon_ink = 0.01
    eps_eff = Matrix{Float64}(undef,round(Integer,1.0/epsilon_ink),2*size(case,1))
    for i = 1:size(case,1)
        m = model_isosceles_triangle(3.0, case[i,2], 0.01)
        vfmat = calc_vfmat(m)
        eps_eff_t = get_epsilon_effective(m, vfmat, 3, epsilon_ink = epsilon_ink)
        eps_eff[:,(2*i-1)] = eps_eff_t[:,1]
        eps_eff[:,(2*i)] = eps_eff_t[:,2]
    end
    plot_epseff(case, eps_eff, "fig_epseff_v_groove")
end

function epsilon_effective_square_groove()
    # square groove
    h = 2.0
    case = Matrix{Float64}(undef,10,2) # L/h, L
    for i = 1:10
        case[i,1] = i
        case[i,2] = i * h
    end
    epsilon_ink = 0.01
    eps_eff = Matrix{Float64}(undef,round(Integer,1.0/epsilon_ink),2*size(case,1))
    for i = 1:size(case,1)
        m = model_rectangle(h, case[i,2], 0.05)
        vfmat = calc_vfmat(m)
        eps_eff_t = get_epsilon_effective(m, vfmat, 3, epsilon_ink = epsilon_ink)
        eps_eff[:,(2*i-1)] = eps_eff_t[:,1]
        eps_eff[:,(2*i)] = eps_eff_t[:,2]
    end
    plot_epseff(case, eps_eff, "fig_epseff_square_groove")
end

function epsilon_effective_elbow_piece()
    # elbow piece
    case = [0.707 1.0 1.0;
            0.850 1.0 0.2]
    epsilon_ink = 0.01
    eps_eff = Matrix{Float64}(undef,round(Integer,1.0/epsilon_ink),2*size(case,1))
    for i = 1:(size(case,1))
        m = model_right_triangle(case[i,2], case[i,3], 0.01)
        vfmat = calc_vfmat(m)
        eps_eff_t = get_epsilon_effective(m, vfmat, 3, epsilon_ink = epsilon_ink)
        eps_eff[:,(2*i-1)] = eps_eff_t[:,1]
        eps_eff[:,(2*i)] = eps_eff_t[:,2]
    end
    plot_epseff(case, eps_eff, "fig_epseff_elbow_piece")
end

function epsilon_effective_trapez()
    # trapez
    a = 2.0
    b = 3.0
    case = collect(1:10)
    epsilon_ink = 0.01
    eps_eff = Matrix{Float64}(undef,round(Integer,1.0/epsilon_ink),2*size(case,1))
    for i = 1:(size(case,1))
        m = model_trapez(a, b, case[i,1], 0.05)
        vfmat = calc_vfmat(m)
        eps_eff_t = get_epsilon_effective(m, vfmat, 4, epsilon_ink = epsilon_ink)
        eps_eff[:,(2*i-1)] = eps_eff_t[:,1]
        eps_eff[:,(2*i)] = eps_eff_t[:,2]
    end
    plot_epseff(case, eps_eff, "fig_epseff_trapez")
end

function epsilon_effective_cosinus()
    # cosinus wave
    case = [0.5 1.0;
            1.0 1.0;
            1.5 1.0;
            2.0 1.0;
            2.5 1.0;
            3.0 1.0]
    epsilon_ink = 0.01
    eps_eff = Matrix{Float64}(undef,round(Integer,1.0/epsilon_ink),2*size(case,1))
    for i = 1:size(case,1)
        m = model_cosinus(case[i,1], case[i,2], 100)
        vfmat = calc_vfmat(m)
        eps_eff_t = get_epsilon_effective(m, vfmat, 2, epsilon_ink = epsilon_ink)
        eps_eff[:,(2*i-1)] = eps_eff_t[:,1]
        eps_eff[:,(2*i)] = eps_eff_t[:,2]
    end
    plot_epseff(case, eps_eff, "fig_epseff_cosinus")
end

#= ######################################################
cases plot
###################################################### =# 

function plot_epseff_round_groove()
    # infinetely long round groove
    # plot all cases in one figure
    case = [0.01 3.56;
            0.05 17.2;
            0.1 33.15;
            0.2 62.55;
            0.5 142.79;
            0.9999 357.19] # cant calculate phi = 1.0
    m = Vector{Model}(undef,size(case,1))
    for i = 1:size(case,1)
        m[i] = model_circle_with_opening_line_centered(2.0, 360, round(Integer,case[i,2]))
    end
    arr = set_arrangement(size(case,1), col = 2)
    plot_cases(case, m, arr, "fig_plot_round_groove")
end

function plot_epseff_v_groove()
    # infinetely long V groove
    # plot all cases in one figure
    case = [1E-3 0.1; # cant calc with phi = 0
            0.224 25.9;
            0.316 36.8;
            0.447 53.1;
            0.592 72.6;
            0.707 90;
            0.9999 178.4] # cant calculate phi = 1.0
    m = Vector{Model}(undef,size(case,1))
    for i = 1:size(case,1)
        m[i] = model_isosceles_triangle(3.0, case[i,2], 0.01)
    end
    arr = set_arrangement(size(case,1), col = 2)
    plot_cases(case, m, arr, "fig_plot_v_groove", refcase = 6)
end

function plot_epseff_square_groove()
    # square groove
    # plot all cases in one figure
    h = 2.0
    case = Matrix{Float64}(undef,10,2) # L/h, L
    for i = 1:10
        case[i,1] = i
        case[i,2] = i * h
    end
    m = Vector{Model}(undef,size(case,1))
    for i = 1:size(case,1)
        m[i] = model_rectangle(h, case[i,2], 0.05)
    end
    arr = set_arrangement(size(case,1), col = 5)
    plot_cases(case, m, arr, "fig_plot_square_groove", refcase = 8)
end

function plot_epseff_elbow_piece()
    # elbow piece
    # plot all cases in one figure
    case = [0.707 1.0 1.0;
            0.850 1.0 0.2]
    m = Vector{Model}(undef,size(case,1))
    for i = 1:size(case,1)
        m[i] = model_right_triangle(case[i,2], case[i,3], 0.01)
    end
    arr = set_arrangement(size(case,1), col = 2)
    plot_cases(case, m, arr, "fig_plot_elbow_piece", refcase = 1)
end

function plot_epseff_trapez()
    # trapez
    # plot all cases in one figure
    a = 2.0
    b = 3.0
    case = collect(1:10)
    m = Vector{Model}(undef,size(case,1))
    for i = 1:size(case,1)
        m[i] = model_trapez(a, b, case[i,1], 0.05)
    end
    arr = set_arrangement(size(case,1), col = 2)
    plot_cases(case, m, arr, "fig_plot_trapez", refcase = 1)
end

function plot_epseff_cosinus()
    # cosinus wave
    # plot all cases in one figure
    case = [0.5 1.0;
            1.0 1.0;
            1.5 1.0;
            2.0 1.0;
            2.5 1.0;
            3.0 1.0]
    m = Vector{Model}(undef,size(case,1))
    for i = 1:size(case,1)
        m[i] = model_cosinus(case[i,1], case[i,2], 100)
    end
    arr = set_arrangement(size(case,1), col = 2)
    plot_cases(case, m, arr, "fig_plot_cosinus", refcase = 1)
end