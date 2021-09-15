#= ######################################################
includes
###################################################### =# 

using DelimitedFiles

# using GLMakie
using CairoMakie
# CairoMakie.activate!(type = "svg")
# filetype = ".svg"
CairoMakie.activate!(type = "png")
filetype = ".png"

include("./epseff2D.jl")
include("./models2D.jl")
include("./plot2D.jl")

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

function setup_axis_epseff(ax)
    linewidth = 2.0
    setup_axis(ax, linewidth)
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
    linewidth = setup_axis_epseff(ax)
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
    vfmat = zeros(Float64, m.nelem, m.nelem)
    existing_vf!(m, vfmat)
    dx, dy = create_tiles(m, n)
    @time t_occ = check_tile_occupation(m, dx, dy, n)
    @time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
    calculating_vf!(m, vfmat, normit = true)
    return vfmat
end


#= ######################################################
cases
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
        m = model_circle_with_opening_line(2.0, 360, round(Integer,case[i,2]))
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
    # rechteckige Nut
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
    # Cosinus Wave
    case = [0.0 1.0;
            0.5 1.0;
            1.0 1.0;
            1.5 1.0;
            2.0 1.0;
            2.5 1.0;
            3.0 1.0]
    epsilon_ink = 0.01
    eps_eff = Matrix{Float64}(undef,round(Integer,1.0/epsilon_ink),2*size(case,1))
    for i = 2:size(case,1)
        m = model_cosinus(case[i,1], case[i,2], 100)
        vfmat = calc_vfmat(m)
        eps_eff_t = get_epsilon_effective(m, vfmat, 2, epsilon_ink = epsilon_ink)
        eps_eff[:,(2*i-1)] = eps_eff_t[:,1]
        eps_eff[:,(2*i)] = eps_eff_t[:,2]
    end
    plot_epseff(case, eps_eff, "fig_epseff_cosinus")
end