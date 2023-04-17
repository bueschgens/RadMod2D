using Pkg
Pkg.activate(".")
Pkg.instantiate()

using RadMod2D

include("./models2D.jl")


function epsilon_effective_groove()
    # infinetely long round groove
    # calculations
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
        m = model_circle_with_opening_line_centered(2.0, 360, round(Integer,case[i,2]))
        # vfmat calculation
        vfmat = zeros(Float64, m.no_elements, m.no_elements)
        existing_vf!(m, vfmat)
        n = 30
        dx, dy = get_tile_dimensions(m, n)
        @time t_occ = check_tile_occupation(m, dx, dy, n)
        @time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
        calculating_vf!(m, vfmat, normit = true)
        # epsilon effective
        eps_eff_t = get_epsilon_effective(m, vfmat, 2, epsilon_ink = epsilon_ink)
        eps_eff[:,(2*i-1)] = eps_eff_t[:,1]
        eps_eff[:,(2*i)] = eps_eff_t[:,2]
    end
    # 2D plot
    # fig = Figure(resolution = (900, 700))
    # ax = fig[1, 1] = Axis(fig)
    # ax.xlabel = "epsilon"
    # ax.ylabel = "epsilon effective"
    # ax.aspect = DataAspect()
    # for i = 1:(size(case,1))
    #     plt_lin[i] = lines!(ax, eps_eff[:,(2*i-1)], eps_eff[:,(2*i)], linewidth = 3, color = colors[i])
    # end
    # xlims!(ax, (0, 1))
    # ylims!(ax, (0, 1))
    # ax.xticks = 0:0.2:1.0
    # ax.yticks = 0:0.2:1.0
    # display(fig)
    # leg_string = ["phi = $(case[i,1])" for i = 1:size(case,1)]
    # leg = fig[1, end+1] = Legend(fig, plt_lin, leg_string)
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

function plot_epseff_groove_allinone()
    # infinetely long round groove
    # plot all cases in one figure
    case = [0.01 3.56;
            0.05 17.2;
            0.1 33.15;
            0.2 62.55;
            0.5 142.79;
            0.9999 357.19] 
    # cant calculate last position with phi = 1.0 therefore using 0.999
    m = Vector{Model}(undef,size(case,1))
    for i = 1:size(case,1)
        m[i] = model_circle_with_opening_line_centered(2.0, 360, round(Integer,case[i,2]))
    end
    # arrangement in figure
    col = 2
    arr = Matrix{Int64}(undef,size(case,1), 2)
    c_col = 0
    c_row = 1
    for i = 1:size(case,1)
        if c_col >= col
            c_col = 1
            c_row += 1
        else
            c_col += 1
        end
        arr[i,1] = c_row
        arr[i,2] = c_col
    end
    # 2D plot
    fig = Figure(resolution = (300, 300), font = "Arial", fontsize = 10)
    for i = 1:size(case,1)
        ax = fig[arr[i,1], arr[i,2]] = Axis(fig)
        linewidth = 1.0
        setup_axis!(ax, linewidth)
        ax.xlabel = "X in m"
        ax.ylabel = "Y in m"
        ax.aspect = DataAspect()
        colors = Vector{Any}(undef,18)
        colors[:] .= :black
        plot_model(fig, ax, m[i], shownvec = false, shownodes = false, showcom = false, showleg = false, colors = colors, linewidth = linewidth)
    end
    # display(fig)
    save("fig_epseff_rundrille_allinone"*filetype, fig)
end

function plot_epseff_groove()
    # infinetely long round groove
    # plot chosen case in figure
    case = [0.01 3.56;
            0.05 17.2;
            0.1 33.15;
            0.2 62.55;
            0.5 142.79;
            0.9999 357.19] 
    # cant calculate last position with phi = 1.0 therefore using 0.999
    i = 4
    m = model_circle_with_opening_line_centered(2.0, 360, round(Integer,case[i,2]))
    @show m.elem2par
    # 2D plot
    fig = Figure(resolution = (270, 270), font = "Arial", fontsize = 15)
    ax = fig[1, 1] = Axis(fig)
    linewidth = 1.0
    setup_axis!(ax, linewidth)
    ax.xlabel = "X in m"
    ax.ylabel = "Y in m"
    ax.aspect = DataAspect()
    colors = Vector{Any}(undef,18)
    colors[:] .= :black
    # plot_model(fig, ax, m, shownvec = false, shownodes = false, showcom = false, showleg = false, colors = colors, linewidth = linewidth)
    plot_model_elements(fig, ax, m, 1:297, shownvec = false, shownodes = false, showcom = false, color = :black, linewidth = linewidth)
    plot_model_elements(fig, ax, m, 298:361, shownvec = false, shownodes = false, showcom = false, color = :black, linewidth = linewidth, linestyle = :dot)
    # display(fig)
    save("fig_epseff_rundrille_neu"*filetype, fig)
end

# epsilon_effective_groove()
# plot_epseff_groove_allinone()
# plot_epseff_groove()

# include("./test/run_epseff_groove.jl")