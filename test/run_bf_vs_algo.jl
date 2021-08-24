using Revise

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

function calc_view2D_algo_n_study(;elemsize = 0.05)
    # calculating blocking with tiles 
    # for different tile numbers from 1 to 200 tiles per dimension
    # choose model
    # m = model_circles_in_circle_rand_full(0.1, 1.8, elemsize)
    # m = model_circles_in_circle_rand_half(0.1, 1.8, elemsize)
    # m = model_circles_in_circle_rand_quarter(0.1, 1.8, elemsize)
    m = model_circles_in_circle_cross(0.1, 1.8, 9, elemsize)
    n_study = [1,2,3,4,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,120,140,160,180,200]
    times = zeros(size(n_study,1),2)
    analy_tiles = zeros(size(n_study,1),4)
    for i = 1:size(n_study,1)
        n = n_study[i]
        println("my algorithm running with ", n)
        vfmat = zeros(Float64, m.nelem, m.nelem)
        existing_vf!(m, vfmat)
        dx, dy = create_tiles(m, n)
        stats_tiles = @timed t_occ = check_tile_occupation(m, dx, dy, n)
        stats_calc = @timed blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
        # calculating_vf!(m, vfmat, normit = false)
        times[i,1] = stats_tiles.time
        times[i,2] = stats_calc.time
        t_max, t_min, t_mean, ratio = tile_occ_analysis(t_occ, printit = false)
        analy_tiles[i,1] = t_max
        analy_tiles[i,2] = t_min
        analy_tiles[i,3] = t_mean
        analy_tiles[i,4] = ratio
    end
    # save results 
    header = ["n" "time_tiles" "time_calc" "t_max" "t_min" "t_mean" "ratio"]
    open("file_algo_n_study.txt", "w") do io
        writedlm(io, header)
        writedlm(io, [n_study times analy_tiles])
    end
end

function fig_n_study_mesh2D_allinone()
    m = Matrix{Model}(undef,2,2)
    m[1,1] = model_circles_in_circle_rand_full(0.1, 1.8, 0.005)
    m[1,2] = model_circles_in_circle_rand_half(0.1, 1.8, 0.005)
    m[2,1] = model_circles_in_circle_rand_quarter(0.1, 1.8, 0.005)
    m[2,2] = model_circles_in_circle_cross(0.1, 1.8, 9, 0.005)
    # 2D plot
    fig = Figure(resolution = (300, 300), font = "Arial", fontsize = 10)
    for i = 1:2
        for j = 1:2
            ax = fig[i, j] = Axis(fig)
            linewidth = 1.0
            setup_axis(ax, linewidth)
            ax.xlabel = "X in m"
            ax.ylabel = "Y in m"
            ax.aspect = DataAspect()
            colors = Vector{Any}(undef,18)
            colors[:] .= :black
            plot_model(fig, ax, m[i,j], shownvec = false, shownodes = false, showcom = false, showleg = false, colors = colors, linewidth = linewidth)
        end
    end
    # display(fig)
    save("fig_cyl_in_cyl_n_study_allinone"*filetype, fig)
end

# calc_view2D_algo_n_study(elemsize = 0.0005)
# fig_n_study_mesh2D_allinone()

# include("./test/run_bf_vs_algo.jl")