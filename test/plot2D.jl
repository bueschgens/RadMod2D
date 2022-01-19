#= ######################################################
model2D - plotting functions
###################################################### =# 

function plot_model(fig, ax, m::Model; shownvec = false, shownodes = false, showcom = false, showleg = true, colors = to_colormap(:rainbow, m.npar), linewidth = 3.0, legpos = "bottom")
    # plot model
    elem_lines = Vector(undef,m.npar)
    plt_elem = Vector(undef,m.npar)
    plt_nodes = Vector(undef,m.npar)
    plt_com = Vector(undef,m.npar)
    plt_nvec = Vector(undef,m.npar)
    # colors = [:red, :green, :blue, :orange]
    # colors = to_colormap(:rainbow, m.npar)
    for p = 1:m.npar
        erange = m.elem2par[p].first:m.elem2par[p].last
        # plot elements
        elem_lines[p] = [Point2f0(m.nodes[m.elem[i].node1].x, m.nodes[m.elem[i].node1].y) => 
                    Point2f0(m.nodes[m.elem[i].node2].x, m.nodes[m.elem[i].node2].y) for i in erange]
        plt_elem[p] = linesegments!(ax, elem_lines[p], color = colors[p], linewidth = linewidth) 
        # plot nodes
        if shownodes == true
            local nodes_points = [Point2f0(m.nodes[m.elem[i].node1].x, m.nodes[m.elem[i].node1].y) for i in erange]
            append!(nodes_points, [Point2f0(m.nodes[m.elem[i].node2].x, m.nodes[m.elem[i].node2].y) for i in erange])
            plt_nodes[p] = scatter!(ax, nodes_points, color = colors[p], markersize = 10)
        end
        # plot com
        if showcom == true
            local com_points = [Point2f0(m.elem[i].com.x, m.elem[i].com.y) for i in erange]
            plt_com[p] = scatter!(ax, com_points, color = :grey, markersize = 10)
        end
        # plot nvec as linesegments
        if shownvec == true
            local com_points = [Point2f0(m.elem[i].com.x, m.elem[i].com.y) for i in erange]
            # local nvec_lines = [Point2f0(m.elem[i].com.x, m.elem[i].com.y) => 
            #             Point2f0(m.elem[i].com.x + m.elem[i].nvec.x, m.elem[i].com.y + m.elem[i].nvec.y) for i in erange]
            # plt_nvec[p] = linesegments!(ax, nvec_lines, color = colors[p], linewidth = linewidth)
            # plot nvec as arrows
            scale = 0.1 # default 0.3
            arrowsize = 0.01 # default 0.05
            local nvec_directions = [Point2f0(m.elem[i].nvec.x * scale, m.elem[i].nvec.y * scale) for i in erange]
            lw_arrow = linewidth * 0.5
            arrows!(ax, com_points, nvec_directions, arrowcolor = colors[p], arrowsize = arrowsize, linecolor = colors[p], linewidth = lw_arrow)
        end
    end
    # legend
    if showleg == true
        leg_string = ["part $i" for i = 1:m.npar]
        if legpos == "right"
            leg = fig[1, end+1] = Legend(fig, plt_elem, leg_string)
        elseif legpos == "bottom"
            leg = fig[end+1, 1] = Legend(fig, plt_elem, leg_string, orientation = :horizontal, tellwidth = false, tellheight = true)
        end
        # leg.padding = (5.0f0, 5.0f0, 4.0f0, 4.0f0)
        # leg.linewidth = 1.0
        # leg.patchsize = (5, 3)
        # leg.rowgap = 0.1
        # leg.colgap = 8
        # axislegend(ax, plt_elem, leg_string, position = :rt, orientation = :horizontal)
    end
end

function plot_model_elements(fig, ax, m::Model, erange; shownvec = false, shownodes = false, showcom = false, color = :cyan, linewidth = 3.0, linestyle = nothing)
    # plot only given elements
    # color = :cyan
    # plot elements
    elem_lines = [Point2f0(m.nodes[m.elem[i].node1].x, m.nodes[m.elem[i].node1].y) => 
                Point2f0(m.nodes[m.elem[i].node2].x, m.nodes[m.elem[i].node2].y) for i in erange]
    plt_elem = linesegments!(ax, elem_lines, color = color, linewidth = linewidth, linestyle = linestyle)
    # plot nodes
    if shownodes == true
        local nodes_points = [Point2f0(m.nodes[m.elem[i].node1].x, m.nodes[m.elem[i].node1].y) for i in erange]
        append!(nodes_points, [Point2f0(m.nodes[m.elem[i].node2].x, m.nodes[m.elem[i].node2].y) for i in erange])
        plt_nodes = scatter!(ax, nodes_points, color = color, markersize = 10)
    end
    # plot com
    if showcom == true
        local com_points = [Point2f0(m.elem[i].com.x, m.elem[i].com.y) for i in erange]
        plt_com = scatter!(ax, com_points, color = :grey, markersize = 10)
    end
    # plot nvec as linesegments
    if shownvec == true
        local com_points = [Point2f0(m.elem[i].com.x, m.elem[i].com.y) for i in erange]
        # plot nvec as arrows
        scale = 0.1
        arrowsize = 0.01 # default 0.05
        local nvec_directions = [Point2f0(m.elem[i].nvec.x * scale, m.elem[i].nvec.y * scale) for i in erange]
        lw_arrow = linewidth * 0.5
        arrows!(ax, com_points, nvec_directions, arrowcolor = color, arrowsize = arrowsize, linecolor = color, linewidth = lw_arrow)
    end
end

#= ######################################################
raycast2D - plotting functions
###################################################### =# 

function plot_empty_tiles(fig, ax, dx::T1, dy::T1, n::T2; linewidth = 3.0) where {T1<:AbstractFloat, T2<:Integer}
    # plot empty tiles
    # plot x lines
    t_lines = [Point2f0(i*dx, 0) => Point2f0(i*dx, n*dy) for i = 0:1:n]
    plt_elem = linesegments!(ax, t_lines, color = :black, linewidth = linewidth)
    # plot y lines
    t_lines = [Point2f0(0, i*dy) => Point2f0(n*dx, i*dy) for i = 0:1:n]
    linesegments!(ax, t_lines, color = :black, linewidth = linewidth)
end

function plot_occupied_tiles(fig, ax, dx::T1, dy::T1, n::T2, t_occ::Matrix{Union{Vector{T2},Missing}}) where {T1<:AbstractFloat, T2<:Integer}
    # plot occ tiles
    for i = 1:n, j = 1:n
        if !ismissing(t_occ[i,j])
            poly!(ax, Rect((i-1)*dx, (j-1)*dy, dx, dy), color = (:red, 0.3))
        end
    end
end

function plot_points_and_connection(fig, ax, p1::Point2D{T}, p2::Point2D{T}) where T<:AbstractFloat
    linesegments!(ax, [Point2f0(p1.x, p1.y) => Point2f0(p2.x, p2.y)], color = :blue, linewidth = 3)
    scatter!(ax, Point2f0(p1.x, p1.y), color = :green, markersize = 15)
    scatter!(ax, Point2f0(p2.x, p2.y), color = :red, markersize = 15)
end

#= ######################################################
view2D - plotting functions
###################################################### =#

function plot_existing(fig, ax, m::Model, e::T, mat; color = :yellow, linewidth = 2) where T<:Integer
    # plot a single element e with all its existing connections
    exist = Vector{T}()
    for i = 1:m.nelem
        if mat[e,i] != 0
            push!(exist, i)
        end
    end
    # plot rays
    # color = :yellow
    ray_lines = [Point2f0(m.elem[e].com.x, m.elem[e].com.y) => 
                Point2f0(m.elem[i].com.x, m.elem[i].com.y) for i in exist]
    plt_ray = linesegments!(ax, ray_lines, color = color, linewidth = linewidth)
end

function plot_model_shadow_1to1(fig, ax, m::Model, mat, source, target)
    # plot specified parts of model with their values in global colormap
    # plot source
    color_source = :yellow
    ns = (m.elem2par[source].last - m.elem2par[source].first + 1)
    source_lines = Vector{Pair{Point{2,Float32},Point{2,Float32}}}(undef,ns)
    erange_source = m.elem2par[source].first:m.elem2par[source].last
    source_lines = [Point2f0(m.nodes[m.elem[i].node1].x, m.nodes[m.elem[i].node1].y) => 
                Point2f0(m.nodes[m.elem[i].node2].x, m.nodes[m.elem[i].node2].y) for i=erange_source]
    # plot target
    nt = (m.elem2par[target].last - m.elem2par[target].first + 1)
    target_lines = Vector{Pair{Point{2,Float32},Point{2,Float32}}}(undef,nt)
    erange_target = m.elem2par[target].first:m.elem2par[target].last
    target_lines = [Point2f0(m.nodes[m.elem[i].node1].x, m.nodes[m.elem[i].node1].y) => 
                Point2f0(m.nodes[m.elem[i].node2].x, m.nodes[m.elem[i].node2].y) for i=erange_target]
    hitcount = zeros(Int64,nt)
    for i = 1:nt
        etarget = m.elem2par[target].first - 1 + i
        hitcount[i] = sum(mat[erange_source,etarget])
    end
    # all parts except source/target
    rest = collect(1:m.npar)
    filter!(x->x != target, rest)
    filter!(x->x != source, rest)
    color_rest = :blue
    nr = 0
    for p in rest
        nr += (m.elem2par[p].last - m.elem2par[p].first + 1)
    end
    rest_lines = Vector{Pair{Point{2,Float32},Point{2,Float32}}}(undef,nr)
    s1 = 0
    s2 = 0
    for p in rest
        s1 = (s2 + 1)
        s2 = (s1 + (m.elem2par[p].last - m.elem2par[p].first))
        erange_rest = m.elem2par[p].first:m.elem2par[p].last
        rest_lines[s1:s2] = [Point2f0(m.nodes[m.elem[i].node1].x, m.nodes[m.elem[i].node1].y) => 
                Point2f0(m.nodes[m.elem[i].node2].x, m.nodes[m.elem[i].node2].y) for i=erange_rest]
    end
    # plot it
    plt_source = linesegments!(ax, source_lines, color = color_source, linewidth = 3)
    plt_target = linesegments!(ax, target_lines, color = hitcount, colormap = Reverse(:Greys_3), linewidth = 3)
    plt_rest = linesegments!(ax, rest_lines, color = color_rest, linewidth = 3)
end

function plot_model_shadow_1toAll(fig, ax, m::Model, mat, source)
    # plot specified parts of model with their values in global colormap
    # plot source
    color_source = :yellow
    ns = (m.elem2par[source].last - m.elem2par[source].first + 1)
    source_lines = Vector{Pair{Point{2,Float32},Point{2,Float32}}}(undef,ns)
    erange_source = m.elem2par[source].first:m.elem2par[source].last
    source_lines = [Point2f0(m.nodes[m.elem[i].node1].x, m.nodes[m.elem[i].node1].y) => 
                Point2f0(m.nodes[m.elem[i].node2].x, m.nodes[m.elem[i].node2].y) for i=erange_source]
    # plot targets
    nt = m.nelem - ns
    target = collect(1:m.npar)
    filter!(x->x != source, target)
    target_lines = Vector{Pair{Point{2,Float32},Point{2,Float32}}}(undef,nt)
    hitcount = zeros(Int64,nt)
    s1 = 0
    s2 = 0
    for p in target
        s1 = (s2 + 1)
        s2 = (s1 + (m.elem2par[p].last - m.elem2par[p].first))
        erange_target = m.elem2par[p].first:m.elem2par[p].last
        target_lines[s1:s2] = [Point2f0(m.nodes[m.elem[i].node1].x, m.nodes[m.elem[i].node1].y) => 
                Point2f0(m.nodes[m.elem[i].node2].x, m.nodes[m.elem[i].node2].y) for i=erange_target]
        k = 0
        for j = s1:s2
            k += 1
            etarget = erange_target[k]
            hitcount[j] = sum(mat[erange_source,etarget])
        end
    end
    # plot it
    plt_source = linesegments!(ax, source_lines, color = color_source, linewidth = 5)
    plt_target = linesegments!(ax, target_lines, color = hitcount, colormap = :coolwarm, linewidth = 5)
    # colorbar
    cbar = fig[1, end+1] = Colorbar(fig, plt_target, label = "hit-intensity")
    cbar.width = 30
    cbar.height = Relative(2/3)
    valmin = minimum(hitcount)
    valmax = maximum(hitcount)
    valdelta = (valmax - valmin) / 10
    cbar.ticks = valmin:valdelta:valmax
end

#= ######################################################
therm2D - plotting functions
###################################################### =#

function plot_model_with_value(fig, ax, m::Model, val::Vector{T}, val_label; parts = 1:m.npar, cmap = :jet, linewidth = 5.0, showcbar = true, cbarpos = "right") where T<:AbstractFloat
    # plot specified parts of model with their values in global colormap
    n = 0
    for p in parts
        n += (m.elem2par[p].last - m.elem2par[p].first + 1)
    end
    pval = Vector{T}(undef,n)
    elem_lines = Vector{Pair{Point{2,Float32},Point{2,Float32}}}(undef,n)
    s1 = 0
    s2 = 0
    for p in parts
        s1 = (s2 + 1)
        s2 = (s1 + (m.elem2par[p].last - m.elem2par[p].first))
        # println(s1, " --> ", s2)
        erange = m.elem2par[p].first:m.elem2par[p].last
        pval[s1:s2] = val[erange]
        elem_lines[s1:s2] = [Point2f0(m.nodes[m.elem[i].node1].x, m.nodes[m.elem[i].node1].y) => 
                Point2f0(m.nodes[m.elem[i].node2].x, m.nodes[m.elem[i].node2].y) for i=erange]
    end
    plt_elem = linesegments!(ax, elem_lines, color = pval, colormap = cmap, linewidth = linewidth)
    # all at once
    # elem_lines = [Point2f0(m.nodes[m.elem[i].node1].x, m.nodes[m.elem[i].node1].y) => 
    #             Point2f0(m.nodes[m.elem[i].node2].x, m.nodes[m.elem[i].node2].y) for i=1:m.nelem]
    # plt_elem = linesegments!(ax, elem_lines, color = val, colormap = :viridis, linewidth = 3)
    # colorbar
    if showcbar == true
        if cbarpos == "right"
            cbar = fig[1, end+1] = Colorbar(fig, plt_elem, label = val_label)
        elseif cbarpos == "bottom"
            cbar = fig[end+1,1] = Colorbar(fig, plt_elem, label = val_label)
            cbar.vertical = false
        end
        # cbar.width = 12
        # cbar.height = Relative(2/3)
        # cbar.labelpadding = 1.0 # default 5.0
        # cbar.ticklabelpad = 1.0 # default 3.0
        # cbar.ticksize = 2.0 # default 6.0
    else
        cbar = 0
    end
    valmin = minimum(pval)
    valmax = maximum(pval)
    valdelta = (valmax - valmin) / 10
    # cbar.ticks = valmin:valdelta:valmax
    println("values between: ", valmin, " and ", valmax)
    return cbar
end

#= ######################################################
others - plotting functions
###################################################### =#

function setup_axis(ax, linewidth)
    # for springer heat and mass transfer
    ax.spinewidth = linewidth
    ax.xtickwidth = linewidth
    ax.ytickwidth = linewidth
    ax.xgridwidth = linewidth
    ax.ygridwidth = linewidth
    ax.xlabelpadding = 0.0 # default 3.0
    ax.ylabelpadding = 0.0
    ax.xticklabelpad = 0.0 # default 2.0
    ax.yticklabelpad = 1.5
    ax.xticksize = 2.0 # default 6.0
    ax.yticksize = 2.0
end