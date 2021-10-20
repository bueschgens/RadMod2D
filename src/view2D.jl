#= ######################################################
fct for struct
###################################################### =# 

function dot(v1::Point2D{T}, v2::Point2D{T}) where T<:AbstractFloat
    d = v1.x * v2.x + v1.y * v2.y
    return d
end

#= ######################################################
inside functions
###################################################### =# 

function get_angle(v1::Point2D{T},v2::Point2D{T}) where T<:AbstractFloat
    # calculate angle between two vectors
    EPSILON = 1E-8
    cosval = dot(v1,v2) / (norm(v1) * norm(v2))
    if cosval < -1 || cosval > 1
        if cosval > (-1 - EPSILON) && cosval < -1
            cosval = -1
        elseif cosval < (1 + EPSILON) && cosval > 1
            cosval = 1
        else
            error("custom error message: cosval in fct: get_angle")
        end
    end
    phir = acos(cosval)
    #phid = phir * 180 / pi
    return phir
end

function line_segment_intersection(p1::Point2D{T}, p2::Point2D{T}, p3::Point2D{T}, p4::Point2D{T}) where T<:AbstractFloat
    # check if two 2D line segments are intersecting
    hit = false
    x1 = p1.x
    x2 = p2.x
    x3 = p3.x
    x4 = p4.x
    y1 = p1.y
    y2 = p2.y
    y3 = p3.y
    y4 = p4.y
    denom = ((y2-y1)*(x3-x4)-(x2-x1)*(y3-y4))
    if denom != 0
        t = ((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/denom
        u = ((x2-x1)*(y1-y3)-(y2-y1)*(x1-x3))/denom
        if (t >= 0.0 && t <= 1.0) && (u >= 0.0 && u <= 1.0)
            hit = true
        end
    end
    return hit
end

function line_segment_intersection_test(p1::Point2D{T}, p2::Point2D{T}, p3::Point2D{T}, p4::Point2D{T}) where T<:AbstractFloat
    # check if two 2D line segments are intersecting
    #print("intersection check")
    x1 = p1.x
    x2 = p2.x
    x3 = p3.x
    x4 = p4.x
    y1 = p1.y
    y2 = p2.y
    y3 = p3.y
    y4 = p4.y
    denom = ((y2-y1)*(x3-x4)-(x2-x1)*(y3-y4))
    #print("denominator = ",denom)
    if denom == 0
        #print("no intersection -> parallel vectors")
        return false
    else
        t = ((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/denom
        u = ((x2-x1)*(y1-y3)-(y2-y1)*(x1-x3))/denom
        #print("t = ",t,"; u = ",u)
        if t >= 0.0 && t <= 1.0
            # pointX = x1+t*(x2-x1)
            # pointY = y1+t*(y2-y1)
            if u >= 0.0 && u <= 1.0
                # pointX = x3+u*(x4-x3)
                # pointY = y3+u*(y4-y3)
                # print("intersection at ",point)
                #plt.plot((x1,x2), (y1,y2))
                #plt.plot((x3,x4), (y3,y4))
                #plt.plot(point[0],point[1],'r+')
                return true
            else
                return false
            end
        else
            return false
        end
    end
end

function is_inside_tile_BB(n1::Point2D{T}, n2::Point2D{T}, tmin_x::T, tmax_x::T, tmin_y::T, tmax_y::T) where T<:AbstractFloat
    # check if element e is inside tile with bounding box
    box1 = @SMatrix [min(n1.x, n2.x) max(n1.x, n2.x); min(n1.y, n2.y) max(n1.y, n2.y)]
    box2 = @SMatrix [tmin_x tmax_x; tmin_y tmax_y]
    return isOverlapping2D(box1, box2)
end

function isOverlapping2D(box1,box2)
    # box1 = (x:(xmin1,xmax1),y:(ymin1,ymax1))
    # box2 = (x:(xmin2,xmax2),y:(ymin2,ymax2))
    # isOverlapping2D(box1,box2)
    # = isOverlapping1D(box1.x, box2.x) and 
    #   isOverlapping1D(box1.y, box2.y)
    box1x = @view box1[1,:]
    box1y = @view box1[2,:]
    box2x = @view box2[1,:]
    box2y = @view box2[2,:]
    inside = false
    if isOverlapping1D(box1x, box2x) && isOverlapping1D(box1y, box2y)
        inside = true
    end
    return inside
end

function isOverlapping1D(box1,box2)
    # box1 = (xmin1, xmax1)
    # box2 = (xmin2, xmax2)
    # isOverlapping1D(box1,box2) 
    # = xmax1 >= xmin2 and xmax2 >= xmin1
    xmin1 = box1[1]
    xmax1 = box1[2]
    xmin2 = box2[1]
    xmax2 = box2[2]
    inside = false
    if xmax1 >= xmin2 && xmax2 >= xmin1
        inside = true
    end
    return inside
end

function get_area_of_part(m, p)
    # returns the area of part p
    e1 = m.elem2par[p].first
    e2 = m.elem2par[p].last
    # println("calculate area between elements ", e1, " and ", e2)
    a = 0.0
    for i = e1:e2
        a += m.elem[i].area
    end
    return a
end

#= ######################################################
final functions
###################################################### =# 

function existing_vf!(m::Model, mat)
    # check if elements are faced towards each other
    # angle between nvec and connector < 90°
    for i1 = 1:m.nelem, i2 = (i1+1):m.nelem
        con = m.elem[i2].com - m.elem[i1].com
        phir1 = get_angle(m.elem[i1].nvec,con)
        phir2 = get_angle(m.elem[i2].nvec,con*(-1))
        tol = 1E-8
        if phir1 < (pi/2 - tol) && phir2 < (pi/2 - tol)
            mat[i1,i2] = 1
            mat[i2,i1] = 1
        end
    end
end

function blocking_vf!(m::Model, mat)
    # check if element pairs are blocked by others with brute force method
    for i1 = 1:m.nelem, i2 = (i1+1):m.nelem
        if mat[i1,i2] == 1
            for ib in 1:m.nelem
                if ib != i1 && ib != i2
                    p1 = m.elem[i1].com
                    p2 = m.elem[i2].com
                    p3 = m.nodes[m.elem[ib].node1]
                    p4 = m.nodes[m.elem[ib].node2]
                    blocked = line_segment_intersection(p1, p2, p3, p4)
                    if blocked
                        # println(i1," and ",i2," intersected by ",ib)
                        mat[i1,i2] = 0
                        mat[i2,i1] = 0
                        break
                    end
                end
            end
        end
    end
end

function create_tiles(m::Model, n::Integer)::Tuple{Float64,Float64}
    # create tiles based on model nodes
    # get min and max x and y of nodes
    nmin, nmax = get_nodes_min_and_max(m)
    dx = (nmax.x - nmin.x) / n
    dy = (nmax.y - nmin.y) / n
    return dx, dy
end

function check_tile_occupation(m::Model, dx::T1, dy::T1, n::T2)::Matrix{Union{Vector{T2},Missing}} where {T1<:AbstractFloat, T2<:Integer}
    # tile sorting
    t_occ = Matrix{Union{Vector{T2},Missing}}(missing,n,n)
    e_in_tile = Vector{T2}(undef,m.nelem) # empty vec for appending
    for t1 = 1:n
        for t2 = 1:n
            tmin_x = dx*(t1-1)
            tmax_x = dx*t1
            tmin_y = dy*(t2-1)
            tmax_y = dy*t2
            hit = 0
            for i = 1:m.nelem
                node1 = m.nodes[m.elem[i].node1]
                node2 = m.nodes[m.elem[i].node2]
                if is_inside_tile_BB(node1, node2, tmin_x, tmax_x, tmin_y, tmax_y)
                    # println("tile (",t1,",",t2,") --> ",i)
                    hit += 1
                    e_in_tile[hit] = i
                end
            end
            if hit > 0
                t_occ[t1,t2] = e_in_tile[1:hit]
            end
        end
    end
    return t_occ
end

function check_tile_occupation_parts(m::Model, dx::T1, dy::T1, n::T2; blockparts = 1:m.npar)::Matrix{Union{Vector{T2},Missing}} where {T1<:AbstractFloat, T2<:Integer}
    # tile sorting
    t_occ = Matrix{Union{Vector{T2},Missing}}(missing,n,n)
    e_in_tile = Vector{T2}(undef,m.nelem) # empty vec for appending
    for t1 = 1:n
        for t2 = 1:n
            tmin_x = dx*(t1-1)
            tmax_x = dx*t1
            tmin_y = dy*(t2-1)
            tmax_y = dy*t2
            hit = 0
            for p in blockparts
                e1 = m.elem2par[p].first
                e2 = m.elem2par[p].last
                for i = e1:e2
                    node1 = m.nodes[m.elem[i].node1]
                    node2 = m.nodes[m.elem[i].node2]
                    if is_inside_tile_BB(node1, node2, tmin_x, tmax_x, tmin_y, tmax_y)
                        # println("tile (",t1,",",t2,") --> ",i)
                        hit += 1
                        e_in_tile[hit] = i
                    end
                end
            end
            if hit > 0
                t_occ[t1,t2] = e_in_tile[1:hit]
            end
        end
    end
    return t_occ
end

function blocking_vf_with_tiles!(m::Model, mat, dx::T1, dy::T1, n::T2, t_occ::Matrix{Union{Vector{T2},Missing}}) where {T1<:AbstractFloat, T2<:Integer}
    # check if element pairs are blocked by others with tiles
    max_steps = get_max_steps(n)
    tiles = Vector{Index2D{T2}}(undef,max_steps)
    @inbounds for i1 = 1:m.nelem, i2 = (i1+1):m.nelem
        # println("--------> checking ",i1," and ",i2)
        p1 = m.elem[i1].com
        p2 = m.elem[i2].com
        if mat[i1,i2] == 1
            # first tilewalk and get all tiles between i1 and i2
            ntiles = tilewalk_with_return!(tiles, p1, p2, dx, dy, n)
            hitten = false
            @inbounds for i = 1:ntiles
                t = tiles[i]
                # println("next tile number: ", t.x, ", ", t.y)
                if !ismissing(t_occ[t.x,t.y])
                    # println("------> occupied ")
                    for it in t_occ[t.x,t.y]
                        if it != i1 && it != i2
                            # println("------> ",i1," and ",i2, " checking with ",it)
                            p3 = m.nodes[m.elem[it].node1]
                            p4 = m.nodes[m.elem[it].node2]
                            blocked = line_segment_intersection(p1, p2, p3, p4)
                            if blocked
                                # println("----------------> ",i1," and ",i2," intersected by ",it)
                                mat[i1,i2] = 0
                                mat[i2,i1] = 0
                                hitten = true
                                break
                            end
                        end
                    end
                end
                if hitten
                    break
                end
            end
        end
    end
end

function blocking_vf_with_tiles_simplified!(m::Model, mat, dx::T1, dy::T1, n::T2, t_occ::Matrix{Union{Vector{T2},Missing}}; elem_in_t = 12, skipped_t = 2) where {T1<:AbstractFloat, T2<:Integer}
    # check if element pairs are blocked by others with tiles
    max_steps = get_max_steps(n)
    tiles = Vector{Index2D{T2}}(undef,max_steps)
    count_simp = 0
    count_all = 0
    @inbounds for i1 = 1:m.nelem, i2 = (i1+1):m.nelem
        # println("--------> checking ",i1," and ",i2)
        p1 = m.elem[i1].com
        p2 = m.elem[i2].com
        if mat[i1,i2] == 1
            # first tilewalk and get all tiles between i1 and i2
            ntiles = tilewalk_with_return!(tiles, p1, p2, dx, dy, n)
            hitten = false
            @inbounds for i = 1:ntiles
                t = tiles[i]
                # println("next tile number: ", t.x, ", ", t.y)
                if !ismissing(t_occ[t.x,t.y])
                    # println("------> occupied ")
                    count_all += 1
                    if size(t_occ[t.x,t.y],1) <= elem_in_t || i <= 1+skipped_t || i >= ntiles-skipped_t
                        for it in t_occ[t.x,t.y]
                            if it != i1 && it != i2
                                # println("------> ",i1," and ",i2, " checking with ",it)
                                p3 = m.nodes[m.elem[it].node1]
                                p4 = m.nodes[m.elem[it].node2]
                                blocked = line_segment_intersection(p1, p2, p3, p4)
                                if blocked
                                    # println("----------------> ",i1," and ",i2," intersected by ",it)
                                    mat[i1,i2] = 0
                                    mat[i2,i1] = 0
                                    hitten = true
                                    break
                                end
                            end
                        end
                    else
                        # println("----------------> SIMPLIFIED -> ",i1," and ",i2," intersected by ",it)
                        mat[i1,i2] = 0
                        mat[i2,i1] = 0
                        count_simp += 1
                    end
                end
                if hitten
                    break
                end
            end
        end
    end
    # println("counter simp/all: ", count_simp, " / ", count_all)
    println("simplification rate: ", round(count_simp/count_all*100, digits = 2), " %")
end

function calculating_vf!(m, mat; normit = false)
    # calculate viewfactors if vf is existing
    for i1 = 1:m.nelem, i2 = (i1+1):m.nelem
        if mat[i1,i2] == 1
            # viewfactor existing
            a = m.nodes[m.elem[i1].node1]
            b = m.nodes[m.elem[i1].node2]
            c = m.nodes[m.elem[i2].node1]
            d = m.nodes[m.elem[i2].node2]
            # punkte müssen richtig sortiert sein
            # nur aus python rauskopiert - nicht mehr kontrolliert
            # hier umgangen mit betrag
            # 4 Strecken berechnen -> Abstand zwischen 2 Punkten
            ac = norm(c - a)
            bd = norm(d - b)
            bc = norm(c - b)
            ad = norm(d - a)
            # println("AC: ", ac," BD: ", bd," BC: ", bc," AD: ", ad)
            # cross strings formula - Fläche berechnen und vf bestimmen
            A1 = norm(b - a)
            A2 = norm(d - c)
            mat[i1,i2] = abs((0.5 * ((ac+bd)-(bc+ad))) / A1)
            mat[i2,i1] = abs((0.5 * ((ac+bd)-(bc+ad))) / A2)
            #print(i1," -> ",i2,":",vf12,vf21)
        end
    end
    if normit
        control = sum(mat, dims=2)
        mat .= mat ./ control
    end
end

#= ######################################################
blocking test with 2 elements
###################################################### =# 

function existing_vf_2elem(m, i1, i2)
    # check if elements are faced towards each other
    # angle between nvec and connector < 90°
    exists = false
    con = m.elem[i2].com - m.elem[i1].com
    phir1 = get_angle(m.elem[i1].nvec,con)
    phir2 = get_angle(m.elem[i2].nvec,con*(-1))
    println("Angles: ", phir1, " and ", phir2)
    tol = 1E-8
    if phir1 < (pi/2 - tol) && phir2 < (pi/2 - tol)
        println("is existing between ", i1, " and ", i2)
        exists = true
    else
        println("is not existing between ", i1, " and ", i2)
    end
    return exists
end

function blocking_vf_with_tiles_2elem(fig, ax, i1, i2, m::Model, dx, dy, n::T2, t_occ) where T2<:Integer
    # check if element pair is blocked by others with tiles
    p1 = m.elem[i1].com
    p2 = m.elem[i2].com
    max_steps = get_max_steps(n)
    tiles = Vector{Index2D{T2}}(undef,max_steps)
    # first tilewalk and get all tiles between i1 and i2
    ntiles = tilewalk_with_return!(tiles, p1, p2, dx, dy, n)
    println(tiles[1:ntiles])
    hitten = false
    for i = 1:ntiles
        t = tiles[i]
        println("next tile number: ", t.x, ", ", t.y)
        # print tile in color
        poly!(ax, Rect((t.x-1) * dx, (t.y-1) * dy, dx, dy), color = (:grey, 0.6))
        if !ismissing(t_occ[t.x,t.y])
            println("------> occupied ")
            for it in t_occ[t.x,t.y]
                if it != i1 && it != i2
                    println("------> ",i1," and ",i2, " checking with ",it)
                    p3 = m.nodes[m.elem[it].node1]
                    p4 = m.nodes[m.elem[it].node2]
                    blocked = line_segment_intersection(p1, p2, p3, p4)
                    if blocked
                        println("----------------> ",i1," and ",i2," intersected by ",it)
                        hitten = true
                        break
                    end
                end
            end
        end
        if hitten
            break
        end
    end
end

function blocking_vf_with_tiles_2elem_tiles(fig, ax, i1, i2, m::Model, dx, dy, n::T2, t_occ) where T2<:Integer
    # check if element pair is blocked by others with tiles
    p1 = m.elem[i1].com
    p2 = m.elem[i2].com
    max_steps = get_max_steps(n)
    tiles = Vector{Index2D{T2}}(undef,max_steps)
    # first tilewalk and get all tiles between i1 and i2
    ntiles = tilewalk_with_return!(tiles, p1, p2, dx, dy, n)
    # println(tiles[1:ntiles])
    for i = 1:ntiles
        t = tiles[i]
        println("next tile number: ", t.x, ", ", t.y)
        # print tile in color
        poly!(ax, Rect((t.x-1) * dx, (t.y-1) * dy, dx, dy), color = (:grey, 0.6))
    end
end

#= ######################################################
post processing of viewfactors
###################################################### =# 

function vfmat_to_parts(m, mat; normit = false)
    # get viewfactor matrix parts out of viewfactor matrix elements
    matp = zeros(Float64, m.npar, m.npar)
    areap = [get_area_of_part(m, i) for i = 1:m.npar]
    for p1 = 1:m.npar, p2 = 1:m.npar
        vfsplit = mat[m.elem2par[p1].first:m.elem2par[p1].last, m.elem2par[p2].first:m.elem2par[p2].last]
        acol = [m.elem[i].area for i = m.elem2par[p1].first:m.elem2par[p1].last]
        asplit = repeat(acol, 1, (m.elem2par[p2].last-m.elem2par[p2].first+1))
        vfsplit = vfsplit .* asplit
        matp[p1,p2] = sum(vfsplit) / areap[p1]
    end
    if normit
        controlp = sum(matp, dims=2)
        matp .= matp ./ controlp
    end
    return matp
end

function save_vfmat_part_to_csv(mat, path, name; control = false)
    # save vfmat parts to parts as csv
    # option of saving with control sum
    if control == true
        control = sum(mat, dims=2)
        mat = hcat(mat, control)
        nparts = size(mat,1)
        header = Matrix{String}(undef, 1, nparts+1)
        for i = 1:nparts
            header[1,i] = "part"*string(i)
        end
        header[1,end] = "control"
    else
        header = []
    end
    filename = path*"/"*name*".txt"
    open(filename, "w") do io
        writedlm(io, header)
        writedlm(io, mat)
    end
end

function get_number_of_elements_in_tiles(t_occ)
    nx = size(t_occ,1)
    ny = size(t_occ,2)
    t_num = zeros(Int64, nx, ny)
    for i = 1:nx, j = 1:ny
        if ismissing(t_occ[i,j])
            t_num[i,j] = 0
        else
            t_num[i,j] = size(t_occ[i,j],1)
        end
    end
    return t_num
end

function tile_occ_analysis(t_occ; printit = true)
    # analysis of tile occupation matrix
    t_num = get_number_of_elements_in_tiles(t_occ)
    t_num_v = vec(t_num)
    t_num_v_occ = t_num_v[t_num_v[:] .> 0]
    t_min = minimum(t_num_v_occ)
    t_max = maximum(t_num_v_occ)
    n_all = size(t_num_v,1)
    n_occ = size(t_num_v_occ,1)
    t_mean = sum(t_num_v_occ) / n_all
    t_mean = round(t_mean, digits=1)
    t_mean_occ = sum(t_num_v_occ) / n_occ
    t_mean_occ = round(t_mean_occ, digits=1)
    if printit
        println("Tile occupation with elements:")
        println("    Tiles occupied: ", n_occ, " / ", n_all)
        println("    Max: ", t_max, " // Min: ", t_min)
        println("    Mean: ", t_mean, " // Mean_occ: ", t_mean_occ)
    else
        return t_max, t_min, t_mean_occ, n_occ/n_all
    end 
end