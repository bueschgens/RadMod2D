#TODO:
# -clean up code (remove unnecessary comments)
# -add docstrings
# -add more tests

"""
    get_tilewalk(p1::Point2D{T1}, p2::Point2D{T1}, n::T2, dx::T1, dy::T1)::Vector{Index2D{T2}}
    where {T1<:AbstractFloat, T2<:Integer}

Calculates the tile walk from coordinate points p1 to p2. The tile walk is a vector of
tile indices. The tile walk is calculated by the Bresenham algorithm. Returns the index
for the next tilewalk step. As (x_ind, y_ind) where x_ind and y_ind are the indices for
the next tile in x and y direction. 
"""
function get_vectors_for_tilewalk(p1::Point2D{T1}, 
            p2::Point2D{T1})::Tuple{Index2D, T1, Point2D{T1}} where T1<:AbstractFloat
    # calculate necessary data for tilewalk
    vec = p2 - p1
    vec_l = norm(vec)
    vec_norm = normit(vec)
    tol = 1E-8
    dirx::typeof(1) = 0

    if vec_norm.x > tol
        dirx = 1
    elseif vec_norm.x < ((-1) * tol)
        dirx = -1
    end

    diry::typeof(1) = 0
    if vec_norm.y > tol
        diry = 1
    elseif vec_norm.y < ((-1) * tol)
        diry = -1
    end
    return Index2D(dirx, diry), vec_l, vec_norm
end



function get_start_tile(p1::Point2D{T1}, dir::Index2D{T2}, n::T2, dx::T1, 
        dy::T1)::Index2D{T2} where {T1<:AbstractFloat, T2 <: Integer}
    # get tile of starting point
    tx_step = 1
    ty_step = 1
    tol = 1E-8
    if dir.x > 0
        tx_step += floor(Integer,(p1.x + dir.x * tol) / dx)
    elseif dir.x < 0
        tx_step += floor(Integer,(p1.x + dir.x * tol) / dx)
    else
        # wenn exakt auf einer tile Linie -> 1D walk
        # entscheidet Vorzeichen vor tol wo der Bucket liegt
        # wichtig wenn linie am Rand
        if p1.x > (dx * n - tol)
            tx_step += floor(Integer,(p1.x - tol) / dx)
        else
            tx_step += floor(Integer,(p1.x + tol) / dx)
        end
    end
    if dir.y > 0
        ty_step += floor(Integer,(p1.y + dir.y * tol) / dy)
    elseif dir.y < 0
        ty_step += floor(Integer,(p1.y + dir.y * tol) / dy)
    else
        # wenn exakt auf einer tile Linie -> 1D walk
        # entscheidet Vorzeichen vor tol wo der Bucket liegt
        # wichtig wenn linie am Rand
        if p1.y > (dy * n - tol)
            ty_step += floor(Integer,(p1.y - tol) / dy)
        else
            ty_step += floor(Integer,(p1.y + tol) / dy)
        end
    end
    return Index2D(tx_step, ty_step)
end



@inline function Sx(m::T)::T where T<:AbstractFloat
    val = sqrt(1 + m^2)
    return val
end



@inline function Sy(m::T)::T where T<:AbstractFloat
    val = sqrt(1 + (m^(-1))^2)
    return val
end



function get_lengths_to_next_tiles(p1::Point2D{T1}, p2::Point2D{T1}, dir::Index2D{T2}, 
        tile::Index2D{T2}, dx::T1, dy::T1)::Tuple{T1,T1} where {T1<:AbstractFloat, T2<:Integer}
    # get lengths in x and y direction to next tiles
    if dir.x > 0
        dx_step = dx * (tile.x)
    elseif dir.x < 0
        dx_step = dx * (tile.x - 1)
    else
        dx_step = NaN
    end
    if dir.y > 0
        dy_step = dy * (tile.y)
    elseif dir.y < 0
        dy_step = dy * (tile.y - 1)
    else
        dy_step = NaN
    end
    # println("dx: ", dx_step, "  dy: ", dy_step)
    lx = dir.x * (dx_step - p1.x) * Sx((p2.y-p1.y) / (p2.x-p1.x))
    ly = dir.y * (dy_step - p1.y) * Sy((p2.y-p1.y) / (p2.x-p1.x))
    # correction for 1D case and changing NaN into Inf
    if isnan(lx)
        lx = Inf
    end

    if isnan(ly)
        ly = Inf
    end
    # println("lx: ",lx,"  ly: ",ly)
    return lx, ly
end



function get_next_tile(tile::Index2D{T2}, dir::Index2D{T2}, lx::T1, 
        ly::T1)::Tuple{Index2D{T2},T1,Symbol} where {T1<:AbstractFloat, T2<:Integer}
    # get next tile
    # bei einem 2D Fall exakt durch die Knoten entscheidet 
    # hier < und <= welche Seite das next tile liegt
    # default von mir: < geht immer in nächste y richtung
    if abs(lx) < abs(ly)
        # println("moving in x direction")
        tile_new = tile + Index2D(dir.x, 0)
        lcurrent = lx
        color = :cyan
    else
        # println("moving in y direction")
        tile_new = tile + Index2D(0, dir.y)
        lcurrent = ly
        color = :coral
    end
    return tile_new, lcurrent, color
end


"""
    create_randomly_occupied_tiles(n::T; density = 0.1)::Matrix{Union{Vector{T},Missing}} where T<:Integer

Debug/Test funciton which randomly creates occupied tiles.
"""
function create_randomly_occupied_tiles(n::T; density = 0.1)::Matrix{Union{Vector{T},Missing}} where T<:Integer
    # create random occ tiles
    # only for testing
    t_occ = rand(n,n)
    t_occ[t_occ .> (1-density)] .= 1
    t_occ[t_occ .<= (1-density)] .= 0
    t_occ = convert.(Integer,t_occ)
    t_occ_n = Matrix{Union{Vector{T},Missing}}(missing,n,n)
    # t_occ_n[t_occ .== 1] .= Vector{T}(undef,2)
    # t_occ = rand([0, 1], n, n)
    # t_occ = zeros(Integer, n, n)
    for i = 1:n, j = 1:n
        if isapprox(t_occ[i,j], 1.0)
            t_occ_n[i,j] = [1,1]
        end
    end
    return t_occ_n
end



function get_max_steps(n::T)::T where T<:Integer
    m1 = round(Int64,sqrt(n*n+n*n)*2)
    # m2 = convert(typeof(n),m1)
    return m1
end



function tilewalk_with_check_for_occ_tiles(fig, ax, p1::Point2D{T1}, p2::Point2D{T1}, 
        dx::T1, dy::T1, n::T2, t_occ::Matrix{Union{Vector{T2},Missing}})::Nothing where {T1<:AbstractFloat, T2<:Integer}
    # do tilewalk between two points and check for occupied tiles
    dir, ltot, lvec = get_vectors_for_tilewalk(p1, p2)
    tile = get_start_tile(p1, dir, n, dx, dy)
    println("start tile number: ", tile.x, ", ", tile.y)
    max_steps = get_max_steps(n)
    if !ismissing(t_occ[tile.x, tile.y])
        println("------> occupied ")
        start_walk = false
    else
        # print start tile in color
        poly!(ax, Rect((tile.x-1) * dx, (tile.y-1) * dy, dx, dy), color = (:blue, 0.3))
        start_walk = true
    end
    if start_walk
        # do walk
        for step = 1:max_steps
            lx, ly = get_lengths_to_next_tiles(p1, p2, dir, tile, dx, dy)
            # check if end reached
            tol2 = 1E-8 # necessary because of not exact numbers
            if lx >= (ltot - tol2) && ly >= (ltot - tol2)
                break
            end
            tile, lcurrent, pcolor = get_next_tile(tile, dir, lx, ly)
            # println("lcurrent: ",lcurrent)
            # plot next point
            p = p1 + lvec * lcurrent
            scatter!(ax, Point2f0(p.x, p.y), color = pcolor, markersize = 15)
            println("next tile number: ", tile.x, ", ", tile.y)
            # check this tile
            if !ismissing(t_occ[tile.x, tile.y])
                println("------> occupied ")
                # plot hitting point
                p = p1 + lvec * lcurrent
                scatter!(ax, Point2f0(p.x, p.y), color = :blue, markersize = 15)
                break
            end
            # print tile in color
            poly!(ax, Rect((tile.x-1) * dx, (tile.y-1) * dy, dx, dy), color = (:blue, 0.3))
        end
    end
end



function tilewalk_with_return(fig, ax, p1::Point2D{T1}, p2::Point2D{T1}, dx::T1, dy::T1, 
        n::T2)::Vector{Index2D{T2}} where {T1<:AbstractFloat, T2<:Integer}
    # do tilewalk between two points and check for occupied tiles
    max_steps = get_max_steps(n)
    tile_list = Vector{Index2D{T2}}(undef,max_steps)
    dir, ltot, lvec = get_vectors_for_tilewalk(p1, p2)
    tile = get_start_tile(p1, dir, n, dx, dy)
    println("start tile number: ", tile.x, ", ", tile.y)
    hit = 1
    tile_list[hit] = Index2D(tile.x, tile.y)
    # print start tile in color
    poly!(ax, Rect((tile.x-1) * dx, (tile.y-1) * dy, dx, dy), color = (:blue, 0.3))
    # do walk
    for step = 1:max_steps
        lx, ly = get_lengths_to_next_tiles(p1, p2, dir, tile, dx, dy)
        # check if end reached
        tol2 = 1E-8 # necessary because of not exact numbers
        if lx >= (ltot - tol2) && ly >= (ltot - tol2)
            break
        end
        tile, lcurrent, pcolor = get_next_tile(tile, dir, lx, ly)
        println("next tile number: ", tile.x, ", ", tile.y)
        hit += 1
        tile_list[hit] = Index2D(tile.x, tile.y)
        # print tile in color
        poly!(ax, Rect((tile.x-1) * dx, (tile.y-1) * dy, dx, dy), color = (:blue, 0.3))
    end
    return tile_list[1:hit]
end



function tilewalk_with_return!(tile_list::Vector{Index2D{T2}}, p1::Point2D{T1}, 
        p2::Point2D{T1}, dx::T1, dy::T1, n::T2)::T2 where {T1<:AbstractFloat, T2<:Integer}
    # do tilewalk between two points and check for occupied tiles
    max_steps = get_max_steps(n)
    dir, ltot, lvec = get_vectors_for_tilewalk(p1, p2)
    tile = get_start_tile(p1, dir, n, dx, dy)
    hit = 1
    tile_list[hit] = tile
    for step = 1:max_steps # do walk
        if step == max_steps
            @warn "reached max allowed tiles in tilewalk"
        end
        lx, ly = get_lengths_to_next_tiles(p1, p2, dir, tile, dx, dy)
        # check if end reached
        tol2 = 1E-8 # necessary because of not exact numbers
        if lx >= (ltot - tol2) && ly >= (ltot - tol2)
            break
        end
        tile, nothing, nothing = get_next_tile(tile, dir, lx, ly)
        hit += 1
        tile_list[hit] = tile
    end
    return hit
end