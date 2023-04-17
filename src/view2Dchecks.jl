
function blocking_vf_with_tiles_2elem(fig, ax, i1, i2, m::ConstModel, dx, dy, n::T2, t_occ) where T2<:Integer
    # check if element pair is blocked by others with tiles
    p1 = m.elements[i1].com
    p2 = m.elements[i2].com
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
                    p3 = m.nodes[m.elements[it].no_node1]
                    p4 = m.nodes[m.elements[it].no_node2]
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



function blocking_vf_with_tiles_2elem_tiles(fig, ax, i1, i2, m::ConstModel, dx, dy, n::T2, t_occ) where T2<:Integer
    # check if element pair is blocked by others with tiles
    p1 = m.elements[i1].com
    p2 = m.elements[i2].com
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




"""
    get_number_of_elements_in_tiles(t_occ)

"""
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