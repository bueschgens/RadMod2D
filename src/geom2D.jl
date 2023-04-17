
"""
    edge(p1::Point2D, p2::Point2D; seed::T = 10, dir::Symbol=:pos) where T<:Integer

Creating edge between two points. The Edge is defined by the number of nodes (seed) and
the direction of the edge (dir).

#Arguments
- `p1::Point2D`: First point of the edge
- `p2::Point2D`: Second point of the edge

#Keyword Arguments
- `seed::T = 10`: Number of nodes of the edge
- `dir::Symbol=:pos`: Direction of the edge. "pos" means that the edge is defined by the
    vector p2 - p1. "neg" means that the edge is defined by the vector p1 - p2.
- `name::String="Edge"`: Name of the edge / resulting part

#Returns
- Part defining the edge
"""
function edge(p1::Point2D{T2}, p2::Point2D{T2}; seed::T1 = 10, dir::Symbol=:pos, 
        name::String="Edge") where {T1<:Integer, T2<:AbstractFloat}

    dx = (p2.x - p1.x) / seed
    dy = (p2.y - p1.y) / seed
    no_nodes = seed+1
    nodes = Vector{Point2D{T2}}(undef,no_nodes)
    for i = 1:no_nodes
        nodes[i] = Point2D(p1.x + (i-1) * dx, p1.y + (i-1) * dy)
    end
    # nodes = [Point2D(p1.x + (i-1) * dx, p1.y + (i-1) * dy) for i = 1:no_nodes]
    no_elements = seed
    elements = Vector{MutableLineElement{T1,T2}}(undef,no_elements)
    for i = 1:no_elements
        i1 = i
        i2 = i + 1
        c = get_com(nodes[i1], nodes[i2])
        nv = get_norm_vec(nodes[i1], nodes[i2], dir)
        a = get_length(nodes[i1], nodes[i2])
        elements[i] = MutableLineElement(i1, i2, c, nv, a)
    end
    return Part(nodes, elements, no_nodes, no_elements, name)
end



"""
    rectangle(x::T1, y::T1, c::Point2D; seedx::T2 = 10, seedy::T2 = 10, 
            dir::Symbol=:pos, name="Rect") where {T1<:Real, T2<:Integer}

Creating rectangle with center. The rectangle is defined by the length in x and y and the center point c.

#Arguments
- `x::T1`: Length in x direction
- `y::T1`: Length in y direction
- `c::Point2D`: Center point of the rectangle

#Keyword Arguments
- `seedx::T2 = 10`: Number of nodes in x direction
- `seedy::T2 = 10`: Number of nodes in y direction
- `dir::Symbol=:pos`: Direction of the edge. "pos" means that the edge is defined by the
    vector p2 - p1. "neg" means that the edge is defined by the vector p1 - p2.
- `name::String="Rect"`: Name of the rectangle / resulting part

#Returns
- Part defining the rectangle
"""
function rectangle(x::T1, y::T1, c::Point2D; seedx::T2 = 10, seedy::T2 = 10, 
        dir::Symbol = :pos, name = "Rectangle") where {T1<:Real, T2<:Integer}

    # creating rectangle with center
    p1 = c + Point2D(-0.5*x, -0.5*y)
    p2 = c + Point2D(-0.5*x, 0.5*y)
    p3 = c + Point2D(0.5*x, 0.5*y)
    p4 = c + Point2D(0.5*x, -0.5*y)
    dx = (p3.x - p2.x) / seedx
    dy = (p2.y - p1.y) / seedy
    no_nodes = 2*seedx + 2*seedy
    nodes = Vector{Point2D{Float64}}(undef, no_nodes)
    for i = 1:seedy
        nodes[i] = Point2D(p1.x, p1.y + (i-1) * dy)
    end
    offset = seedy
    for i = 1:seedx
        nodes[offset + i] = Point2D(p2.x + (i-1) * dx, p2.y)
    end
    offset = seedy + seedx
    for i = 1:seedy
        nodes[offset + i] = Point2D(p3.x, p3.y - (i-1) * dy)
    end
    offset = seedy + seedx + seedy
    for i = 1:seedx
        nodes[offset + i] = Point2D(p4.x - (i-1) * dx, p4.y)
    end
    no_elements = no_nodes
    elements = Vector{MutableLineElement{Int64,Float64}}(undef, no_elements)
    for i = 1:no_elements-1
        i1 = i
        i2 = i + 1
        c = get_com(nodes[i1], nodes[i2])
        nv = get_norm_vec(nodes[i1], nodes[i2], dir)
        a = get_length(nodes[i1], nodes[i2])
        elements[i] = MutableLineElement(i1, i2, c, nv, a)
    end
    i1 = no_nodes
    i2 = 1
    c = get_com(nodes[i1], nodes[i2])
    nv = get_norm_vec(nodes[i1], nodes[i2], dir)
    a = get_length(nodes[i1], nodes[i2])
    elements[end] = MutableLineElement(i1, i2, c, nv, a)
    return Part(nodes, elements, no_nodes, no_elements, name)
end



"""
    circle(d::T1, c::Point2D; seed::T2 = 12, dir::Symbol=:pos, 
        name::String = "Circle") where {T1<:Real, T2<:Integer}

Creating circle based on diameter `d` and center `c`.

#Arguments
- `d::T1`: Diameter of the circle
- `c::Point2D`: Center point of the circle

#Keyword Arguments
- `seed = 12`: Number of nodes on the circle
- `dir::Symbol=:pos`: Direction of the edge. "pos" means that the edge is defined by the
    vector p2 - p1. "neg" means that the edge is defined by the vector p1 - p2.
- `name::String="Circle"`: Name of the circle / resulting part

#Returns 
- Part defining the circle
"""
function circle(d::T1, c::Point2D; seed::T2 = 12, dir::Symbol = :pos,
         name::String = "Circle") where {T1<:Real, T2<:Integer}
    # creating circle based on diameter and center
    r = d/2
    phi = 2 * pi / seed
    no_nodes = seed
    nodes = Vector{Point2D{Float64}}(undef, no_nodes)
    for i = 1:no_nodes
        nodes[i] = Point2D(c.x + r * cos(i*phi), c.y + r * sin(i*phi))
    end
    no_elements = seed
    elements = Vector{MutableLineElement{Int64,Float64}}(undef, no_elements)
    for i = 1:no_elements-1
        i1 = i
        i2 = i + 1
        c = get_com(nodes[i1], nodes[i2])
        nv = get_norm_vec(nodes[i1], nodes[i2], dir)
        a = get_length(nodes[i1], nodes[i2])
        elements[i] = MutableLineElement(i1, i2, c, nv, a)
    end
    i1 = no_nodes
    i2 = 1
    c = get_com(nodes[i1], nodes[i2])
    nv = get_norm_vec(nodes[i1], nodes[i2], dir)
    a = get_length(nodes[i1], nodes[i2])
    elements[end] = MutableLineElement(i1, i2, c, nv, a)
    return Part(nodes, elements, no_nodes, no_elements, name)
end



"""
    circle_open(d::T1, c::Point2D; seed::T2 = 12, leftout::T2 = 2, open::String = "right", dir::Symbol=:pos) where {T1<:Real, T2 <: Integer}

Creating circle based on diameter `d` and center `c` which is open at one side. The side is defined by `open` and the number of elements left out is defined by `leftout`.

#Arguments
- `d::T1`: Diameter of the circle
- `c::Point2D`: Center point of the circle

#Keyword Arguments
- `seed = 12`: Number of nodes on the circle
- `leftout = 2`: Number of elements left out
- `open = "right"`: Side which is open. "right" or "left"
- `dir::Symbol=:pos`: Direction of the edge. "pos" means that the edge is defined by the
    vector p2 - p1. "neg" means that the edge is defined by the vector p1 - p2.
- `name::String="Circle"`: Name of the circle / resulting part


#Returns
- Part defining the open circle
"""
function circle_open(d::T1, c::Point2D; seed::T2 = 12, leftout::T2 = 2, 
        open::String = "right", dir::Symbol = :pos, 
        name::String = "Open Circle") where {T1<:Real, T2 <: Integer}
    # creating circle based on diameter and center which is open at one side
    r = d/2
    phi = 2 * pi / seed
    # use half of leftout as start
    phio = 2 * pi / seed * (leftout/2) # phi offset
    if open == "left"
        phio = phio + 1 * pi
    end
    no_nodes = seed-leftout+1
    nodes = Vector{Point2D{Float64}}(undef,no_nodes)
    for i = 1:no_nodes
        nodes[i] = Point2D(c.x + r * cos((i-1) * phi + phio), c.y + r * sin((i-1) * phi + phio))
    end
    no_elements = seed - leftout
    elements = Vector{MutableLineElement{Int64,Float64}}(undef,no_elements)
    for i = 1:no_elements
        i1 = i
        i2 = i + 1
        c = get_com(nodes[i1], nodes[i2])
        nv = get_norm_vec(nodes[i1], nodes[i2], dir)
        a = get_length(nodes[i1], nodes[i2])
        elements[i] = MutableLineElement(i1, i2, c, nv, a)
    end
    return Part(nodes, elements, no_nodes, no_elements, name)
end



"""
    cosinus(a::T1, b::T1, p1::Point2D; seed::T2 = 30, dir::Symbol=:pos) where {T1<:Real, T2<:Integer}

Creating cosinus wave based on amplitude `a`, wavelength `b` and start point `p1`.

#Arguments
- `a::T1`: Amplitude of the cosinus wave
- `b::T1`: Wavelength of the cosinus wave
- `p1::Point2D`: Start point of the cosinus wave

#Keyword Arguments
- `seed = 30`: Number of nodes on the cosinus wave
- `dir::Symbol=:pos`: Direction of the edge. "pos" means that the edge is defined by the
    vector p2 - p1. "neg" means that the edge is defined by the vector p1 - p2.
- `name::String="Cosinus"`: Name of the cosinus wave / resulting part

#Returns
- Part defining the cosinus wave
"""
function cosinus(a::T1, b::T1, p1::Point2D; seed::T2 = 30, dir::Symbol = :pos, 
        name::String = "Cosinus Wave") where {T1<:Real, T2<:Integer}
    # creating cosinus wave
    no_nodes = seed+1
    nodes = Vector{Point2D{Float64}}(undef, no_nodes)
    for i = 1:no_nodes
        nodes[i] = Point2D(p1.x + b*(i-1)*2*pi/(seed), p1.y + a*cos((i-1)*2*pi/(seed)))
    end
    no_elements = seed
    elements = Vector{MutableLineElement{Int64,Float64}}(undef, no_elements)
    for i = 1:no_elements
        i1 = i
        i2 = i + 1
        c = get_com(nodes[i1], nodes[i2])
        nv = get_norm_vec(nodes[i1], nodes[i2], dir)
        a = get_length(nodes[i1], nodes[i2])
        elements[i] = MutableLineElement(i1, i2, c, nv, a)
    end
    return Part(nodes, elements, no_nodes, no_elements, name)
end
