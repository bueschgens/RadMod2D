#= ######################################################
model structs
###################################################### =# 

mutable struct ElementRaw{T1<:Integer, T2<:AbstractFloat}
    node1::T1
    node2::T1
    com::Point2D{T2}
    nvec::Point2D{T2}
    area::T2
end

struct Element{T1<:Integer, T2<:AbstractFloat}
    node1::T1
    node2::T1
    com::Point2D{T2}
    nvec::Point2D{T2}
    area::T2
end

struct Part{T1<:Integer, T2<:AbstractFloat}
    nodes::Vector{Point2D{T2}}
    elem::Vector{ElementRaw{T1,T2}}
    nnodes::T1
    nelem::T1
end

struct ElementAssign{T<:Integer}
    first::T
    last::T
end

mutable struct ModelRaw{T1<:Integer, T2<:AbstractFloat}
    nodes::Vector{Point2D{T2}}
    elem::Vector{ElementRaw{T1,T2}}
    elem2par::Vector{ElementAssign{T1}}
    nnodes::T1
    nelem::T1
    npar::T1
end

struct Model{T1<:Integer, T2<:AbstractFloat}
    nodes::Vector{Point2D{T2}}
    elem::Vector{Element{T1,T2}}
    elem2par::Vector{ElementAssign{T1}}
    nnodes::T1
    nelem::T1
    npar::T1
end

#= ######################################################
Functions for working with structs
###################################################### =# 

function create_empty_model()
    return ModelRaw(Vector{Point2D{Float64}}(), Vector{ElementRaw{Int64,Float64}}(), Vector{ElementAssign{Int64}}(), 0, 0, 0)
end

function norm(v::Point2D)
    l = sqrt(v.x^2 + v.y^2)
    return l
end

function normit(v::Point2D)::Point2D
    vn = v / norm(v)
    return vn
end

function add_offset!(p::Part, offset::Integer)
    for i = 1:p.nelem
        p.elem[i].node1 += offset
        p.elem[i].node2 += offset
    end
end

function get_com(p1::Point2D{T}, p2::Point2D{T})::Point2D{T} where T<:AbstractFloat
    # calculate center of line (center of mass)
    c = p1 + ((p2 - p1) / 2)
    return c
end

function get_nvec(p1::Point2D{T}, p2::Point2D{T}, dir::String)::Point2D{T} where T<:AbstractFloat
    # calculate normal vector
    n = Point2D((-1) * (p2.y - p1.y), (p2.x - p1.x))
    (dir == "neg") && (n = n * (-1.0))
    n_norm = normit(n)
    return n_norm
end

function get_length(p1::Point2D, p2::Point2D)
    # calculate length of line element
    # is used as area here
    l = sqrt((p1.x-p2.x)^2 + (p1.y-p2.y)^2)
    return l
end

function add!(m::ModelRaw, p::Part)
    append!(m.nodes, p.nodes)
    add_offset!(p, m.nnodes)
    append!(m.elem, p.elem)
    push!(m.elem2par, ElementAssign(m.nelem+1, m.nelem+p.nelem))
    m.nnodes += p.nnodes
    m.nelem += p.nelem
    m.npar += 1
end

function get_nodes_min_and_max(m)
    # get minimum and maximum value of nodes
    # as Vector of Point2D
    xmin = minimum(m.nodes[i].x for i = 1:m.nnodes)
    xmax = maximum(m.nodes[i].x for i = 1:m.nnodes)
    ymin = minimum(m.nodes[i].y for i = 1:m.nnodes)
    ymax = maximum(m.nodes[i].y for i = 1:m.nnodes)
    return Point2D(xmin, ymin), Point2D(xmax, ymax)
end

function offset_model!(m::ModelRaw)
    # offset all nodes to only positiv coords
    nmin, nmax = get_nodes_min_and_max(m)
    for i = 1:m.nnodes
        m.nodes[i] = m.nodes[i] - nmin
    end
    for i = 1:m.nelem
        m.elem[i].com = m.elem[i].com - nmin
    end
end

function make_model_immutable(m::ModelRaw)::Model
    nodes = m.nodes
    elem = Vector{Element{Int64,Float64}}(undef,m.nelem)
    for i = 1:m.nelem
        i1 = m.elem[i].node1
        i2 = m.elem[i].node2
        c = m.elem[i].com
        nv = m.elem[i].nvec
        a = m.elem[i].area
        elem[i] = Element(i1, i2, c, nv, a)
    end
    elem2par = m.elem2par
    nnodes = m.nnodes
    nelem = m.nelem
    npar = m.npar
    return Model(nodes, elem, elem2par, nnodes, nelem, npar)
end

#= ######################################################
Define discretized geometries
###################################################### =# 

function edge(p1::Point2D, p2::Point2D; seed::T = 10, dir::String = "pos") where T<:Integer
    """
    creating edge between two points
    # println("edge between ", p1, " and ",p2)
    """
    dx = (p2.x - p1.x) / seed
    dy = (p2.y - p1.y) / seed
    nnodes = seed+1
    nodes = Vector{Point2D{Float64}}(undef,nnodes)
    for i = 1:nnodes
        nodes[i] = Point2D(p1.x + (i-1) * dx, p1.y + (i-1) * dy)
    end
    # nodes = [Point2D(p1.x + (i-1) * dx, p1.y + (i-1) * dy) for i = 1:nnodes]
    nelements = seed
    elements = Vector{ElementRaw{Int64,Float64}}(undef,nelements)
    for i = 1:nelements
        i1 = i
        i2 = i + 1
        c = get_com(nodes[i1], nodes[i2])
        nv = get_nvec(nodes[i1], nodes[i2], dir)
        a = get_length(nodes[i1], nodes[i2])
        elements[i] = ElementRaw(i1, i2, c, nv, a)
    end
    return Part(nodes, elements, nnodes, nelements)
end

function rectangle(x::T1, y::T1, c::Point2D; seedx::T2 = 10, seedy::T2 = 10, dir::String = "pos") where {T1<:Real, T2<:Integer}
    # creating rectangle with center
    p1 = c + Point2D(-0.5*x, -0.5*y)
    p2 = c + Point2D(-0.5*x, 0.5*y)
    p3 = c + Point2D(0.5*x, 0.5*y)
    p4 = c + Point2D(0.5*x, -0.5*y)
    dx = (p3.x - p2.x) / seedx
    dy = (p2.y - p1.y) / seedy
    nnodes = 2*seedx + 2*seedy
    nodes = Vector{Point2D{Float64}}(undef,nnodes)
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
    nelements = nnodes
    elements = Vector{ElementRaw{Int64,Float64}}(undef,nelements)
    for i = 1:nelements-1
        i1 = i
        i2 = i + 1
        c = get_com(nodes[i1], nodes[i2])
        nv = get_nvec(nodes[i1], nodes[i2], dir)
        a = get_length(nodes[i1], nodes[i2])
        elements[i] = ElementRaw(i1, i2, c, nv, a)
    end
    i1 = nnodes
    i2 = 1
    c = get_com(nodes[i1], nodes[i2])
    nv = get_nvec(nodes[i1], nodes[i2], dir)
    a = get_length(nodes[i1], nodes[i2])
    elements[end] = ElementRaw(i1, i2, c, nv, a)
    return Part(nodes, elements, nnodes, nelements)
end

function circle(d::T1, c::Point2D; seed::T2 = 12, dir::String = "pos") where {T1<:Real, T2<:Integer}
    # creating circle based on diameter and center
    r = d/2
    phi = 2 * pi / seed
    nnodes = seed
    nodes = Vector{Point2D{Float64}}(undef,nnodes)
    for i = 1:nnodes
        nodes[i] = Point2D(c.x + r * cos(i*phi), c.y + r * sin(i*phi))
    end
    nelements = seed
    elements = Vector{ElementRaw{Int64,Float64}}(undef,nelements)
    for i = 1:nelements-1
        i1 = i
        i2 = i + 1
        c = get_com(nodes[i1], nodes[i2])
        nv = get_nvec(nodes[i1], nodes[i2], dir)
        a = get_length(nodes[i1], nodes[i2])
        elements[i] = ElementRaw(i1, i2, c, nv, a)
    end
    i1 = nnodes
    i2 = 1
    c = get_com(nodes[i1], nodes[i2])
    nv = get_nvec(nodes[i1], nodes[i2], dir)
    a = get_length(nodes[i1], nodes[i2])
    elements[end] = ElementRaw(i1, i2, c, nv, a)
    return Part(nodes, elements, nnodes, nelements)
end

function circle_open(d::T1, c::Point2D; seed::T2 = 12, leftout::T2 = 2, open::String = "right", dir::String = "pos") where {T1<:Real, T2 <: Integer}
    # creating circle based on diameter and center which is open at one side
    r = d/2
    phi = 2 * pi / seed
    # use half of leftout as start
    phio = 2 * pi / seed * (leftout/2) # phi offset
    if open == "left"
        phio = phio + 1 * pi
    end
    nnodes = seed-leftout+1
    nodes = Vector{Point2D{Float64}}(undef,nnodes)
    for i = 1:nnodes
        nodes[i] = Point2D(c.x + r * cos((i-1) * phi + phio), c.y + r * sin((i-1) * phi + phio))
    end
    nelements = seed - leftout
    elements = Vector{ElementRaw{Int64,Float64}}(undef,nelements)
    for i = 1:nelements
        i1 = i
        i2 = i + 1
        c = get_com(nodes[i1], nodes[i2])
        nv = get_nvec(nodes[i1], nodes[i2], dir)
        a = get_length(nodes[i1], nodes[i2])
        elements[i] = ElementRaw(i1, i2, c, nv, a)
    end
    return Part(nodes, elements, nnodes, nelements)
end

function cosinus(a::T1, b::T1, p1::Point2D; seed::T2 = 30, dir::String = "pos") where {T1<:Real, T2<:Integer}
    # creating cosinus wave
    nnodes = seed+1
    nodes = Vector{Point2D{Float64}}(undef,nnodes)
    for i = 1:nnodes
        nodes[i] = Point2D(p1.x + b*(i-1)*2*pi/(seed), p1.y + a*cos((i-1)*2*pi/(seed)))
    end
    nelements = seed
    elements = Vector{ElementRaw{Int64,Float64}}(undef,nelements)
    for i = 1:nelements
        i1 = i
        i2 = i + 1
        c = get_com(nodes[i1], nodes[i2])
        nv = get_nvec(nodes[i1], nodes[i2], dir)
        a = get_length(nodes[i1], nodes[i2])
        elements[i] = ElementRaw(i1, i2, c, nv, a)
    end
    return Part(nodes, elements, nnodes, nelements)
end

#= ######################################################
Other functions
###################################################### =# 

function element_analysis(m::Model; printit = true)
    # analysis of elements
    area = [m.elem[i].area for i = m.elem2par[1].first:m.elem2par[end].last]
    e_length_mean = sum(area) / m.nelem
    # e_length_mean = round(e_length_mean, digits=2)
    e_length_min = minimum(area)
    e_length_max = maximum(area)
    if printit
        println("Element length:")
        println("    Mean: ", e_length_mean)
        println("    Min: ", e_length_min)
        println("    Max: ", e_length_max)
    else
        return t_max, t_min, t_mean_occ, n_occ/n_all
    end 
end