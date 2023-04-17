#= ###################################################### Geoemtry representation for 2D

Models are represented by a list of nodes and a list of line elements.

LineElements only refer to model (part) node list for computation efficiency.
LineElement normal vectors (norm_vec), center of mass (com) and length are precomputed
at geometry generation for efficiency.

For geometry generation, the model is mutable. After geometry generation, the model is
made immutable for efficiency within the vief factor and blocking calculations.
Therefore, two different structs are used.
###################################################### =# 

abstract type AbstractLineElement end
abstract type AbstractModel end


mutable struct MutableLineElement{T1<:Integer, T2<:AbstractFloat} <: AbstractLineElement
    no_node1::T1
    no_node2::T1
    com::Point2D{T2}
    norm_vec::Point2D{T2}
    length::T2

    # function MutableLineElement(no_node1::T1, no_node2::T1, com::Point2D{T2}, norm_vec::Point2D{T2}, length::T2) where {T1<:Integer, T2<:AbstractFloat}
    #     new{T1, T2}(no_node1, no_node2, com, norm_vec, length)
    # end
end


struct ConstLineElement{T1<:Integer, T2<:AbstractFloat} <: AbstractLineElement
    no_node1::T1
    no_node2::T1
    com::Point2D{T2}
    norm_vec::Point2D{T2}
    length::T2

    # function ConstLineElement(no_node1::T1, no_node2::T1, com::Point2D{T2}, norm_vec::Point2D{T2}, length::T2) where {T1<:Integer, T2<:AbstractFloat}
    #     new{T1,T2}(no_node1, no_node2, com, norm_vec, length)
    # end
end


mutable struct Part{T1<:Integer, T2<:AbstractFloat}
    nodes::Vector{Point2D{T2}}
    elements::Vector{MutableLineElement{T1,T2}}
    no_nodes::T1
    no_elements::T1
    name::String

    # function Part(nodes::Vector{Point2D{T2}}, elements::Vector{MutableLineElement{T1,T2}}, no_nodes::T1, no_elements::T1, name::String) where {T1<:Integer, T2<:AbstractFloat}
    #     new{T1,T2}(nodes, elements, no_nodes, no_elements, name)
    # end
end


struct ElementToPartAssign{T<:Integer}
    first::T
    last::T

    # function ElementToPartAssign(first::T, last::T) where {T<:Integer}
    #     new{T}(first, last)
    # end
end


mutable struct MutableModel{T1<:Integer, T2<:AbstractFloat} <: AbstractModel
    nodes::Vector{Point2D{T2}}
    elements::Vector{MutableLineElement{T1,T2}}
    elem2par::Vector{ElementToPartAssign{T1}}
    no_nodes::T1
    no_elements::T1
    no_parts::T1
    part_names::Vector{String}

    function MutableModel(int_val::T1,float_val::T2)::MutableModel{T1, T2} where {T1<:Integer, T2<:AbstractFloat}
        return new{T1,T2}(Vector{Point2D{T2}}(), Vector{MutableLineElement{T1,T2}}(), 
                Vector{ElementToPartAssign{T1}}(), 0, 0, 0, [""])
    end
end


struct ConstModel{T1<:Integer, T2<:AbstractFloat} <: AbstractModel
    nodes::Vector{Point2D{T2}}
    elements::Vector{ConstLineElement{T1,T2}}
    elem2par::Vector{ElementToPartAssign{T1}}
    no_nodes::T1
    no_elements::T1
    no_parts::T1
    part_names::Vector{String}

    function ConstModel(m::MutableModel{T1, T2})::ConstModel{T1, T2} where {T1<:Integer, T2<:AbstractFloat}
        nodes = m.nodes
        elements = Vector{ConstLineElement{Int64,Float64}}(undef,m.no_elements)
        for i = 1:m.no_elements
            i1 = m.elements[i].no_node1
            i2 = m.elements[i].no_node2
            com = m.elements[i].com
            norm_vec = m.elements[i].norm_vec
            len = m.elements[i].length
            elements[i] = ConstLineElement(i1, i2, com, norm_vec, len)
        end
        elem2par = m.elem2par
        no_nodes = m.no_nodes
        no_elements = m.no_elements
        no_parts = m.no_parts
        part_names = m.part_names

        return new{T1, T2}(nodes, elements, elem2par, no_nodes, no_elements, no_parts, 
                part_names)
    end
end



"""
    create_empty_model()

Initializes an empty model.
"""
function create_empty_model()
    return MutableModel(1, 1.0)
end



"""
    add_offset!(p::Part, offset::Integer)

Adds offset `offset` to all node numbers in part `p`.
"""
function add_offset!(p::Part, offset::Integer)
    for i = 1:p.no_elements
        p.elements[i].no_node1 += offset
        p.elements[i].no_node2 += offset
    end
end



"""
    add!(m::MutableModel, p::Part)

Adds part ``p` to mutable model `m`.
"""
function add!(m::MutableModel, p::Part)
    if m.no_parts == 0
        m.nodes = p.nodes
        m.elements = p.elements
        m.elem2par = [ElementToPartAssign(1, p.no_elements)]
        m.no_nodes = p.no_nodes
        m.no_elements = p.no_elements
        m.no_parts = 1
        m.part_names = [p.name]
    else
        append!(m.nodes, p.nodes)
        add_offset!(p, m.no_nodes)
        append!(m.elements, p.elements)
        push!(m.elem2par, ElementToPartAssign(m.no_elements+1, m.no_elements+p.no_elements))
        m.no_nodes += p.no_nodes
        m.no_elements += p.no_elements
        m.no_parts += 1
         
        if p.name == nothing || p.name == "" 
            push!(m.part_names, "part $(m.no_parts)")
        else
            push!(m.part_names, p.name)
        end
    end
end



"""
add!(p1::Part, p2::Part)

Adds part ``p1` to mutable part `p2`.
"""
function add!(p1::Part, p2::Part)
        append!(p1.nodes, p2.nodes)
        add_offset!(p2, p1.no_nodes)
        append!(p1.elements, p2.elements)
        p1.no_nodes += p2.no_nodes
        p1.no_elements += p2.no_elements
         
        if p1.name == nothing || p1.name == "" 
            p1.name = p2.name
        else
            # not change part name!
        end
end



"""
    get_min_max_coordinates(m::MutableModel)

#Arguments
- `m::Model`: Model

#Returns
- `Tuple{Point2D,Point2D}` with minimum and maximum coordinates of model nodes
"""
function get_min_max_coordinates(m::AbstractModel)::Tuple{Point2D,Point2D}
    # get minimum and maximum value of nodes
    # as Vector of Point2D
    xmin = minimum(m.nodes[i].x for i = 1:m.no_nodes)
    xmax = maximum(m.nodes[i].x for i = 1:m.no_nodes)
    ymin = minimum(m.nodes[i].y for i = 1:m.no_nodes)
    ymax = maximum(m.nodes[i].y for i = 1:m.no_nodes)
    return Point2D(xmin, ymin), Point2D(xmax, ymax)
end



"""
    offset_model!(m::MutableModel)

Modifies m by offsetting all nodes and elements to only positiv coordinates.

#Arguments
- `m::MutableModel`
"""
function offset_model!(m::MutableModel)
    # offset all nodes to only positiv coords
    p_xy_min, p_xy_max = get_min_max_coordinates(m)
    for i = 1:m.no_nodes
        m.nodes[i] = m.nodes[i] - p_xy_min
    end
    for i = 1:m.no_elements
        m.elements[i].com = m.elements[i].com - p_xy_min
    end
end


"""
    element_analysis(m::AbstractModel; printit = true)

Quick analysis of model elements.

#Arguments
- `m::MutableModel`: Model

#Keyword Arguments
- `printit::Bool`: If true, print results to console

#Returns
- `e_length_max::Float64`: Maximum element length
- `e_length_min::Float64`: Minimum element length
- `e_length_mean::Float64`: Mean element length
"""
function element_analysis(m::AbstractModel; printit = true)::Tuple{Float64,Float64,Float64}
    # analysis of elements
    length = [m.elements[i].length for i = m.elem2par[1].first:m.elem2par[end].last]
    e_length_mean = sum(length) / m.no_elements
    # e_length_mean = round(e_length_mean, digits=2)
    e_length_min = minimum(length)
    e_length_max = maximum(length)
    if printit
        println("Element length:")
        println("    Mean: ", e_length_mean)
        println("    Min: ", e_length_min)
        println("    Max: ", e_length_max)
    end
    
    return e_length_max, e_length_min, e_length_mean
end