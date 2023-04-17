################################################################
# Code part of FARADS2D                                        #
# by dominik.bueschgens                                        #
# March 2021                                                   #
################################################################


abstract type Pair2D <: AbstractVector{Real} end
abstract type Vector2D <: Pair2D end

struct Point2D{T<:AbstractFloat} <: Vector2D
    x::T
    y::T
end

Point2D(x::Real, y::Real) = Point2D(promote(x,y)...)
Point2D(x::Integer, y::Integer) = Point2D(convert(typeof(1.0),x), convert(typeof(1.0),y))

struct Index2D{T<:Integer} <: Pair2D
    x::T
    y::T
end

import Base: size, getindex, unsafe_getindex, setindex!, unsafe_setindex!, +, -, *, /

size(v::Pair2D) = (2, )
# looks like using perf if used as vector
getindex(v::Pair2D, i::Int) = i < 3 ? (i == 1 ? v.x : v.y) : error("max index = 2")

unsafe_getindex(v::Pair2D, i::Int) = (i == 1 ? v.x : v.y)


function setindex!(v::Pair2D, a::Real, i::Integer) 
    if i < 3
        if i == 1
            v.x = a
        else
            v.y = a
        end
    else
        error("max index = 2")
    end
    return v
end



function unsafe_setindex!(v::Pair2D, a::Real, i::Integer) 
    if i == 1
        v.x = a
    else
        v.y = a
    end

    return v
end

+(a::T, b::T) where T <: Pair2D = T(a.x+b.x,a.y+b.y)
-(a::T, b::T) where T <: Pair2D = T(a.x-b.x,a.y-b.y)
*(a::T1, b::T2) where {T1 <: Pair2D, T2 <: Real} = T1(a.x*b,a.y*b)
/(a::T1, b::T2) where {T1 <: Pair2D, T2 <: Real} = T1(a.x/b,a.y/b)



"""
    norm(v::Point2D)

l-2 norm of vector

#Arguments
- `v::Point2D`: Point2D vector
#Returns
- `l::Float64`: length of vector
"""
function norm(v::Point2D{T})::T where T<:AbstractFloat
    l = LinearAlgebra.norm([v.x, v.y], 2.0)
    # l = sqrt(v.x^2 + v.y^2)
    return l
end



"""
    normit(v::Point2D)::Point2D

Return normalized vector coordinates as Point2D
"""
function normit(v::Point2D)::Point2D
    vn = v / norm(v)
    return vn
end



"""
    get_com(p1::Point2D{T}, p2::Point2D{T})::Point2D{T} where T<:AbstractFloat

#Arguments
- `p1::Point2D{T}`: first point of line
- `p2::Point2D{T}`: second point of line

#Returns
- `c::Point2D{T}`:  The center of line ((center of mass)) defined by two points p1 and p2 
"""
function get_com(p1::Vector2D, p2::Vector2D)::Point2D
    c = p1 + ((p2 - p1) / 2)
    return c
end



function get_norm_vec(p1::Vector2D, p2::Vector2D, dir::Symbol)::Point2D
    # calculate normal vector
    n = Point2D((-1) * (p2.y - p1.y), (p2.x - p1.x))
    (dir == :neg) && (n = n * (-1.0))
    n_norm = normit(n)
    return n_norm
end



@inline function get_length(p1::Vector2D, p2::Vector2D)
    return norm(p1-p2)
end



function get_angle(v1::Vector2D,v2::Vector2D)
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