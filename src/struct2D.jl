################################################################
# Code part of FARADS2D                                        #
# by dominik.bueschgens                                        #
# March 2021                                                   #
################################################################


abstract type Pair2D <: AbstractVector{Real} end

struct Point2D{T<:AbstractFloat} <: Pair2D
    x::T
    y::T
end

Point2D(x::Real, y::Real) = Point2D(promote(x,y)...)
Point2D(x::Integer, y::Integer) = Point2D(convert(Float64,x), convert(Float64,y))

struct Index2D{T<:Integer} <: Pair2D
    x::T
    y::T
end

Pair2D(x::T,y::T) where T <: AbstractFloat = Point2D(x,y)
Pair2D(x::Integer,y::T) where T <: AbstractFloat = Point2D(x,y)
Pair2D(x::T,y::Integer) where T <: AbstractFloat = Point2D(x,y)
Pair2D(x::T,y::T) where T <: Integer = Index2D(x,y)

import Base: size, getindex, setindex!, +, -, *, /

size(v::Pair2D) = (2, )
getindex(v::Pair2D, i::Int) = i < 3 ? (i == 1 ? v.x : v.y) : error("max index = 2")

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

+(a::T, b::T) where T <: Pair2D = T(a.x+b.x,a.y+b.y)
-(a::T, b::T) where T <: Pair2D = T(a.x-b.x,a.y-b.y)
*(a::T1, b::T2) where {T1 <: Pair2D, T2 <: Real} = T1(a.x*b,a.y*b)
/(a::T1, b::T2) where {T1 <: Pair2D, T2 <: Real} = T1(a.x/b,a.y/b)


