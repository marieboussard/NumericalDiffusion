abstract type Domain{T<:Real} end

mutable struct Interval{T<:Real} <: Domain{T}
    Nx::Int
    xmin::T
    xmax::T
    t0::T
    Tf::T
    dx::T
    x::Vector{T}
    interfaces::Vector{T}
    sourceVec::Union{Nothing, Matrix{T}, Vector{T}}
    DSourceVec::Union{Nothing, Array{T}}
end

#Interval(Nx::Int, xmin::Real, xmax::Real, t0::Real, Tf::Real, dx::Real, x, interfaces, sourceVec, DSourceVec) = Interval(Nx, promote(xmin, xmax, t0, Tf, dx)..., x, interfaces, sourceVec, DSourceVec)

# # A general constructor for Interval handling different argument types
# function Interval(Nx::Int, xmin::Real, xmax::Real, t0::Real, Tf::Real, dx::Real, x::Vector{Real}, interfaces::Vector{Real}, sourceVec, DSourceVec)

# end

createInterval(Nx::Int, xmin::Real, xmax::Real, t0::Real, Tf::Real) = createInterval(Nx, promote(xmin, xmax, t0, Tf)...)

function createInterval(Nx::Int, xmin::T, xmax::T, t0::T, Tf::T) where {T<:Real}
    dx = (xmax - xmin) / Nx
    x = collect(LinRange(xmin+dx/2, xmax-dx/2, Nx))
    #interfaces = LinRange(xmin - dx / 2, xmax + dx / 2, Nx + 1)
    #interfaces = collect(LinRange(xmin - dx, xmax, Nx + 1))
    interfaces = collect(LinRange(xmin, xmax, Nx + 1))
    Interval(Nx, xmin, xmax, t0, Tf, dx, x, interfaces, nothing, nothing)
end

createUnitInterval(Nx::Int, t0::Real, Tf::Real) = createInterval(Nx, 0.0, 1.0, t0, Tf)
createInterval(interval::Interval) = Interval(interval.Nx, interval.xmin, interval.xmax, interval.t0, interval.Tf, interval.dx, interval.x, interval.interfaces, interval.sourceVec, interval.DSourceVec)

mutable struct CartesianMesh{T<:Real} <: Domain{T}
    Nx::Int
    Ny::Int
    xmin::T
    xmax::T
    ymin::T
    ymax::T
    t0::T
    Tf::T
    dx::T
    dy::T
    x::Vector{T}
    y::Vector{T}
    # interfaces::Vector{T}
    sourceVec::Union{Nothing, Matrix{T}, Vector{T}}
    DSourceVec::Union{Nothing, Array{T}}
end

function CartesianMesh(Nx::Int, Ny::Int, xmin::T, xmax::T, ymin::T, ymax::T, t0::T, Tf::T) where {T<:Real}
    dx = (xmax - xmin) / Nx
    dy = (ymax - ymin) / Ny 
    x = collect(LinRange(xmin+dx/2, xmax-dx/2, Nx))
    y = collect(LinRange(ymin+dy/2, ymax-dy/2, Ny))
    CartesianMesh(Nx, Ny, xmin, xmax, ymin, ymax, t0, Tf, dx, dy, x, y, nothing, nothing)
end

CartesianMesh(Nx::Int, Ny::Int, xmin, xmax, ymin, ymax, t0, Tf) = CartesianMesh(Nx::Int, Ny::Int, promote(xmin, xmax, ymin, ymax, t0, Tf)...)