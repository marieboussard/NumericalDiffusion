abstract type Domain end

mutable struct Interval <: Domain

    xmin::Real
    xmax::Real
    Nx::Int
    t0::Real
    Tf::Real
    dx::Real
    x::Vector{Real}
    interfaces::Vector{Real}
    sourceVec
    DSourceVec
end

function createInterval(xmin::Real, xmax::Real, Nx::Int, t0::Real, Tf::Real)
    dx = (xmax - xmin) / Nx
    x = LinRange(xmin, xmax, Nx)
    interfaces = LinRange(xmin - dx / 2, xmax + dx / 2, Nx + 1)
    Interval(xmin, xmax, Nx, t0, Tf, dx, x, interfaces, nothing, nothing)
end

createUnitInterval(Nx::Int, t0::Real, Tf::Real) = createInterval(0.0, 1.0, Nx, t0, Tf)
createInterval(interval::Interval) = Interval(interval.xmin, interval.xmax, interval.Nx, interval.t0, interval.Tf, interval.dx, interval.x, interval.interfaces, interval.sourceVec, interval.DSourceVec)


