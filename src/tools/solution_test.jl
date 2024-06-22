include("domain.jl")

struct EntropyOptSolution
    domain::Domain
end


abstract type A end

abstract type B <: A end

struct C{T,T2}
    x::T
    y::T2
    C(z) = new{typeof(z)}(z * 2)
end

Base.eltype(::C{T}) where {T} = T

f(x) = x
f(::A) = print("zouzou A")

f(::B) = print("zouzou B")

f(c::C) = c.x^2


abstract type Method end

struct Method1 end
struct Method2 end

flux(::Method1, args; kwarsg...)
flux(::Method2, args; kwarsg...)