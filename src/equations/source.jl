# abstract type AbstractSource{eqFunType<:AbstractEquationFun} end
# struct NoSource <: AbstractSource{EquationFun} end
abstract type AbstractSource end
struct NoSource <: AbstractSource end

# has_source(::AbstractSource{T}) where T = true
has_source(::AbstractSource) = true
has_source(::NoSource) = false

# s(::NoSource) = nothing
# Ds(::NoSource) = nothing

# struct Source{F1<:Base.Callable, F2<:Base.Callable} <: AbstractSource 
#     s::F1
#     Ds::F2
# end

# s(source::Source, x) = source.s(x)
# Ds(source::Source, x) = source.Ds(x)