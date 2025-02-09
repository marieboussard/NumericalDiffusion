
# abstract type ScalarEquation <: AbstractEquation end
# abstract type SystemEquation <: AbstractEquation end

abstract type EquationType end
struct Scalar <: EquationType end
struct System <: EquationType end

abstract type AbstractEquation{eqType<:EquationType} end

# struct Equation{eqType<:EquationType} <: AbstractEquation{eqType}
#     p::Int
#     function Equation(p::Int)
#         if p==1
#             new{Scalar}(p)
#         else
#             new{System}(p)
#         end
#     end
# end

struct Equation{T<:EquationType, F1} <: AbstractEquation{T}
    p::Int
    eqType::T
    fun::F1
end

struct Equation2{T<:EquationType} <: AbstractEquation{T}
    p::Int
    eqType::T
    fun
end