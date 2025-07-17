# abstract type AbstractEquation end
# abstract type ScalarEquation <: AbstractEquation end

abstract type EquationType end
struct Scalar <: EquationType end
struct System <: EquationType end

abstract type AbstractEquation{eqType<:EquationType} end