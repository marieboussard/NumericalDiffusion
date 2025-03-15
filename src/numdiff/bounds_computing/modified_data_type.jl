abstract type ModifiedDataType end
abstract type SymmetricMD <: ModifiedDataType end
struct LinearCombMD <: SymmetricMD
    weights
end
struct MeanMD <: SymmetricMD end
struct MaxMD <: SymmetricMD end
struct MinMD <: SymmetricMD end
struct AsymmetricMD <: ModifiedDataType end

computeK(::MeanMD, u::AbstractVector) = mean(u)
computeK!(::MeanMD, u::AbstractVector, res::AbstractVector) = res .= mean(u)
computeK!(::MeanMD, u::AbstractMatrix, res::AbstractVector) = res .= mean(u, dims=1)'

computeK(::MaxMD, u::AbstractVector) = maximum(u)
computeK!(::MaxMD, u::AbstractVector, res::AbstractVector) = res .= maximum(u)
computeK!(::MaxMD, u::AbstractMatrix, res::AbstractVector) = @. res = maximum(u, dims=1)