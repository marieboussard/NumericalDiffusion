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

computeK(::MaxMD, u::AbstractVector) = maximum(u)
computeK(::MaxMD, u::AbstractArray) = maximum(u, dims=1)
computeK!(::MaxMD, u::AbstractArray, res::AbstractArray) = @. res = maximum(u, dims=1)