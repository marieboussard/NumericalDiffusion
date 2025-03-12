abstract type ModifiedDataType end
abstract type SymmetricMD <: ModifiedDataType end
struct LinearCombMD <: SymmetricMD
    weights
end
struct MaxMD <: SymmetricMD end
struct MinMD <: SymmetricMD end
struct AsymmetricMD <: ModifiedDataType end

computeK(::MaxMD, u::AbstractArray) = maximum(u, dims=1)
computeK!(::MaxMD, u::AbstractArray, res::AbstractArray) = @. res = maximum(u, dims=1)