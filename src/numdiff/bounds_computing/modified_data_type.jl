abstract type ModifiedDataType end
abstract type SymmetricMD <: ModifiedDataType end
struct LinearCombMD <: SymmetricMD
    weights
end
struct MeanMD <: SymmetricMD end
struct MaxMD <: SymmetricMD end
struct MinMD <: SymmetricMD end
struct AsymmetricMD <: ModifiedDataType end


get_name(::SymmetricMD) = "symmetric"
get_name(::AsymmetricMD) = "asymmetric"

computeK(::MeanMD, u::AbstractVector) = mean(u)
computeK!(::MeanMD, u::AbstractVector, res::AbstractVector) = res .= mean(u)
# computeK!(::MeanMD, u::AbstractMatrix, res::AbstractVector) = res .= mean(u, dims=1)'
# computeK!(::MeanMD, u::AbstractMatrix, res::AbstractVector) = res .= vec(mean(u, dims=1))
function computeK!(::MeanMD, u::AbstractMatrix, res::AbstractVector)
    fill!(res, zero(eltype(res)))
    for I in CartesianIndices(u)
        res[I[2]] += u[I]
    end
    res ./= size(u, 1)
end


computeK(::MaxMD, u::AbstractVector) = maximum(u)
computeK!(::MaxMD, u::AbstractVector, res::AbstractVector) = res .= maximum(u)
computeK!(::MaxMD, u::AbstractMatrix, res::AbstractVector) = @. res = maximum(u, dims=1)

# CACHE FOR MODIFIED DATA

abstract type ModifiedDataCache <: Cache end