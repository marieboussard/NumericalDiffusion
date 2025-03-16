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
# computeK!(::MeanMD, u::AbstractMatrix, res::AbstractVector) = res .= mean(u, dims=1)'
# computeK!(::MeanMD, u::AbstractMatrix, res::AbstractVector) = res .= vec(mean(u, dims=1))
function computeK!(::MeanMD, u::AbstractMatrix, res::AbstractVector)
    fill!(res, 0)  # Met res à zéro sans allocation
    for i in axes(u, 1)  # Parcours les lignes
        res .+= u[i, :]  # Ajoute chaque ligne sans allocation
    end
    res ./= size(u, 1)  # Normalise par le nombre de lignes
end


computeK(::MaxMD, u::AbstractVector) = maximum(u)
computeK!(::MaxMD, u::AbstractVector, res::AbstractVector) = res .= maximum(u)
computeK!(::MaxMD, u::AbstractMatrix, res::AbstractVector) = @. res = maximum(u, dims=1)

# CACHE FOR MODIFIED DATA

abstract type ModifiedDataCache <: Cache end