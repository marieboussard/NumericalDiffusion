abstract type ModifiedDataType end
abstract type SymmetricModifiedData <: ModifiedDataType end
struct CLModifiedData <: SymmetricModifiedData
    weights
end
struct MaxModifiedData <: SymmetricModifiedData end
struct MinModifiedData <: SymmetricModifiedData end
struct AsymmetricModifiedData <: ModifiedDataType end

function extractLocalData(u::AbstractArray{T}, j, sL, sR) where T

    # To keep only cells that K takes as arguments
    #Nx = length(u)
    Nx, p = size(u)
    u_short = zeros(T, sL + sR, p)
    i = 1

    for k in j-sL+1:j+sR
        u_short[i,:] = u[mod1(k, Nx),:]
        i += 1
    end

    return u_short
end

function extractExtendedLocalData(u, j, sL, sR)

    # To keep only cells for which the modified flux is not G(K)
    Nx, p = size(u)

    u_short = zeros(2*(sR+sL),p)
    i=1
    for k in j-sL-sR+1:j+sR+sL
        u_short[i,:] = u[mod1(k, Nx), :]
        i+=1
    end

    return u_short
end

extractExtendedLocalData(u::Nothing, j, sL, sR) = nothing

computeK(clModifiedData::CLModifiedData, u) = sum(clModifiedData.weights .* u) / length(u)
function computeK(clModifiedData::CLModifiedData, u::AbstractArray{T}) where T
    Nx, p = size(u)
    K = zeros(T, p)
    for k in 1:p
        K[k] = sum(clModifiedData.weights .* u[:,k]) / sum(clModifiedData.weights)
    end
    K
end

#computeK(::MaxModifiedData, u) = max(u...)
computeK(::MaxModifiedData, u) = transpose(maximum(u, dims=1))
computeK(::MinModifiedData, u) = transpose(minimum(u, dims=1))

meanK(sL, sR) = CLModifiedData(ones(sL + sR))
maxK() = MaxModifiedData()
minK() = MinModifiedData()
midLeftK(sL, sR) = CLModifiedData(vcat(zeros(sL-1), 1, zeros(sR)))
midRightK(sL, sR) = CLModifiedData(vcat(zeros(sL), 1, zeros(sR-1)))

computeZ(KFun::ModifiedDataType, z, j, sL, sR) = computeK(KFun, extractLocalData(z, j, sL, sR))
computeZ(::ModifiedDataType, ::Nothing, j, sL, sR) = nothing