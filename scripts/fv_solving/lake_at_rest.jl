using FinitVolumes

function init_lake_at_rest(x::T, znum::T; c=one(eltype(x))) where T<:AbstractVector
    v = zeros(eltype(x), (length(x), 2))
    for i in eachindex(x)
        v[i,1] = max(zero(eltype(x)), c - znum[i])
        v[i,2] = zero(eltype(x))
    end
    return v
end