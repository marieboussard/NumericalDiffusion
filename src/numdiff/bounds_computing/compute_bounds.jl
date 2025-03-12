function utilde!(::SymmetricMD, estimator::Estimator, j::Int)
    # return ̃uᵢʲ⁺¹/²
    @unpack u = indices
    @unpack Nx = estimator.params.mesh
    @unpack sL, sR, indices, utilde = estimator.cache

    indices = mod1.(j-sL-sR+1 : j+sL+sR, Nx)
    @view ushort = u[indices]

    computeK!(estimator.md, ushort, estimator.cache.K)

    for i in 1:2*(sL+sR)
        utilde[i] = 0
    end

end

function uhat!(estimator::Estimator, j::Int)
    @unpack utilde = estimator.cache
    @unpack uhat = estimator.cache
end

function compute_G_bounds!(estimator::Estimator)
    @unpack Nx = estimator.params.mesh
    for j in 1:Nx
        utilde!(estimator.mdtype, estimator, j)
        uhat!(estimator, j)
    end
end