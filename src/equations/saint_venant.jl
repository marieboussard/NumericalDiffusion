const g = 9.8
const treshold = 1e-6

struct SaintVenant <: AbstractEquationFun end

function flux(::SaintVenant, v)
    res = similar(v)
    g_half = g * 0.5 
    for i in eachindex(view(v, :, 1))
        h = v[i, 1]
        if h > treshold
            hu = v[i, 2]
            res[i, 1] = hu
            res[i, 2] = hu^2 / h + g_half * h^2
        else
            res[i, 1] = 0
            res[i, 2] = 0
        end
    end
    res
end

function Dflux(::SaintVenant, u)
    [0]
end


function CFL_cond(::SaintVenant, v)
    max = 0.0
    lamb = 0.0
    for i in eachindex(view(v, :, 1))
        h = v[i,1]
        hu = v[i,2]
        if h > treshold
            lamb = abs(hu / h) + sqrt(g * abs(h))
            if lamb > max
                max = lamb
            end
        end
    end
    return max
end

function CFL_cond!(::SaintVenant, integrator::Integrator)

    @unpack uprev = integrator
    integrator.cfl = 0.0
    lamb = 0.0

    for i in eachindex(view(uprev, :, 1))
        h = uprev[i,1]
        hu = uprev[i,2]
        if h > treshold
            lamb = abs(hu / h) + sqrt(g * abs(h))
            if lamb > integrator.cfl
                integrator.cfl = lamb
            end
        end
    end
end

function CFL_local!(::SaintVenant, integrator::Integrator)

    @unpack uprev, cache = integrator
    @unpack stencil = cache
    cache.cfl_loc = 0.0
    lamb = 0.0

    for i in eachindex(stencil)
        h = uprev[i,1]
        hu = uprev[i,2]
        if h > treshold
            lamb = abs(hu / h) + sqrt(g * abs(h))
            if lamb > cache.cfl_loc
                cache.cfl_loc = lamb
            end
        end
    end
end

function init_lake_at_rest(x::AbstractVector)
    zeros(eltype(x), (length(x), 2))
end

SaintVenantAtRest = Equation(2, System(), SaintVenant(), init_lake_at_rest)