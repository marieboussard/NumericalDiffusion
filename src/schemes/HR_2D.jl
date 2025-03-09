struct HR2D <: SpaceScheme
    subscheme::SpaceScheme
end

function fillcache!(::HR2D, integrator::Integrator)
    @unpack space_cache, source_cache, params = integrator
    @unpack znum = source_cache
    @unpack Nx, Ny = params.mesh
    for j in 1:Nx
        for k in 1:Ny
            space_cache.zinter[j,k,1] = max(znum[j,k], znum[mod1(j+1, Nx), k])
            space_cache.zinter[j,k,2] = max(znum[j,k], znum[j, mod1(k+1, Ny)])
        end
    end
end

get_sL(scheme::HR2D) = get_sL(scheme.subscheme)
get_sR(scheme::HR2D) = get_sR(scheme.subscheme)

function hplus2!(integrator::Integrator, j::Int, k::Int)
    @unpack params, space_cache, source_cache, uprev = integrator
    @unpack Nx, Ny = params.mesh
    @unpack hplus, zinter = space_cache
    @unpack znum = source_cache
    hplus[j,k,1] = max(zero(eltype(znum)),  uprev[mod1(j+1, Nx), k, 1] + znum[mod1(j+1, Nx), k] - zinter[j, k, 1])
    hplus[j,k,2] = max(zero(eltype(znum)),  uprev[j, mod1(k+1, Ny), 1] + znum[j, mod1(k+1, Ny)] - zinter[j, k, 2])
end

function hminus2!(integrator::Integrator, j::Int, k::Int)
    @unpack space_cache, source_cache, uprev = integrator
    @unpack hminus, zinter = space_cache
    @unpack znum = source_cache
    hminus[j,k,1] = max(zero(eltype(znum)),  uprev[j, k, 1] + znum[j,k] - zinter[j,k,1])
    hminus[j,k,2] = max(zero(eltype(znum)),  uprev[j, k, 1] + znum[j,k] - zinter[j,k,2])
end

function numflux!(scheme::HR2D, integrator::Integrator, j::Int, k::Int, args...)
    @unpack equation, params, space_cache, uprev = integrator
    @unpack Nx, Ny = params.mesh
    hplus2!(integrator, j, k)
    hminus2!(integrator, j, k)
    @unpack hplus, hminus, uplus, uminus, fplus, fminus = space_cache

    uplus[1,1], uplus[2,1], uplus[3,1] = hplus[j,k,1] .* (one(eltype(hplus)), uprev[mod1(j+1,Nx),k,2], uprev[mod1(j+1,Nx),k,3])
    uplus[1,2], uplus[2,2], uplus[3,2] = hplus[j,k,2] .*(one(eltype(hplus)), uprev[j,mod1(k+1,Ny),2], uprev[j,mod1(k+1,Ny),3])

    uminus[1,1], uminus[2,1], uminus[3,1] = hminus[j,k,1] .* (one(eltype(hminus)), uprev[j,k,2], uprev[j,k,3])
    uminus[1,2], uminus[2,2], uminus[3,2] = hminus[j,k,2] .*(one(eltype(hminus)), uprev[j,k,2], uprev[j,k,3])

    flux_f!(equation.funcs, view(uplus,:,1), view(fplus,:,1))
    flux_h!(equation.funcs, view(uplus,:,2), view(fplus,:,2))
    flux_f!(equation.funcs, view(uminus,:,1), view(fminus,:,1))
    flux_h!(equation.funcs, view(uminus,:,2), view(fminus,:,2))

    numflux!(scheme.subscheme, view(uminus,:,1), view(uplus,:,1), view(fplus,:,1), view(fminus,:,1), view(integrator.fnum.fnum,j,k,:), space_cache.subcache, equation)
    numflux!(scheme.subscheme, view(uminus,:,2), view(uplus,:,2), view(fplus,:,2), view(fminus,:,2), view(integrator.fnum.hnum,j,k,:), space_cache.subcache, equation)
end

# SOURCE DISCRETIZATION

function discretize_sourceterm!(::TwoD, ::HRDisc2, integrator::Integrator)
    @unpack params, cache, source_cache = integrator
    @unpack hplus, hminus = integrator.space_cache
    @unpack Nx, Ny, dx, dy = params.mesh
    for j in 1:Nx
        for k in 1:Ny
            cache.sourceterm[j,k,1] = zero(eltype(cache.sourceterm))
            cache.sourceterm[j,k,2] = 0.5*g*(hminus[j,k,1]^2-hplus[mod1(j-1, Nx),k,1]^2)/dx
            cache.sourceterm[j,k,3] = 0.5*g*(hminus[j,k,2]^2-hplus[j,mod1(k-1, Ny),2]^2)/dy
        end
    end
end