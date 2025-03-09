struct HR <: SpaceScheme
    subscheme::SpaceScheme
end

mutable struct HRCache{dtype<: AbstractArray, utype<:AbstractArray, subcache_type<:scacheType} <: scacheType
    hplus::dtype    # hi+1/2,+
    hminus::dtype   # hi+1/2,-
    zinter::dtype   # zi+1/2
    uplus::utype    # Uj+1/2,+
    uminus::utype   # Uj+1/2,-
    fplus::utype    # f(Uj+1/2,+)
    fminus::utype   # f(Uj+1/2,-)
    subcache::subcache_type
    function HRCache(scheme::SpaceScheme, uinit, eqdim::EquationDim, args...)
        # hplus = zero(selectdim(uinit, ndims(uinit), 1))
        hplus = init_hplus(eqdim, uinit)
        hminus = copy(hplus)
        zinter = copy(hplus)
        # uplus = zeros(eltype(uinit), size(uinit)[end])
        uplus = init_uplus(eqdim, uinit)
        uminus = copy(uplus)
        fplus = copy(uplus)
        fminus = copy(uplus)
        subcache = init_cache(scheme.subscheme, args...)
        new{typeof(hplus), typeof(uplus), typeof(subcache)}(hplus, hminus, zinter, uplus, uminus, fplus, fminus, subcache)
    end
end

init_hplus(::OneD, uinit::AbstractArray) = zero(selectdim(uinit, ndims(uinit), 1))
init_hplus(::TwoD, uinit::AbstractArray) = zeros(eltype(uinit), size(uinit)[1:2]..., 2)

init_uplus(::OneD, uinit::AbstractArray) = zeros(eltype(uinit), size(uinit)[end])
init_uplus(::TwoD, uinit::AbstractArray) = zeros(eltype(uinit), size(uinit)[end], 2)

function fillcache!(::HR, integrator::Integrator)
    @unpack space_cache, source_cache, params = integrator
    @unpack znum = source_cache
    @unpack Nx = params.mesh
    for j in eachindex(znum)
        space_cache.zinter[j] = max(znum[j], znum[mod1(j+1, Nx)])
    end
end

get_sL(scheme::HR) = get_sL(scheme.subscheme)
get_sR(scheme::HR) = get_sR(scheme.subscheme)

function hplus!(integrator::Integrator, j::Int)
    @unpack params, space_cache, source_cache, uprev = integrator
    @unpack Nx = params.mesh
    @unpack hplus, zinter = space_cache
    @unpack znum = source_cache
    hplus[j] = max(zero(eltype(znum)),  uprev[mod1(j+1, Nx), 1] + znum[mod1(j+1, Nx)] - zinter[j])
end

function hminus!(integrator::Integrator, j::Int)
    @unpack space_cache, source_cache, uprev = integrator
    @unpack hminus, zinter = space_cache
    @unpack znum = source_cache
    hminus[j] = max(zero(eltype(znum)),  uprev[j, 1] + znum[j] - zinter[j])
end

function numflux!(scheme::HR, integrator::Integrator, j::Int, args...)
    @unpack equation, params, space_cache, uprev = integrator
    @unpack Nx = params.mesh
    hplus!(integrator, j)
    hminus!(integrator, j)
    @unpack hplus, hminus, uplus, uminus, fplus, fminus = space_cache
    uplus[1], uplus[2] = hplus[j], hplus[j]*uprev[mod1(j+1,Nx),2]
    uminus[1], uminus[2] = hminus[j], hminus[j]*uprev[j,2]
    flux!(equation.funcs, uplus, fplus)
    flux!(equation.funcs, uminus, fminus)
    numflux!(scheme.subscheme, uminus, uplus, fplus, fminus, space_cache.subcache, integrator, j)
end

# SOURCE DISCRETIZATION

function discretize_sourceterm!(::OneD, ::HRDisc, integrator::Integrator)
    @unpack params, cache, source_cache = integrator
    @unpack hplus, hminus = integrator.space_cache
    @unpack Nx, dx = params.mesh
    for j in eachindex(hminus)
        cache.sourceterm[j,1] = zero(eltype(cache.sourceterm))
        cache.sourceterm[j,2] = 0.5*g*(hminus[j]^2-hplus[mod1(j-1, Nx)]^2)/dx
    end
end