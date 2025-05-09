struct RK2 <: TimeScheme end

struct RK2Cache{utype<:AbstractArray, ctype<:Cache} <: TimeCache 
    fnum_euler::utype
    u_bar::utype
    fcont_bar::utype
    fnum_bar::utype
    subcache::ctype
    function RK2Cache(space_scheme::SpaceScheme, equation::Equation, uinit::AbstractVector)
        fnum_euler = zero(uinit)
        u_bar = zero(uinit)
        fcont_bar = zero(uinit)
        fnum_bar = zero(uinit)
        subcache = init_subcache(space_scheme, equation.dim, equation.eqtype, equation.funcs, uinit)
        new{typeof(u_bar), typeof(subcache)}(fnum_euler, u_bar, fcont_bar, fnum_bar, subcache)
    end
end


get_sL(::RK2, scheme::SpaceScheme) = 2*get_sL(scheme)
get_sR(::RK2, scheme::SpaceScheme) = 2*get_sR(scheme)

### GLOBAL VERSION 

function numflux!(::RK2, integrator::Integrator,)
    @unpack space_scheme, time_cache, uprev, equation, dt = integrator
    @unpack Nx, dx = integrator.params.mesh
    @unpack fnum_euler, u_bar, fcont_bar, fnum_bar, subcache = time_cache

    # First, compute u_bar by applying an iteration of the space scheme
    # numflux!(Euler(), integrator)
    numflux!(Euler(), space_scheme, integrator.params, fnum_euler, integrator.space_cache, integrator.cache.cfl_cache, integrator.equation, integrator.uprev, integrator.fcont)
    # @unpack fnum = integrator
    for i in 1:Nx
        for j in 1:equation.p
            u_bar[i,j] = uprev[i,j] - dt / dx * (fnum_euler[i,j] - fnum_euler[mod1(i-1, Nx),j])
        end
    end

    # Then, compute auxiliary values like f(u_bar) and others depending on the space scheme
    flux!(equation.funcs, u_bar, fcont_bar)
    update_subcache!(space_scheme, equation.dim, equation.eqtype, subcache, equation, u_bar)

    # Apply the flux to u_bar
    numflux!(Euler(), integrator.space_scheme, integrator.params, fnum_bar, integrator.space_cache, subcache, integrator.equation, u_bar, fcont_bar)

    # Final form of the flux
    @. integrator.fnum = 0.5*(fnum_euler + fnum_bar)

end

### LOCAL VERSION 

function numflux!(::RK2, space_scheme::SpaceScheme, j::Int, i::Int, params::Parameters, equation::Equation, cache::Cache, space_cache::SpaceCache, fnum::AbstractArray, fcont::AbstractArray, u::AbstractArray)

    # In this version, it is assumed that update_timecache! has been called before, for the computations that need to be done in a global way
    for r in 1:equation.p
        fnum[i,r] = 0.5*(fnum[j,r] + fnum_bar[j,r])
    end

end