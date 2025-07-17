struct RK2 <: TimeScheme end

struct RK2Cache <: TimeCache 
    u_bar
    fcont_bar
    Dfcont_bar
    fnum_bar
end

get_sL(::RK2, scheme::SpaceScheme) = 2*get_sL(scheme)
get_sR(::RK2, scheme::SpaceScheme) = 2*get_sR(scheme)

function numflux!(::RK2, integrator::Integrator, args...)
    @unpack time_cache, uprev = integrator
    @unpack u_bar, fcont_bar, fnum_bar = time_cache

    # First, compute u_bar by applying an iteration of the space scheme
    numflux!(Euler(), integrator, args...)
    @unpack fnum = integrator
    for i in 1:Nx
        for j in 1:equation.p
            u_bar[i,j] = uprev[i,j] - dt / dx * (fnum[i,j] - fnum[mod1(i-1, Nx),j])
        end
    end

    # Then, compute auxiliary values like f(u_bar) and others depending on the space scheme
    flux!(equation.funcs, u_bar, fcont_bar)

    # Apply the flux to u_bar
    numflux!(Euler(), integrator.space_scheme, integrator.params, integrator.equation, integrator.cache, integrator.space_cache, fnum_bar, fcont_bar, u_bar)

    # Final form of the flux
    @. fnum = 0.5*(fnum + fnum_bar)

end