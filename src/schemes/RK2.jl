struct RK2 <: TimeScheme end

# CACHE

mutable struct RK2Cache{utype<:AbstractArray, tctype<:TimeCache} <: TimeCache 

    dx::Float64
    dt::Float64

    fnum_euler::utype
    u_bar::utype
    fnum_bar::utype

    time_cache_euler::tctype
    time_cache_bar::tctype


    function RK2Cache(space_scheme::SpaceScheme, u::AbstractArray, dim::EquationDim, params::Parameters, dt::Float64)

        @unpack Nx = params.mesh
        dx = params.mesh.dx

        # If we are working on the whole domain
        if length(u) == Nx

            fnum_euler = zero(u)
            u_bar = zero(u)
            fnum_bar = zero(u)
            
            time_cache_euler = EulerCache(init_cache(space_scheme, u, dim))
            time_cache_bar = EulerCache(init_cache(space_scheme, u, dim))

        else # If we are working on a subdomain (we cannot use periodic BCs)

            sL_space, sR_space = get_sL(space_scheme), get_sR(space_scheme)
            fnum_euler = zeros(eltype(u), length(u) - sL_space - sR_space + 1)
            u_bar = zeros(eltype(u), length(fnum_euler) - 1)
            fnum_bar = zeros(eltype(u), length(u_bar) - sL_space - sR_space + 1)

            time_cache_euler = EulerCache(init_cache(space_scheme, u, dim))
            time_cache_bar = EulerCache(init_cache(space_scheme, u_bar, dim))

        end
        new{typeof(u_bar), typeof(time_cache_euler)}(dx, dt, fnum_euler, u_bar, fnum_bar, time_cache_euler, time_cache_bar)
    end
end

# AUXILIARY FUNCTIONS

get_sL(::RK2, scheme::SpaceScheme) = 2*get_sL(scheme)
get_sR(::RK2, scheme::SpaceScheme) = 2*get_sR(scheme)

# NUMERICAL FLUX COMPUTATION

function global_numflux!(::RK2, space_scheme::SpaceScheme, time_cache::TimeCache, equation::Equation, u::AbstractArray, fnum::AbstractArray, jstart::Int=1, jend::Int=length(u), Nx=length(u), shift::Int=0)

    @unpack dx, dt, fnum_euler, u_bar, fnum_bar, time_cache_euler, time_cache_bar = time_cache

    if length(fnum_bar) == length(fnum_euler)

        
    
        global_numflux!(Euler(), space_scheme, time_cache_euler, equation, u, fnum_euler)
        for i in 1:Nx
            for j in 1:equation.p
                u_bar[i,j] = u[i,j] - dt / dx * (fnum_euler[i,j] - fnum_euler[mod1(i-1, Nx),j])
            end
        end

        # Then, apply the flux to u_bar
        global_numflux!(Euler(), space_scheme, time_cache_bar, equation, u_bar, fnum_bar)

        # The final flux is the mean of the two temporary fluxes
        @. fnum = 0.5*(fnum_euler + fnum_bar)

    else

        sL_space, sR_space = get_sL(space_scheme), get_sR(space_scheme)

        # First, compute u_bar by applying an iteration of the space scheme
        global_numflux!(Euler(), space_scheme, time_cache_euler, equation, u, fnum_euler, sL_space, length(u)-sR_space, Nx, -sL_space+1)
        for i in 1:length(u_bar)
            for r in 1:equation.p
                u_bar[i,r] = u[i+sL_space, r] - dt / dx * (fnum_euler[i+1, r] - fnum_euler[i,r])
            end
        end

        # Then, apply the flux to u_bar
        global_numflux!(Euler(), space_scheme, time_cache_bar, equation, u_bar, fnum_bar, sL_space, length(u_bar)-sR_space, Nx, -sL_space+1)

        # The final flux is the mean of the two temporary fluxes
        short_fnum_euler = view(fnum_euler, 1+sL_space:length(fnum_euler)-sR_space)
        #short_fnum_euler = view(fnum_euler, sL_space:length(fnum_euler)-sR_space-1)
        @. fnum = 0.5*(short_fnum_euler + fnum_bar)

    end

end
    

function global_numflux!(ts::RK2, integrator::Integrator)
    integrator.time_cache.dt = integrator.dt
    global_numflux!(ts, integrator.space_scheme, integrator.time_cache, integrator.equation, integrator.uprev, integrator.fnum)
end