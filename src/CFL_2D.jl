# CFL CONDITION FOR TWO DIMENSIONNAL EQUATIONS 

mutable struct CFLCacheScalar2D{dataType<:AbstractArray} <: CFLCacheType
    cflx::Float64
    cfly::Float64
    Dfcont::dataType
    # Dhcont::dataType
    function CFLCacheScalar2D(equation::Equation, uinit)
        Dfcont = zeros(eltype(uinit), size(uinit)...,2)
        # Dhcont = zero(uinit)
        # Dfcont .= Dflux_f(equation.funcs, uinit)
        # Dhcont .= Dflux_h(equation.funcs, uinit)
        Dflux!(equation.funcs, uinit, selectdim(Dfcont, ndims(Dfcont),1), selectdim(Dfcont, ndims(Dfcont),2))
        new{typeof(Dfcont)}(zero(Float64), zero(Float64), Dfcont)
    end
end

function CFL_cond2D!(::Scalar, integrator::Integrator)
    @unpack equation, cache = integrator
    @unpack cfl_cache = cache
    # @unpack Dfcont, Dhcont = cfl_cache
    @unpack Dfcont = cfl_cache
    cfl_cache.cflx = 0.0
    cfl_cache.cfly = 0.0
    # for k in eachindex(Dfcont)
    #     cfl_cache.cflx = max(cfl_cache.cflx, abs(Dfcont[k]))
    #     cfl_cache.cfly = max(cfl_cache.cfly, abs(Dhcont[k]))
    # end
    for I in CartesianIndices(view(Dfcont,:,:,1))
        cfl_cache.cflx = max(cfl_cache.cflx, abs(Dfcont[I,1]))
        cfl_cache.cfly = max(cfl_cache.cfly, abs(Dfcont[I,2]))
    end
end

function CFL_local!(::TwoD, ::Scalar, integrator::Integrator, j::Int, k::Int)
    @unpack cache, space_cache = integrator
    @unpack cfl_cache = cache
    @unpack Nx, Ny = integrator.params.mesh
    @unpack Dfcont = cfl_cache
    # @unpack Dfcont, Dhcont = cfl_cache
    # space_cache.cflx_loc = max(abs(Dfcont[j,k]), abs(Dfcont[mod1(j+1,Nx),k]))
    # space_cache.cfly_loc = max(abs(Dhcont[j,k]), abs(Dhcont[j, mod1(k+1,Ny)]))
    space_cache.cflx_loc = max(abs(Dfcont[j,k,1]), abs(Dfcont[mod1(j+1,Nx),k,1]))
    space_cache.cfly_loc = max(abs(Dfcont[j,k,2]), abs(Dfcont[j, mod1(k+1,Ny),2]))
end

CFL_cond2D!(::System, integrator::Integrator, args...) = CFL_cond2D!(integrator.equation.funcs, integrator, args...)
# CFL_local2D!(::System, integrator::Integrator, args...) = CFL_local2D!(integrator.equation.funcs, integrator, args...)

function dt_CFL!(::TwoD, integrator::Integrator)
    @unpack params, t, cache = integrator
    @unpack mesh, tf, CFL_factor = params
    @unpack cfl_cache = cache
    CFL_cond2D!(integrator.equation.eqtype, integrator)
    integrator.dt = min(tf - t, CFL_factor/(cfl_cache.cflx/mesh.dx + cfl_cache.cfly/mesh.dy))
    nothing
end