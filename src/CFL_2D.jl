# CFL CONDITION FOR TWO DIMENSIONNAL EQUATIONS 

mutable struct CFLCacheScalar2D <: CFLCacheType
    clfx::Float64
    cfly::Float64
    Dfcont::Matrix{Float64}
    Dhcont::Matrix{Float64}
    function CFLCacheScalar2D(equation::Equation, uinit)
        Dfcont = zero(uinit)
        Dhcont = zero(uinit)
        Dfcont .= Dflux_f(equation.funcs, uinit)
        Dhcont .= Dflux_h(equation.funcs, uinit)
        new(zero(Float64), zero(Float64), Dfcont, Dhcont)
    end
end

function CFL_cond2D!(::Scalar, integrator::Integrator)
    @unpack equation, cache = integrator
    @unpack cfl_cache = cache
    @unpack Dfcont, Dhcont = cfl_cache
    cfl_cache.cflx = 0.0
    cfl_cache.cfly = 0.0
    for k in eachindex(Dfcont)
        cfl_cache.cflx = max(cfl_cache.cflx, abs(Dfcont[k]))
        cfl_cache.cfly = max(cfl_cache.cfly, abs(Dhcont[k]))
    end
end

# function CFL_local2D!(::Scalar, integrator::Integrator)
#     @unpack cache, space_cache = integrator
#     @unpack stencil, cfl_cache = cache
#     @unpack Dfcont, Dhcont = cfl_cache
#     space_cache.cfl_loc = 0.0
#     for k in eachindex(stencil)
#         space_cache.cfl_loc = max(space_cache.cfl_loc, abs(Dfcont[stencil[k]]))
#     end
# end


function dt_CFL2D!(integrator::Integrator)
    @unpack params, t, cache = integrator
    @unpack mesh, tf, CFL_factor = params
    @unpack cfl_cache = cache
    CFL_cond2D!(integrator.equation.eqtype, integrator)
    integrator.dt = min(tf - t, CFL_factor/(cfl_cache.cflx/mesh.dx + cfl_cache.cfly/mesh.dy))
end