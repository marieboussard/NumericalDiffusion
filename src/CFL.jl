mutable struct CFLCacheScalar <: CFLCacheType
    cfl::Float64
    Dfcont::Vector{Float64}
    function CFLCacheScalar(equation, uinit)
        Dfcont= zero(uinit)
        Dfcont .= Dflux(equation.funcs, uinit)
        new(zero(Float64), Dfcont)
    end
end

function CFL_cond(u, equation::Equation)
    maximum(abs.(equation.funcs.Dflux(u)))
end

function CFL_cond!(::Scalar, integrator::Integrator)
    @unpack equation, cache = integrator
    @unpack cfl_cache = cache
    @unpack Dfcont = cfl_cache
    # integrator.cfl = maximum(abs.(Dflux(equation.funcs, uprev)))
    cfl_cache.cfl = 0.0
    for k in eachindex(Dfcont)
        cfl_cache.cfl = max(cfl_cache.cfl, abs(Dfcont[k]))
    end
end

function CFL_local!(::Scalar, integrator::Integrator)
    @unpack cache, space_cache = integrator
    @unpack stencil, cfl_cache = cache
    @unpack Dfcont = cfl_cache
    # cache.cfl_loc = max(abs(Dfcont[stencil[1],j]), abs(Dfcont[stencil[2],j]))
    space_cache.cfl_loc = 0.0
    for k in eachindex(stencil)
        space_cache.cfl_loc = max(space_cache.cfl_loc, abs(Dfcont[stencil[k]]))
    end
end

CFL_cond!(::System, integrator::Integrator) = CFL_cond!(integrator.equation.funcs, integrator)
CFL_local!(::System, integrator::Integrator) = CFL_local!(integrator.equation.funcs, integrator)



# function CFL_local!(integrator::Integrator, u)
#     @unpack equation, cache = integrator
#     cache.cfl_loc = maximum(abs.(equation.Dflux(u)))
# end


function dt_CFL!(integrator::Integrator)
    @unpack params, t, cache = integrator
    @unpack mesh, tf, CFL_factor = params
    @unpack cfl_cache = cache
    CFL_cond!(integrator.equation.eqtype, integrator)
    integrator.dt = min(CFL_factor * mesh.dx / cfl_cache.cfl, tf - t)
end

# CFL CONDITION FOR TWO DIMENSIONNAL EQUATIONS 

function CFL_cond2D(u, equation)

end

function dt_CFL2D!(integrator::Integrator)
    @unpack params, t, cache = integrator
    @unpack mesh, tf, CFL_factor = params
    @unpack cfl_cache = cache
    CFL_cond!(integrator.equation.eqtype, integrator)
    integrator.dt = min(tf - t, CFL_factor/(cfl_cache.cfl))
    
    
    min(CFL_factor * mesh.dx / cfl_cache.cfl, tf - t)
end