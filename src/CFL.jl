function CFL_cond(u, equation::Equation)
    maximum(abs.(equation.funcs.Dflux(u)))
end

function CFL_cond!(::Scalar, integrator::Integrator)
    @unpack equation, cache = integrator
    @unpack Dfcont = cache
    # integrator.cfl = maximum(abs.(Dflux(equation.funcs, uprev)))
    integrator.cfl = 0.0
    for k in eachindex(Dfcont)
        integrator.cfl = max(integrator.cfl, abs(Dfcont[k]))
    end
end

function CFL_local!(::Scalar, integrator::Integrator)
    @unpack cache = integrator
    @unpack stencil, Dfcont = cache
    # cache.cfl_loc = max(abs(Dfcont[stencil[1],j]), abs(Dfcont[stencil[2],j]))
    cache.cfl_loc = 0.0
    for k in eachindex(stencil)
        cache.cfl_loc = max(cache.cfl_loc, abs(Dfcont[stencil[k]]))
    end
end

CFL_cond!(::System, integrator::Integrator) = CFL_cond!(integrator.equation.funcs, integrator)
CFL_local!(::System, integrator::Integrator) = CFL_local!(integrator.equation.funcs, integrator)



# function CFL_local!(integrator::Integrator, u)
#     @unpack equation, cache = integrator
#     cache.cfl_loc = maximum(abs.(equation.Dflux(u)))
# end


function dt_CFL!(integrator::Integrator)
    @unpack params, t, = integrator
    @unpack mesh, tf, CFL_factor = params
    CFL_cond!(integrator.equation.eqtype, integrator)
    integrator.dt = min(CFL_factor * mesh.dx / integrator.cfl, tf - t)
end
