function CFL_cond(u, equation::Equation)
    maximum(abs.(equation.funcs.Dflux(u)))
end

function CFL_cond!(::Scalar, integrator::Integrator)
    @unpack equation, uprev = integrator
    integrator.cfl = maximum(abs.(equation.funcs.Dflux(uprev)))
end

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
