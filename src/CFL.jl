mutable struct CFLCacheScalar <: CFLCacheType
    cfl::Float64
    absDfcont::Vector{Float64}
    function CFLCacheScalar(equation, uinit)
        absDfcont= zero(uinit)
        absDfcont .= abs.(Dflux(equation.funcs, uinit))
        new(zero(Float64), absDfcont)
    end
end

function CFL_cond(u, equation::Equation)
    maximum(abs.(equation.funcs.Dflux(u)))
end

function CFL_cond!(::Scalar, integrator::Integrator)
    @unpack equation, cache = integrator
    @unpack cfl_cache = cache
    @unpack absDfcont = cfl_cache
    # integrator.cfl = maximum(abs.(Dflux(equation.funcs, uprev)))
    cfl_cache.cfl = 0.0
    for k in eachindex(absDfcont)
        # cfl_cache.cfl = max(cfl_cache.cfl, abs(Dfcont[k]))
        cfl_cache.cfl = max(cfl_cache.cfl, absDfcont[k])
    end
end

function CFL_local!(::OneD, ::Scalar, integrator::Integrator, j::Int)
    @unpack cache, space_cache = integrator
    @unpack cfl_cache = cache
    @unpack absDfcont = cfl_cache
    @unpack Nx = integrator.params.mesh

    space_cache.cfl_loc = absDfcont[j]
    space_cache.cfl_loc = max(space_cache.cfl_loc, absDfcont[mod1(j+1,Nx)])

    # D1, D2 = absDfcont[j], absDfcont[mod1(j+1,Nx)]
    # space_cache.cfl_loc = max(D1,D2)

    # cache.cfl_loc = max(abs(Dfcont[stencil[1],j]), abs(Dfcont[stencil[2],j]))
    # D1, D2 = Dfcont[j], Dfcont[mod1(j+1,Nx)]
    # space_cache.cfl_loc = max(abs(D1), abs(D2))
    # space_cache.cfl_loc = max(abs(Dfcont[j]), abs(Dfcont[mod1(j+1,Nx)]))
    # space_cache.cfl_loc = 0.0
    # for k in eachindex(stencil)
    #     space_cache.cfl_loc = max(space_cache.cfl_loc, abs(Dfcont[stencil[k]]))
    # end
end

CFL_cond!(::System, integrator::Integrator) = CFL_cond!(integrator.equation.funcs, integrator)
CFL_local!(dim::EquationDim, ::System, integrator::Integrator, args...) = CFL_local!(dim, integrator.equation.funcs, integrator, args...)



# function CFL_local!(integrator::Integrator, u)
#     @unpack equation, cache = integrator
#     cache.cfl_loc = maximum(abs.(equation.Dflux(u)))
# end


function dt_CFL!(::OneD, integrator::Integrator)
    @unpack params, t, cache = integrator
    @unpack mesh, tf, CFL_factor = params
    @unpack cfl_cache = cache
    CFL_cond!(integrator.equation.eqtype, integrator)
    integrator.dt = min(CFL_factor * mesh.dx / cfl_cache.cfl, tf - t)
end

