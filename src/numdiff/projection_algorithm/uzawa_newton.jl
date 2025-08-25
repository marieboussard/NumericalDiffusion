"
output of function `compute_entropic_G`
See also [`compute_entropic_G`](@ref).
"
struct UzawaNewtonSolution{soltype<:Solution, wtype<:AbstractMatrix, atype<:AbstractMatrix, gtype<:AbstractVector, Mtype<:AbstractMatrix, ptype<:AbstractVector, bound_mode_type<:BoundMode, weights_type<:AbstractNormWeights}

    # Data info
    sol::soltype

    # Primal formulation
    #> QP components
    W::wtype # Weights matrix
    A::atype # Constraints matrix
    b::ptype # Constraint vector 
    Gc::gtype # Consistent flux 
    Dc::gtype # Diffusion associated to consistent flux 

    #> Primal variables
    G_uz::gtype
    D_uz::gtype
    G_newt::gtype
    D_newt::gtype

    # Dual formulation
    #> LCP components
    M::Mtype 
    q::ptype
    #> Dual variables 
    p_uz::ptype
    w_uz::ptype
    p_newt::ptype
    w_newt::ptype

    # Settings
    bound_mode::bound_mode_type
    weights::weights_type

    # Convergence info
    niter_uz::Int
    niter_newt::Int

end

"""
    compute_entropic_G(params::Parameters, equation::Equation; <keywords arguments>)

Compute an entropic numerical entropy flux for the finite volume resolution of `equation` with parameters `params`.

The method used is Uzawa algorithm for the projection, and a Newton step can be added for speed-up.

By default, the scheme is Euler + Rusanov, but it can also be specified.
"""
function compute_entropic_G(params::Parameters, equation::Equation; bound_mode::BoundMode=SingleBound(), weights::AbstractNormWeights=AbsWeights(), maxiter_uzawa::Int=1000, maxiter_newton::Int=1000, time_scheme::TimeScheme=Euler(), space_scheme::SpaceScheme=Rusanov(), use_newton::Bool=true, ent_numflux::AbstractEntNumFlux=CenteredG(), kwargs...)

    # Finite volumes resolution
    sol = solve(equation, params, time_scheme, space_scheme; log_config=LogConfig(true, false, true, false, false), kwargs...);

    # Multidimensional bounds for ΔG
    estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
    @unpack uinit, u, l, L = estimate

    Gc, A, b, W = init_optim_components(bound_mode, estimate, weights; ent_numflux=ent_numflux)

    optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=maxiter_uzawa, eps=1e-12);

    @show optsol.constraint_residual

    # Associated flux and diffusion
    Gopt = optsol.gamma_opt
    Dopt = zero(L)
    @unpack etacont_init, etacont = estimate
    diffusion!(Posteriori(), Gopt, etacont_init, etacont, estimate.dt, params.mesh, Dopt)

    # Diffusion associated with consistent flux 
    Dc = zero(L)
    diffusion!(Posteriori(), Gc, etacont_init, etacont, estimate.dt, params.mesh, Dc)

    if use_newton

        # Now, solve the dual LCP
        M = A*inv(W)*A'
        q = b - A*Gc
        w0 = b - A*optsol.gamma_opt
        pend, wend, niter = newton_lcp(M, -q, optsol.popt, w0; maxiter=maxiter_newton)

        # We get the flux back 
        Gopt_newt = Gc - inv(W)*A'*pend
        newton_residual = norm(max.(0.0, A*Gopt_newt .- b))
        @show newton_residual
        Dopt_newt = zero(L)
        diffusion!(Posteriori(), Gopt_newt, etacont_init, etacont, estimate.dt, params.mesh, Dopt_newt)

        return UzawaNewtonSolution(sol, W, A, b, Gc, Dc, Gopt, Dopt, Gopt_newt, Dopt_newt, M, q, optsol.popt, w0, pend, wend, bound_mode, weights, optsol.niter, niter)

    end
    
    UzawaNewtonSolution(sol, W, A, b, Gc, Dc, Gopt, Dopt, Gopt, Dopt, zero(A'), zero(optsol.popt), optsol.popt, zero(optsol.popt), optsol.popt, zero(optsol.popt), bound_mode, weights, optsol.niter, 0)

end