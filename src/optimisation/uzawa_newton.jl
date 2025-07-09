include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")
include("../../src/optimisation/newton_for_LCP.jl")

struct UzawaNewtonSolution{wtype<:AbstractMatrix, atype<:AbstractMatrix, gtype<:AbstractVector, Mtype<:AbstractMatrix, ptype<:AbstractVector, bound_mode_type<:BoundMode, weights_type<:AbstractNormWeights}

    # Primal formulation
    #> QP components
    W::wtype # Weights matrix
    A::atype # Constraints matrix
    b::ptype # Constraint vector 
    Gc::gtype
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

function compute_entropic_G(params::Parameters, equation::Equation; bound_mode=SingleBound(), weights=AbsWeights())

    # Finite volumes resolution
    sol = solve(equation, params, Euler(), Rusanov(); log_config=LogConfig(true, false, true, false, false));

    # Multidimensional bounds for ΔG
    estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
    @unpack uinit, u, l, L = estimate

    Gc, A, b, W = init_optim_components(bound_mode, estimate, weights)

    optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=1000, eps=1e-12);

    @show optsol.constraint_residual

    # Associated flux and diffusion
    Gopt = optsol.gamma_opt
    Dopt = zero(L)
    @unpack etacont_init, etacont = estimate
    diffusion!(Posteriori(), Gopt, etacont_init, etacont, estimate.dt, mesh, Dopt)

    # Now, initialize the corresponding LCP
    M = A*inv(W)*A'
    q = b - A*Gc
    w0 = b - A*optsol.gamma_opt
    pend, wend, niter = newton_lcp(M, -q, optsol.popt, w0; maxiter=10000)

    # We get the flux back 
    Gopt_newt = Gc - inv(W)*A'*pend
    Dopt_newt = zero(L)
    diffusion!(Posteriori(), Gopt_newt, etacont_init, etacont, estimate.dt, mesh, Dopt_newt)
    
    UzawaNewtonSolution(W, A, b, Gc, Gopt, Dopt, Gopt_newt, Dopt_newt, M, q, optsol.popt, w0, pend, wend, bound_mode, weights, optsol.niter, niter)

end