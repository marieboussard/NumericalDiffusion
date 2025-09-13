using NumericalDiffusion
using UnPack
using LinearAlgebra
using GLM, DataFrames
using CairoMakie


function compute_beta_slope(Nxvec::Vector{Int}, equation::Equation; xmin::Real=-2, xmax::Real=2, t0::Real=0.0, tf::Real=0.4, CFL_factor::Real=0.5, weights=AbsWeights(), time_scheme::TimeScheme=Euler(), space_scheme::SpaceScheme=Rusanov(), ent_numflux::AbstractEntNumFlux=CenteredG(), maxiter_uzawa::Int=1000, maxiter_newton::Int=1000, bound_mode=SingleBound())

    N = length(Nxvec)
    dxvec = zeros(N)
    Hvec = zeros(N)
    nitervec = zeros(Int, N)
    maxweightvec = zeros(N)

    for k in 1:N

        Nx = Nxvec[k]
        @show Nx
        mesh = OneDMesh(Nx, xmin, xmax)
        params = Parameters(mesh, t0, tf, CFL_factor)

        dxvec[k] = mesh.dx

        sol = NumericalDiffusion.solve(equation, params, time_scheme, space_scheme, log_config=LogConfig(true, true, true, false, false))
        @unpack niter = sol

        nitervec[k] = niter

        H = zero(eltype(sol.u))

        for n in 1:niter

            # Multidimensional bounds for ΔG
            estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()), n)
            @unpack uinit, u, l, L = estimate
            Gc, A, b, W = init_optim_components(bound_mode, estimate, weights; ent_numflux=ent_numflux)

            # UZAWA
            optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=maxiter_uzawa, eps=1e-12)
            @show optsol.constraint_residual

            # NEWTON
            # Now, solve the dual LCP
            M = A * inv(W) * A'
            q = b - A * Gc
            w0 = b - A * optsol.gamma_opt
            pend, wend, niter = newton_lcp(M, -q, optsol.popt, w0; maxiter=maxiter_newton)
            # We get the flux back 
            G_newt = Gc - inv(W) * A' * pend
            newton_residual = norm(max.(0.0, A * G_newt .- b))
            @show newton_residual

            for i in 1:Nx
                maxweightvec[k] = max(maxweightvec[k], 1.0/abs(W[i,i]))
                H += W[i, i] * (G_newt[i] - Gc[i])^2
            end
        end

        Hvec[k] = H

    end

    dxvec, Hvec, nitervec, maxweightvec

end


function add_slope!(axes::AbstractVector, time_scheme_vec::AbstractVector, space_scheme_vec::AbstractVector, Nxvec::Vector{Int}, equation::Equation; kwargs...)

    ax, ax2, ax3, ax4 = axes

    N = length(time_scheme_vec)
    nitermat = zeros(N, length(Nxvec))

    for k in 1:N

        time_scheme = time_scheme_vec[k]
        space_scheme = space_scheme_vec[k]

        dxvec, Hvec, nitervec, maxweightvec = compute_beta_slope(Nxvec, equation; time_scheme=time_scheme, space_scheme=space_scheme, kwargs...)

        view(nitermat, k, :) .= nitervec

        df = DataFrame(x=log10.(dxvec[end-3:end]), y=log10.(abs.(Hvec[end-3:end])))  # Création d'un DataFrame
        # Ajustement du modèle de régression linéaire
        model = lm(@formula(y ~ x), df)
        # Affichage des résultats
        println(model)
        coeffs = round.(coef(model); digits=2)
        dx_be = [log10(dxvec[1]), log10(dxvec[end])]

        label = get_name(time_scheme) * " + " * get_name(space_scheme)*", "*get_name(bound_mode)*", alpha="*string(weights.alpha)

        scatter!(ax, log10.(dxvec), log10.(abs.(Hvec)), label=label * ", beta=" * string(-round(coeffs[2], digits=2)))
        lines!(ax, dx_be, dx_be * coeffs[2] .+ coeffs[1])

        scatter!(ax2, Nxvec, abs.(Hvec), label=label * ", beta=" * string(round(coeffs[2], digits=2)))

        barplot!(ax3, log10.(dxvec), nitervec, label=label, width=0.1)

        scatter!(ax4, log10.(dxvec), maxweightvec, label=label)
        lines!(ax4, log10.(dxvec), maxweightvec)
    end

    nitermat

end

function plot_beta_slopes(Nxvec::AbstractVector, equation, time_scheme_vec::AbstractVector, space_scheme_vec::AbstractVector; xmin::Real=-2, xmax::Real=2, t0::Real=0.0, tf::Real=0.03, weights::AbstractNormWeights=AbsWeights(), bound_mode=SingleBound(), kwargs...)

    fig = Figure(size = (800, 1500))
    ax = Axis(fig[1, 1], title="Beta slopes", xlabel="log(dx)", ylabel="log|H|")
    ax3 = Axis(fig[2,1], title="Number of time iterations", xlabel="log(dx)", ylabel="Niter")
    ax4 = Axis(fig[3,1], title="Maximal inverse weights", xlabel="log(dx)", ylabel="max(1/w)")

    fig2 = Figure()
    ax2 = Axis(fig2[1, 1], title="H", xlabel="Nx", ylabel="H")

    axes = [ax, ax2, ax3, ax4]

    nitermat = add_slope!(axes, time_scheme_vec, space_scheme_vec, Nxvec, equation; xmin=xmin, xmax=xmax, t0=t0, tf=tf, weights=weights, bound_mode=bound_mode, kwargs...)
    axislegend(ax, position=:lb)
    axislegend(ax2, position=:lt)
    axislegend(ax3, position=:rt)
    axislegend(ax4, position=:rt)
    fig

end

# Nxlog = LinRange(1, 2, 6)
# #Nxlog = LinRange(1,1.5,4)
# Nxvec = Int.(floor.(10 .^ (Nxlog)))
# xmin, xmax = -2, 2
# t0, tf = 0.0, 0.03
# equation = BurgersArticle
# alpha = 1
# weights = AbsWeights(alpha)
# bound_mode = SingleBound()

# # time_scheme_vec = [Euler()]
# # space_scheme_vec = [Rusanov()]

# # time_scheme_vec = [Euler(), Euler()]
# # space_scheme_vec = [Roe(), Rusanov()]

# time_scheme_vec = [RK2(), RK2()]
# space_scheme_vec = [MUSCL(Rusanov(), Superbee()), MUSCL(Rusanov(), Minmod())]

# fig = plot_beta_slopes(Nxvec, equation, time_scheme_vec, space_scheme_vec; xmin=xmin, xmax=xmax, t0=t0, tf=tf, weights=weights, bound_mode=bound_mode)