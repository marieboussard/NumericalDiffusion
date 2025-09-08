using NumericalDiffusion
using UnPack
using LinearAlgebra
using GLM, DataFrames


function compute_beta_slope(Nxvec::Vector{Int}, equation::Equation; xmin::Real=-2, xmax::Real=2, t0::Real=0.0, tf::Real=0.4, CFL_factor::Real=0.5, weights=AbsWeights(), time_scheme::TimeScheme=Euler(), space_scheme::SpaceScheme=Rusanov(), ent_numflux::AbstractEntNumFlux=CenteredG(), maxiter_uzawa::Int=1000, maxiter_newton::Int=1000, kwargs...)

    N = length(Nxvec)
    dxvec = zeros(N)
    Hvec = zeros(N)
    nitervec = zeros(Int, N)

    for k in 1:N

        Nx = Nxvec[k]
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
                H += W[i, i] * (G_newt[i] - Gc[i])^2
            end
        end

        Hvec[k] = H

    end

    dxvec, Hvec, nitervec

end


function add_slope!(ax::Axis, ax2::Axis, time_scheme_vec::AbstractVector, space_scheme_vec::AbstractVector, Nxvec::Vector{Int}, equation::Equation; kwargs...)
    N = length(time_scheme_vec)

    nitermat = zeros(N, length(Nxvec))

    for k in 1:N

        time_scheme = time_scheme_vec[k]
        space_scheme = space_scheme_vec[k]

        dxvec, Hvec, nitervec = compute_beta_slope(Nxvec, equation; time_scheme=time_scheme, space_scheme=space_scheme, kwargs...)

        view(nitermat, k, :) .= nitervec

        df = DataFrame(x=log10.(dxvec), y=log10.(abs.(Hvec)))  # Création d'un DataFrame
        # Ajustement du modèle de régression linéaire
        model = lm(@formula(y ~ x), df)
        # Affichage des résultats
        println(model)
        coeffs = round.(coef(model); digits=2)
        dx_be = [log10(dxvec[1]), log10(dxvec[end])]

        scatter!(ax, log10.(dxvec), log10.(abs.(Hvec)), label=get_name(time_scheme) * " + " * get_name(space_scheme) * ": " * string(-round(coeffs[2], digits=2)))
        lines!(ax, dx_be, dx_be * coeffs[2] .+ coeffs[1])

        scatter!(ax2, Nxvec, abs.(Hvec), label=get_name(time_scheme) * " + " * get_name(space_scheme) * ": " * string(round(coeffs[2], digits=2)))
    end

    nitermat

end

Nxlog = LinRange(1, 3, 5)
Nxvec = Int.(floor.(10 .^ (Nxlog)))
xmin, xmax = -2, 2
t0, tf = 0.0, 0.05
equation = BurgersArticle
weights = AbsWeights(1.5)

time_scheme_vec = [Euler(), RK2()]
space_scheme_vec = [Rusanov(), Rusanov()]

using CairoMakie

fig = Figure()
ax = Axis(fig[1, 1], title="Beta slopes", xlabel="log(dx)", ylabel="log|H|")

fig2 = Figure()
ax2 = Axis(fig2[1, 1], title="H", xlabel="Nx", ylabel="H")

nitermat = add_slope!(ax, ax2, time_scheme_vec, space_scheme_vec, Nxvec, equation; xmin=xmin, xmax=xmax, t0=t0, tf=tf, weights=weights)
axislegend(ax, position=:lb)
axislegend(ax2, position=:lt)
fig2
