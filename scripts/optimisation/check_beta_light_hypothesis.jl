using NumericalDiffusion
using UnPack

function check_beta_hypothesis(Nxvec::AbstractVector, beta::Real, equation::Equation; xmin::Real=-4, xmax::Real=4, t0::Real=0.0, tf::Real=0.2, CFL_factor::Real=0.5,bound_mode=SingleBound(), weights::AbstractNormWeights=AbsWeights(), maxiter_uzawa::Int=1000, time_scheme::TimeScheme=Euler(), space_scheme::SpaceScheme=Rusanov(), ent_numflux=CenteredG())
    N = length(Nxvec)
    hvec = zeros(N)
    for k in 1:N
        Nx = Nxvec[k]
        mesh = OneDMesh(Nx, xmin, xmax)
        params = Parameters(mesh, t0, tf, CFL_factor)
        sol = solve(equation, params, time_scheme, space_scheme; log_config=LogConfig(true, true, true, false, false));
        for n in 1:sol.niter

            # Multidimensional bounds for ΔG
            estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()), n);
            @unpack uinit, u, l, L = estimate

            Gc, A, b, W = init_optim_components(bound_mode, estimate, weights; ent_numflux=ent_numflux)

            optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=maxiter_uzawa, eps=1e-12);

            @show optsol.constraint_residual

            # Associated flux and diffusion
            Gopt = optsol.gamma_opt

            for j in 1:Nx
                hvec[k] += W[j,j]*(Gopt[j]-Gc[j])^2
            end
        end
        hvec[k] *= (mesh.dx)^beta
    end
    hvec
end

#Nxvec = [10, 20, 50, 100, 200]
Nxvec = [2^i for i in 3:10]

beta = 1
equation = BurgersArticle

hvec = check_beta_hypothesis(Nxvec, beta, equation; maxiter_uzawa=10000)
@show Nxvec


using Plots
scatter(log.(Nxvec), log.(hvec))