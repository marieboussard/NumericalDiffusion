using FiniteVolumes
using BenchmarkTools
using UnPack
using Plots
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")

# Domain definition
Nx = 100
xmin, xmax = -4, 4
t0, tf = 0.0, 0.2
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)


# A sequence of function which converges to the constant function equal to 1, as n tends to infinity
u0n(n::Int, x::Real) =  exp(-x^2)/n

function quantify_consistency(nvec::Vector{Int}, params::Parameters, eqfun::AbstractEquationFun, alpha::Real=0)
    N = length(nvec)
    gapvec = zeros(N)
    for i in 1:N
        n = nvec[i]
        u0(x::Real) = u0n(n, x)
        # function u0(x::Real)
        #     if x <= -1
        #         return zero(x)
        #     elseif x >=1
        #         return zero(x)
        #     else
        #         return one(x)
        #     end
        # end
        u0(x::AbstractVector) = u0.(x)
        # equation = Equation(OneD(), 1, Scalar(), Burgers(), u0)
        # aa = 3
        # equation = Equation(OneD(), 1, Scalar(), Advection(aa), u0)
        equation = Equation(OneD(), 1, Scalar(), eqfun, u0)

        sol = solve(equation, params, Euler(), Rusanov(); maxiter=1, log_config=LogConfig(true, false, true, false, false));

        @show sol.niter
        plot(mesh.x, sol.uinit, label="uinit")
        display(plot!(mesh.x, sol.u, label="u"))

        # estimate_post = quantify_diffusion(sol, Posteriori(AsymmetricMD()));
        

        # Multidimensional bounds for ΔG
        estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
        @unpack uinit, u, l, L, etacont, etacont_init = estimate

        # Entropy flux from theory
        Gexact = G_from_theory(estimate)

        # DG from theory
        DGexact = zeros(Nx)
        for i in eachindex(DGexact)
            DGexact[i] = (Gexact[i] - Gexact[mod1(i-1,Nx)])*estimate.dt/mesh.dx
        end

        # # display(plot!(mesh.x, uinit, label="estimate"))

        # Defining optimization components
        Gc = zeros(eltype(u), Nx)
        A = zeros(eltype(u), 2*Nx, Nx)
        b = zeros(eltype(u), 2*Nx)
        W = zeros(eltype(u), Nx, Nx)
        Gflux!(CenteredG(), Gc, estimate)
        fill_A!(A, estimate)
        fill_b!(b, estimate)
        fill_W!(W, estimate, alpha)

        # Uzawa solving
        optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=100000, eps=1e-12);
        gapvec[i] = norm(optsol.Gcgap)

        # DG and diffusion for optimised flux
        @unpack gamma_opt = optsol
        DG = zeros(Nx)
        for i in eachindex(DG)
            DG[i] = (gamma_opt[i] - gamma_opt[mod1(i-1,Nx)])*estimate.dt/mesh.dx
        end
        D = zero(L)
        @unpack etacont_init, etacont = estimate
        diffusion!(Posteriori(), gamma_opt, etacont_init, etacont, estimate.dt, mesh, D)

        # Diffusion for consistent flux
        Dc = zero(L)
        diffusion!(Posteriori(), Gc, etacont_init, etacont, estimate.dt, mesh, Dc)

        # Diffusion for exact flux
        Dexact = zero(L)
        diffusion!(Posteriori(), Gexact, etacont_init, etacont, estimate.dt, mesh, Dexact)

        # Visualisation
        pltA = []

        plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
            legendfontsize=15,
            titlefontsize=21,
            guidefontsize=21,
            tickfontsize=18)

        plot!(mesh.x, Gc, label="consistent", lw=2)
        plot!(mesh.x, Gexact, label="exact", lw=2)
        plot!(mesh.x, optsol.gamma_opt, label="uzawa", lw=2)

        xlabel!("x")
        ylabel!("Numerical Entropy Flux")
        title!("Nx="*string(Nx)*", alpha="*string(alpha))

        push!(pltA, plt1)


        plt2 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomleft,
        legendfontsize=15,
        titlefontsize=21,
        guidefontsize=21,
        tickfontsize=18)

        # scatter!(mesh.x, estimate.D, label="priori multidim")
        plot!(mesh.x, D, label="uzawa", lw=2)
        plot!(mesh.x, Dexact, label="from theory", lw=2)
        plot!(mesh.x, Dc, label="Gc", lw=2)
        xlabel!("x")
        ylabel!("Numerical Diffusion")
        title!("Max Diff:"*string(maximum(D)))
        push!(pltA, plt2)


        plt3 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottom,
        legendfontsize=15,
        titlefontsize=21,
        guidefontsize=21,
        tickfontsize=18)

        plot!(mesh.x, estimate.l, label="l", lw=2)
        plot!(mesh.x, estimate.L, label="L", lw=2)
        plot!(mesh.x, DGexact, label="exact", lw=2)
        plot!(mesh.x, DG, label="DG", lw=2)
        xlabel!("x")
        ylabel!("Bounds for numerical Entropy Flux")

        push!(pltA, plt3)


        plt4 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topright,
        legendfontsize=15,
        titlefontsize=21,
        guidefontsize=21,
        tickfontsize=18)

        plot!(mesh.x, optsol.popt[1:Nx], label="popt for l", lw=2)
        plot!(mesh.x, optsol.popt[Nx+1:end], label="popt for L", lw=2)
        xlabel!("x")
        ylabel!("Lagrange multiplier")

        push!(pltA, plt4)

        display(plot(pltA..., layout=(2,2), size=(1600, 1200)))

    end
    gapvec
end

nvec = [1, 10, 100, 1000]
# nvec = [1]
# gapvec = quantify_consistency(nvec, params, Advection(3.0))
gapvec = quantify_consistency(nvec, params, Burgers())
scatter(log10.(nvec), log10.(gapvec))