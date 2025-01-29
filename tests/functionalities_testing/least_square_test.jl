include("../../src/include_file.jl")

# Parameters
xmin, xmax, Nx, t0, Tf = -2, 2, 20, 0, 0.4
CFL_factor = 0.5
equation = burgers()
scheme = FVScheme(Euler(), Rusanov(CFL_factor))
# scheme = FVScheme(Euler(), Roe(CFL_factor))
testcase = ArticleTestcase()
#testcase = SimpleShock()
#domain, u0 = createOneTimestepInterval(Nx, t0, xmin, xmax, equation, testcase, CFL_factor)
domain = createInterval(Nx, xmin, xmax, t0, Tf)
u0 = initialData(domain, testcase)

# Finite volume resolution
fv_sol = fv_solve(domain, u0, equation, scheme)
display(plot_fv_sol(fv_sol, uexact_burgers_article))

# Consistent numerical entropy Flux
Gc = CenteredG()
G_centered = vecNumFlux(equation.source, Gc, equation, fv_sol.u_approx[end-1])
D_centered = diffusion(fv_sol.u_approx[end-1], fv_sol.u_approx[end], G_centered, domain.dx, fv_sol.dt_vec[end], equation, domain)

sol = lsEntropicConsistentG(domain, equation, scheme, u0, Gc);
D_ls = diffusion(u0, fv_sol.u_approx[end], sol.Gopt, domain.dx, domain.Tf, equation, domain)

solEnt = optimize_for_entropy(u0, domain, equation, scheme)

plot(domain.interfaces, sol.Gopt, label="Least Squares")
#scatter!(domain.interfaces, solEnt.Gopt, label="A posteriori")
#plot!(domain.interfaces, G_centered, label="Centered")
display(title!("Numerical Entropy Flux"))

plot(domain.x, D_ls, label="Least Squares")
#scatter!(domain.x, solEnt.Dopt, label="A posteriori")
#plot!(domain.x, D_centered, label="Centered")
display(title!("Numerical Diffusion"))