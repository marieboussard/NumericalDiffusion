include("../../src/include_file.jl")

# Parameters
xmin, xmax, Nx, t0, Tf = -2, 2, 20, 0, 0.4
CFL_factor = 0.5
equation = burgers()
scheme = FVScheme(Euler(), Rusanov(CFL_factor))
# scheme = FVScheme(Euler(), Roe(CFL_factor))
testcase = ArticleTestcase()
#testcase = SimpleShock()
domain, u0 = createOneTimestepInterval(Nx, t0, xmin, xmax, equation, testcase, CFL_factor)
#u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_fun(testcase, domain.x[i])] end; res)

# Finite volume resolution
fv_sol = fv_solve(domain, u0, equation, scheme)
display(plot_fv_sol(fv_sol, uexact_burgers_article))

# Consistent numerical entropy Flux
Gc = CenteredG()

sol = lsEntropicConsistentG(domain, equation, scheme, u0, Gc);
D_ls = diffusion(u0, fv_sol.u_approx[end], sol.Gopt, domain.dx, domain.Tf, equation, domain)

solEnt = optimize_for_entropy(u0, domain, equation, scheme)

plot(domain.interfaces, sol.Gopt, label="Least Squares")
plot!(domain.interfaces, solEnt.Gopt, label="A posteriori")
display(title!("Numerical Entropy Flux"))

plot(domain.x, D_ls, label="Least Squares")
plot!(domain.x, solEnt.Dopt, label="A posteriori")
display(title!("Numerical Diffusion"))