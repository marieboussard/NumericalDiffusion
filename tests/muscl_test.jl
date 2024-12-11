include("../src/include_file.jl")

# Solving Burgers equation with MUSCL (constructed with Rusanov)

xmin, xmax, Nx, t0, Tf = -2, 2, 50, 0, 0.4
CFL_factor = 0.5
domain = createInterval(Nx, xmin, xmax, t0, Tf)
testcase = ArticleTestcase()
u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_fun(testcase, domain.x[i])] end; res)
method = MUSCL(0.5, Roe(0.5), Minmod(), domain)

solBurgers_muscl = fv_solve(domain, u0, burgers(), method)
solBurgers_Rusanov = fv_solve(domain, u0, burgers(), Rusanov(0.5))
#plot_fv_sol(solBurgers, nb_plots=6)
plot_fv_sol(solBurgers, uexact_burgers_article)
plot!(domain.x, solBurgers_Rusanov.u_approx[end], label="Rusanov")
display(title!("MUSCL"))