include("../src/include_file.jl")

xmin, xmax, Nx, t0, Tf = -2, 2, 100, 0, 0.2
CFL_factor = 0.5
domain = createInterval(Nx, xmin, xmax, t0, Tf)
testcase = SimpleShock()
#testcase = ArticleTestcase()
u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_fun(testcase, domain.x[i])] end; res)

method = Roe(CFL_factor)
#method = Rusanov(CFL_factor)

solBurgers = fv_solve(domain, u0, burgers(), method)
plot_fv_sol(solBurgers, nb_plots=6)
#plot_fv_sol(solBurgers, uexact_burgers_article)
title!(get_name(method))