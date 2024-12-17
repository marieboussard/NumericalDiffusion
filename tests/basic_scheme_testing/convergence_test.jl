include("../src/include_file.jl")

xmin, xmax = -2, 2
Nx_vec = [2^i for i in 3:15]

equation = burgers()
#testcase = ArticleTestcase()
testcase = SmoothTestcase()
CFL_factor = 0.5
#scheme = FVScheme(Euler(), Rusanov(CFL_factor))
#scheme = FVScheme(Euler(), MUSCL(CFL_factor, Rusanov(CFL_factor), Minmod()))
scheme = FVScheme(RK2(), MUSCL(CFL_factor, Rusanov(CFL_factor), Minmod()))

error_vec = convergence_error(scheme, equation, testcase, xmin, xmax, CFL_factor, Nx_vec; show=false)
display(plot_convergence_order(Nx_vec, error_vec, scheme, 2))