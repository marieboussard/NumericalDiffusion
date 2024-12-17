include("../src/include_file.jl");

Nx = 10
xmin, xmax = 0.0, 1.0
CFL_factor = 0.5
equation = SaintVenant(flat_zb(height=0.0), 1e-10)
method = createHydrostatic(CFL_factor, Rusanov)
boxBounds=[1.0 10;-5.0 5.0]
sourceBounds=[-5.0, 5.0]

sol = iterate_WID(xmin, xmax, Nx, equation, method; nb_it=1, boxBounds=boxBounds, sourceBounds=sourceBounds)
@show sol.worstLowDiffVec

u, domain = correct_extend_initial_data(sol)
dx = domain.dx
dt = method.CFL_factor * dx / CFL_cond(equation, u)
@show domain.Tf
domain.Tf = dt

fv_sol_test = fv_solve(domain, u, equation, method)
display(plot_fv_sol(fv_sol_test, equation, nb_plots=2))
solEnt_test = optimize_for_entropy(u, domain, equation, method)
plot_solution(solEnt_test)
plot(solEnt_test.domain.interfaces, solEnt_test.Copt, label="Copt")