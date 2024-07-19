include("../src/include_file.jl")

# 2 # Saint Venant

Nx, t0, Tf = 10, 0, 0.4
CFL_factor = 0.5
domain = createUnitInterval(Nx, t0, Tf)
height_bump = 0.7
width_bump=0.5
source_term = bump_zb(height=height_bump, width=width_bump)
eq = SaintVenant(source_term, 1e-10)
addSource!(eq.source, domain)

#method = Rusanov(CFL_factor)
method = createHydrostatic(CFL_factor, Rusanov)

#v0 = v0_lake_at_rest(domain.x, source_term)
v0 = v0_lake_at_rest_perturbated(domain.x, source_term)
#v0 = v0_lake_at_rest(domain.x, NullSource())

solSV = fv_solve(domain, v0, eq, method)
nb_plots = 10
display(plot_fv_sol(solSV, solSV.equation; nb_plots=nb_plots))

solEnt = optimize_for_entropy(v0, domain, eq, method)
solEntAbs = optimize_for_entropy(v0, domain, eq, method, optimFunctional=AbsMinFun())
solEntSqrt = optimize_for_entropy(v0, domain, eq, method, optimFunctional=SqrtMinFun())
solEntAbs.label = "Abs"
solEnt.label = "Square"
solEntSqrt.label = "Sqrt"
display(plot_solutions([solEnt, solEntAbs, solEntSqrt]))