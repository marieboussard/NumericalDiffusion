include("../src/include_file.jl")

Nx, t0, Tf = 10, 0.0, 0.4
CFL_factor = 0.5
domain = create_unit_interval(Nx, t0, Tf)
topography = FlatTopo()
equation = SaintVenant(topography, 1e-10)
testcase = LakeAtRest(topography)

u_init = initialData(domain, testcase)
#uexact = exactData(domain, testcase)
equation = burgers()
saveLog = true

scheme = FVScheme(Euler(), HR(CFL_factor, Rusanov(CFL_factor)))

problem = FVProblem(domain, equation, testcase, scheme, saveLog, u_init)

fv_sol = solve(problem)
#plot_fv_sol(fv_sol, 4)