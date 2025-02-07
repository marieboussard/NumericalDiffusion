# A serie of tests that need to be checked at each major modification of the code

include("../src/include_file.jl")

using ProfileView

# 1.1 # Solving Burgers equation with Rusanov

xmin, xmax, Nx, t0, Tf = -2, 2, 10000, 0.0, 0.4
CFL_factor = 0.5
domain = Interval(Nx, xmin, xmax, t0, Tf)
equation = burgers()
testcase = ArticleTestcase()
scheme = FVScheme(Euler(), Rusanov(CFL_factor))
saveLog = false

u_init = initialData(domain, testcase)
uexact = exactData(domain, testcase)

problem = FVProblem(domain, equation, testcase, scheme, saveLog, u_init)
#VSCodeServer.@profview 

fv_sol = solve(problem);