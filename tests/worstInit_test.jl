include("../src/include_file.jl")

# 1 # With Burgers

#xmin, xmax, Nx, t0, Tf = -2, 2, 5, 0, 0.4
CFL_factor = 0.5
#omega = createInterval(xmin, xmax, Nx, t0, Tf)

Nx = 100

@time sol = iterate_WID(Nx, burgers(), Rusanov(CFL_factor))
#@time solRoe = iterate_WID(Nx, burgers(), Roe(CFL_factor))

#worstData, worstLowDiff = find_worst_initial_data(Nx, burgers(), Rusanov(CFL_factor))