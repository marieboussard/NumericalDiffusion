include("../src/include_file.jl")

# # 1 # With Burgers

# #xmin, xmax, Nx, t0, Tf = -2, 2, 5, 0, 0.4
# CFL_factor = 0.5
# #omega = createInterval(xmin, xmax, Nx, t0, Tf)

# Nx = 100

# #@time sol = iterate_WID(Nx, burgers(), Rusanov(CFL_factor), nb_it=10, boxBounds=[-3 3;])
# @time sol = iterate_WID(Nx, burgers(), Roe(CFL_factor), nb_it=10, boxBounds=[-3 3;])
# plotWorstWD(sol)

# 2 # Saint-Venant
Nx = 100
CFL_factor = 0.5
topoHeight = 2.0
eq = SaintVenant(Bump_zb(topoHeight), 1e-10)

#method = createHydrostatic(CFL_factor, Rusanov)
method = Rusanov(CFL_factor)

sol = iterate_WID(Nx, eq, method; nb_it=2, boxBounds=[0.1 3;-1.0 1.0], sourceBounds=[-1.0, 1.0])

plotWorstWD(sol, eq)