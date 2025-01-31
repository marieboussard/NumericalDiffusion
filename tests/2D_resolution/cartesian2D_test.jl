include("../../src/include_file.jl")

Nx, Ny, xmin, xmax, ymin, ymax, t0, Tf = 3, 3, -1, 1, -1, 1, 0, 0.5
domain = CartesianMesh(Nx, Ny, xmin, xmax, ymin, ymax, t0, Tf)