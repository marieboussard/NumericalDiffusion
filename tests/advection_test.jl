include("../src/include_file.jl")

Nx = 50
t0, Tf = 0.0, 1.0
domain = createUnitInterval(Nx, t0, Tf)
eq = advection(0.2)
CFL_factor = 0.99
u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_gaussian(domain.x[i])] end; res)
display(plot(domain.x, u0))

# With classic Rusanov
method = Rusanov(CFL_factor)
fv_sol = fv_solve(domain, u0, eq, method)
plot_fv_sol(fv_sol; nb_plots=5)
display(title!(get_name(method)))

@show CFL_cond(eq, u0)

# With modified Rusanov
A_vec = LinRange(0, 0.75, 10)
for A in A_vec
    method = ModifiedRusanov(CFL_factor, A)
    fv_sol = fv_solve(domain, u0, eq, method)
    plot_fv_sol(fv_sol; nb_plots=5)
    display(title!(get_name(method)))
end