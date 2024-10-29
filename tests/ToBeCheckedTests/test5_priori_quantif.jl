# A serie of tests that need to be checked at each major modification of the code

include("../../src/include_file.jl")

# 4 # Quantifying numerical diffusion with an a priori method for Burgers

xmin, xmax, Nx, t0, Tf = -2, 2, 50, 0, 0.4
CFL_factor = 0.5
domain = createInterval(xmin, xmax, Nx, t0, Tf)
eq = burgers()
method = Rusanov(CFL_factor)

u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(domain.x[i])] end; res)

solEnt = optimize_for_entropy(u0, domain, eq, method)
D_priori = diffusion_a_priori(u0, domain, eq, method)

D_low = D_priori.D_low_norm
D_up = D_priori.D_up_norm

plot(domain.x, D_low, label="D priori")
plot!(domain.x, solEnt.Dopt, label="Dopt")
xlabel!("x")
ylabel!("Numerical Diffusion")