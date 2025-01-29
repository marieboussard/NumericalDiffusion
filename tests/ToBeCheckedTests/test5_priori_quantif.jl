# A serie of tests that need to be checked at each major modification of the code

include("../../src/include_file.jl")

# 4 # Quantifying numerical diffusion with an a priori method for Burgers

xmin, xmax, Nx, t0, Tf = -2, 2, 50, 0, 0.4
CFL_factor = 0.5
domain = createInterval(Nx, xmin, xmax, t0, Tf)
equation = burgers()
scheme = FVScheme(Euler(), Rusanov(CFL_factor))
testcase = ArticleTestcase()

u0 = initialData(domain, testcase)

solEnt = optimize_for_entropy(u0, domain, equation, scheme)
D_priori = diffusion_a_priori(u0, domain, equation, scheme)

D_low = D_priori.D_low_norm
D_up = D_priori.D_up_norm

plot(domain.x, D_low, label="D priori")
plot!(domain.x, solEnt.Dopt, label="Dopt")
xlabel!("x")
ylabel!("Numerical Diffusion")