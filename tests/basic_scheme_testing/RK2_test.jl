include("../src/include_file.jl")

# Settings

xmin, xmax, Nx, t0, Tf = -2, 2, 50, 0, 0.4
CFL_factor = 0.5
testcase = ArticleTestcase()
domain = createInterval(Nx, xmin, xmax, t0, Tf)
u0 = initialData(domain, testcase)
#domain, u0 = createOneTimestepInterval(Nx, t0, xmin, xmax, equation, testcase, CFL_factor)

equation = burgers()
scheme = FVScheme(RK2(), Rusanov(CFL_factor))
#scheme2 = FVScheme(Euler(), Rusanov(CFL_factor))

# A posteriori quantification of numerical diffusion
modifiedDataType = AsymmetricModifiedData()
#modifiedDataType = meanK(2,2)
solEnt = optimize_for_entropy(u0, domain, equation, scheme; modifiedDataType = modifiedDataType)
# println("Euler")
# solEnt2 = optimize_for_entropy(u0, domain, equation, scheme2; modifiedDataType = modifiedDataType)

# Exact solution
u_exact = [uexact_fun(testcase, xi, domain.Tf) for xi in domain.x]

# Entropic numerical flux
Gexact = exactG(scheme, equation, solEnt.u_approx[end-1], solEnt.dt_vec[end], domain)

# If we give the exact flux as initial guess to the optimization procedure
#solEnt = optimize_for_entropy(u0, domain, equation, scheme; modifiedDataType = modifiedDataType, vectorInitGuess=Gexact)

Dexact = diffusion(solEnt, Gexact)

plot(domain.interfaces, Gexact, label = "Entropic G")
plot!(domain.interfaces, solEnt.Gopt, label = "Gopt")
plot!(domain.interfaces, solEnt.m_vec, label="m")
plot!(domain.interfaces, solEnt.M_vec, label="M")
# scatter!(domain.interfaces, solEnt2.m_vec, label="m 2")
# scatter!(domain.interfaces, solEnt2.M_vec, label="M 2")
display(title!("Numerical entropy flux"))

plot(domain.x, solEnt.u_approx[end-1], label=get_name(scheme))
#plot!(domain.x, solEnt2.u_approx[end-1], label=get_name(scheme2))
plot!(domain.x, u_exact, label="exact")
display(title!("u"))

plot(domain.x, Dexact, label="Dexact")
plot!(domain.x, solEnt.Dopt, label="Dopt")
#plot!(domain.x, solEnt2.Dopt, label="Dopt 2")
display(title!("Numerical Diffusion"))

@show solEnt.Jopt