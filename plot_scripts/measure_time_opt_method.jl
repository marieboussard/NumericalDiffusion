include("../src/include_file.jl")

# For bounds not very well constrained
Nx = 100
CFL_factor = 0.5
eq = SaintVenant(bump_zb(height=0.5, width=0.4), 1e-10)
domain = createUnitInterval(Nx, 0.0, 0.1)
method = createHydrostatic(CFL_factor, Rusanov)
addSource!(eq.source, domain)
u_init = v0_lake_at_rest(domain.x, eq.source)
to = TimerOutput()

nb_its = 10

for k in 1:nb_its

    #println("=======NELDER MEAD======")
    @timeit to "Nelder Mead" solNelderMead = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, modifiedDataType=maxK(), g_tol=1e-10)
    #@show solNelderMead.summary

    #println("======LBFGS========")
    @timeit to "LBFGS" solLBFGS = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=LBFGS(), autodiff=:forward, modifiedDataType=maxK(), g_tol=1e-10)
    #@show solLBFGS.summary

    #println("======NEWTON========")
    @timeit to "Newton" solNewton = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=Newton(), autodiff=:forward, modifiedDataType=maxK(), g_tol=1e-10)
    #@show solNewton.summary

    #println("======A PRIORI======")
    @timeit to "Priori" solPriori = diffusion_a_priori(u_init, domain, eq, method)

end

show(to)

# # For bounds well constrained
# Nx = 100
# CFL_factor = 0.5
# eq = SaintVenant(bump_zb(height=0.5, width=0.4), 1e-10)
# domain = createUnitInterval(Nx, 0.0, 0.1)
# method = createHydrostatic(CFL_factor, Rusanov)
# addSource!(eq.source, domain)
# u_init = v0_lake_at_rest_perturbated(domain.x, eq.source)
# to = TimerOutput()

# nb_its = 10

# for k in 1:nb_its
#     @timeit to "Nelder Mead" solNelderMead = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, g_tol=1e-10)
#     @timeit to "LBFGS" solLBFGS = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=LBFGS(), autodiff=:forward, g_tol=1e-10)
#     @timeit to "Newton" solNewton = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=Newton(), autodiff=:forward, g_tol=1e-10)
#     @timeit to "Priori" solPriori = diffusion_a_priori(u_init, domain, eq, method)
# end

# show(to)