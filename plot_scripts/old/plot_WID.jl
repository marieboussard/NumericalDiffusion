include("../src/include_file.jl")

nb_repetitions = 3
worstDiffFirstStep = zeros(nb_repetitions)
worstDiffLastStep = zeros(nb_repetitions)

Nx = 50
t0, Tf = 0.0, 0.5
CFL_factor = 0.5
eq = SaintVenant(flat_zb(height=0), 1e-10)
method = createHydrostatic(CFL_factor, Rusanov)
boxBounds=2*[0.1 2.0;-1.0 1.0] # We stay away from the dry states
sourceBounds=2*[-1.0, 1.0]

for k in 1:nb_repetitions

    sol = iterate_WID(Nx, eq, method; nb_it=nb_repetitions, boxBounds=boxBounds, sourceBounds=sourceBounds)

    # Reconstruction of the initial data from the optimization results
    u_init, z = extendInitialDataToLinear(sol, Nx, boxBounds=boxBounds, sourceBounds=sourceBounds)

    # Optimization at last timestep
    domain = createUnitInterval(Nx, t0, Tf)
    domain.sourceVec = z
    solEnt = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, g_tol=1e-10)#, modifiedDataType=maxK())
    #display(plot_solution(solEnt))
    worstDiffLastStep[k] = maximum(solEnt.Dopt)

    # Optimization at First timestep
    dt = sol.method.CFL_factor * domain.dx / CFL_cond(sol.equation, u_init) # Timestep given by CFL condition
    domain = createUnitInterval(Nx, t0, dt) # Redefining the domain with dt as final time
    domain.sourceVec = z
    solEnt = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, g_tol=1e-10)#, modifiedDataType=maxK())
    worstDiffFirstStep[k] = maximum(solEnt.Dopt)

end

plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
scatter!(worstDiffFirstStep, label="First step")
scatter!(worstDiffLastStep, label="Last step")

xlabel!("Repetitions")
#title!("Maximum")
display(ylabel!("Maximum Diffusion value"))