# A serie of tests that need to be checked at each major modification of the code

include("../../src/include_file.jl")

# 2.1 # Quantifying numerical diffusion for solution of Burgers with Rusanov

Nx = 100
CFL_number = 0.5
domain = createInterval(Nx, -2, 2, 0, 0.4)
u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(domain.x[i])] end; res)
modifiedDataType = maxK()
#modifiedDataType = AsymmetricModifiedData()

@time sol = optimize_for_entropy(u0, domain, burgers(), Rusanov(CFL_number), modifiedDataType=modifiedDataType)

display(plot_solution(sol))
@show sol.Jopt

# # 2.2 # Quantifying numerical diffusion for solution of Burgers with Roe

# @time sol = optimize_for_entropy(u0, domain, burgers(), Roe(CFL_number), modifiedDataType=maxK())

# display(plot_solution(sol))
# @show sol.Jopt