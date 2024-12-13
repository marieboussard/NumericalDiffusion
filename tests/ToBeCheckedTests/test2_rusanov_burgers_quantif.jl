# A serie of tests that need to be checked at each major modification of the code

include("../../src/include_file.jl")

# 2.1 # Quantifying numerical diffusion for solution of Burgers with Rusanov

Nx = 100
CFL_number = 0.5
domain = createInterval(Nx, -2, 2, 0, 0.4)
testcase = ArticleTestcase()
u0 = initialData(domain, testcase)
scheme = FVScheme(Euler(), Rusanov(CFL_number))
#modifiedDataType = maxK()
modifiedDataType = AsymmetricModifiedData()

@time sol = optimize_for_entropy(u0, domain, burgers(), scheme, modifiedDataType=modifiedDataType)
@show sol.Jopt
display(plot_solution(sol))


# # 2.2 # Quantifying numerical diffusion for solution of Burgers with Roe

# @time sol = optimize_for_entropy(u0, domain, burgers(), Roe(CFL_number), modifiedDataType=maxK())

# display(plot_solution(sol))
# @show sol.Jopt