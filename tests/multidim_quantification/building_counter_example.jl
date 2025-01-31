include("../../src/include_file.jl")

Nx = 10
CFL_factor = 0.5
a1, b1, a2, b2 = 1.0, 2.0, 2.0, 4.5
testcase = PiecewiseLinear(a1, b1, a2, b2)
xmin, xmax = spaceBounds(testcase)
@show Nx = integerNx(Nx, testcase)
equation = burgers()
scheme = FVScheme(Euler(), Roe(CFL_factor))
boundsType = NormalBounds()
#boundsType = LightBounds()

domain, u0 = createOneTimestepInterval(Nx, 0.0, xmin, xmax, equation, testcase, CFL_factor)
# Tf = 0.2
# domain = createInterval(Nx, xmin, xmax, 0.0, Tf)
# u0 = initialData(domain, testcase)

sol = fv_solve(domain, u0, equation, scheme)
display(plot_fv_sol(sol, testcase))

## Part One : We check that the scheme is non entropic for this initial data

modifiedDataType = AsymmetricModifiedData()
m_vec, M_vec = compute_G_bounds(sol.u_approx[end-1], Nx, domain.dx, sol.dt_vec[end], equation, domain, scheme, modifiedDataType, boundsType)

display(plot(domain.interfaces, M_vec .- m_vec))

@show minimum(M_vec .- m_vec)
@show domain.Tf

## Part Two : We check that it satisfies however the three criteria

#modifiedDataType = meanK_multidim(1, 1)
modifiedDataType = AsymmetricModifiedData()

#D_priori = diffusion_a_priori(u0, domain, equation, scheme)
D_priori_multidim = diffusion_a_priori_multidim(u0, domain, equation, scheme; modifiedDataType=modifiedDataType, boundsType=boundsType)

# m_vec, M_vec = D_priori.m_vec, D_priori.M_vec
# ll_vec = domain.Tf/domain.dx*(m_vec[begin+1:end] .- M_vec[begin:end-1])
# LL_vec = domain.Tf/domain.dx*(M_vec[begin+1:end] .- m_vec[begin:end-1])

l_vec, L_vec = D_priori_multidim.l_vec, D_priori_multidim.L_vec

@show minimum(L_vec .- l_vec)
println("Maximal value of D low: " * string(maximum(D_priori_multidim.D_low)))
println("Maximal value of D low norm: " * string(maximum(D_priori_multidim.D_priori)))
@show sum(L_vec[begin+1:end])
@show sum(l_vec[begin+1:end]);