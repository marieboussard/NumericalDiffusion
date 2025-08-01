using Optimization
using OptimizationBBO
# using Optim
# using OptimizationMOI, Ipopt
# using ForwardDiff, Zygote
# using ModelingToolkit
# using JuMP, OSQP
using BenchmarkTools
include("../../../src/numdiff/include_file.jl")

struct CounterExample_WID
    nb_it::Int
    params
    eps
end

# const N = 50
# const T0 = 0.0
# const Tf = 0.2
# const CFLconst = 0.5

# mutable struct OptimizerBuffer
#     N::Int
#     t0::Float64
#     tf::Float64
#     CFL::Float64
#     mesh::meshType
#     params::paramsType
#     equation::eqtype
#     sol::soltype
#     estimate_priori::estimatePrioType
#     estimate_posteriori::estimatePostType
# end

function epsilon(a::AbstractVector, Nx::Int, CFL_factor::Float64, scheme::SpaceScheme, mdtype_prio::ModifiedDataType, mdtype_post::ModifiedDataType)
    a1, b1, a2, b2 = a
    xmin, xmax = -b1/a1, b2/a2
    mesh = OneDMesh(Nx, xmin, xmax)
    params = Parameters(mesh, 0.0, 0.2, CFL_factor)
    # INITIAL CONDITION
    u0(x::Real, α1::Real, β1::Real, α2::Real, β2::Real) = x<=0 ? -α1*x-β1 : -α2*x+β2
    u0(x::AbstractArray, α1::Real, β1::Real, α2::Real, β2::Real) = u0.(x, α1, β1, α2, β2)
    equation = Equation(OneD(), 1, Scalar(), Burgers(), x -> u0(x, a1, b1, a2, b2))
    sol = FiniteVolumes.solve(equation, params, Euler(), scheme; maxiter=1, log_config=LogConfig(true, false, true, false));
    estimate_prioristd = quantify_diffusion(sol, Priori(mdtype_post))
    estimate_priori = quantify_diffusion(sol, PrioriMultidim(mdtype_prio))
    eps1 = minimum(estimate_prioristd.M .- estimate_prioristd.m)
    eps2 = maximum(estimate_priori.l .- estimate_priori.L)
    eps3 = -sum(estimate_priori.L)
    eps4 = sum(estimate_priori.l)
    eps3 == 0.0 ? max(eps1, eps2, eps4) : max(eps1, eps2, eps3, eps4)
end

function find_counterex(Nx::Int, CFL_factor::Float64, scheme::SpaceScheme; mdtype_prio::ModifiedDataType=AsymmetricMD(), mdtype_post::ModifiedDataType=AsymmetricMD(), nb_it::Int=1, paramBound=10)
    minimizers_vec = zeros(nb_it, 4)
    minima_vec = zeros(nb_it)

    for k in 1:nb_it

    # Generating a set of random parameters

        @show init_params = [rand()*paramBound for it in 1:4]
        lower = [0.0 for it in 1:4]
        upper = [paramBound for it in 1:4]

        # Minimizing the cost epsilon

        prob = OptimizationProblem((x, p) -> epsilon(x, Nx, CFL_factor, scheme, mdtype_prio, mdtype_post), init_params, p=0.0, lb = lower, ub = upper)

        println("Starting resolution...")
        @time sol = Optimization.solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxtime=60, maxiters=100)#, abstol=1e-8
        minimizers_vec[k,:], minima_vec[k] = sol.minimizer, sol.minimum
        println("Resolution ended!")

    end

    i = argmin(minima_vec)

    CounterExample_WID(nb_it, minimizers_vec[i,:], minima_vec[i])
end

Nx = 102
CFL_factor = 0.5
scheme = Roe()
nb_it = 1
paramBound = 10

counter_ex = find_counterex(Nx, CFL_factor, scheme; mdtype_prio=MeanMD(), mdtype_post=MeanMD(), nb_it=nb_it, paramBound=paramBound)
# a = [1.0, 2.0, 2.0, 4.5]
# a = [8.85, 10.0, 8.85, 10.0]
# a = [0.005, 8.0, 1.3, 10.0]
# a = [0.0052025404873337015, 8.071319267808494, 1.3084796721409913, 9.999999999999982]
# @show epsilon(a, Nx, CFL_factor, Roe(), MeanMD(), AsymmetricMD())