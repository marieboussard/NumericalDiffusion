using BenchmarkTools
include("../../../src/numdiff/include_file.jl")

# PARAMETRIZED INITIAL DATA

const a1 = 1.0
const b1 = 2.0
const a2 = 2.0
const b2 = 4.5

xmin, xmax = -b1/a1, b2/a2
t0, tf = 0.0, 0.2
CFL_factor = 0.5

# INITIAL CONDITION
u0(x::Real) = x<=0 ? -a1*x-b1 : -a2*x+b2
u0(x::AbstractArray) = u0.(x)

# EXACT SOLUTION
function uexact_fun(x::Real, t::Real)
    if t >= min(1/a1, 1/a2)
        @warn "Warning: This solution is not valid for t ≥ min(1/a1, 1/a2)"
    end
    if x <= -b1*t
        return (a1*x+b1) / (a1*t-1)
    elseif x >= b2*t
        return (a2*x-b2) / (a2*t-1)
    else
        return x/t
    end
end

equation = Equation(OneD(), 1, Scalar(), Burgers(), u0)

function refine_counterex(M_vec, xmin::Float64, xmax::Float64, t0::Float64, tf::Float64, CFL_factor::Float64, equation::Equation; mdtype_post=AsymmetricMD(), mdtype_prio=MaxMD())
    N = length(M_vec)
    eps1 = zeros(N)
    eps2 = zeros(N)
    eps3 = zeros(N)
    eps4 = zeros(N)
    for k in 1:N
        M = M_vec[k]
        Nx = Int(floor((1 + floor(M/(a1*b2+a2*b1)))*(a2*b1+a1*b2)))
        @show Nx
        mesh =  OneDMesh(Nx, xmin, xmax)
        params = Parameters(mesh, t0, tf, CFL_factor)
        sol = FiniteVolumes.solve(equation, params, Euler(), Roe();maxiter=1, log_config=LogConfig(true, true, true, false));
        estimate_prioristd = quantify_diffusion(sol, Priori(mdtype_post))
        estimate_priori = quantify_diffusion(sol, PrioriMultidim(mdtype_prio))
        eps1[k] = minimum(estimate_prioristd.M .- estimate_prioristd.m)
        eps2[k] = -maximum(estimate_priori.l .- estimate_priori.L)
        eps3[k] = sum(estimate_priori.L)
        eps4[k] = sum(estimate_priori.l)
        if eps1[k] < 0 && eps2[k]>=0 && eps3[k] >=0 && eps4[k]<=0
            println("For Nx = "*string(Nx)*", all counter example condition are satisfied")
            println(eps1[k])
            println(eps2[k])
            println(eps3[k])
            println(eps4[k])
        else
            println("For Nx = "*string(Nx)*", at least one counter example condition is not satisfied")
            println(eps1[k])
            println(eps2[k])
            println(eps3[k])
            println(eps4[k])
        end
    end
    eps1, eps2, eps3, eps4
end

# M_vec = [10, 20, 40, 80, 160, 320, 640, 1280]
M_vec = [10*(2^k) for k in 0:15]
#M_vec = [100]
eps1, eps2, eps3, eps4 = refine_counterex(M_vec, xmin, xmax, t0, tf, CFL_factor, equation, mdtype_post=MaxMD())

# flog(x::Real) = x>0 ? log(x) : -log(x)

display(scatter(log.(M_vec)./log(10), eps1, xlabel="log10(Nx)", title="min(M-m)"))
# savefig("images/counter_example/refinement/Mm.png")
display(scatter(log.(M_vec)./log(10), log.(eps2)./log(10), xlabel="log10(Nx)", title="log10(min(L-l))"))
# savefig("images/counter_example/refinement/Ll.png")
display(scatter(log.(M_vec)./log(10), eps3, xlabel="log10(Nx)", title="sum(L)"))
# savefig("images/counter_example/refinement/sumL.png")
display(scatter(log.(M_vec)./log(10), eps4, xlabel="log10(Nx)", title="sum(l)"))
# savefig("images/counter_example/refinement/sum_l.png")