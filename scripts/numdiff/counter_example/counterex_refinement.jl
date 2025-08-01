using BenchmarkTools
include("../../../src/numdiff/include_file.jl")

# PARAMETRIZED INITIAL DATA

const a1 = 1.0
const b1 = 2.0
const a2 = 2.0
const b2 = 4.5

xmin, xmax = -b1/a1, b2/a2
M = 60
@show Nx = Int(floor((1 + floor(M/(a1*b2+a2*b1)))*(a2*b1+a1*b2)))
t0, tf = 0.0, 0.2
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

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
uexact = [uexact_fun(xi, tf) for xi in mesh.x]



# FINITE VOLUMES RESOLUTION
equation = Equation(OneD(), 1, Scalar(), Burgers(), u0)
sol = FiniteVolumes.solve(equation, params, Euler(), Roe(); log_config=LogConfig(true, true, true, false, false));

@show sol.niter

# using Plots
# plot(mesh.x, sol.uinit, label="uinit")
# plot!(mesh.x, uexact, label="exact")
# display(plot!(mesh.x, sol.u, label="t = "*string(sol.t)))

# QUANTIFICATION OF NUMERICAL DIFFUSION

function quantify_each_step(sol::Solution)
    umat = zeros(mesh.Nx, sol.niter)
    Dprio_mat = zero(umat)
    Dpost_mat = zero(umat)
    # plot(mesh.x, sol.uinit)
    plt1=plot(title="u")
    plt2=plot(title="D priori multidim")
    plt3=plot(title="D posteriori")
    for k in 1:sol.niter
        println("Timestep "*string(k)*" reaching t = "*string(sol.log.tlog[k]))
        priori_multidim = quantify_diffusion(sol, PrioriMultidim(MaxMD()), k)
        # posteriori = quantify_diffusion(sol, Posteriori(AsymmetricMD()), k)
        # posteriori = quantify_diffusion(sol, Posteriori(MaxMD()), k)
        posteriori = quantify_diffusion(sol, Priori(MaxMD()), k)

        println("min(M-m):")
        println(minimum(posteriori.M .- posteriori.m))
        println("min(L-l)")
        println(minimum(priori_multidim.L .- priori_multidim.l))
        println("sum(L)")
        println(sum(priori_multidim.L))
        println("sum(l)")
        println(sum(priori_multidim.l))
        umat[:,k] .= sol.log.ulog[k]
        Dprio_mat[:,k] .= priori_multidim.D
        Dpost_mat[:,k] .= posteriori.D
        plot!(plt1, mesh.x, sol.log.ulog[k], label="t="*string(round(sol.log.tlog[k],digits=2)))
        plot!(plt2, mesh.x, priori_multidim.D, label="t="*string(round(sol.log.tlog[k],digits=2)))
        plot!(plt3, mesh.x, posteriori.D, label="t="*string(round(sol.log.tlog[k],digits=2)))
    end
    # plot!(plt1, mesh.x, sol.log.ulog[end])
    display(plot!(plt1, mesh.x, uexact, label="exact"))
    display(plt2)
    display(plt3)
    umat, Dprio_mat, Dpost_mat
end

umat, Dprio_mat, Dpost_mat = quantify_each_step(sol);

animation(umat, sol.log.tlog, mesh.x, "u_"*string(Nx);d=5,filename="test1", fps=20, ylabel="u")
# animation(Dprio_mat, sol.log.tlog, mesh.x, "Dprio_"*string(Nx);filename="test2", d=5,fps=20, ylabel="D priori multidim")
# animation(Dpost_mat, sol.log.tlog, mesh.x, "Dprio_std_max_"*string(Nx);d=5, filename="test3",fps=20, ylabel="D priori std");


#estimate_priori = quantify_diffusion(sol, PrioriMultidim(MaxMD()), i=3; name="priori multidim");
# estimate_posteriori = quantify_diffusion(sol, Posteriori(AsymmetricMD()); name="asymmetric posteriori");
# estimate_prioristd = quantify_diffusion(sol, Priori(AsymmetricMD()); name="priori standard");

# CHECKING COUNTER EXAMPLE CONDITIONS

# println("min(M-m):")
# println(minimum(estimate_posteriori.M .- estimate_posteriori.m))
# println("min(L-l)")
# println(minimum(estimate_priori.L .- estimate_priori.l))
# println("sum(L)")
# println(sum(estimate_priori.L))
# println("sum(l)")
# println(sum(estimate_priori.l))

# plot(mesh.x, estimate_priori.l, label="l")
# display(plot!(mesh.x, estimate_priori.L, label="L"))

# plot(mesh.x, estimate_posteriori.m, label ="m")
# display(plot!(mesh.x, estimate_posteriori.M, label ="M"))

# plot(mesh.x, estimate_posteriori.D, label="posteriori")
# plot!(mesh.x, estimate_priori.D, label="priori multidim")
# # plot!(mesh.x, estimate_prioristd.D, label = "priori standard")
# plot!(mesh.x, estimate_priori.Dlow, label="Dlow")
# display(title!("Numerical Diffusion"))