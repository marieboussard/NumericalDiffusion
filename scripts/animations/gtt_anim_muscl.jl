using NumericalDiffusion
using UnPack
#include("numdiff/animation.jl")
include("numdiff/animation2.jl")

# Domain definition
Nx = 50
xmin, xmax = -2, 2
t0, tf = 0.0, 0.2
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
equation = BurgersArticle

# Exact solution 
function uexact(x::Real, t::Real)
    if x <= -2*t
        return -(x+2)/(1-t)
    elseif x <= 3*t 
        return x/t 
    else
        return -3*(x-2)/(2-3*t)
    end
end

thinmesh = OneDMesh(2000, xmin, xmax)

time_scheme = RK2()
space_scheme = MUSCL(Rusanov(), Minmod())
sol = solve(equation, params, time_scheme, space_scheme; log_config=LogConfig(true, true, true, false, false));
sol2 = solve(equation, params, time_scheme, MUSCL(Rusanov(), Superbee()); log_config=LogConfig(true, true, true, false, false))

# QUANTIFICATION OF NUMERICAL DIFFUSION

function quantify_each_step(sol::Solution, sol2::Solution)

    @unpack mesh = sol.params

    umat = zeros(mesh.Nx, sol.niter)
    umat2 = zeros(mesh.Nx, sol.niter)
    uexmat = zeros(thinmesh.Nx, sol.niter)
    #Dprio_mat = zero(umat)
    Dpost_mat = zero(umat)
    Dpost_mat2 = zero(umat)
    # plot(mesh.x, sol.uinit)
    #plt1=plot(title="u")
    #plt2=plot(title="D priori multidim")
    #plt3=plot(title="D posteriori")
    for k in 1:sol.niter
        println("Timestep "*string(k)*" reaching t = "*string(sol.log.tlog[k]))
        #priori_multidim = quantify_diffusion(sol, PrioriMultidim(MaxMD()), k)
        posteriori = quantify_diffusion(sol, Posteriori(AsymmetricMD()), k)
        posteriori2 = quantify_diffusion(sol2, Posteriori(AsymmetricMD()), k)
        # posteriori = quantify_diffusion(sol, Posteriori(MaxMD()), k)
        # posteriori = quantify_diffusion(sol, Priori(MaxMD()), k)

        # println("min(M-m):")
        # println(minimum(posteriori.M .- posteriori.m))
        # println("min(L-l)")
        # println(minimum(priori_multidim.L .- priori_multidim.l))
        # println("sum(L)")
        # println(sum(priori_multidim.L))
        # println("sum(l)")
        # println(sum(priori_multidim.l))
        umat[:,k] .= sol.log.ulog[k]
        umat2[:,k] .= sol2.log.ulog[k]
        for j in 1:thinmesh.Nx
            uexmat[j,k] = uexact(thinmesh.x[j], sol.log.tlog[k])
        end
        #Dprio_mat[:,k] .= priori_multidim.D
        Dpost_mat[:,k] .= posteriori.D
        Dpost_mat2[:,k] .= posteriori2.D
        # plot!(plt1, mesh.x, sol.log.ulog[k], label="t="*string(round(sol.log.tlog[k],digits=2)))
        # #plot!(plt2, mesh.x, priori_multidim.D, label="t="*string(round(sol.log.tlog[k],digits=2)))
        # plot!(plt3, mesh.x, posteriori.D, label="t="*string(round(sol.log.tlog[k],digits=2)))
    end
    # plot!(plt1, mesh.x, sol.log.ulog[end])
    # #display(plot!(plt1, mesh.x, uexact, label="exact"))
    # #display(plt2)
    # display(plt3)
    umat, umat2, uexmat, Dpost_mat, Dpost_mat2
end

umat, umat2, uexmat, Dpost_mat, Dpost_mat2 = quantify_each_step(sol, sol2);

# pltA = []

# # plt1 = plot(size=(900, 600), margin=1Plots.cm, legend=:topleft,
# # legendfontsize=30,
# # titlefontsize=42,
# # guidefontsize=42,
# # tickfontsize=36)

# plt1 = plot()

# animation([umat, umat2, uexmat], sol.log.tlog, [mesh.x, mesh.x, thinmesh.x] ;d=5,filename="MUSCL_u", fps=20, ylabel="u", labels=["u_"*string(Nx)*" Minmod", "u_"*string(Nx)*" Superbee", "exact"], colors=["blue", "red", "green"])
# # animation(Dprio_mat, sol.log.tlog, mesh.x, "Dprio_"*string(Nx);filename="test2", d=5,fps=20, ylabel="D priori multidim")
# animation([Dpost_mat, Dpost_mat2], sol.log.tlog, [mesh.x, mesh.x] ;d=5, filename="MUSCL_D",fps=20, ylabel="D posteriori", labels= ["Dpost_"*string(Nx)*" Minmod", "Dpost_"*string(Nx)*" Superbee"], colors=["blue", "red"]);

animation_makie([umat, umat2, uexmat], sol.log.tlog, [mesh.x, mesh.x, thinmesh.x] ;d=5,filename="mMUSCL_u", fps=20, ylabel="u", labels=["u_"*string(Nx)*" Minmod", "u_"*string(Nx)*" Superbee", "exact"], colors=["blue", "red", "green"])

animation_makie([Dpost_mat, Dpost_mat2], sol.log.tlog, [mesh.x, mesh.x] ;d=5, filename="mMUSCL_D",fps=20, ylabel="D posteriori", labels= ["Dpost_"*string(Nx)*" Minmod", "Dpost_"*string(Nx)*" Superbee"], colors=["blue", "red"]);


#gif(anim, "./images/GTT/anim_"*filename*".gif"; fps = fps)


