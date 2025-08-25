using NumericalDiffusion
using UnPack
#include("../../src/optimisation/optimisation_process.jl")
using CairoMakie

"""
    testhyp(Nxvec::Vector{Int}, equation::Equation; xmin::Real=-2, xmax::Real=2, t0::Real=0.0, tf::Real=0.5, CFL_factor::Real=0.5, weights_type=AbsWeights(), beta::Real=2, display=false, kwargs...)

Test strong and relaxed hypothesis over a single timestep
"""
function testhyp(Nxvec::Vector{Int}, equation::Equation; xmin::Real=-2, xmax::Real=2, t0::Real=0.0, tf::Real=0.5, CFL_factor::Real=0.5, weights_type=AbsWeights(), beta::Real=2, display=false, kwargs...)

    N = length(Nxvec)

    strong_res = zeros(N)
    relaxed_res = zeros(N)

    scheme_name=""

    for k in 1:N

        Nx = Nxvec[k]
        mesh = OneDMesh(Nx, xmin, xmax)
        params = Parameters(mesh, t0, tf, CFL_factor)

        uzanewtsol = compute_entropic_G(params, equation; weights=weights_type, maxiter=1, kwargs...)

        if k==1
            scheme_name = get_name(uzanewtsol.sol.time_scheme)*" + "*get_name(uzanewtsol.sol.space_scheme) # getting numerical scheme
            if display
                
                f = Figure()
                ax = Axis(f[1, 1], title = "FV resolution", xlabel = "x", ylabel = "u",)
                lines!(ax, mesh.x, uzanewtsol.sol.u, color = :tomato, label = "u",)
                axislegend(position = :rb)
                save("images/relaxed_hypothesis/fv_test.png", f)
            end
        end

        @unpack W, Gc, G_newt = uzanewtsol
        
        H = zero(eltype(Gc))
        for i in 1:Nx
            H += W[i,i]*(G_newt[i] - Gc[i])
        end

        strong_res[k] = H
        relaxed_res[k] = (mesh.dx)^beta * H

    end

    strong_res, relaxed_res, scheme_name
end

#Nxvec = [10, 100, 300, 500]
#Nxvec = [10, 40, 70, 100, 130, 160]
Nxvec = [10*i for i in 2:20]
equation = BurgersArticle

beta = 1.5
time_scheme = Euler()
#space_scheme = Rusanov()
space_scheme = MUSCL(Rusanov(), Minmod())
#space_scheme = MUSCL(Rusanov(), Superbee())

bound_mode = SingleBound()

strong_res, relaxed_res, scheme_name = testhyp(Nxvec, equation; beta=beta, maxiter_uzawa=10000, maxiter_newton=10000, time_scheme=time_scheme, space_scheme = space_scheme, use_newton=true, display=false, bound_mode=bound_mode)


f = Figure()
ax = Axis(f[1,1], title = scheme_name, xlabel="Nx", ylabel="H(Gopt, dx)")
scatter!(ax, Nxvec, strong_res, color = :tomato, label="strong")
scatter!(ax, Nxvec, relaxed_res, color = :navy, label="relaxed, beta="*string(beta))
axislegend(position = :lt)
f