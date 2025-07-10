include("../../src/numdiff/include_file.jl")
include("../../src/optimisation/optimisation_process.jl")
using CairoMakie


function testhyp(Nxvec::Vector{Int}, equation::Equation; xmin::Real=-2, xmax::Real=2, t0::Real=0.0, tf::Real=0.5, CFL_factor::Real=0.5, weights_type=AbsWeights(), beta::Real=2, kwargs...)

    N = length(Nxvec)

    res = zeros(N)

    scheme_name=""

    for k in 1:N

        Nx = Nxvec[k]
        mesh = OneDMesh(Nx, xmin, xmax)
        params = Parameters(mesh, t0, tf, CFL_factor)

        uzanewtsol = compute_entropic_G(params, equation; weights=weights_type, kwargs...)

        if k==1
            scheme_name = get_name(uzanewtsol.sol.space_scheme) # getting numerical scheme
        end

        @unpack W, Gc, G_newt = uzanewtsol
        
        H = zero(eltype(Gc))
        for i in 1:N
            H += W[i,i]*(G_newt[i] - Gc[i])
        end

        res[k] = (mesh.dx)^beta * H

    end

    res, scheme_name
end

Nxvec = [10, 100, 300, 500]
equation = BurgersArticle

beta = 0
space_scheme = Roe()

res, scheme_name = testhyp(Nxvec, equation; beta=beta, maxiter_uzawa=10000, maxiter_newton=10000, space_scheme = space_scheme, use_newton=false)

f = Figure()
ax = Axis(f[1,1], title = scheme_name, xlabel="Nx", ylabel="H(Gopt, dx)")
scatter!(ax, Nxvec, res, color = :tomato)
f