# Objective : find (numerically) bounds on A^T p by varying input parameters Gc and W.

using NumericalDiffusion


function ATp_on_grid(theta::AbstractVector, lamb::AbstractVector, params::Parameters, equation::Equation)
    NGc = length(theta)
    # NW = length(lamb)

    # res = zeros(NGc, NW, Nx)
    res = zeros(NGc, Nx)
    uinit, u = zeros(Nx), zeros(Nx)

    for i in 1:NGc
        sol = compute_entropic_G(params, equation; ent_numflux=ConvexCombinationG(theta[i]))
        res[i,:] .= (sol.A)'*sol.p_newt
        if i==1
            uinit .= sol.sol.uinit
            u .= sol.sol.u
        end
    end

    # for j in 1:NW
    # end
    res, uinit, u
end


# DOMAIN AND EQUATION

Nx = 50
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
equation = BurgersArticle


# GRID SPECIFICATIONS

# Gc is a convex combination: Gc(uj, uj+1) = θ Gc(uj) + (1 - θ) Gc(uj+1)
NGc = 100
theta = LinRange(0.0, 1.0, NGc)

# The weights are guided by a multiplicative factor between 1 and 2
NW = 5
lamb = LinRange(1.0, 2.0, NW)


# RESOLUTION

ATp, uinit, u = ATp_on_grid(theta, lamb, params, equation)


# VISUALIZATION

using Plots

plot(mesh.x, uinit, label="uinit")
display(plot!(mesh.x, u, label="u"))
heatmap(mesh.x, theta, ATp)