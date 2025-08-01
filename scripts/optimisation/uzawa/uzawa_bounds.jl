using FiniteVolumes
using BenchmarkTools
using UnPack
using Plots
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")

function uzawa_bounds(estimate::DiffEstimate, optsol::UzawaSol)
    @unpack dt, l, L = estimate
    @unpack Nx, dx = estimate.params.mesh
    @unpack W, Gc = optsol
    m, M = zeros(Nx), zeros(Nx)
    ilambda = dx / dt
    CW, wtot = 0.0, 0.0
    for i in 1:Nx
        CW += Gc[i]*W[i,i]
        wtot += W[i,i]
    end
    m[Nx] = L[1]*W[1,1]
    M[Nx] = l[1]*W[1,1]
    # for i in 1:Nx-1
    #     for j in i+1:Nx-1
    #         m[Nx] += W[j,j]*L[i]
    #         M[Nx] += W[j,j]*l[i]
    #     end
    # end
    for i in 2:Nx
        for j in 1:i-1
            m[Nx] += W[j,j]*l[i]
            M[Nx] += W[j,j]*L[i]
        end
    end
    # # m[Nx] = 1.0/wtot * (CW - ilambda*m[Nx])
    # # M[Nx] = 1.0/wtot * (CW - ilambda*M[Nx])
    m[Nx] = 1.0/wtot * (CW + ilambda*m[Nx])
    M[Nx] = 1.0/wtot * (CW + ilambda*M[Nx])
    # # for k in 1:Nx-1
    # #     m[k] = m[mod1(k-1, Nx)] + ilambda*l[k]
    # #     M[k] = M[mod1(k-1, Nx)] + ilambda*L[k]
    # # end
    for k in Nx-1:-1:1
        m[k] = m[k+1] - ilambda*L[k+1]
        M[k] = M[k+1] - ilambda*l[k+1]
    end

    for k in 1:Nx-1
        m[k] = max(m[Nx] - ilambda*sum([L[i] for i in k+1:Nx]), m[Nx] + ilambda*sum([l[i] for i in 1:k]))
        M[k] = min(M[Nx] - ilambda*sum([l[i] for i in k+1:Nx]), M[Nx] + ilambda*sum([L[i] for i in 1:k]) )
    end

    # for k in 1:Nx
    #     m[k] = Nx * (Gc[k] * W[k,k]/wtot + ilambda * W[k,k]/wtot*l[k])*10000
    #     M[k] = Nx*(Gc[k] * W[k,k]/wtot + ilambda * W[k,k]/wtot*L[k])*10000
    # end
    m, M, CW, wtot
end

# Domain definition
Nx = 100
xmin, xmax = -4, 4
t0, tf = 0.0, 0.5
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
# equation = BurgersConstant
equation = BurgersArticle

# Finite volumes resolution
sol = solve(equation, params, Euler(), Rusanov(); maxiter=1, log_config=LogConfig(true, false, true, false, false));

# Multidimensional bounds for ΔG
estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
@unpack uinit, u, l, L = estimate

# Defining optimization components
Gc = zeros(eltype(u), Nx)
A = zeros(eltype(u), 2*Nx, Nx)
b = zeros(eltype(u), 2*Nx)
W = zeros(eltype(u), Nx, Nx)
alpha=1.0

Gflux!(CenteredG(), Gc, estimate)
fill_A!(A, estimate)
fill_b!(b, estimate)
fill_W!(AlphaWeights(alpha), W, estimate)

# Uzawa algorithm
# optsol = optimize_uzawa(Gexact, A, b; W=W, maxiter=100000, eps=1e-12);
optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=500000, eps=1e-12, eps_cons=1e-12);

m, M, CW, wtot = uzawa_bounds(estimate, optsol)

plot(mesh.x, m, label="m")
plot!(mesh.x, optsol.gamma_opt, label="G uzawa")
display(plot!(mesh.x, M, label="M"))

# plot(mesh.x, M.-optsol.gamma_opt, label="M-G")
du = [abs(sol.uinit[mod1(j+1, Nx)] - sol.uinit[j]) for j in 1:Nx]
# plot(mesh.x, du, label="|uj+1-uj|")

plot(mesh.x, M.-optsol.gamma_opt .- du)