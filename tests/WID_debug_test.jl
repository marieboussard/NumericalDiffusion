include("../src/include_file.jl");

Nx = 10
xmin, xmax = 0.0, 1.0
CFL_factor = 0.5
equation = SaintVenant(flat_zb(height=0.0), 1e-10)
method = createHydrostatic(CFL_factor, Rusanov)
boxBounds=[5.0 10;-5.0 5.0]
sourceBounds=[-5.0, 5.0]
nb_it = 1

sol = iterate_WID(xmin, xmax, Nx, equation, method; nb_it=1, boxBounds=boxBounds, sourceBounds=sourceBounds)

rw = sol.reducedWD
print(rw.eps_uz)
domain = rw.domain
uz_unk = rw.uz
equation = sol.equation
method = sol.method
modifiedDataType = sol.modifiedDataType
boundsType = NormalBounds()

p = get_unknowns_number(equation)
Nx = domain.Nx
sL, sR = get_sL(method), get_sR(method)
j = sL + sR + 2

# Reconstructing u from the given data
u = zeros(Nx,p)
u_unk = reshape(uz_unk[1:p*(sL+sR+1)], (sL+sR+1,p))
# for k in sR+1:sL+2*sR+1
#     u[k,:] .= u_unk[k-sR,:]
# end
K = computeK(modifiedDataType, extractLocalData(u_unk, sL+1, sL, sR))
for k in 1:Nx
    if sR+2 ≤ k ≤ sL+2*sR+2
        u[k,:] .= u_unk[k-sR-1,:]
    else
        u[k,:] = K
    end
end

# Reconstructing the source from the given data
if length(uz_unk) > p*(sL+sR+1)
    zVec = zeros(Nx,1)
    Z = computeK(modifiedDataType, extractLocalData(reshape(uz_unk[p*(sL+sR+1)+1:end,:], (sL+sR+1,1)), sL+1, sL, sR))
    for k in 1:Nx
        if sR+2 ≤ k ≤ sL+2*sR+2
            zVec[k,:] .= uz_unk[p*(sL+sR+1) + k-sR-1]
        else
            zVec[k,:] = Z
        end
    end
    domain.sourceVec = zVec
end

# dx = domain.dx
# dt = method.CFL_factor * dx / CFL_cond(equation, u)
# domain.Tf = 2*dt

fv_sol = fv_solve(domain, u, equation, method)
#display(plot_fv_sol(fv_sol, equation, nb_plots=5))
#@show fv_sol.dt_vec
@show tt = domain.t0
domain.Tf = domain.t0
for i in 1:length(fv_sol.dt_vec)
    domain.Tf += fv_sol.dt_vec[i]
    solEnt = optimize_for_entropy(u, domain, equation, method)
    plot_solution(solEnt)
    plot(solEnt.domain.interfaces, solEnt.Copt, label="Copt")
    global tt += fv_sol.dt_vec[i]
    @show tt
    @show maximum(solEnt.Dopt)
    @show maximum(solEnt.Copt)
end