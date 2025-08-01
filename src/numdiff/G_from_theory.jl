# ENTROPIC FLUX G (KNOWN FROM THEORY) FOR RUSANOV SCHEME
function G_from_theory(::Rusanov, equation::Equation, params::Parameters, uinit::AbstractVector)
    @unpack mesh = params
    @unpack Nx = mesh
    entfun = entropy(equation.funcs)
    etavec = zero(mesh.x)
    Gvec = zero(mesh.x)
    eta!(entfun,  uinit, etavec)
    G!(entfun,  uinit, Gvec)
    Gexact = zero(mesh.x)
    cfl_loc = CFL_local(equation.dim, equation.eqtype, equation.funcs, params, uinit)
    for i in 1:Nx
        Gexact[i] = 0.5*(Gvec[i] + Gvec[mod1(i+1,Nx)]) - cfl_loc[i]*0.5*(etavec[mod1(i+1,Nx)] - etavec[i])
    end
    Gexact
end

G_from_theory(estimate::DiffEstimate) = G_from_theory(estimate.space_scheme, estimate.equation, estimate.params, estimate.uinit)