struct PrioriDiffSol
    domain::Domain
    equation::Equation
    scheme::FVScheme
    modifiedDataType::ModifiedDataType
    boundsType::BoundsType
    u_approx
    dt_vec
    m_vec
    M_vec
    D_low
    D_up
    D_low_norm
    D_up_norm
    D_CL
    alpha
end

function compute_D_CL(D_low, D_up, D_base)
    alpha = sum(D_up .- D_base)/sum(D_up .- D_low)
    alpha * D_low + (1-alpha) * D_up, alpha
end

function diffusion_a_priori(u_init, domain::Domain, equation::Equation, scheme::FVScheme; modifiedDataType::ModifiedDataType=meanK(get_sL(scheme), get_sR(scheme)), boundsType::BoundsType=NormalBounds())
    
    Nx, dx = domain.Nx, domain.dx

    FVsol = fv_solve(domain, u_init, equation, scheme)
    u_approx, dt_vec = FVsol.u_approx, FVsol.dt_vec

    # source
    z = isnothing(domain.sourceVec) ? zeros((Nx, 1)) : reshape(domain.sourceVec, (domain.Nx,1))

    m_vec, M_vec = compute_G_bounds(u_approx[end-1], Nx, dx, dt_vec[end], equation, domain, scheme, modifiedDataType, boundsType)

    # D_low = [eta(equation, u_approx[end][i,:]; z=domain.sourceVec[i])[1] - eta(equation, u_approx[end-1][i,:]; z=domain.sourceVec[i])[1] for i in 1:length(u_approx[end-1][:,1])].+ dt_vec[end]/dx*(m_vec[begin+1:end] .- M_vec[begin:end-1])
    # D_up = [eta(equation, u_approx[end][i,:]; z=domain.sourceVec[i])[1] - eta(equation, u_approx[end-1][i,:]; z=domain.sourceVec[i])[1] for i in 1:length(u_approx[end-1][:,1])].+ dt_vec[end]/dx*(M_vec[begin+1:end] .- m_vec[begin:end-1])

    D_base = [eta(equation, u_approx[end][i,:], z[i])[1] - eta(equation, u_approx[end-1][i,:], z[i])[1] for i in 1:length(u_approx[end-1][:,1])]

    D_low = D_base .+ dt_vec[end]/dx*(m_vec[begin+1:end] .- M_vec[begin:end-1])
    D_up = D_base.+ dt_vec[end]/dx*(M_vec[begin+1:end] .- m_vec[begin:end-1])

    # D_low_norm = D_low * sum([eta(equation, u_approx[end][i,:]; z=domain.sourceVec[i])[1] - eta(equation, u_approx[end-1][i,:]; z=domain.sourceVec[i])[1] for i in 1:length(u_approx[end-1][:,1])]) / sum(D_low)
    # D_up_norm = D_up * sum([eta(equation, u_approx[end][i,:]; z=domain.sourceVec[i])[1] - eta(equation, u_approx[end-1][i,:]; z=domain.sourceVec[i])[1] for i in 1:length(u_approx[end-1][:,1])]) / sum(D_up)

    D_low_norm = D_low * sum(D_base) / sum(D_low)
    D_up_norm = D_up * sum(D_base) / sum(D_up)

    PrioriDiffSol(domain, equation, scheme, modifiedDataType, boundsType, u_approx, dt_vec, m_vec, M_vec, D_low, D_up, D_low_norm, D_up_norm, compute_D_CL(D_low, D_up, D_base)...)
end