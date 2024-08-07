struct PrioriDiffSol
    domain::Domain
    equation::Equation
    method::FVMethod
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
end

function diffusion_a_priori(u_init, domain::Domain, equation::Equation, method::FVMethod; modifiedDataType::ModifiedDataType=meanK(get_sL(method), get_sR(method)), boundsType::BoundsType=NormalBounds())
    
    Nx, dx = domain.Nx, domain.dx

    FVsol = fv_solve(domain, u_init, equation, method)
    u_approx, dt_vec = FVsol.u_approx, FVsol.dt_vec

    m_vec, M_vec = compute_G_bounds(u_approx[end-1], Nx, dx, dt_vec[end], equation, domain, method, modifiedDataType, boundsType)

    D_low = [get_eta(equation, u_approx[end][i,:]; z=domain.sourceVec[i])[1] - get_eta(equation, u_approx[end-1][i,:]; z=domain.sourceVec[i])[1] for i in 1:length(u_approx[end-1][:,1])].+ dt_vec[end]/dx*(m_vec[begin+1:end] .- M_vec[begin:end-1])
    D_up = [get_eta(equation, u_approx[end][i,:]; z=domain.sourceVec[i])[1] - get_eta(equation, u_approx[end-1][i,:]; z=domain.sourceVec[i])[1] for i in 1:length(u_approx[end-1][:,1])].+ dt_vec[end]/dx*(M_vec[begin+1:end] .- m_vec[begin:end-1])

    D_low_norm = D_low * sum([get_eta(equation, u_approx[end][i,:]; z=domain.sourceVec[i])[1] - get_eta(equation, u_approx[end-1][i,:]; z=domain.sourceVec[i])[1] for i in 1:length(u_approx[end-1][:,1])]) / sum(D_low)
    D_up_norm = D_up * sum([get_eta(equation, u_approx[end][i,:]; z=domain.sourceVec[i])[1] - get_eta(equation, u_approx[end-1][i,:]; z=domain.sourceVec[i])[1] for i in 1:length(u_approx[end-1][:,1])]) / sum(D_up)

    PrioriDiffSol(domain, equation, method, modifiedDataType, boundsType, u_approx, dt_vec, m_vec, M_vec, D_low, D_up, D_low_norm, D_up_norm)
end