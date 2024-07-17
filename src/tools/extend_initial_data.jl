function extendInitialDataToConstant(wd::WorstData, Nx::Int)
    i = argmin(wd.worstLowDiffVec)
    p = get_unknowns_number(wd.equation)
    sL, sR = get_sL(wd.method), get_sR(wd.method)
    u_init, z = zeros(Nx, p), zeros(Nx, 1)
    j = Int(round(Nx/2))

    for k in 1:Nx
        if k<j-sL 
            u_init[k,:] = sol.worstDataMat[i,begin,:]
            z[k] = sol.worstSource[i][begin]
        elseif j-sL <= k && k <= j+sR 
            u_init[k,:] = sol.worstDataMat[i,k-(j-sL)+1,:]
            z[k] = sol.worstSource[i][k-(j-sL)+1]
        else
            u_init[k,:] = sol.worstDataMat[i,end,:]
            z[k] = sol.worstSource[i][end]
        end
    end

    u_init, z
end

function extendInitialDataToLinear(wd::WorstData, Nx::Int; boxBounds=nothing, sourceBounds=[-1.0 1.0])
    i = argmin(wd.worstLowDiffVec)
    p = get_unknowns_number(wd.equation)
    sL, sR = get_sL(wd.method), get_sR(wd.method)
    s_init, u_init, z = zeros(Nx, p), zeros(Nx, p), zeros(Nx, 1)
    j = Int(round(Nx/2))
    dx = wd.domain.dx

    if isnothing(boxBounds)
        boxBounds = zeros((p,2))
        boxBounds[:,2] = 1.0
    end

    # Reconstruction on the water height
    left_slope = (wd.worstDataMat[i,begin+1,:] .- wd.worstDataMat[i,begin,:])./dx
    right_slope = (wd.worstDataMat[i,end,:] .- wd.worstDataMat[i,end-1,:])./dx

    # # Reconstruction on the surface height
    # left_slope = (wd.worstDataMat[i,begin+1,:] .+ wd.worstSource[i][begin+1] .- wd.worstDataMat[i,begin,:] .- wd.worstSource[i][begin])./dx
    # right_slope = (wd.worstDataMat[i,end,:] .+ wd.worstSource[i][end] .- wd.worstDataMat[i,end-1,:] .- wd.worstSource[i][end-1])./dx

    # Reconstruction of the topography
    source_left_slope = (wd.worstSource[i][begin+1] .- wd.worstSource[i][begin])./dx
    source_right_slope = (wd.worstSource[i][end] .- wd.worstSource[i][end-1])./dx

    for k in 1:Nx
        if k<j-sL 
            u_init[k,:] = min.(max.(wd.worstDataMat[i,begin,:] - left_slope*dx*(j-sL-k), boxBounds[:,1]), boxBounds[:,2])
            #s_init[k,:] = min.(max.(wd.worstDataMat[i,begin,:] - left_slope*dx*(j-sL-k), boxBounds[:,1]), boxBounds[:,2])
            z[k] = min.(max.(wd.worstSource[i][begin] - source_left_slope*dx*(j-sL-k), sourceBounds[1]), sourceBounds[2])
        elseif j-sL <= k && k <= j+sR 
            u_init[k,:] = wd.worstDataMat[i,k-(j-sL)+1,:]
            z[k] = wd.worstSource[i][k-(j-sL)+1]
        else
            u_init[k,:] = min.(max.(wd.worstDataMat[i,end,:] + right_slope*dx*(k-j-sR), boxBounds[:,1]), boxBounds[:,2])
            z[k] = min.(max.(wd.worstSource[i][end] + source_right_slope*dx*(k-j-sR), sourceBounds[1]), sourceBounds[2])
        end
    end

    u_init, z
end

function extendInitialDataToK(wd::WorstData, Nx::Int)
    i = argmin(wd.worstLowDiffVec)
    p = get_unknowns_number(wd.equation)
    sL, sR = get_sL(wd.method), get_sR(wd.method)
    u_init, z = zeros(Nx, p), zeros(Nx, 1)
    j = Int(round(Nx/2))

    @show wd.worstDataMat[i,:,:]
    @show extractLocalData(wd.worstDataMat[i,:,:], sL+1, sL, sR)

    @show wd.worstSource[i]

    @show K = computeK(wd.modifiedDataType, extractLocalData(wd.worstDataMat[i,:,:], sL+1, sL, sR))
    Z = computeK(wd.modifiedDataType, extractLocalData(reshape(wd.worstSource[i], (sL+sR+1,1)), sL+1, sL, sR))

    for k in 1:Nx
        if j-sL ≤ k ≤ j+sR
            u_init[k,:] = wd.worstDataMat[i,k-j+sL+1,:]
            z[k,:] .= wd.worstSource[i][k-j+sL+1]
        else
            u_init[k,:] = K
            z[k,:] = Z 
        end
    end

    u_init, z
end