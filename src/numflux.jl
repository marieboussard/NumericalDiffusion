function view_stencil(u, j, sL, sR)
    indices = mod1.(j - sL + 1 : j + sR, size(u, 1))  
    if ndims(u) == 1 
        return view(u,indices)
    else
        return view(u,indices, :)
    end
end

view_stencil!(integrator::Integrator, j::Int) = view_stencil(integrator.uprev, j, compute_sL(integrator.time_scheme, integrator.space_scheme), compute_sR(integrator.time_scheme, integrator.space_scheme))

function numflux!(integrator::Integrator)
    for i ∈ 2:integrator.params.mesh.Nx+1
        numflux!(integrator.time_scheme, integrator, view_stencil!(integrator, i-1), i)
    end
    integrator.flux[1,:] .= integrator.flux[end,:]
end