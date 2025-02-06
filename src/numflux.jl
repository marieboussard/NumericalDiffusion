function view_stencil(u, j, sL, sR)
    indices = mod1.(j - sL + 1 : j + sR, size(u, 1))  
    if ndims(u) == 1 
        return view(u,indices)
    else
        return view(u,indices, :)
    end
end

view_stencil(integrator::Integrator, j::Int) = view_stencil(integrator.uprev, j, integrator.time_cache.sL, integrator.time_cache.sR)

view_stencil(integrator::Integrator, j::Int, scheme::Scheme) = view_stencil(integrator.uprev, j, scheme.cache.sL, scheme.cache.sR)

function numflux(integrator::Integrator)
    for i ∈ 2:Nx+1
        integrator.flux[i,:] .= numflux(time_scheme, integrator, view_stencil(integrator, i-1))
    end
end