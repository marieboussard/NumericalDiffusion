function view_stencil(u, j, sL, sR)
    indices = mod1.(j - sL + 1 : j + sR, size(u, 1))  
    if ndims(u) == 1 
        return view(u,indices)
    else
        return view(u,indices, :)
    end
end

function view_stencil!(integrator::Integrator, j::Int)
    @unpack cache, params, uprev, fcont = integrator
    @unpack sL, sR, stencil = cache
    stencil .= mod1.(j - sL + 1 : j + sR, params.mesh.Nx) 
    nothing
    #=
    if ndims(uprev) == 1 
        return view(uprev,stencil), view(fcont, stencil)
    else
        return view(uprev,stencil, :), view(fcont,stencil, :)
    end
    =#
end


#view_stencil!(integrator::Integrator, j::Int) = view_stencil(integrator.uprev, j, compute_sL(integrator.time_scheme, integrator.space_scheme), compute_sR(integrator.time_scheme, integrator.space_scheme))

function numflux!(integrator::Integrator)
    for i ∈ 1:integrator.params.mesh.Nx
        view_stencil!(integrator, i)
        numflux!(integrator.time_scheme, integrator, i+1)
    end
    # for j in 1:integrator.equation.p
    #     integrator.fnum[end,j] = integrator.fnum[1,j]
    # end
    integrator.fnum[end,:] .= integrator.fnum[1,:]
    nothing
end
