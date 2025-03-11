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


#view_stencil!(integrator::Integrator, j::Int) = view_stencil(integrator.uprev, j, get_sL(integrator.time_scheme, integrator.space_scheme), get_sR(integrator.time_scheme, integrator.space_scheme))

# function numflux!(integrator::Integrator)
#     for i ∈ 1:integrator.params.mesh.Nx
#         view_stencil!(integrator, i)
#         numflux!(integrator.time_scheme, integrator, i+1)
#     end
#     for j in 1:integrator.equation.p
#         integrator.fnum[1,j] = integrator.fnum[end,j]
#     end
#     # integrator.fnum[end,:] .= integrator.fnum[1,:]
#     nothing
# end

function numflux!(integrator::Integrator)
    for j ∈ 1:integrator.params.mesh.Nx
        numflux!(integrator.time_scheme, integrator, j)
    end
    nothing
end

function numflux2D!(integrator::Integrator)
    for j ∈ 1:integrator.params.mesh.Nx
        for k ∈ 1:integrator.params.mesh.Ny
            numflux!(integrator.time_scheme, integrator, j, k)
        end
    end
    nothing
end

flux!(::OneD, eqfun::AbstractEquationFun, integrator::Integrator)= flux!(eqfun, integrator.u, integrator.fcont)
# flux!(::TwoD, eqfun::AbstractEquationFun, integrator::Integrator)= flux!(eqfun, integrator.u, integrator.fcont.fcont, integrator.fcont.hcont)
# flux!(::TwoD, eqfun::AbstractEquationFun, integrator::Integrator)= flux!(eqfun, integrator.u, selectdim(integrator.fcont, ndims(integrator.fcont), 1), selectdim(integrator.fcont, ndims(integrator.fcont), 2))

function flux!(::TwoD, eqfun::AbstractEquationFun, integrator::Integrator)
    # @show @allocated fcont = selectdim(integrator.fcont, ndims(integrator.fcont), 1)
    # @show @allocated hcont = selectdim(integrator.fcont, ndims(integrator.fcont), 2)
    fcont = view(integrator.fcont, :,:,:,1)
    hcont = view(integrator.fcont, :,:,:, 2)
    flux!(eqfun, integrator.u, fcont, hcont)
end