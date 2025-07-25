"""
    global_numflux!(integrator::Integrator)

Return the numerical flux for all interfaces of the mesh.
Called by `performstep!`.

See also [`performstep!`](@ref).
"""
# function global_numflux!(integrator::Integrator)
#     for j ∈ 1:integrator.params.mesh.Nx
#         numflux!(integrator.time_scheme, integrator, j)
#     end
#     nothing
# end

# function global_numflux!(time_scheme::TimeScheme, space_scheme::SpaceScheme, params::Parameters, args...)
#     for j ∈ 1:params.mesh.Nx
#         numflux!(time_scheme, space_scheme, j, params, args...)
#     end
# nothing
# end

#global_numflux(integrator::Integrator) = global_numflux!()

# numflux!(integrator::Integrator) = numflux!(integrator.time_scheme, integrator.space_scheme, integrator.params, integrator.equation, integrator.cache, integrator.space_cache, integrator.fnum, integrator.fcont, integrator.uprev)

function numflux2D!(integrator::Integrator)
    for j ∈ 1:integrator.params.mesh.Nx
        for k ∈ 1:integrator.params.mesh.Ny
            numflux!(integrator.time_scheme, integrator, j, k)
        end
    end
    nothing
end


"""
    flux!(::OneD, eqfun::AbstractEquationFun, integrator::Integrator)

Return an array of size (size(`integrator.u`) * equation numbers) where the cell [I,j] contains the jth component of f(u[I]) (note that I can be a couple of indices).
Here, f is the continuous flux function.
"""
flux!(::OneD, eqfun::AbstractEquationFun, integrator::Integrator)= flux!(eqfun, integrator.u, integrator.fcont)
# flux!(::TwoD, eqfun::AbstractEquationFun, integrator::Integrator)= flux!(eqfun, integrator.u, integrator.fcont.fcont, integrator.fcont.hcont)
# flux!(::TwoD, eqfun::AbstractEquationFun, integrator::Integrator)= flux!(eqfun, integrator.u, selectdim(integrator.fcont, ndims(integrator.fcont), 1), selectdim(integrator.fcont, ndims(integrator.fcont), 2))

function flux!(::TwoD, eqfun::AbstractEquationFun, integrator::Integrator)
    # @show @allocated fcont = selectdim(integrator.fcont, ndims(integrator.fcont), 1)
    # @show @allocated hcont = selectdim(integrator.fcont, ndims(integrator.fcont), 2)
    fcont = view(integrator.fcont, :, :, :, 1)
    hcont = view(integrator.fcont, :, :, :, 2)
    flux!(eqfun, integrator.u, fcont, hcont)
end