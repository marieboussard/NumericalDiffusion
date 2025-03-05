# FUNCTIONS TO SIMPLIFY THE NOTATIONS IN SOME SPECIFIC CASES

# SAINT-VENANT WITH HYDROSTATIC RECONSTRUCTION
function hrsolve(params::Parameters, time_scheme::TimeScheme, subscheme::SpaceScheme, z::Base.Callable, Dz::Base.Callable; kwargs...)
    equation = saintvenant_with_topo(z, Dz; sourcedisc = HRDisc())
    space_scheme = HR(subscheme)
    solve(equation, params, time_scheme, space_scheme; kwargs...)
end