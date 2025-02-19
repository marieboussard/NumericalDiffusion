initial_data(domain::Domain, testcase::Testcase) = (res=zeros(domain.Nx, 1); for i in 1:domain.Nx res[i,:]=[u0_fun(testcase, domain.x[i])] end; res)

exact_data(domain::Domain, testcase::Testcase) = (res=zeros(domain.Nx, 1); for i in 1:domain.Nx res[i,:]=[uexact_fun(testcase, domain.x[i], domain.Tf)] end; res)

# Creating an interval when one wants exactly one timestep to be performed
function single_timestep_interval(Nx::Int, t0::Real, xmin::Real, xmax::Real, equation::Equation, testcase::Testcase, CFL_factor::Real)
    # First creating a domain with arbitrary final time
    domain = Interval(Nx, xmin, xmax, t0, 1.0)
    u0 = initialData(domain, testcase)
    dt = CFL_factor * domain.dx / CFL_cond(equation, u0)
    domain.Tf = dt
    domain, u0
end

function compute_dt_with_CFL(scheme, equation, u, dx)
    scheme.spaceScheme.CFL_factor * dx / CFL_cond(equation, u)
end