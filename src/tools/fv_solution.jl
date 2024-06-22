include("utils.jl")
include("domain.jl")
include("method.jl")

struct FVSolution
    domain::Domain
    u_approx
    Nt::Int
    dt_vec::Vector{Real}
    t_vec::Vector{Real}
end

function scheme_step(v, dt, dx, equation::Equation, method::FVMethod)
    Nx = length(v)
    numericalFluxVec = vcat(numFlux(method, equation, v[end], v[1]), [numFlux(method, equation, v[j], v[j+1]) for j in 1:Nx-1], numFlux(method, equation, v[end], v[1]))
    v - dt / dx * (numericalFluxVec[2:end] - numericalFluxVec[1:end-1])
end


function fv_solve(domain::Domain, u_init, equation::Equation, method::FVMethod)

    t0, Tf, dx = domain.t0, domain.Tf, domain.dx
    t = t0
    Nt::Int = 0

    u_approx = Vector{Real}[]
    dt_vec = Real[]
    t_vec = Real[]
    push!(u_approx, apply(u_init, domain.x))

    while t < Tf

        dt = min(method.CFL_factor * dx / max(D_flux(equation, u_approx[end])...), Tf - t)

        push!(dt_vec, dt)
        push!(t_vec, t + dt)
        push!(u_approx, scheme_step(u_approx[end], dt, dx, equation, method))

        t += dt
        Nt += 1

    end

    FVSolution(domain, u_approx, Nt, dt_vec, t_vec)

end


function plot_fv_sol(sol::FVSolution)

    plot!(sol.domain.x, sol.u_approx[end])
    xlabel!("x")
    ylabel!("u")

end

function plot_fv_sol(sol::FVSolution, exact_sol::Base.Callable)

    x = sol.domain.x

    u_exact = [exact_sol(xi, sol.domain.Tf) for xi in x]

    plot(x, sol.u_approx[end], label="Approx")
    plot!(x, u_exact, label="Exact")
    xlabel!("x")
    ylabel!("u")


end