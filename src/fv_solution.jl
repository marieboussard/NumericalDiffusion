# include("utils.jl")
# include("domain.jl")
# include("method.jl")
# include("source_term.jl")


struct FVSolution
    domain::Domain
    equation::Equation
    method::FVMethod
    u_approx
    Nt::Int
    dt_vec
    t_vec
end

function scheme_step(::NullSource, v, dt, domain::Domain, equation::Equation, method::FVMethod)
    Nx, p = length(v[:,1]), get_unknowns_number(equation)
    numericalFluxMat = zeros(Nx+1, p)
    #Nx = length(v)
    #numericalFluxVec = Vector{eltype(v)}(undef, Nx + 1)#zeros(eltype(v), Nx + 1)
    #vcat(numFlux(method, equation, v[end], v[1]), [numFlux(method, equation, v[j], v[j+1]) for j in 1:Nx-1], numFlux(method, equation, v[end], v[1]))
    for i ∈ 2:Nx
        numericalFluxMat[i,:] = giveNumFlux(method, equation, v[i-1,:], v[i,:])
    end
    numericalFluxMat[1,:] = giveNumFlux(method, equation, v[end,:], v[1,:])
    numericalFluxMat[end,:] = numericalFluxMat[1,:]
    
    v - dt / domain.dx * (numericalFluxMat[2:end,:] - numericalFluxMat[1:end-1,:])
end

function scheme_step(zbSource::ZbSource, v, dt, domain::Domain, equation::Equation, method::FVMethod)
    
    Nx, p = length(v[:,1]), get_unknowns_number(equation)
    numericalFluxMat = zeros(Nx+1, p)
    #numericalFluxVec = Vector{eltype(v)}(undef, Nx + 1)
    for i ∈ 2:Nx
        v[i-1,:], v[i,:]
        numericalFluxMat[i,:] = giveNumFlux(method, equation, v[i-1,:], v[i,:]; zL=domain.sourceVec[i-1], zR=domain.sourceVec[i])
    end
    numericalFluxMat[1,:] = giveNumFlux(method, equation, v[end,:], v[1,:]; zL=domain.sourceVec[end], zR=domain.sourceVec[1])
    numericalFluxMat[end,:] = numericalFluxMat[1,:]
    numericalFluxMat
    v - dt / domain.dx * (numericalFluxMat[2:end,:] - numericalFluxMat[1:end-1,:]) + dt * sourceTerm(method, zbSource, domain, v)
end


function fv_solve(domain::Domain, u_init, equation::Equation, method::FVMethod)

    t0, Tf, dx = domain.t0, domain.Tf, domain.dx
    t = t0
    Nt = 0

    u_approx = [u_init]
    dt_vec = Float64[]
    t_vec = Float64[0.0]

    while t < Tf

        # Find the next time step with a CFL condition
        dt = min(method.CFL_factor * dx / CFL_cond(equation, u_approx[end]), Tf - t)

        push!(dt_vec, dt)
        push!(t_vec, t + dt)
        
        #@show size(scheme_step(equation.source, u_approx[end], dt, domain, equation, method))
        push!(u_approx, scheme_step(equation.source, u_approx[end], dt, domain, equation, method))

        t += dt
        Nt += 1

    end

    FVSolution(domain, equation, method, u_approx, Nt, dt_vec, t_vec)

end


function plot_fv_sol(sol::FVSolution; nb_plots::Int64=2)

    p = div(sol.Nt, nb_plots)

    plt = plot()

    for k in 0:nb_plots-2
        plot!(sol.domain.x, sol.u_approx[k*p+1], label="t = " * string(round(sol.t_vec[k*p+1], sigdigits=2)))
    end
    plot!(sol.domain.x, sol.u_approx[end], label="t = " * string(round(sol.t_vec[end], sigdigits=2)))
    xlabel!("x")
    ylabel!("u")

end

function plot_fv_sol(sol::FVSolution, equation::SaintVenant; nb_plots::Int64=2)

    p = div(sol.Nt, nb_plots)

    plt = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
    legendfontsize=14,
    titlefontsize=14,
    guidefontsize=14,
    tickfontsize=14)

    for k in 0:nb_plots-2
        #plot!(domain.x, [solSV.u_approx[k*p+1][i][1] for i in 1:Nx] .+ zb(Bump_zb(), domain.x), label="t = " * string(round(solSV.t_vec[k*p+1], sigdigits=2)))
        plot!(sol.domain.x, sol.u_approx[k*p+1][:,1] .+ zb(sol.equation.source, sol.domain.x), label="t = " * string(round(sol.t_vec[k*p+1], sigdigits=2)))
    end
    plot!(sol.domain.x, sol.u_approx[end][:,1] .+ zb(sol.equation.source, sol.domain.x), label="t = " * string(round(sol.t_vec[end], sigdigits=2)))
    plot!(sol.domain.x, zb(sol.equation.source, sol.domain.x), label="zb")
    xlabel!("x")
    title!(get_name(sol.method))
    display(ylabel!("Surface of the lake"))

end

function plot_fv_sol(sol::FVSolution, exact_sol::Base.Callable)

    x = sol.domain.x

    u_exact = [exact_sol(xi, sol.domain.Tf) for xi in x]

    plot(x, sol.u_approx[end], label="Approx")
    plot!(x, u_exact, label="Exact")
    xlabel!("x")
    ylabel!("u")


end