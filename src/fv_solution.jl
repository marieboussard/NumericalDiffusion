struct FVSolution
    domain::Domain
    equation::Equation
    scheme::FVScheme
    u_approx
    Nt::Int
    dt_vec
    t_vec
end

#=
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

function scheme_step(::ZbSource, v, dt, domain::Domain, equation::Equation, method::FVMethod)
    
    Nx, p = length(v[:,1]), get_unknowns_number(equation)
    numericalFluxMat = zeros(Nx+1, p)
    #numericalFluxVec = Vector{eltype(v)}(undef, Nx + 1)
    for i ∈ 2:Nx
        # @show v[i-1,:], v[i,:]
        # @show domain.sourceVec
        numericalFluxMat[i,:] = giveNumFlux(method, equation, v[i-1,:], v[i,:]; zL=domain.sourceVec[i-1], zR=domain.sourceVec[i])
    end
    numericalFluxMat[1,:] = giveNumFlux(method, equation, v[end,:], v[1,:]; zL=domain.sourceVec[end], zR=domain.sourceVec[1])
    numericalFluxMat[end,:] = numericalFluxMat[1,:]
    numericalFluxMat
    v - dt / domain.dx * (numericalFluxMat[2:end,:] - numericalFluxMat[1:end-1,:]) + dt * sourceTerm(equation, method, domain, v)
end

=#

########## TO DO : MODIFY IN ORDER TO UNCOMMENT ##########################################################""

# function localized_scheme_step(j::Int, ns::NullSource, v, dt, domain::Domain, equation::Equation, method::MixedMethod{T,U}) where T where U
#     #@show T, U
#     sL, sR = get_sL(method), get_sR(method)
#     v_local = extractExtendedLocalData(v, j, sL, sR)
#     numericalFluxMat_local = giveNumFlux(ns, method, equation, v_local)
#     numericalFluxMat = zeros(U, (domain.Nx+1, get_unknowns_number(equation)))
#     for k in j-sL-sR+1:j+sR+sL
#         # @show typeof(numericalFluxMat[mod1(k,Nx),:])
#         # @show typeof(numericalFluxMat_local[k-j+sL+sR,:])
#         numericalFluxMat[mod1(k,Nx),:] =  numericalFluxMat_local[k-j+sL+sR,:]
#         # numericalFluxMat[k,:] = giveNumFlux(method, equation, v_local[k-,:], v[i,:])
#         # numericalFluxMat[k] = giveNumFlux(ns, method, equation, v_local)
#     end
#     v - dt / domain.dx * (numericalFluxMat[2:end,:] - numericalFluxMat[1:end-1,:])
# end

function scheme_step(ns::NullSource, v, dt, domain::Domain, equation::Equation, scheme::FVScheme)
    #@code_warntype giveNumFlux(ns, method, equation, v)
    # numericalFluxMat = giveNumFlux(ns, method, equation, v)
    # v - dt / domain.dx * (numericalFluxMat[2:end,:] - numericalFluxMat[1:end-1,:])
    #next_timestep(scheme.timeScheme, v, dt, domain, equation, scheme.spaceScheme)
    #next_timestep(scheme, v, dt, domain, equation)
    numericalFluxMat = vecNumFlux(ns, scheme, equation, collect(v); dt=dt, domain=domain)
    v .- dt / domain.dx * (numericalFluxMat[2:end,:] .- numericalFluxMat[1:end-1,:])
end

function scheme_step(zb::ZbSource, v, dt, domain::Domain, equation::Equation, scheme::FVScheme;kwargs...)
    #@show sourceTerm(equation, scheme, domain, v)
    numericalFluxMat = vecNumFlux(zb, scheme, equation, v; dt=dt, domain=domain, kwargs...)
    v - dt / domain.dx * (numericalFluxMat[2:end,:] - numericalFluxMat[1:end-1,:]) + dt * sourceTerm(equation, scheme.spaceScheme, domain, v; kwargs...)
end


function fv_solve(domain::Domain, u_init, equation::Equation, scheme::FVScheme)

    t0, Tf, dx = domain.t0, domain.Tf, domain.dx
    t = t0
    Nt = 0

    u_approx = [u_init]
    dt_vec = Float64[]
    t_vec = Float64[0.0]

    while t < Tf

        # Find the next time step with a CFL condition
        dt = min(scheme.spaceScheme.CFL_factor * dx / CFL_cond(equation, u_approx[end]), Tf - t)

        push!(dt_vec, dt)
        push!(t_vec, t + dt)
        
        #@show size(scheme_step(equation.source, u_approx[end], dt, domain, equation, method))
        push!(u_approx, scheme_step(equation.source, u_approx[end], dt, domain, equation, scheme))

        t += dt
        Nt += 1

    end

    FVSolution(domain, equation, scheme, u_approx, Nt, dt_vec, t_vec)

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

function plot_fv_sol(sol::FVSolution, ::SaintVenant; nb_plots::Int64=2, plotMode::PlottingMode=DisplayMode())

    p = div(sol.Nt, nb_plots)

    pltA = []

    plt = createPlot(plotMode)
    
    # plt = plot(size=(750, 600), margin=0.5Plots.cm, legend=:bottomright,
    # legendfontsize=15,
    # titlefontsize=21,
    # guidefontsize=21,
    # tickfontsize=18)

    for k in 0:nb_plots-2
        #plot!(domain.x, [solSV.u_approx[k*p+1][i][1] for i in 1:Nx] .+ zb(Bump_zb(), domain.x), label="t = " * string(round(solSV.t_vec[k*p+1], sigdigits=2)))
        plot!(sol.domain.x, sol.u_approx[k*p+1][:,1] .+ sol.domain.sourceVec, label="t = " * string(round(sol.t_vec[k*p+1], sigdigits=4)), lw=2)
    end
    plot!(sol.domain.x, sol.u_approx[end][:,1] .+ sol.domain.sourceVec, label="t = " * string(round(sol.t_vec[end], sigdigits=4)), lw=3)
    plot!(sol.domain.x, sol.domain.sourceVec, label="zb", lw=2)
    xlabel!("x")
    title!(get_name(sol.scheme)*", Nx = "*string(sol.domain.Nx))
    ylabel!("Surface of the lake")

    push!(pltA, plt)

    plt2 = createPlot(plotMode)

    # plt2 = plot(size=(750, 600), margin=0.5Plots.cm, legend=:bottomright,
    # legendfontsize=15,
    # titlefontsize=21,
    # guidefontsize=21,
    # tickfontsize=18)

    for k in 0:nb_plots-2
        plot!(sol.domain.x, sol.u_approx[k*p+1][:,2], label="t = " * string(round(sol.t_vec[k*p+1], sigdigits=2)), lw=2)
    end
    plot!(sol.domain.x, sol.u_approx[end][:,2], label="t = " * string(round(sol.t_vec[end], sigdigits=2)), lw=2)
    xlabel!("x")
    #title!(get_name(sol.method))
    ylabel!("Water flow")

    push!(pltA, plt2)

    display(assemblePlot(plotMode, pltA))

    #display(plot(pltA..., layout=(2, 1), size=(1500, 1200)))
    #display(plot(pltA..., layout=(2, 1), size=(750, 600)))

end

function plot_fv_sol(sol::FVSolution, exact_sol::Base.Callable)

    x = sol.domain.x

    u_exact = [exact_sol(xi, sol.domain.Tf) for xi in x]

    plot(x, sol.u_approx[end], label=get_name(sol.scheme))
    plot!(x, u_exact, label="Exact")
    xlabel!("x")
    ylabel!("u")


end

function plot_fv_sol(sol::FVSolution, testcase::Testcase)

    x = sol.domain.x

    u_exact = exactData(sol.domain, testcase)

    plot(x, sol.u_approx[end], label=get_name(sol.scheme))
    plot!(x, u_exact, label="Exact")
    xlabel!("x")
    ylabel!("u")


end