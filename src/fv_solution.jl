struct FVProblem
    domain::Domain
    equation::Equation
    testcase::Testcase
    scheme::FVScheme
    saveLog::Bool
    u_init
end

FVProblem(domain::Domain, equation::Equation, testcase::Testcase, scheme::FVScheme, u_init) = FVProblem(domain, equation, testcase, scheme, false, u_init)

mutable struct FVLog
    u_log
    dt_log
    t_log
end

mutable struct FVSolution
    problem::FVProblem
    u_approx
    Nt::Int
    t
    dt
    log::Union{Nothing, FVLog}
end

# FVSolution(problem, domain, equation, scheme, u_init) = FVSolution(problem, domain, equation, scheme, u_init, copy(u_init), 0, domain.t0, 0.0, false, nothing)

############### TO DO ####################################
# FIX THE PROBLEM OF TYPE IN FVLOG

function initialize_FV(problem::FVProblem)
    @unpack domain, equation, scheme, saveLog, u_init = problem
    if saveLog
        FVSolution(problem, copy(u_init), 0, domain.t0, 0.0, FVLog([copy(u_init)], typeof(t0)[], [t0]))
    else
        FVSolution(problem, copy(u_init), 0, domain.t0, 0.0, nothing)
    end
end

function performstep!(fv_sol::FVSolution)
    @unpack problem, u_approx, Nt, t = fv_sol
    @unpack domain, equation, scheme = problem
    @unpack dx, Tf = domain
    # Find the next time step with a CFL condition
    dt = min(compute_dt_with_CFL(scheme, equation, u_approx, dx), Tf - t)
    numericalFluxMat = vecNumFlux(equation.source, scheme, equation, collect(u_approx); dt=dt, domain=domain)
    u_approx .= u_approx .- dt / domain.dx * (numericalFluxMat[2:end,:] .- numericalFluxMat[1:end-1,:])
    fv_sol.u_approx = u_approx
    fv_sol.dt = dt
    fv_sol.t += dt
    fv_sol.Nt += 1
end

function saveStep!(fv_sol::FVSolution)
    @unpack log = fv_sol
    if isnothing(log)
        @error "Trying to update log, but it has not been initialized."
    end
    @unpack u_log, dt_log, t_log = log
    @unpack u_approx, dt, t = fv_sol
    push!(u_log, copy(u_approx))
    push!(dt_log, dt)
    push!(t_log, t)
end

function solve(problem::FVProblem)
    fv_sol = initialize_FV(problem)
    while fv_sol.t < Tf
        performstep!(fv_sol)
        if fv_sol.problem.saveLog
            saveStep!(fv_sol)
        end
    end
    fv_sol
end



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

# function scheme_step(ns::NullSource, v, dt, domain::Domain, equation::Equation, scheme::FVScheme)
#     #@code_warntype giveNumFlux(ns, method, equation, v)
#     # numericalFluxMat = giveNumFlux(ns, method, equation, v)
#     # v - dt / domain.dx * (numericalFluxMat[2:end,:] - numericalFluxMat[1:end-1,:])
#     #next_timestep(scheme.timeScheme, v, dt, domain, equation, scheme.spaceScheme)
#     #next_timestep(scheme, v, dt, domain, equation)
#     numericalFluxMat = vecNumFlux(ns, scheme, equation, collect(v); dt=dt, domain=domain)
#     v .- dt / domain.dx * (numericalFluxMat[2:end,:] .- numericalFluxMat[1:end-1,:])
# end

# function scheme_step(zb::ZbSource, v, dt, domain::Domain, equation::Equation, scheme::FVScheme;kwargs...)
#     #@show sourceTerm(equation, scheme, domain, v)
#     numericalFluxMat = vecNumFlux(zb, scheme, equation, v; dt=dt, domain=domain, kwargs...)
#     v - dt / domain.dx * (numericalFluxMat[2:end,:] - numericalFluxMat[1:end-1,:]) + dt * sourceTerm(equation, scheme.spaceScheme, domain, v; kwargs...)
# end


# function fv_solve(domain::Domain, u_init, equation::Equation, scheme::FVScheme)

#     t0, Tf, dx = domain.t0, domain.Tf, domain.dx
#     t = t0
#     Nt = 0

#     u_approx = [u_init]
#     dt_vec = Float64[]
#     t_vec = Float64[0.0]

#     while t < Tf

#         # Find the next time step with a CFL condition
#         dt = min(scheme.spaceScheme.CFL_factor * dx / CFL_cond(equation, u_approx[end]), Tf - t)

#         push!(dt_vec, dt)
#         push!(t_vec, t + dt)
        
#         #@show size(scheme_step(equation.source, u_approx[end], dt, domain, equation, method))
#         push!(u_approx, scheme_step(equation.source, u_approx[end], dt, domain, equation, scheme))

#         t += dt
#         Nt += 1

#     end

#     FVSolution(domain, equation, scheme, u_approx, Nt, dt_vec, t_vec)

# end


function plot_fv_sol(sol::FVSolution, nb_plots::Int64=2)

    @unpack problem, Nt, log = sol
    @unpack domain = problem
    @unpack x = domain

    @show Nt

    if isnothing(log)
        @error "Trying to plot hisory, but it has not been saved during solving. Put saveLog=true to allow saving."
    end

    @unpack u_log, dt_log, t_log = log
    p = div(Nt, nb_plots)

    plt = plot()

    for k in 0:nb_plots-2
        plot!(x, u_log[k*p+1], label="t = " * string(round(t_log[k*p+1], sigdigits=2)))
    end
    plot!(x, u_log[end], label="t = " * string(round(t_log[end], sigdigits=2)))
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

# function plot_fv_sol(sol::FVSolution, exact_sol::Base.Callable)

#     x = sol.domain.x

#     u_exact = [exact_sol(xi, sol.domain.Tf) for xi in x]

#     plot(x, sol.u_approx[end], label=get_name(sol.scheme))
#     plot!(x, u_exact, label="Exact")
#     xlabel!("x")
#     ylabel!("u")


# end

function plot_fv_sol(sol::FVSolution)

    @unpack u_approx = sol
    @unpack domain, testcase, scheme = sol.problem
    @unpack x = domain

    u_exact = exactData(domain, testcase)

    plot(x, u_approx, label=get_name(scheme))
    plot!(x, u_exact, label="Exact")
    xlabel!("x")
    ylabel!("u")

end