#using Plots

mutable struct OptForEntropySol{T<:Real, S<:Domain{T}}
    domain::S
    equation::Equation
    scheme::FVScheme
    modifiedDataType::ModifiedDataType
    boundsType::BoundsType
    u_approx::Vector{Matrix{T}}
    dt_vec::Vector{T}
    Gopt
    Jopt::Real
    Dopt
    Copt
    m_vec
    M_vec
    summary
    label
end

function plot_solution(sol::OptForEntropySol, plotMode::PlottingMode=DisplayMode())

    plt = []

    plt1 = createPlot(plotMode)

    # plt1 = plot(size=(750, 600), margin=0.5Plots.cm, legend=:bottomright,
    # legendfontsize=15,
    # titlefontsize=21,
    # guidefontsize=21,
    # tickfontsize=18)
    xlabel!("x")
    ylabel!("Numerical Entropy Flux")
    title!(get_name(sol.scheme)*", Nx = "*string(sol.domain.Nx))
    plot!(sol.domain.interfaces, sol.m_vec, label="m", lw=2)
    plot!(sol.domain.interfaces, sol.Gopt, label="Gopt", lw=2)
    plot!(sol.domain.interfaces, sol.M_vec, label="M", lw=2)
    #xlims!((0.45, 0.55))

    push!(plt, plt1)

    plt2 = createPlot(plotMode)

    # plt2 = plot(sol.domain.x, sol.Dopt, label="Dopt", size=(750, 600), margin=0.5Plots.cm, legend=:bottomright,
    # legendfontsize=15,
    # titlefontsize=21,
    # guidefontsize=21,
    # tickfontsize=18, lw=2)
    plot!(sol.domain.x, sol.Dopt, label="Dopt")
    xlabel!("x")
    ylabel!("Numerical Diffusion")
    title!("Max Diff : "*string(maximum(sol.Dopt)))
    #ylims!((-0.001, 0.015))

    push!(plt, plt2)

    display(assemblePlot(plotMode, plt))

    #display(plot(plt..., layout=(2, 1), size=(1500, 1200)))
    println("Maximal diffusion value : ", maximum(sol.Dopt))
end

function plot_solution_EntBreak(sol::OptForEntropySol, plotMode::PlottingMode=DisplayMode())
    plt = []

    plt1 = createPlot(plotMode)

    xlabel!("x")
    ylabel!("Numerical Entropy Flux")
    title!(get_name(sol.scheme)*", Nx = "*string(sol.domain.Nx))
    plot!(sol.domain.interfaces, sol.m_vec, label="m", lw=2)
    plot!(sol.domain.interfaces, sol.Gopt, label="Gopt", lw=2)
    plot!(sol.domain.interfaces, sol.M_vec, label="M", lw=2)

    push!(plt, plt1)

    plt2 = createPlot(plotMode)

    plot!(sol.domain.x, sol.Dopt, label="Dopt", lw=2)
    xlabel!("x")
    ylabel!("Numerical Diffusion")
    title!("Max Diff : "*string(maximum(sol.Dopt)))

    push!(plt, plt2)

    plt3 = createPlot(plotMode)

    plot!(sol.domain.x, sol.Dopt, label="Dopt", lw=2)
    xlabel!("x")
    ylabel!("Numerical Diffusion")
    #title!("Max Diff : "*string(maximum(sol.Dopt)))
    title!("Zoom on positive diffusions")
    ylims!((-0.005, maximum(sol.Dopt)*1.2))

    push!(plt, plt3)

    plt4 = createPlot(plotMode)

    plot!(sol.domain.interfaces, sol.Copt, label="Copt", lw=2)
    xlabel!("x")
    ylabel!("Consistency term at optimality")

    push!(plt, plt4)

    #display(assemblePlot(plotMode, plt))
    display(plot(plt..., layout=(4, 1), size=(800, 1200)))
    println("Maximal diffusion value : ", maximum(sol.Dopt))
end

function plot_solutions(solVec, plotMode::PlottingMode=DisplayMode())

    plt = []

    # Numerical Entropy Flux
    plt1 = createPlot(plotMode)
    xlabel!("x")
    ylabel!("Numerical Entropy Flux")
    title!(get_name(solVec[1].scheme)*", Nx = "*string(solVec[1].domain.Nx))
    for sol in solVec
        plot!(sol.domain.interfaces, sol.Gopt, label=sol.label, lw=2)
    end
    push!(plt, plt1)

    # Numerical Diffusion
    plt2 = createPlot(plotMode)
    xlabel!("x")
    ylabel!("Numerical Diffusion")
    for sol in solVec
        plot!(sol.domain.x, sol.Dopt, label=sol.label)
    end
    push!(plt, plt2)

    display(assemblePlot(plotMode, plt))
end

function plot_bounds(sol::OptForEntropySol; exactG_known=false)

    if exactG_known

        Gexact = exactG(sol.scheme, sol.equation, sol.u_approx[end-1])

        Gexact_test = zeros(sol.domain.Nx + 1)
        Gexact_test[2:end] = Gexact[1:end-1]
        Gexact_test[1] = Gexact[end]
    end

    plot(sol.domain.interfaces, sol.m_vec, label="m")
    plot!(sol.domain.interfaces, sol.Gopt, label="Gopt")
    if exactG_known
        plot!(sol.domain.interfaces, Gexact_test, label="Gexact")
    end
    plot!(sol.domain.interfaces, sol.M_vec, label="M")
    xlabel!("x")
    display(ylabel!("Numerical Entropy Flux"))

    # print("m_vec", "\n")
    # print(sol.m_vec, "\n")
    # if exactG_known
    #     print("Gexact", "\n")
    #     print(Gexact, "\n")
    # end
    
    # print("M_vec", "\n")
    # print(sol.M_vec, "\n")

end

function compare_exact_flux(sol::OptForEntropySol)
    Gexact = exactG(sol.scheme, sol.equation, sol.u_approx[end-1])
    Dexact = diffusion(sol.u_approx[end-1], sol.u_approx[end], Gexact, sol.domain.dx, sol.dt_vec[end], sol.equation, sol.domain)

    plot(sol.domain.interfaces, sol.Gopt, label="Gopt")
    plot!(sol.domain.interfaces, Gexact, label="Gexact")
    xlabel!("x")
    #ylabel!("Numerical Entropy Flux")
    display(ylabel!("Numerical Entropy Flux"))

    plot(sol.domain.x, sol.Dopt, label="Dopt")
    plot!(sol.domain.x, Dexact, label="Dexact")
    xlabel!("x")
    ylabel!("Numerical Diffusion")
end