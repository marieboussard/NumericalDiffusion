#using Plots

struct OptForEntropySol
    domain::Domain
    equation::Equation
    method::FVMethod
    u_approx
    dt_vec
    Gopt
    Jopt::Float64
    Dopt
    m_vec
    M_vec
end

function plot_solution(sol::OptForEntropySol)

    plt = []

    plt1 = plot(sol.domain.interfaces, sol.Gopt, label="Gopt")
    xlabel!("x")
    ylabel!("Numerical Entropy Flux")

    push!(plt, plt1)

    plt2 = plot(sol.domain.x, sol.Dopt, label="Dopt")
    xlabel!("x")
    ylabel!("Numerical Diffusion")

    push!(plt, plt2)

    pltfin = plot(plt..., layout=(1, 2), size=(700, 350))

    display(title!(get_name(sol.method)))

    println("Maximal diffusion value : ", maximum(sol.Dopt))
end

function plot_bounds(sol::OptForEntropySol; exactG_known=false)

    if exactG_known

        Gexact = exactG(sol.method, sol.equation, sol.u_approx[end-1])

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
    Gexact = exactG(sol.method, sol.equation, sol.u_approx[end-1])
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