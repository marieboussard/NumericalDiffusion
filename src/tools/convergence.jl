function convergence_error(scheme::FVScheme, equation::Equation, testcase::Testcase, xmin, xmax, CFL_factor, Nx_vec; show=false)
    error_vec = zeros(length(Nx_vec))
    i = 1
    for Nx in Nx_vec
        @show Nx
        domain, u0 = createOneTimestepInterval(Nx, 0.0, xmin, xmax, equation, testcase, CFL_factor)
        fv_sol = fv_solve(domain, u0, equation, scheme)
        uexact = exactData(domain, testcase)
        error_vec[i] = 1/sqrt(Nx) * sqrt(sum((fv_sol.u_approx[end] .- uexact).^2)) #norm(fv_sol.u_approx[end] .- uexact, 2)
        i+=1
        if show
            display(plot_fv_sol(fv_sol, testcase))
        end
    end
    error_vec
end

function plot_convergence_order(Nx_vec, error_vec, scheme, orderEstimate=1)
    scatter(log.(Nx_vec)/log(2), log.(error_vec)/log(2))
    plot!(log.(Nx_vec)/log(2), log(error_vec[begin])/log(2) .+ (log.(Nx_vec) .- log(Nx_vec[begin]))/log(2)*(-1*orderEstimate), label="slope "*string(orderEstimate))
    xlabel!("log(Nx)")
    ylabel!("log(norm 2 error)")
    title!(get_name(scheme))
end