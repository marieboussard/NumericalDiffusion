using ForwardDiff
using UnPack

# Functions to solve exactly Poisson equation -z0''(x) = f(x) depending on the nature of f 

abstract type BoundaryConditions end
struct Periodic <: BoundaryConditions end
struct Dirichlet <: BoundaryConditions end

abstract type SourceType end

mutable struct CPMSource{ftype<:AbstractVector, xtype<:AbstractVector} <: SourceType
    fvec::ftype
    xvec::xtype
    alpha::ftype
    beta::ftype
    function CPMSource(fvec::AbstractVector, xvec::AbstractVector)
        new{typeof(fvec), typeof(xvec)}(fvec, xvec, zero(fvec), zero(fvec))
    end
end

mutable struct PoissonProblem{stype<:SourceType, bcstype<:BoundaryConditions}
    source::stype
    bcs::bcstype
end

function eval_f(source::CPMSource, x::Real)
    @unpack fvec, xvec = source
    N = length(fvec)
    for k in 1:N
        if (xvec[k] <= x) && (x <= xvec[k+1])
            return fvec[k]
        end
    end
    fvec[N]
end


function z0(source::CPMSource, x::Real)
    @unpack fvec, xvec, alpha, beta = source
    N = length(alpha)
    for k in 1:N
        if (xvec[k] <= x) && (x <= xvec[k+1])
            return -fvec[k]*0.5*x^2 + alpha[k]*x + beta[k]
        end
    end
    -fvec[N]*0.5*x^2 + alpha[N]*x + beta[N]
end

z0p(source::SourceType, x::Real) = ForwardDiff.derivative(t -> z0(source, t), x)

z(source::SourceType, x::Real) = max(zero(x), z0(source, x))
zp(source::SourceType, x::Real) = ForwardDiff.derivative(t -> z(source, t), x)
zpp(source::SourceType, x::Real) = ForwardDiff.derivative(t -> zp(source, t), x)

w(source::SourceType, x::Real) = -zpp(source, x) - eval_f(source, x)

function coeffs_poisson_cpm!(::Periodic, source::CPMSource, betaone=0.0)

    @unpack fvec, xvec, alpha, beta = source
    N = length(fvec)

    beta[1] = betaone
    for k in 2:N
        beta[k] = beta[k-1] +0.5*xvec[k]^2*(fvec[k-1]-fvec[k])
    end
    alpha[1] = 0.5*(fvec[1]*xvec[1]^2-fvec[N]*xvec[N+1]^2)
    for j in 1:N-1
        alpha[1] += (fvec[j+1]-fvec[j])*xvec[j+1]*(xvec[N+1]-0.5*xvec[j+1])
    end
    alpha[1] *= 1.0/(xvec[1]-xvec[N+1])
    for j in 1:N-1
        alpha[j+1] = alpha[j] + (fvec[j+1]-fvec[j])*xvec[j+1]
    end
    # alpha[N] = fvec[N]/2 + beta[1] - beta[N]
    # for k in N-1:-1:1
    #     alpha[k] = alpha[k+1] - (fvec[k+1]-fvec[k])*xvec[k+1]
    # end


end

function coeffs_poisson_cpm!(::Dirichlet, source::CPMSource, args...)

    @unpack fvec, xvec, alpha, beta = source
    N = length(fvec)

    alpha[1] = 0.5*(fvec[1]*xvec[1]^2-fvec[N]*xvec[N+1]^2)
    for j in 1:N-1
        alpha[1] += (fvec[j+1]-fvec[j])*xvec[j+1]*(xvec[N+1]-0.5*xvec[j+1])
    end
    alpha[1] *= 1.0/(xvec[1]-xvec[N+1])
    for j in 1:N-1
        alpha[j+1] = alpha[j] + (fvec[j+1]-fvec[j])*xvec[j+1]
    end
    beta[1] = 0.5*fvec[1]*xvec[1]^2 - alpha[1]*xvec[1]
    for j in 1:N-1
        beta[j+1] = beta[j] - 0.5*xvec[j+1]^2*(fvec[j+1]-fvec[j])
    end

end


"""
Solve Poisson with a piecewise continuous source term
"""
function solve_poisson_cpm!(zvec::AbstractVector, wvec::AbstractVector, x::AbstractVector, prob::PoissonProblem; betaone=0.0)

    coeffs_poisson_cpm!(prob.bcs, prob.source, betaone)

    for i in eachindex(x)
        zvec[i] = z(prob.source, x[i])
        wvec[i] = w(prob.source, x[i])
    end

    # [z(source, xi) for xi in x], [w(source, xi) for xi in x]
    # #z.(source, x), w.(source, x)

end



# N = 5
# fvec = (2*rand(N) .- 1)*10
# #fvec = zeros(N)
# xvec = LinRange(-1, 1, N+1)
# source = CPMSource(fvec, xvec)
# problem = PoissonProblem(source, Dirichlet())

# x = LinRange(-1.0, 1.0, 1000)

# zvec, wvec = zero(x), zero(x)

# solve_poisson_cpm!(zvec, wvec, x, problem)

# z0pvec = [zp(source, xi) for xi in x]

# using CairoMakie 
# fig = Figure(size=(1000, 1500))
# ax = Axis(fig[1,1], xlabel="x", title="z(x)")
# scatter!(ax, x, zvec)
# ax2 = Axis(fig[2,1], xlabel="x", title="w(x)")
# scatter!(ax2, x, wvec)
# ax3 = Axis(fig[3,1], xlabel="x", title="z0'(x)")
# scatter!(ax3, x, z0pvec)
# fig