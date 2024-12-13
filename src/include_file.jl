""" Including all necessary files """

using Plots
using LaTeXStrings
using LinearAlgebra
using Optim
using Optimization
using OptimizationBBO
using OptimizationMOI, Ipopt
using ForwardDiff, Zygote
using ModelingToolkit
#using Debugger
using TimerOutputs
using Test
using JuMP, OSQP

include("tools/plotting_tools.jl")
include("domain.jl")
include("equations/abstract_equations.jl")
include("schemes/abstract_schemes.jl")
include("equations/source_terms.jl")

# Equations
include("equations/eq_burgers.jl")
include("equations/eq_advection.jl")
include("equations/eq_saint_venant.jl")
include("equations/new_eq.jl")

# Finite Volumes schemes
include("schemes/finite_volume_schemes.jl")
include("schemes/numerical_entropy_flux.jl")

    # Schemes in space
    include("schemes/space_schemes/rusanov.jl")
    include("schemes/space_schemes/roe.jl")
    include("schemes/space_schemes/centered.jl")
    include("schemes/space_schemes/modified_rusanov.jl")
    include("schemes/space_schemes/hydrostatic.jl")
    # include("schemes/space_schemes/mixed_method.jl")
    include("schemes/space_schemes/MUSCL.jl")

    # Schemes in time
    include("schemes/time_schemes/Euler.jl")
    include("schemes/time_schemes/RungeKutta.jl")


# Various tools
include("tools/utils.jl")

# Finite Volumes resolution
include("tools/num_flux_computation.jl")
include("fv_solution.jl")

# # Optimization process
include("tools/modified_data_type.jl")
include("tools/bounds_type.jl")
include("tools/solution.jl")
include("opt_diffusion.jl")
# include("worst_initial_data.jl")
# include("tools/extend_initial_data.jl")
# include("quantification_a_priori.jl")
# include("diffusion_compensation/diffusion_compensation.jl")
# # include("diffusion_compensation/diffusion_compensation_v2.jl")
# include("diffusion_compensation/diffusion_compensation_v4.jl")
# include("diffusion_compensation/diffusion_compensation_v5.jl")
# include("diffusion_compensation/diffusion_compensation_v6.jl")

# # Two dimensionnal extensions
# include("multidims/multidim_priori.jl")
# include("least_squares_minimization.jl")
# include("new_inequalities.jl")