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

include("tools/plotting_tools.jl")
include("domain.jl")
include("equations/abstract_equations.jl")
include("methods/abstract_methods.jl")
include("equations/source_terms.jl")

# Equations
include("equations/eq_burgers.jl")
include("equations/eq_advection.jl")
include("equations/eq_saint_venant.jl")
include("equations/new_eq.jl")

# Finite Volumes methods
include("methods/rusanov.jl")
include("methods/centered.jl")
include("methods/modified_rusanov.jl")
include("methods/roe.jl")
include("methods/hydrostatic.jl")
include("methods/mixed_method.jl")

# Finite Volumes resolution
include("tools/num_flux_computation.jl")
include("fv_solution.jl")

# Optimization process
include("tools/modified_data_type.jl")
include("tools/bounds_type.jl")
include("tools/solution.jl")
include("opt_diffusion.jl")
include("worst_initial_data.jl")
include("tools/extend_initial_data.jl")
include("quantification_a_priori.jl")
include("diffusion_compensation/diffusion_compensation.jl")
# include("diffusion_compensation/diffusion_compensation_v2.jl")
include("diffusion_compensation/diffusion_compensation_v4.jl")
include("diffusion_compensation/diffusion_compensation_v5.jl")