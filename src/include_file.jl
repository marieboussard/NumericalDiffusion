""" Including all necessary files """

using Plots
using LinearAlgebra
using Optim
using Optimization
using OptimizationBBO

include("tools/plotting_tools.jl")
include("domain.jl")
include("equations/abstract_equations.jl")
include("methods/abstract_methods.jl")
include("equations/source_terms.jl")

# Equations
include("equations/eq_burgers.jl")
include("equations/eq_saint_venant.jl")

# Finite Volumes methods
include("methods/rusanov.jl")
include("methods/roe.jl")
include("methods/hydrostatic.jl")

# Finite Volumes resolution
include("fv_solution.jl")

# Optimization process
include("tools/solution.jl")
include("opt_diffusion.jl")
include("worst_initial_data.jl")
include("tools/extend_initial_data.jl")
