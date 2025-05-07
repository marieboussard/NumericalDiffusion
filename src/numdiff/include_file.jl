using LinearAlgebra
using UnPack
using Optim
using FiniteVolumes

using FiniteVolumes:Cache
using FiniteVolumes:SpaceScheme
using FiniteVolumes:TimeScheme
using FiniteVolumes:SpaceCache
using FiniteVolumes:TimeCache
using FiniteVolumes:CFLCache
using FiniteVolumes:SourceCache
using FiniteVolumes:EulerCache
using FiniteVolumes:RK2Cache
using FiniteVolumes:init_cfl_cache
using FiniteVolumes:init_subcache
using FiniteVolumes:update_subcache!
import FiniteVolumes:init_cache
import FiniteVolumes:init_sourceterm

using FiniteVolumes:get_sL
using FiniteVolumes:get_sR

abstract type QuantifMethod end
abstract type MethodCache <: FiniteVolumes.Cache end

abstract type BoundsType end
struct DefaultBounds <: BoundsType end
struct MultiBounds <: BoundsType end

include("bounds_computing/modified_data_type.jl")
include("entropy.jl")

include("estimator.jl")
include("subcaches.jl")
include("diff_estimate.jl")

# BOUNDS COMPUTING
include("bounds_computing/symmetric_md.jl")
include("bounds_computing/asymmetric_md.jl")
include("bounds_computing/multisymmetric_md.jl")
include("bounds_computing/multiasymmetric_md.jl")

# METHODS OF QUANTIFICATION
include("methods/posteriori.jl")
include("methods/priori.jl")
include("methods/priori_multidim.jl")

# MAIN FILE
include("quantify_diffusion.jl")

# AUXILIARY FUNCTIONS
include("G_from_theory.jl")
include("../uzawa/consistent_fluxes.jl")
include("../uzawa/fill_matrices.jl")