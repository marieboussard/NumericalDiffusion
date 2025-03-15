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
using FiniteVolumes:init_cfl_cache
import FiniteVolumes:init_cache
import FiniteVolumes:init_sourceterm

using FiniteVolumes:get_sL
using FiniteVolumes:get_sR

abstract type QuantifMethod end
abstract type MethodCache <: FiniteVolumes.Cache end

include("bounds_computing/modified_data_type.jl")
include("entropy.jl")

include("estimator.jl")
include("diff_estimate.jl")

# BOUNDS COMPUTING
include("bounds_computing/symmetric_md.jl")

# METHODS OF QUANTIFICATION
include("methods/posteriori.jl")

# MAIN FILE
include("quantify_diffusion.jl")