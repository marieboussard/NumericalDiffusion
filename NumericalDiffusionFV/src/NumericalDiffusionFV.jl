module NumericalDiffusionFV

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

    abstract type BoundsType end
    struct DefaultBounds <: BoundsType end
    struct MultiBounds <: BoundsType end


    include("bounds_computing/modified_data_type.jl")
    include("entropy.jl")

    include("estimator.jl")
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


    export quantify_diffusion

    # QUANTIFICATION METHODS
    export Priori
    export PrioriMultidim
    export Posteriori

end