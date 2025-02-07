module FiniteVolumes

    using LinearAlgebra
    using UnPack

    include("parameters.jl")

    abstract type Cache end


    include("schemes/abstract_schemes.jl")
    include("equation.jl")

    export DefaultLogConfig
    include("log.jl")

    include("integrator.jl")
    include("problem.jl")
    include("solution.jl")

    include("CFL.jl")

    include("schemes/euler.jl")
    include("schemes/rusanov.jl")
    include("numflux.jl")
    include("cache.jl")

    include("solve.jl")

    export OneDMesh
    export Parameters

    export Equation

    export Euler
    export Rusanov

    export CFL_cond
    export CFL_cond!
    export CFL_local!

    export Integrator
    export numflux!
    export solve

end