module FiniteVolumes


    using LinearAlgebra
    using Plots
    using UnPack

    include("parameters.jl")
    include("abstract_cache.jl")
    include("schemes/abstract_schemes.jl")
    include("equation.jl")

    include("log.jl")

    include("integrator.jl")
    include("problem.jl")
    include("solution.jl")

    include("schemes/Euler.jl")
    include("schemes/Rusanov.jl")
    include("numflux.jl")
    include("cache.jl")

    include("solve.jl")


    export solve

end