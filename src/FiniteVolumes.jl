module FiniteVolumes

    using LinearAlgebra
    using UnPack

    include("parameters.jl")

    abstract type Cache end
    abstract type tcacheType <: Cache end
    abstract type scacheType <: Cache end


    include("schemes/abstract_schemes.jl")
    include("equations/abstract_equation.jl")
    include("equations/equation_fun.jl")
    include("equations/equation.jl")

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

    export loopfooter!
    export loopheader!
    export performstep!
    export update_log!
    export dt_CFL!

    export Scalar
    export System
    export EquationFun

    export LogConfig

end