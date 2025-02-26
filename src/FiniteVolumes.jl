module FiniteVolumes

    using LinearAlgebra
    using UnPack

    include("parameters.jl")

    abstract type AbstractSourceDiscretize end

    abstract type Cache end
    abstract type tcacheType <: Cache end
    abstract type scacheType <: Cache end
    abstract type sourceCacheType <: Cache end
    abstract type CFLCacheType <: Cache end


    include("schemes/abstract_schemes.jl")
    include("equations/abstract_equation.jl")
    include("equations/equation_fun.jl")
    include("equations/source.jl")
    include("equations/equation.jl")

    include("two_dim_tools.jl")

    abstract type SourceDiscretize{eqFunType<:AbstractEquationFun} end

    export DefaultLogConfig
    include("log.jl")

    include("integrator.jl")
    include("problem.jl")
    include("solution.jl")

    include("CFL.jl")
    include("CFL_2D.jl")
    
    include("equations/burgers.jl")
    include("equations/saint_venant.jl")
    include("equations/advection.jl")

    include("schemes/euler.jl")
    include("schemes/rusanov.jl")
    include("schemes/rusanov_2D.jl")
    include("schemes/roe.jl")
    include("numflux.jl")
    include("cache.jl")

    include("solve.jl")

    export OneDMesh
    export TwoDCartesian
    export Parameters

    export Equation
    export OneD
    export TwoD

    export Euler
    export Rusanov
    export Rusanov2D
    export Roe

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

    export Burgers
    export BurgersArticle

    export SaintVenant
    export SaintVenantAtRest
    export SaintVenantFlat

    export Advection
    export AdvectionExample
    export Advection2Example
    export advection2_vecfield
    export exact_advection_sol
    export u0_gauss2

    export discretize_sourceterm!
    export Pointwise

    export flux 
    export Dflux
    export has_source

    export z
    export Dz

end