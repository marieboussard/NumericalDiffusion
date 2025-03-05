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

    # abstract type SourceDiscretize{eqFunType<:AbstractEquationFun} end
    abstract type SourceDiscretize end

    export DefaultLogConfig
    include("log.jl")

    include("integrator.jl")
    include("problem.jl")
    include("solution.jl")

    include("CFL.jl")
    include("CFL_2D.jl")
    
    include("equations/burgers.jl")
    include("equations/saint_venant.jl")
    include("equations/saint_venant_2D.jl")
    include("equations/advection.jl")

    include("schemes/euler.jl")
    include("schemes/rusanov.jl")
    include("schemes/rusanov_2D.jl")
    include("schemes/roe.jl")
    include("schemes/HR.jl")
    include("numflux.jl")
    include("cache.jl")

    include("solve.jl")
    include("solve_functions.jl")

    export Mesh
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
    export HR

    export CFL_cond
    export CFL_cond!
    export CFL_local!

    export Integrator
    export numflux!
    export solve
    export hrsolve

    export loopfooter!
    export loopheader!
    export performstep!
    export update_flux!
    export update_cflcache!
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
    export SaintVenant2D
    export SaintVenant2Flat
    export SaintVenantFlat2
    export SaintVenantAtRest2
    export saintvenant_with_topo
    export saintvenant_2d_with_topo

    export Advection
    export AdvectionExample
    export Advection2Example
    export advection2_vecfield
    export exact_advection_sol
    export u0_gauss2

    export discretize_sourceterm!
    export Pointwise
    export HRDisc
    export TopoSource
    export TopoSourceCache

    export flux 
    export Dflux
    export flux_f
    export flux_h
    export has_source

    # export z
    # export Dz

    export g
    export treshold

end