module FiniteVolumes

    using LinearAlgebra
    using UnPack

    include("tools.jl")
    include("parameters.jl")

    abstract type Cache end
    abstract type TimeCache <: Cache end
    abstract type SpaceCache <: Cache end
    abstract type SourceCache <: Cache end
    abstract type CFLCache <: Cache end
    
    abstract type AbstractSourceDiscretize end
    abstract type SourceDiscretize end

    include("schemes/abstract_schemes.jl")
    include("equations/abstract_equation.jl")
    include("equations/equation_fun.jl")
    include("equations/source.jl")
    include("equations/equation.jl")

    # abstract type SourceDiscretize{eqFunType<:AbstractEquationFun} end

    include("log.jl")

    include("integrator.jl")
    include("problem.jl")
    include("solution.jl")

    include("CFL.jl")
    include("CFL_2D.jl")
    
    # CLASSICAL EQUATIONS
    include("equations/burgers.jl")
    include("equations/saint_venant.jl")
    include("equations/saint_venant_2D.jl")
    include("equations/advection.jl")

    # NUMERICAL SCHEMES
    include("schemes/euler.jl")
    include("schemes/rusanov.jl")
    include("schemes/rusanov_2D.jl")
    include("schemes/roe.jl")
    include("schemes/HR.jl")
    include("schemes/HR_2D.jl")
    include("numflux.jl")
    include("cache.jl")

    # MAIN FUNCTIONS
    include("solve.jl")
    include("solve_functions.jl")


    # SOLVING INTERFACE
    export solve
    export hrsolve
    export hrsolve2
    export DefaultLogConfig
    export Solution

    # BUILDING MESH
    export Mesh
    export OneDMesh
    export TwoDCartesian
    export Parameters

    # BUILDING EQUATIONS
    export Equation
    export OneD
    export TwoD
    export Scalar
    export System
    export EquationFun

    # NUMERICAL SCHEMES
    export TimeScheme
    export SpaceScheme
    export Euler
    export Rusanov
    export Rusanov2D
    export Roe
    export HR
    export HR2D

    export CFL_cond
    export CFL_cond!
    export CFL_local
    export CFL_local!


    # SOME CLASSICAL EQUATIONS

    export AbstractEquationFun
    export EquationFun

    # BURGERS
    export Burgers
    export BurgersArticle
    export BurgersConstant

    # ADVECTION 
    export Advection
    export AdvectionExample
    export Advection2Example
    export advection2_vecfield
    export exact_advection_sol
    export u0_gauss2

    # SAINT-VENANT
    export SaintVenant
    export SaintVenantAtRest
    export SaintVenantFlat
    export saintvenant_with_topo

    # SAINT-VENANT 2D 
    export SaintVenant2D
    export SaintVenant2Flat
    export SaintVenantFlat2
    export SaintVenantAtRest2
    export saintvenant_2d_with_topo
    export g
    export treshold


    # SOURCE DISCRETIZATION
    export discretize_sourceterm!
    export Pointwise
    export HRDisc
    export HRDisc2
    export TopoSource
    export TopoSourceCache


    # INTERNAL COMPONENTS
    export Integrator
    export performstep!
    export loopheader!
    export loopfooter!
    export flux!
    export Dflux!
    export abs!
    export numflux!
    export update_flux!
    export update_cflcache!
    export update_log!
    export dt_CFL!
    export LogConfig
    export has_source


    # SOME FUNCTIONS USEFUL FOR PLOTTING
    export get_name

end