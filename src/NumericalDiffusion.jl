module NumericalDiffusion

    using LinearAlgebra
    using UnPack
    using Optim
    using BenchmarkTools

    #===============================================================#
    #                         FINITE VOLUMES                        #
    #===============================================================#

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
    #include("problem.jl")
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
    include("schemes/RK2.jl")
    include("schemes/rusanov.jl")
    include("schemes/rusanov_2D.jl")
    include("schemes/roe.jl")
    include("schemes/HR.jl")
    include("schemes/HR_2D.jl")
    include("schemes/MUSCL.jl")
    include("numflux.jl")
    include("cache.jl")

    # MAIN FUNCTIONS
    include("solve.jl")
    include("solve_functions.jl")


    #===============================================================#
    #                     NUMERICAL DIFFUSION                       #
    #===============================================================#

    abstract type QuantifMethod end
    abstract type MethodCache <: Cache end

    abstract type BoundsType end
    struct DefaultBounds <: BoundsType end
    struct MultiBounds <: BoundsType end


    include("numdiff/bounds_computing/modified_data_type.jl")
    include("numdiff/entropy.jl")

    include("numdiff/estimator.jl")
    include("numdiff/diff_estimate.jl")

    # BOUNDS COMPUTING
    include("numdiff/bounds_computing/symmetric_md.jl")
    include("numdiff/bounds_computing/asymmetric_md.jl")
    include("numdiff/bounds_computing/multisymmetric_md.jl")
    include("numdiff/bounds_computing/multiasymmetric_md.jl")

    # METHODS OF QUANTIFICATION
    include("numdiff/methods/posteriori.jl")
    include("numdiff/methods/priori.jl")
    include("numdiff/methods/priori_multidim.jl")

    # MAIN FILE
    include("numdiff/quantify_diffusion.jl")

    # AUXILIARY FUNCTIONS
    include("numdiff/G_from_theory.jl")

    # PROJECTION ALGORITHM 
    include("numdiff/projection_algorithm/uzawa/fill_matrices.jl")
    include("numdiff/projection_algorithm/uzawa/optimizer.jl")
    include("numdiff/projection_algorithm/uzawa/uzawa_solution.jl")
    include("numdiff/projection_algorithm/uzawa/optimize.jl")
    include("numdiff/projection_algorithm/uzawa/consistent_fluxes.jl")
    include("numdiff/projection_algorithm/newton_for_LCP.jl")
    include("numdiff/projection_algorithm/uzawa_newton.jl")


    #===============================================================#
    #                           EXPORTATIONS                        #
    #===============================================================#


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
    export RK2
    export Rusanov
    export Rusanov2D
    export Roe
    export HR
    export HR2D
    export MUSCL

    # Limiters
    export Minmod
    export Superbee

    export get_sL
    export get_sR

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



    # NUMERICAL DIFFUSION
    export quantify_diffusion


    # QUANTIFICATION METHODS
    export Priori
    export PrioriMultidim
    export Posteriori

    # MODIFIED DATA TYPE
    export AsymmetricMD
    export SymmetricMD

    # BOUNDS TYPE
    export DefaultBounds

    # ENTROPY
    export eta!
    export G!
    export BurgersEnt

    # INTERNAL COMPONENTS
    export Estimator
    export compute_G_bounds!
    export global_numflux!
    export numflux!
    export RK2Cache
    export EulerCache
    export RusanovCache
    export init_utilde
    export init_uhat
    export init_ftilde
    export init_indices
    export init_cache
    export utilde!
    export uhat!
    export J

    # PROJECTION ALGORITHM
    export BoundMode
    export SingleBound
    export DoubleBound
    export init_optim_components
    export optimize_uzawa
    export compute_entropic_G

end # module NumericalDiffusion
