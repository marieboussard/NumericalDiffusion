struct PrioriMultidim{mdtype_type<:ModifiedDataType, btype_type<:BoundsType} <: QuantifMethod
    mdtype::mdtype_type
    boundstype::btype_type
    PrioriMultidim(mdtype::ModifiedDataType=MeanMD(), boundstype=MultiBounds()) = new{typeof(mdtype), typeof(boundstype)}(mdtype, boundstype)
end

get_name(method::PrioriMultidim) = "Priori Multidim"*" ("*get_name(method.mdtype)*")"

# CACHE

struct PrioriMultidimCache{ltype<:AbstractArray} <: MethodCache
    l::ltype
    L::ltype
    Dlow::ltype
    function PrioriMultidimCache(equation::Equation, u::AbstractArray)
        l = init_G(equation.dim, equation.eqtype, u)
        L = zero(l)
        Dlow = zero(l)
        new{typeof(l)}(l, L, Dlow)
    end
end
init_cache(::PrioriMultidim, args...) = PrioriMultidimCache(args...)

# ESTIMATION

function perform_estimation!(method::PrioriMultidim, estimator::Estimator; kwargs...)
    compute_DG_bounds!(method.mdtype, estimator)
    diffusion!(method, estimator)
end

# BOUNDS

function compute_DG_bounds!(::SymmetricMD, estimator::Estimator)
    @unpack entfun, etacont, etacont_init = estimator
    @unpack Nx = estimator.params.mesh
    @unpack mdtype, boundstype = estimator.method
    @unpack l, L = estimator.method_cache
    @unpack sL, sR, eta_tilde, eta_hat, utilde, uhat = estimator.cache
    for j in 1:Nx
        L[j] = etacont_init[j] - etacont[j]
        l[j] = zero(eltype(l))
        utilde!(mdtype, boundstype, estimator, j)
        uhat!(mdtype, boundstype, estimator)
        @views utilde_short = selectdim(utilde, 1, sL+1:3*sL+2*sR+1)
        eta!(entfun, utilde_short, eta_tilde)
        eta!(entfun, uhat, eta_hat)
        for i in 1:sL+sR
            l[j] += eta_hat[i] - eta_tilde[i]
        end
        for i in sL+sR+2:2*(sL+sR)+1
            l[j] += eta_hat[i] - eta_tilde[i]
        end
    end
end

function diffusion!(::PrioriMultidim, estimator::Estimator)
    @unpack l, L, Dlow = estimator.method_cache
    @unpack D = estimator
    # LOWER DIFFUSION BOUND
    @. Dlow = l - L
    # DIFFUSION ESTIMATE AS THE NORMALIZATION OF LOWER BOUND
    suml, sumL = sum(l), sum(L)
    @. D = -sumL/(sumL - suml)* (L - l)
end

function diffusion!(::PrioriMultidim, l::AbstractVector, L::AbstractVector, Dlow::AbstractVector, D::AbstractVector)
    # LOWER DIFFUSION BOUND
    @. Dlow = l - L
    # DIFFUSION ESTIMATE AS THE NORMALIZATION OF LOWER BOUND
    suml, sumL = sum(l), sum(L)
    @. D = -sumL/(sumL - suml)* (L - l)
end