struct Posteriori{mdtype_type<:ModifiedDataType} <: QuantifMethod
    mdtype::mdtype_type
    Posteriori(mdtype::ModifiedDataType=MeanMD()) = new{typeof(mdtype)}(mdtype)
end


# struct Posteriori{qparamType <: PosterioriParams} <: QuantifMethod
#     qparams::qparamType
#     function Posteriori(mdtype::ModifiedDataType=SymmetricMD())
#         qparams = PosterioriParams(mdtype)
#         new{typeof(qparams)}(qparams)
#     end
# end

# # PARAMETERS

# struct PosterioriParams{mdtype_type<:ModifiedDataType} <: QuantifParams
#     mdtype::mdtype_type
#     PosterioriParams(mdtype::ModifiedDataType) = new{typeof(mdtype)}(mdtype)
# end

# CACHE

mutable struct PosterioriCache{Gtype<:AbstractArray} <: MethodCache
    Ginit::Gtype
    Gopt::Gtype
    Jopt::Float64
    function PosterioriCache(equation::Equation, u::AbstractArray)
        Ginit = init_G(equation.dim, equation.eqtype, u)
        Gopt = zero(Ginit)
        Jopt = zero(Float64)
        new{typeof(Ginit)}(Ginit, Gopt, Jopt)
    end
end

init_G(::OneD, ::Scalar, u::Vector{Float64}) = zero(u)
init_G(::OneD, ::System, u::Matrix{Float64}) = zero(selectdim(u,2,1))
init_cache(::Posteriori, args...) = PosterioriCache(args...)

# COST FUNCTION

function J(estimator::Estimator, gamma::AbstractVector)
    @unpack params, etacont_init, etacont, m, M, dt = estimator
    @unpack Nx, dx = params.mesh
    T = eltype(etacont)
    res = zero(T)
    for j in 1:Nx
        # DIFFUSION TERM
        res += (max(zero(T), etacont[j] - etacont_init[j] + dt / dx *(gamma[j] - gamma[mod1(j-1, Nx)])))^2
        # CONSISTENCY TERM
        res += (dt / dx)^2 *(max(zero(T), gamma[j] - M[j])^2 + max(zero(T), m[j] - gamma[j])^2)
    end
    res += (dt / dx)^2 *(max(zero(T), gamma[end] - M[end])^2 + max(zero(T), m[end] - gamma[end])^2)
    res
end

# NUMERICAL DIFFUSION

function diffusion!(::Posteriori, estimator::Estimator)
    @unpack D, etacont_init, etacont, params, dt, method_cache = estimator
    @unpack Nx, dx = params.mesh
    @unpack Gopt = method_cache
    for j in 1:Nx
        D[j] = etacont[j] - etacont_init[j] + dt / dx * (Gopt[j] - Gopt[mod1(j-1, Nx)])
    end
end