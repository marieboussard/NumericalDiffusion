struct Priori{mdtype_type<:ModifiedDataType, btype_type<:BoundsType} <: QuantifMethod
    mdtype::mdtype_type
    boundstype::btype_type
    Priori(mdtype::ModifiedDataType=MeanMD(), boundstype=DefaultBounds()) = new{typeof(mdtype), typeof(boundstype)}(mdtype, boundstype)
end

get_name(method::Priori) = "Priori"*" ("*get_name(method.mdtype)*")"

struct PrioriCache{mtype<:AbstractArray} <: MethodCache
    m::mtype
    M::mtype
    Dlow::mtype
    Dup::mtype
    function PrioriCache(equation::Equation, u::AbstractArray)
        m = init_G(equation.dim, equation.eqtype, u)
        M = zero(m)
        Dlow = zero(m)
        Dup = zero(m)
        new{typeof(m)}(m, M, Dlow, Dup)
    end
end
init_cache(::Priori, args...) = PrioriCache(args...)

function perform_estimation!(method::Priori, estimator::Estimator; kwargs...)
    compute_G_bounds!(estimator)
    diffusion!(method, estimator)
end

function diffusion!(::Priori, estimator::Estimator)
    @unpack Dlow, Dup, m, M = estimator.method_cache
    @unpack etacont, etacont_init, dt, D = estimator
    @unpack Nx, dx = estimator.params.mesh
    Dtot = zero(eltype(D))
    for j in 1:Nx
        Dlow[j] = etacont[j] - etacont_init[j] + dt / dx *(m[j] - M[mod1(j-1, Nx)])
        Dup[j] = etacont[j] - etacont_init[j] + dt / dx *(M[j] - m[mod1(j-1, Nx)])
        Dtot += etacont[j] - etacont_init[j]
    end
    D .= (Dtot / sum(Dlow)) * Dlow
end