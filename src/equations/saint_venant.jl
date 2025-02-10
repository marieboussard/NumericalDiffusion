const g = 9.8
const treshold = 1e-6 # treshold for dealing with dry water states (h=0)

# DEFINITION OF SAINT VENANT FUNCTIONS 

struct SaintVenant <: AbstractEquationFun end

function flux(::SaintVenant, v)
    res = similar(v)
    g_half = g * 0.5 
    for i in eachindex(view(v, :, 1))
        h = v[i, 1]
        if h > treshold
            hu = v[i, 2]
            res[i, 1] = hu
            res[i, 2] = hu^2 / h + g_half * h^2
        else
            res[i, 1] = 0
            res[i, 2] = 0
        end
    end
    res
end

function Dflux(::SaintVenant, u)
end

# COMPUTING CFL CONDITION FOR SAINT VENANT

function CFL_cond(::SaintVenant, v)
    max = 0.0
    lamb = 0.0
    for i in eachindex(view(v, :, 1))
        h = v[i,1]
        hu = v[i,2]
        if h > treshold
            lamb = abs(hu / h) + sqrt(g * abs(h))
            if lamb > max
                max = lamb
            end
        end
    end
    return max
end

function CFL_cond!(::SaintVenant, integrator::Integrator)

    @unpack uprev = integrator
    integrator.cfl = 0.0
    lamb = 0.0

    for i in eachindex(view(uprev, :, 1))
        h = uprev[i,1]
        hu = uprev[i,2]
        if h > treshold
            lamb = abs(hu / h) + sqrt(g * abs(h))
            if lamb > integrator.cfl
                integrator.cfl = lamb
            end
        end
    end
end

function CFL_local!(::SaintVenant, integrator::Integrator)

    @unpack uprev, cache = integrator
    @unpack stencil = cache
    cache.cfl_loc = 0.0
    lamb = 0.0

    for i in eachindex(stencil)
        h = uprev[i,1]
        hu = uprev[i,2]
        if h > treshold
            lamb = abs(hu / h) + sqrt(g * abs(h))
            if lamb > cache.cfl_loc
                cache.cfl_loc = lamb
            end
        end
    end
end

# SOME INITIAL CONDITIONS

function init_lake_at_rest(x::AbstractVector, z::Base.Callable; c=one(eltype(x)))
    znum = z.(x)
    v = zeros(eltype(x), (length(x), 2))
    for i in eachindex(x)
        v[i,1] = max(zero(eltype(x)), c - znum[i])
        v[i,2] = zero(eltype(x))
    end
    return v
end

# SOURCE TERM FOR SAINT VENANT EQUATION

struct Pointwise <: SourceDiscretize{SaintVenant} end
struct HRDisc <: SourceDiscretize{SaintVenant} end

struct TopoSource{sourceDiscretizeType<:SourceDiscretize{SaintVenant}, F1<:Base.Callable, F2<:Base.Callable} <: AbstractSource{SaintVenant}
    z::F1
    Dz::F2
    source_discretize::sourceDiscretizeType
end

z(ts::TopoSource, x::AbstractVector) = ts.z.(x)
Dz(ts::TopoSource, x::AbstractVector) = ts.Dz.(x)

struct TopoSourceCache{znumType, DznumType} <: sourceCacheType
    znum::znumType
    Dznum::DznumType
    function TopoSourceCache(topo_source::TopoSource, source_discretize::SourceDiscretize{SaintVenant}, x)
        znum = z(topo_source, x)
        Dznum = init_Dznum(source_discretize, topo_source, x)
        new{typeof(znum), typeof{Dznum}}(znum, Dznum)
    end
end

# INIT CACHE CONTENT
init_Dznum(::Pointwise, topo_source::TopoSource, x) = Dz(topo_source, x)
init_Dznum(::HRDisc, args...) = nothing

# BY DEFAULT, POINTWISE DISCRETIZATION

function discretize_sourceterm(::Pointwise, topo_source::TopoSource , v, x)
    res = zero(v)
    for i in eachindex(view(v, :, 1))
        res[i,2] = - v[:,1] .*g.* Dz(topo_source, x)
    end
    res
end
function discretize_sourceterm!(::Pointwise, integrator::Integrator)
    @unpack uprev, sourceterm, source_cache = integrator
    @views h = uprev[:,1]
    for i in eachindex(h)
       cache.sourceterm[i,2] = -h[i]*g*source_cache.Dznum[i]
    end
end


# EXAMPLE OF CONFIGURATIONS FOR SAINT VENANT EQUATION

# 1 # LAKE AT REST WITH SINUSOIDAL TOPOGRAPHY

freq = 1.0
height = 0.5
z(x) = (-cos.(2*pi*freq * x) .+ 1)*height*0.5
Dz(x) = pi*freq*(sin.(2*pi*freq * x))*height
BumpTopo = TopoSource(z, Dz, Pointwise())

SaintVenantAtRest = Equation(2, System(), SaintVenant(), x -> init_lake_at_rest(x,z), BumpTopo)