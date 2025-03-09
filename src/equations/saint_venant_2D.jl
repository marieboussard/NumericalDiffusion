struct SaintVenant2D <: AbstractEquationFun end

function flux_f(::SaintVenant2D, v)
    res = similar(v)
    g_half = g * 0.5 
    # h = view(v, :, :, 1)
    # hu = view(v, :, :, 2)
    # hv = view(v, :, :, 3)
    h = selectdim(v, 3, 1)
    hu = selectdim(v, 3, 2)
    hv = selectdim(v, 3, 3)
    for I in CartesianIndices(h)
        if h[I] > treshold
            res[I, 1] = hu[I]
            res[I, 2] = hu[I]^2 / h[I] + g_half * h[I]^2
            res[I, 3] = hu[I]*hv[I] / h[I]
        else
            res[I, 1] = zero(eltype(v))
            res[I, 2] = zero(eltype(v))
            res[I, 3] = zero(eltype(v))
        end
    end
    res
end

function flux_f!(::SaintVenant2D, integrator)
    @unpack fcont = integrator.fcont
    @unpack u = integrator
    g_half = g * 0.5 
    # h = view(u, :, :, 1)
    # hu = view(u, :, :, 2)
    # hv = view(u, :, :, 3)
    h = selectdim(u, 3, 1)
    hu = selectdim(u, 3, 2)
    hv = selectdim(u, 3, 3)
    for I in CartesianIndices(h)
        if h[I] > treshold
            fcont[I, 1] = hu[I]
            fcont[I, 2] = hu[I]^2 / h[I] + g_half * h[I]^2
            fcont[I, 3] = hu[I]*hv[I] / h[I]
        else
            fcont[I, 1] = zero(eltype(u))
            fcont[I, 2] = zero(eltype(u))
            fcont[I, 3] = zero(eltype(u))
        end
    end
end

function flux_h(::SaintVenant2D, v)
    res = similar(v)
    g_half = g * 0.5 
    h = selectdim(v, 3, 1)
    hu = selectdim(v, 3, 2)
    hv = selectdim(v, 3, 3)
    # h = view(v, :, :, 1)
    # hu = view(v, :, :, 2)
    # hv = view(v, :, :, 3)
    for I in CartesianIndices(h)
        if h[I] > treshold
            # hu = v[i, 2]
            # hv = v[i, 3]
            res[I, 1] = hv[I]
            res[I, 2] = hu[I]*hv[I] / h[I]
            res[I, 3] = hv[I]^2 / h[I] + g_half * h[I]^2
        else
            res[I, 1] = zero(eltype(v))
            res[I, 2] = zero(eltype(v))
            res[I, 3] = zero(eltype(v))
        end
    end
    res
end

function flux_h!(::SaintVenant2D, integrator)
    @unpack hcont = integrator.fcont
    @unpack u = integrator
    g_half = g * 0.5 
    # h = view(u, :, :, 1)
    # hu = view(u, :, :, 2)
    # hv = view(u, :, :, 3)
    h = selectdim(u, 3, 1)
    hu = selectdim(u, 3, 2)
    hv = selectdim(u, 3, 3)
    for I in CartesianIndices(h)
        if h[I] > treshold
            hcont[I, 1] = hv[I]
            hcont[I, 2] = hu[I]*hv[I] / h[I]
            hcont[I, 3] = hv[I]^2 / h[I] + g_half * h[I]^2
        else
            hcont[I, 1] = zero(eltype(u))
            hcont[I, 2] = zero(eltype(u))
            hcont[I, 3] = zero(eltype(u))
        end
    end
end

function flux!(::SaintVenant2D, u::AbstractArray, resf::AbstractArray, resh::AbstractArray)
    g_half = g * 0.5
    h = selectdim(u, 3, 1)
    hu = selectdim(u, 3, 2)
    hv = selectdim(u, 3, 3)
    for I in CartesianIndices(h)
        if h[I] > treshold
            # f(U)
            resf[I, 1] = hu[I]
            resf[I, 2] = hu[I]^2 / h[I] + g_half * h[I]^2
            resf[I, 3] = hu[I]*hv[I] / h[I]
            # h(U)
            resh[I, 1] = hv[I]
            resh[I, 2] = hu[I]*hv[I] / h[I]
            resh[I, 3] = hv[I]^2 / h[I] + g_half * h[I]^2
        else
            fill!(resf[I], zero(eltype(u)))
        end
    end
end

function flux!(::SaintVenant2D, u::Matrix, resf::Matrix, resh::Matrix)
    g_half = g * 0.5
    h = selectdim(u, 2, 1)
    hu = selectdim(u, 2, 2)
    hv = selectdim(u, 2, 3)
    for I in CartesianIndices(h)
        if h[I] > treshold
            # f(U)
            resf[I, 1] = hu[I]
            resf[I, 2] = hu[I]^2 / h[I] + g_half * h[I]^2
            resf[I, 3] = hu[I]*hv[I] / h[I]
            # h(U)
            resh[I, 1] = hv[I]
            resh[I, 2] = hu[I]*hv[I] / h[I]
            resh[I, 3] = hv[I]^2 / h[I] + g_half * h[I]^2
        else
            fill!(resf[I], zero(eltype(u)))
        end
    end
end

function flux_f!(::SaintVenant2D, u::AbstractVector, res::AbstractVector)
    g_half = g * 0.5
    if u[1] > treshold
        # f(U)
        res[1] = u[2]
        res[2] = u[2]^2 / u[1] + g_half * u[1]^2
        res[3] = u[2]*u[3] / u[1]
    else
        fill!(res, zero(eltype(u)))
    end
end

function flux_h!(::SaintVenant2D, u::AbstractVector, res::AbstractVector)
    g_half = g * 0.5
    if u[1] > treshold
        # h(U)
        res[1] = u[3]
        res[2] = u[2]*u[3] / u[1]
        res[3] = u[3]^2 / u[1] + g_half * u[1]^2
    else
        fill!(res, zero(eltype(u)))
    end
end

flux!(eqfun::SaintVenant2D, integrator::Integrator) = flux!(eqfun, integrator.u, integrator.fcont.fcont, integrator.fcont.hcont)
# COMPUTING CFL CONDITION FOR SAINT VENANT 2D

mutable struct CFLCacheSaintVenant2D <: CFLCacheType
    cflx::Float64
    cfly::Float64
    xeigenmax::Matrix{Float64}
    yeigenmax::Matrix{Float64}
    function CFLCacheSaintVenant2D(uinit)
        xeigenmax = zero(uinit[:,:,1])
        yeigenmax = zero(uinit[:,:,1])
        h = view(uinit, :, :, 1)
        hu = view(uinit, :, :, 2)
        hv = view(uinit, :, :, 3)
    for i in eachindex(h)
        xeigenmax[i] = h[i] > treshold ? abs(hu[i] / h[i]) + sqrt(g*h[i]) : zero(eltype(xeigenmax))
        yeigenmax[i] = h[i] > treshold ? abs(hv[i] / h[i]) + sqrt(g*h[i]) : zero(eltype(yeigenmax))
    end
        new(zero(Float64), zero(Float64), xeigenmax, yeigenmax)
    end
end

init_cfl_cache(::TwoD, ::System, ::SaintVenant2D, equation, uinit) = CFLCacheSaintVenant2D(uinit)
function update_cflcache!(::TwoD, ::System, ::SaintVenant2D, integrator::Integrator)
    @unpack u, cache = integrator
    @unpack xeigenmax, yeigenmax = cache.cfl_cache
    h = view(u, :, :, 1)
    hu = view(u, :, :, 2)
    hv = view(u, :, :, 3)
    for i in eachindex(h)
        xeigenmax[i] = h[i] > treshold ? abs(hu[i] / h[i]) + sqrt(g*h[i]) : zero(eltype(xeigenmax))
        yeigenmax[i] = h[i] > treshold ? abs(hv[i] / h[i]) + sqrt(g*h[i]) : zero(eltype(yeigenmax))
    end
end

function CFL_cond2D!(::SaintVenant2D, integrator::Integrator)
    @unpack cfl_cache = integrator.cache
    cfl_cache.cflx = maximum(cfl_cache.xeigenmax)
    cfl_cache.cfly = maximum(cfl_cache.yeigenmax)
end

function CFL_local!(::TwoD, ::SaintVenant2D, integrator::Integrator, j::Int, k::Int)
    @unpack uprev, cache, space_cache = integrator
    @unpack xeigenmax, yeigenmax = cache.cfl_cache
    @unpack Nx, Ny = integrator.params.mesh
    space_cache.cflx_loc = max(xeigenmax[j], xeigenmax[mod1(j+1,Nx)])
    space_cache.cfly_loc = max(yeigenmax[k], yeigenmax[mod1(k+1,Ny)])
end

# SOME INITIAL CONDITIONS

# 1 # NO SOURCE TERM

function init_sv(x, y, args...)
    # return ones(eltype(x), (length(x), length(y), 2))
    return [1.0, 0.0, 0.0]
end

SaintVenant2Flat = Equation(TwoD(), 3, System(), SaintVenant2D(), init_sv)

# SOURCE TERM

struct HRDisc2 <: SourceDiscretize end

z(ts::TopoSource, x::AbstractVector, y::AbstractVector) = [ts.z(xj, yk) for xj in x for yk in y]
# function Dz(ts::TopoSource, x::AbstractVector, y::AbstractVector)
#     Nx, Ny = length(x), length(y)
#     res = zeros(Nx, Ny, 2)
#     for j in 1:Nx
#         for k in 1:Ny
#             dz = 
#             res[j,k] = 
#     [ts.Dz(xj, yk) for xj in x for yk in y]
# end
z(ts::TopoSource, x::Real, y::Real) = ts.z(x, y)
Dz(ts::TopoSource, x::Real, y::Real) = ts.Dz(x, y)

# z(ts::TopoSource, mesh::TwoDMesh) = ts.z(mesh.x, mesh.y)
# Dz(ts::TopoSource, mesh::TwoDMesh) = ts.Dz(mesh.x, mesh.y)

function z(ts::TopoSource, mesh::TwoDMesh)
    @unpack x, y, Nx, Ny = mesh
    res = zeros(eltype(x), Nx, Ny)
    for j in 1:Nx
        for k in 1:Ny
            res[j,k] = ts.z(x[j], y[k])
        end
    end
    res
end
function Dz(ts::TopoSource, mesh::TwoDMesh)
    @unpack x, y, Nx, Ny = mesh
    res = zeros(eltype(x), Nx, Ny, 2)
    for j in 1:Nx
        for k in 1:Ny
            dz = ts.Dz(x[j], y[k])
            res[j,k,1] = dz[1]
            res[j,k,2] = dz[2]
        end
    end
    res
end


init_Dznum(::HRDisc2, args...) = nothing

# BY DEFAULT, POINTWISE DISCRETIZATION

function discretize_sourceterm!(::TwoD, ::Pointwise, integrator::Integrator)
    @unpack u, cache, source_cache = integrator
    @views h = u[:,:,1]
    for I in CartesianIndices(h)
        cache.sourceterm[I,1] = zero(eltype(cache.sourceterm))
        cache.sourceterm[I,2] = -h[I]*g*source_cache.Dznum[I,1]
        cache.sourceterm[I,3] = -h[I]*g*source_cache.Dznum[I,2]
    end
end

# 1 # LAKE AT REST WITH FLAT TOPOGRAPHY

zflat2(x, y) = zero(x)
Dzflat2(x, y) = (zero(x), zero(y))
FlatTopo2 = TopoSource(zflat2, Dzflat2, Pointwise())

SaintVenantFlat2 = Equation(TwoD(), 3, System(), SaintVenant2D(), (x,znum) -> init_lake_at_rest(x,znum), FlatTopo2)

# 2 # LAKE AT REST WITH SINUSOIDAL TOPOGRAPHY

zsinus2(x,y) = (-cos.(2*pi*freq * x) .+ 1)*height*0.5 .+ (-cos.(2*pi*freq * y) .+ 1)*height*0.5
Dzsinus2(x,y) = (pi*freq*(sin(2*pi*freq * x))*height, pi*freq*(sin(2*pi*freq * y))*height)
BumpTopo2 = TopoSource(zsinus2, Dzsinus2, Pointwise())

SaintVenantAtRest2 = Equation(TwoD(), 3, System(), SaintVenant2D(), (x,znum) -> init_lake_at_rest(x,znum), BumpTopo2)

# A FUNCTION TO CREATE ANY TOPOGRAPHY SPECIFIED BY THE USER (THE DERIVATIVES MUST BE KNOWN)

function saintvenant_2d_with_topo(z, Dz; sourcedisc::SourceDiscretize=Pointwise(), init_fun=init_lake_at_rest)
    topo = TopoSource(z, Dz, sourcedisc)
    Equation(TwoD(), 3, System(), SaintVenant2D(), (x,znum) -> init_fun(x, znum), topo)
end