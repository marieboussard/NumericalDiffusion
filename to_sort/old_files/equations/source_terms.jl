abstract type Source end
struct NullSource <: Source end

# Source for Saint-Venant equation

abstract type Topography <: Source end

# Different kinds of source

# Quadratic bump
struct BumpTopo{T <: Real} <: Topography
    height::T
    width::T
end
BumpTopo(;height::Real=0.5, width::Real=1.0) = BumpTopo(height, width)
zb(bt::BumpTopo, x) = max.(0, -1/bt.width*(bt.height/0.5^2)*(x .-0.5).^2 .+ bt.height)
function Dzb(bt::BumpTopo, x)
    res = zeros(size(x))
    for i in eachindex(x)
        if -1/bt.width*(bt.height/0.5^2)*(x[i] .-0.5).^2 .+ bt.height < 0
            res[i] = 0
        else
            res[i] = -2.0/bt.width*(bt.height/0.5^2)*(x[i] .-0.5)
        end
    end
    res
end

# Sinusoidal
struct SinusTopo{T <: Real} <: Topography
    height::T
    freq::T
end
SinusTopo(;height::Real=0.5, freq::Real=1.0) = SinusTopo(height, freq)
zb(st::SinusTopo, x) = (-cos.(2*pi*st.freq * x) .+ 1)*st.height/2
Dzb(st::SinusTopo, x) = pi*st.freq*(sin.(2*pi*st.freq * x))*st.height

# Flat
struct FlatTopo{T <: Real} <: Topography
    height::T
end
FlatTopo(;height::Real=0.5) = FlatTopo(height)
zb(ft::FlatTopo, x) = zero(x) .+ ft.height
Dzb(::FlatTopo, x) = zero(x)

zb(::NullSource, x) = zero(x)

# A counter example for the discrete entropy inequality ?
struct DiscontinuousTopo <: Topography end
# zb(::DiscontinuousTopo, x) = (0.5*x) .+ (x.>0.5)
Dzb(::DiscontinuousTopo, x) = 0.5 .+ zero(x)
zb(::DiscontinuousTopo, x) = 2.0 .+ (10*x.^2) .*(x.<0.5).*(x.>0.48) .- (10*x.^2)  .*(x.<0.52).*(x.>0.5)


zb_tilde(::BumpTopo, ut, c) = c - ut[1]

#sourceTerm(::FVMethod, zbSource::Topography, domain::Domain, v) = [[0.0, -v[i][1] * g * Dzb(zbSource, domain.x)[i]] for i in eachindex(v)]

# function sourceTerm(::FVMethod, zbSource::Topography, domain::Domain, v)
#     S = zero(v)
#     for i in 1:length(v[:,1])
#         S[i,1] = 0.0
#         S[i,2] = -v[i,1]*g*domain.DSourceVec[i]
#     end
#     S
# end

function source_term(::Equation, ::Method, domain::Domain, v; kwargs...)
    return zero(v)
end


function add_source!(zbSource::Topography, domain::Domain)
    domain.sourceVec = zb(zbSource, domain.x)
    domain.DSourceVec = Dzb(zbSource, domain.x)
end

function add_source!(::NullSource, domain::Domain) end

function give_source(mSource::AbstractArray{T}) where T
    mSource
end

function manage_source(domain::Domain{T}) where T
    # @show typeof(zeros(T, (domain.Nx,1)))
    # @show typeof(domain.sourceVec)
    isnothing(domain.sourceVec) ? zeros(T, (domain.Nx,1)) : give_source(domain.sourceVec)
end