abstract type Source end
struct NullSource <: Source end

abstract type ZbSource <: Source end
struct Bump_zb <: ZbSource
    height::Float64
end

bump_zb(height=3) = Bump_zb(height)

zb(::NullSource, x) = zero(x)

# zb(bz::Bump_zb, x) = (-0.5 .* x .* (-1 .+ x)) .* bz.height
# Dzb(bz::Bump_zb, x) = (-x .+ 0.5) .* bz.height

# zb(bz::Bump_zb, x) = x .* bz.height
# Dzb(bz::Bump_zb, x) = (ones(size(x))) .* bz.height

zb(bz::Bump_zb, x) = -0.1 * (cos.(2*pi * x) .- 1)*bz.height
Dzb(bz::Bump_zb, x) = 0.1 * pi * sin.(2*pi * x)*bz.height

# zb(::Bump_zb, x) = zero(x)
# Dzb(::Bump_zb, x) = zero(x)

zb_tilde(::Bump_zb, ut, c) = c - ut[1]

#sourceTerm(::FVMethod, zbSource::ZbSource, domain::Domain, v) = [[0.0, -v[i][1] * g * Dzb(zbSource, domain.x)[i]] for i in eachindex(v)]

# function sourceTerm(::FVMethod, zbSource::ZbSource, domain::Domain, v)
#     S = zero(v)
#     for i in 1:length(v[:,1])
#         S[i,1] = 0.0
#         S[i,2] = -v[i,1]*g*domain.DSourceVec[i]
#     end
#     S
# end

function sourceTerm(::FVMethod, domain::Domain, v; z=domain.sourceVec, Dz=domain.DSourceVec)
    S = zero(v)
    for i in 1:length(v[:,1])
        S[i,1] = 0.0
        S[i,2] = -v[i,1]*g*domain.DSourceVec[i]
    end
    S
end


function addSource!(zbSource::ZbSource, domain::Domain)
    domain.sourceVec = zb(zbSource, domain.x)
    domain.DSourceVec = Dzb(zbSource, domain.x)
end

function addSource!(::NullSource, domain::Domain) end