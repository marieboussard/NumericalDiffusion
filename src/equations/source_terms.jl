abstract type Source end
struct NullSource <: Source end

abstract type ZbSource <: Source end

# Different kinds of source

# Quadratic bump
struct Bump_zb <: ZbSource
    height::Real
    width::Real
end
bump_zb(;height::Real=3, width::Real=1.0) = Bump_zb(height, width)
zb(bz::Bump_zb, x) = max.(0, -1/bz.width*(bz.height/0.5^2)*(x .-0.5).^2 .+ bz.height)
Dzb(bz::Bump_zb, x) = max.(0, 2*-1/bz.width*(bz.height/0.5^2)*(x .-0.5))

# Sinusoidal
struct Sinus_zb <: ZbSource
    height::Real
    freq::Real
end

# Flat
struct Flat_zb <: ZbSource
    height::Real
end

zb(::NullSource, x) = zero(x)



# Linear topography
# zb(bz::Bump_zb, x) = x .* bz.height
# Dzb(bz::Bump_zb, x) = (ones(size(x))) .* bz.height

# # Constant topography
# zb(bz::Bump_zb, x) = zero(x) .+ bz.height
# Dzb(bz::Bump_zb, x) = zero(x)

# # Sinusoidal bump
# zb(bz::Bump_zb, x) = -0.1 * (cos.(2*pi * x) .- 1)*bz.height
# Dzb(bz::Bump_zb, x) = 0.1 * pi * sin.(2*pi * x)*bz.height

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