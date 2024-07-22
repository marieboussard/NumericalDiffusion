abstract type Source end
struct NullSource <: Source end

abstract type ZbSource <: Source end

# Different kinds of source

# Quadratic bump
struct Bump_zb <: ZbSource
    height::Real
    width::Real
end
bump_zb(;height::Real=0.5, width::Real=1.0) = Bump_zb(height, width)
zb(bz::Bump_zb, x) = max.(0, -1/bz.width*(bz.height/0.5^2)*(x .-0.5).^2 .+ bz.height)
function Dzb(bz::Bump_zb, x)
    res = zeros(size(x))
    for i in eachindex(x)
        if -1/bz.width*(bz.height/0.5^2)*(x[i] .-0.5).^2 .+ bz.height < 0
            res[i] = 0
        else
            res[i] = -2.0/bz.width*(bz.height/0.5^2)*(x[i] .-0.5)
        end
    end
    res
end

# Sinusoidal
struct Sinus_zb <: ZbSource
    height::Real
    freq::Real
end
sinus_zb(;height::Real=0.5, freq::Real=1.0) = Sinus_zb(height, freq)
zb(sz::Sinus_zb, x) = (-cos.(2*pi*sz.freq * x) .+ 1)*sz.height/2
Dzb(sz::Sinus_zb, x) = pi*sz.freq*(sin.(2*pi*sz.freq * x))*sz.height

# Flat
struct Flat_zb <: ZbSource
    height::Real
end
flat_zb(;height::Real=0.5) = Flat_zb(height)
zb(fz::Flat_zb, x) = zero(x) .+ fz.height
Dzb(::Flat_zb, x) = zero(x)

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