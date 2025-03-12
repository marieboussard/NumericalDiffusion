struct PosterioriParams <: QuantifParams end

struct PosterioriCache <: MethodCache
    m::mtype
    M::mtype

    Gopt::mtype
end