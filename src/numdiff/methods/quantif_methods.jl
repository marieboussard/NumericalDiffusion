abstract type QuantifMethod end
struct Priori <: QuantifMethod end
struct Posteriori <: QuantifMethod end

abstract type MethodCache end