struct Y
    a::Int
    b::Float64
end

struct X
    y::Y
end

function Base.getproperty(x::X, name::Symbol)
    if name in fieldnames(Y)  # Vérifie si le champ appartient à Inner
        #return getproperty(getfield(x, :y), name)
        return getfield(getfield(x, :y), name)
    else
        return getfield(x, name)  # Accès normal aux autres champs de Outer
    end
end

# Test
y = Y(42, 3.14)
x = X(y)

println(x.a)  # 42 (au lieu de outer.inner.a)
println(x.b)  # 3.14 (au lieu de outer.inner.b)


abstract type A end
abstract type B end

struct A1 <: A end
struct A2 <: A end

struct B1 <: B end
struct B2 <: B end

# Fonction qui associe A à B
type_map(::Type{A1}) = B1
type_map(::Type{A2}) = B2

# Utilisation
println(type_map(A1))  # B1
println(type_map(A2))  # B2

