using Base: @view

"""
    periodic_view(u::AbstractVector, j::Int, l::Int, r::Int)

Retourne une vue sur les éléments de `u` entre `j-l` et `j+r`, avec un accès périodique sans allocations.

# Arguments
- `u`: Un vecteur de taille `N`
- `j`: Indice central
- `l`, `r`: Définissent l'intervalle `[j-l, j+r]`

# Retour
Un `view` sur le vecteur `u` respectant la périodicité.
"""
struct PeriodicView{T,N} <: AbstractVector{T}
    u::AbstractVector{T}
    j::Int
    l::Int
    r::Int
    N::Int

    function PeriodicView(u::AbstractVector{T}, j::Int, l::Int, r::Int) where T
        new{T,length(u)}(u, j, l, r, length(u))
    end
end

Base.size(pv::PeriodicView) = (pv.l + pv.r + 1,)
Base.getindex(pv::PeriodicView, i::Int) = pv.u[mod1(pv.j - pv.l + i - 1, pv.N)]

function periodic_view(u::AbstractVector, j::Int, l::Int, r::Int)
    return PeriodicView(u, j, l, r)
end

# Exemple d'utilisation
u = collect(1:10)  # Vecteur exemple
j, l, r = 3, 2, 3  # Paramètres
v = periodic_view(u, j, l, r)
println(v)
