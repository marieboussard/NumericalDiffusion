using Base: @view

"""
    periodic_view(u::AbstractVector, j::Int, l::Int, r::Int)

Retourne une vue sur les éléments de `u` entre `j-l` et `j+r`, avec un accès périodique.

# Arguments
- `u`: Un vecteur de taille `N`
- `j`: Indice central
- `l`, `r`: Définissent l'intervalle `[j-l, j+r]`

# Retour
Une `view` sur le vecteur `u` respectant la périodicité.
"""
function periodic_view1(u::AbstractVector, j::Int, l::Int, r::Int)
    N = length(u)
    indices = mod1.(j .- l:j + r, N)  # Indices périodiques
    return @view u[indices]  # Vue sans allocation
end

# Exemple d'utilisation
u = collect(1:10)  # Vecteur exemple
j, l, r = 3, 2, 3  # Paramètres
v = periodic_view1(u, j, l, r)
println(v)
