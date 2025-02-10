v = Matrix(transpose([1.0 2 3 4 5 6; 7 8 9 10 11 12]))
const g = 9.8
const treshold = 1e-6

function flux1(v)
    @show h = v[:,1]
    # if h > 1e-6
    #     hu = v[:,2]; return [hu, hu.^2 ./ h + g .* h.^2 / 2]
    # else
    #     return zero(v)
    # end
    @show hu = v[:,2]
    @show hu.^2 ./ h .+ g .* h.^2 / 2
    return h.>1e-6 .* hcat(hu, hu.^2 ./ h .+ g .* h.^2 / 2)
end

function flux2(v)
    h = v[1]
    h > 1e-6 ? (hu = v[2]; [hu, hu^2 / h + 9.8 * h^2 / 2]) : [0.0, 0.0]
end

function flux3(v)
    res = zero(v)
    h = v[:,1]
    for i in 1:length(h)
        if h[i] > 1e-6
            res[i,1] = v[i,2]
            res[i,2] = v[i,2]^2 / h[i] + g * h[i]^2 *0.5
        end
    end
    res
end

function flux4(v)
    res = similar(v)
    g_half = g * 0.5 
    for i in eachindex(view(v, :, 1))
        h = v[i, 1]
        if h > 1e-6
            hu = v[i, 2]
            res[i, 1] = hu
            res[i, 2] = hu^2 / h + g_half * h^2
        else
            res[i, 1] = 0
            res[i, 2] = 0
        end
    end
    res
end

ff = similar(v)

function flux5!(res, v, g)
    g_half = g * 0.5  # Pré-calcul pour optimiser
    h = view(v, :, 1)
    hu = view(v, :, 2)

    mask = h .> 1e-6  # Masque booléen pour filtrer les valeurs
    res[:, 1] .= mask .* hu  # Multiplication par 0 pour éviter l'affectation conditionnelle
    res[:, 2] .= mask .* (hu.^2 ./ h .+ g_half .* h.^2)
    return res
end

function CFL_cond(v::Matrix)
    n = size(v)[1]
    lamb = zeros(n)
    for i in 1:n
        if v[i,1] > treshold
            lamb[i] = abs.(v[i,2] / v[i,1]) + sqrt(g * abs.(v[i,1]))
        else
            lamb[i] = 0.0
        end
    end
    return maximum(lamb)
end

function CFL_cond2(v::Matrix)
    max = 0.0
    lamb = 0.0
    for i in eachindex(view(v, :, 1))
        h = v[i,1]
        hu = v[i,2]
        if h > treshold
            lamb = abs(hu / h) + sqrt(g * abs(h))
            if lamb > max
                max = lamb
            end
        end
    end
    return max
end

function CFL_cond3(v::Matrix)
    h = view(v, :, 1)  # Évite la copie des colonnes
    hu = view(v, :, 2)

    mask = h .> treshold  # Masque booléen
    lambdas = abs.(hu ./ h) .+ sqrt.(g .* abs.(h))  # Opérations vectorielles

    return maximum(lambdas .* mask)  # Assure que les valeurs en dessous du seuil restent à 0
end

function init_lake_at_rest(x, source::Source; c=1)
    v0 = zeros(size(x)..., 2)
    for I in CartesianIndices(x)
        v0[I,1] = max(0, c - zb(source, x[I]))
        v0[I,2] = 0
    end
    v0
end