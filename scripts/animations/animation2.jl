# using CairoMakie
# using Interpolations

# """
#     animation_makie(Mvec::AbstractVector, time::AbstractVector, space::AbstractVector;
#                     filename="", d::Real=10, fps::Int=30, labels=[...],
#                     ylabel="", colors=[:black for i in eachindex(Mvec)],
#                     stack::Bool=false)

# Crée une animation avec CairoMakie et l'enregistre en `.mp4`.

# - `Mvec` : vecteur de matrices (chacune de taille NxT, espace × temps).
# - `time` : vecteur des temps.
# - `space` : vecteur des maillages (un pour chaque matrice de `Mvec`).
# - `stack` : si `true`, deux figures l’une au-dessus de l’autre.
# """
# function animation_makie(Mvec::AbstractVector, time::AbstractVector, space::AbstractVector;
#                          filename="", d::Real=10, fps::Int=30,
#                          labels=[string(i) for i in eachindex(Mvec)],
#                          ylabel="", colors=[:black for i in eachindex(Mvec)],
#                          stack::Bool=false)

#     # total frames
#     nbframe = fps * d
#     Nt = length(time)
#     rescaltime = range(0, 1, Nt)

#     # échelle verticale (min/max global)
#     allvals = reduce(vcat, [vec(M) for M in Mvec])
#     minV, maxV = minimum(allvals), maximum(allvals)
#     minV *= 1.05^sign(-minV)
#     maxV *= 1.05^sign(maxV)

#     # figure
#     fig = Figure(resolution = (900, stack ? 1000 : 600))

#     # un ou deux axes
#     ax1 = Axis(fig[1, 1], xlabel="x", ylabel=ylabel, title="Animation 1")
#     ax2 = stack ? Axis(fig[2, 1], xlabel="x", ylabel=ylabel, title="Animation 2") : nothing

#     # lignes pour chaque matrice
#     linesets = []
#     for (m, M) in enumerate(Mvec)
#         x = space[m]
#         lineobj = lines!(ax1, x, M[:,1], color=colors[m], label=labels[m])
#         push!(linesets, lineobj)
#     end

#     if stack
#         # doublon pour deuxième axe (même données par défaut)
#         for (m, M) in enumerate(Mvec)
#             x = space[m]
#             lineobj = lines!(ax2, x, M[:,1], color=colors[m], label=labels[m])
#             push!(linesets, lineobj)
#         end
#     end

#     # limites
#     ylims!(ax1, minV, maxV)
#     if stack
#         ylims!(ax2, minV, maxV)
#     end

#     # fonction mise à jour
#     function updateframe(n)
#         t = n / (fps*d)
#         for (m, M) in enumerate(Mvec)
#             interpol = [linear_interpolation(rescaltime, M[i, :])(t) for i in 1:size(M,1)]
#             linesets[m][1] = space[m]
#             linesets[m][2] = interpol
#         end
#         if stack
#             offset = length(Mvec)
#             for (m, M) in enumerate(Mvec)
#                 interpol = [linear_interpolation(rescaltime, M[i, :])(t) for i in 1:size(M,1)]
#                 linesets[offset+m][1] = space[m]
#                 linesets[offset+m][2] = interpol
#             end
#         end
#     end

#     # enregistrement vidéo
#     out = isempty(filename) ? "animation.mp4" : filename*".mp4"
#     record(fig, out, 1:nbframe; framerate=fps) do n
#         updateframe(n)
#     end

#     return out
# end


# using CairoMakie
# using Interpolations

# function animation_makie(Mvec::AbstractVector, time::AbstractVector, space::AbstractVector;
#                          filename="", d::Real=10, fps::Int=30,
#                          labels=[string(i) for i in eachindex(Mvec)],
#                          ylabel="", colors=[:black for i in eachindex(Mvec)],
#                          stack::Bool=false)

#     # total frames
#     nbframe = fps * d
#     Nt = length(time)
#     rescaltime = range(0, 1, Nt)

#     # échelle verticale (min/max global)
#     allvals = reduce(vcat, [vec(M) for M in Mvec])
#     minV, maxV = minimum(allvals), maximum(allvals)
#     minV *= 1.05^sign(-minV)
#     maxV *= 1.05^sign(maxV)

#     # figure
#     fig = Figure(resolution = (900, stack ? 1000 : 600))

#     # un ou deux axes
#     ax1 = Axis(fig[1, 1], xlabel="x", ylabel=ylabel)
#     ax2 = stack ? Axis(fig[2, 1], xlabel="x", ylabel=ylabel) : nothing

#     # label pour titre dynamique
#     title_label1 = Label(fig, "t = 0.00 s", fontsize=24)
#     fig[1, 1] = hbox(title_label1, ax1)

#     title_label2 = stack ? Label(fig, "t = 0.00 s", fontsize=24) : nothing
#     if stack
#         fig[2, 1] = hbox(title_label2, ax2)
#     end

#     # lignes pour chaque matrice
#     linesets = []
#     for (m, M) in enumerate(Mvec)
#         x = space[m]
#         lineobj = lines!(ax1, x, M[:,1], color=colors[m], label=labels[m])
#         push!(linesets, lineobj)
#     end

#     if stack
#         # doublon pour deuxième axe
#         for (m, M) in enumerate(Mvec)
#             x = space[m]
#             lineobj = lines!(ax2, x, M[:,1], color=colors[m], label=labels[m])
#             push!(linesets, lineobj)
#         end
#     end

#     # limites
#     ylims!(ax1, minV, maxV)
#     if stack
#         ylims!(ax2, minV, maxV)
#     end

#     # fonction mise à jour
#     function updateframe(n)
#         t = n / (fps*d) * maximum(time)  # temps réel
#         # mettre à jour les lignes
#         for (m, M) in enumerate(Mvec)
#             interpol = [linear_interpolation(rescaltime, M[i, :])(t/maximum(time)) for i in 1:size(M,1)]
#             linesets[m][1] = space[m]
#             linesets[m][2] = interpol
#         end
#         if stack
#             offset = length(Mvec)
#             for (m, M) in enumerate(Mvec)
#                 interpol = [linear_interpolation(rescaltime, M[i, :])(t/maximum(time)) for i in 1:size(M,1)]
#                 linesets[offset+m][1] = space[m]
#                 linesets[offset+m][2] = interpol
#             end
#         end

#         # mettre à jour le titre
#         title_label1.text = "t = $(round(t, digits=2)) s"
#         if stack
#             title_label2.text = "t = $(round(t, digits=2)) s"
#         end
#     end

#     # enregistrement vidéo
#     out = isempty(filename) ? "animation.mp4" : filename*".mp4"
#     record(fig, out, 1:nbframe; framerate=fps) do n
#         updateframe(n)
#     end

#     return out
# end


# using CairoMakie
# using Interpolations

# function animation_makie(Mvec::AbstractVector, time::AbstractVector, space::AbstractVector;
#                          filename="", d::Real=10, fps::Int=30,
#                          labels=[string(i) for i in eachindex(Mvec)],
#                          ylabel="", colors=[:black for i in eachindex(Mvec)],
#                          stack::Bool=false)

#     nbframe = fps * d
#     Nt = length(time)
#     rescaltime = range(0, 1, Nt)

#     # échelle verticale (min/max global)
#     allvals = reduce(vcat, [vec(M) for M in Mvec])
#     minV, maxV = minimum(allvals), maximum(allvals)
#     minV *= 1.05^sign(-minV)
#     maxV *= 1.05^sign(maxV)

#     # figure
#     fig = Figure(resolution = (900, stack ? 1000 : 600))

#     # un ou deux axes
#     ax1 = Axis(fig[1, 1], xlabel="x", ylabel=ylabel, title="t = 0.00 s")
#     ax2 = stack ? Axis(fig[2, 1], xlabel="x", ylabel=ylabel, title="t = 0.00 s") : nothing

#     # lignes pour chaque matrice
#     linesets = []
#     for (m, M) in enumerate(Mvec)
#         x = space[m]
#         lineobj = lines!(ax1, x, M[:,1], color=colors[m], label=labels[m])
#         push!(linesets, lineobj)
#     end

#     if stack
#         # doublon pour deuxième axe
#         for (m, M) in enumerate(Mvec)
#             x = space[m]
#             lineobj = lines!(ax2, x, M[:,1], color=colors[m], label=labels[m])
#             push!(linesets, lineobj)
#         end
#     end

#     # limites
#     ylims!(ax1, minV, maxV)
#     if stack
#         ylims!(ax2, minV, maxV)
#     end

#     # fonction mise à jour
#     function updateframe(n)
#         t = n / (fps*d) * maximum(time)  # temps réel
#         # mettre à jour les lignes
#         for (m, M) in enumerate(Mvec)
#             interpol = [linear_interpolation(rescaltime, M[i, :])(t/maximum(time)) for i in 1:size(M,1)]
#             linesets[m][1] = space[m]
#             linesets[m][2] = interpol
#         end
#         if stack
#             offset = length(Mvec)
#             for (m, M) in enumerate(Mvec)
#                 interpol = [linear_interpolation(rescaltime, M[i, :])(t/maximum(time)) for i in 1:size(M,1)]
#                 linesets[offset+m][1] = space[m]
#                 linesets[offset+m][2] = interpol
#             end
#         end

#         # mettre à jour le titre des axes
#         ax1.title = "t = $(round(t, digits=2)) s"
#         if stack
#             ax2.title = "t = $(round(t, digits=2)) s"
#         end
#     end

#     # enregistrement vidéo
#     out = isempty(filename) ? "animation.mp4" : filename*".mp4"
#     record(fig, out, 1:nbframe; framerate=fps) do n
#         updateframe(n)
#     end

#     return out
# end


# using CairoMakie
# using Interpolations

# function animation_makie(Mvec::AbstractVector, time::AbstractVector, space::AbstractVector;
#                          filename="", d::Real=10, fps::Int=30,
#                          labels=[string(i) for i in eachindex(Mvec)],
#                          ylabel="", colors=[:black for i in eachindex(Mvec)],
#                          stack::Bool=false,
#                          linewidth::Real=3.0,      # épaisseur des courbes
#                          fontsize_title::Int=30,  # taille du titre
#                          fontsize_label::Int=24,  # taille des labels axes
#                          fontsize_legend::Int=20) # taille de la légende

#     nbframe = fps * d
#     Nt = length(time)
#     rescaltime = range(0, 1, Nt)

#     # échelle verticale (min/max global)
#     allvals = reduce(vcat, [vec(M) for M in Mvec])
#     minV, maxV = minimum(allvals), maximum(allvals)
#     minV *= 1.05^sign(-minV)
#     maxV *= 1.05^sign(maxV)

#     # figure
#     fig = Figure(resolution = (900, stack ? 1000 : 600))

#     # un ou deux axes
#     ax1 = Axis(fig[1, 1], xlabel="x", ylabel=ylabel,
#                title="t = 0.00 s", 
#                xlabelsize=fontsize_label, ylabelsize=fontsize_label,
#                titlesize=fontsize_title)
#     ax2 = stack ? Axis(fig[2, 1], xlabel="x", ylabel=ylabel,
#                        title="t = 0.00 s",
#                        xlabelsize=fontsize_label, ylabelsize=fontsize_label,
#                        titlesize=fontsize_title) : nothing

#     # lignes pour chaque matrice
#     linesets = []
#     for (m, M) in enumerate(Mvec)
#         x = space[m]
#         lineobj = lines!(ax1, x, M[:,1], color=colors[m], label=labels[m], linewidth=linewidth)
#         push!(linesets, lineobj)
#     end

#     if stack
#         # doublon pour deuxième axe
#         for (m, M) in enumerate(Mvec)
#             x = space[m]
#             lineobj = lines!(ax2, x, M[:,1], color=colors[m], label=labels[m], linewidth=linewidth)
#             push!(linesets, lineobj)
#         end
#     end

#     # limites
#     ylims!(ax1, minV, maxV)
#     if stack
#         ylims!(ax2, minV, maxV)
#     end

#     # activer la légende avec taille de texte
#     axislegend(ax1, fontsize=fontsize_legend)
#     if stack
#         axislegend(ax2, fontsize=fontsize_legend)
#     end

#     # fonction mise à jour
#     function updateframe(n)
#         t = n / (fps*d) * maximum(time)  # temps réel
#         # mettre à jour les lignes
#         for (m, M) in enumerate(Mvec)
#             interpol = [linear_interpolation(rescaltime, M[i, :])(t/maximum(time)) for i in 1:size(M,1)]
#             linesets[m][1] = space[m]
#             linesets[m][2] = interpol
#         end
#         if stack
#             offset = length(Mvec)
#             for (m, M) in enumerate(Mvec)
#                 interpol = [linear_interpolation(rescaltime, M[i, :])(t/maximum(time)) for i in 1:size(M,1)]
#                 linesets[offset+m][1] = space[m]
#                 linesets[offset+m][2] = interpol
#             end
#         end

#         # mettre à jour le titre des axes
#         ax1.title = "t = $(round(t, digits=2)) s"
#         if stack
#             ax2.title = "t = $(round(t, digits=2)) s"
#         end
#     end

#     # enregistrement vidéo
#     out = isempty(filename) ? "animation.mp4" : filename*".mp4"
#     record(fig, out, 1:nbframe; framerate=fps) do n
#         updateframe(n)
#     end

#     return out
# end


using CairoMakie
using Interpolations

function animation_makie(Mvec::AbstractVector, time::AbstractVector, space::AbstractVector;
                         filename="", d::Real=10, fps::Int=30,
                         labels=[string(i) for i in eachindex(Mvec)],
                         ylabel="", colors=[:black for i in eachindex(Mvec)],
                         stack::Bool=false,
                         linewidth::Real=3.0,
                         fontsize_title::Int=30,
                         fontsize_label::Int=24,
                         fontsize_legend::Int=20,
                         fontsize_ticks::Int=20)

    nbframe = fps * d
    Nt = length(time)
    rescaltime = range(0, 1, Nt)

    # échelle verticale (min/max global)
    allvals = reduce(vcat, [vec(M) for M in Mvec])
    minV, maxV = minimum(allvals), maximum(allvals)
    minV *= 1.05^sign(-minV)
    maxV *= 1.05^sign(maxV)

    # figure
    fig = Figure(resolution = (900, stack ? 1000 : 600))

    # un ou deux axes
    ax1 = Axis(fig[1, 1], xlabel="x", ylabel=ylabel,
               title="t = 0.00 s", 
               xlabelsize=fontsize_label, ylabelsize=fontsize_label,
               titlesize=fontsize_title,
               xticklabelsize=fontsize_ticks, yticklabelsize=fontsize_ticks)
    ax2 = stack ? Axis(fig[2, 1], xlabel="x", ylabel=ylabel,
                       title="t = 0.00 s",
                       xlabelsize=fontsize_label, ylabelsize=fontsize_label,
                       titlesize=fontsize_title,
                       xticklabelsize=fontsize_ticks, yticklabelsize=fontsize_ticks) : nothing

    # lignes pour chaque matrice
    linesets = []
    for (m, M) in enumerate(Mvec)
        x = space[m]
        lineobj = lines!(ax1, x, M[:,1], color=colors[m], label=labels[m], linewidth=linewidth)
        push!(linesets, lineobj)
    end

    if stack
        # doublon pour deuxième axe
        for (m, M) in enumerate(Mvec)
            x = space[m]
            lineobj = lines!(ax2, x, M[:,1], color=colors[m], label=labels[m], linewidth=linewidth)
            push!(linesets, lineobj)
        end
    end

    # limites
    ylims!(ax1, minV, maxV)
    if stack
        ylims!(ax2, minV, maxV)
    end

    # activer la légende
    axislegend(ax1, fontsize=fontsize_legend)
    if stack
        axislegend(ax2, fontsize=fontsize_legend)
    end

    # fonction mise à jour
    function updateframe(n)
        t = n / (fps*d) * maximum(time)
        for (m, M) in enumerate(Mvec)
            interpol = [linear_interpolation(rescaltime, M[i, :])(t/maximum(time)) for i in 1:size(M,1)]
            linesets[m][1] = space[m]
            linesets[m][2] = interpol
        end
        if stack
            offset = length(Mvec)
            for (m, M) in enumerate(Mvec)
                interpol = [linear_interpolation(rescaltime, M[i, :])(t/maximum(time)) for i in 1:size(M,1)]
                linesets[offset+m][1] = space[m]
                linesets[offset+m][2] = interpol
            end
        end

        # titre dynamique
        ax1.title = "t = $(round(t, digits=2)) s"
        if stack
            ax2.title = "t = $(round(t, digits=2)) s"
        end
    end

    # enregistrement vidéo
    out = isempty(filename) ? "animation.mp4" : filename*".mp4"
    record(fig, out, 1:nbframe; framerate=fps) do n
        updateframe(n)
    end

    return out
end
