using Plots
using Interpolations

function animation(Mvec::AbstractVector, time::AbstractVector, space::AbstractVector; filename="", d::Real = 10, fps::Int = 30, labels = ["" for i in eachindex(Mvec)], ylabel="", colors=[:black for i in eachindex(Mvec)])#{<AbstractMatrix{<:Real}}

    # INPUTS
    #   - M         : Matrix to plot
    #   - time      : Time
    #   - space     : Mesh
    #   - d         : Duration wanted for the animation,
    #   - fps       : Number of frame per secondes,

    # set the min and maximum values of the plot, to scale the axis
    function scale_axis(V; _min = min(V...), _max = max(V...))
        minV = min(min(V...), _min)
        maxV = max(max(V...), _max)
        minV *= 1.05^sign(-minV)
        maxV *= 1.05^sign(maxV)
        [minV, maxV]
    end

    axis = scale_axis(Mvec[end])

    # Total frame
    nbframe = fps * d

    # Space and time vector
    Nt = length(time)
    rescaltime = range(0,1,Nt)
    xLim = [minimum([space[i][begin] for i in eachindex(space)]), minimum([space[i][end] for i in eachindex(space)])]

    anim = @animate for n in 1:nbframe
        

        t = n / (fps*d)

        #pltsol = plot()
        pltsol = plot(size=(900, 600), margin=1Plots.cm, legend=:topleft,
legendfontsize=30,
titlefontsize=42,
guidefontsize=42,
tickfontsize=36)

        for m in 1:length(Mvec)

            M = Mvec[m]

            # Computation of interpolations
            interpol = [linear_interpolation(rescaltime, M[i,:])(t) for i in eachindex(M[:,1])]

            # Create plots at rescaled time t
            plot!(space[m], interpol, xlabel="x", ylabel=ylabel, label = labels[m], xlim = xLim, ylim = axis, legend=:bottomleft, lc=colors[m], title = "t="*string(round(t*time[end];digits=4),"s,   step="*string(floor(Nt*t))))
        end



        
    end

    gif(anim, "./images/GTT/anim_"*filename*".gif"; fps = fps)

end