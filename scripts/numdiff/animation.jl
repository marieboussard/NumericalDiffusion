using Plots
using Interpolations

function animation(M::AbstractMatrix{<:Real}, time::AbstractVector, space::AbstractVector, filename=""; d::Real = 10, fps::Int = 30, label = "", ylabel="")

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

    axis = scale_axis(M)

    # Total frame
    nbframe = fps * d

    # Space and time vector
    Nt = length(time)
    rescaltime = range(0,1,Nt)
    xLim = [space[begin], space[end]]

    anim = @animate for n in 1:nbframe

        t = n / (fps*d)

        # Computation of interpolations
        interpol = [linear_interpolation(rescaltime, M[i,:])(t) for i in eachindex(M[:,1])]

        # Create plots at rescaled time t
        pltsol = plot(space, interpol, xlabel="x", ylabel=ylabel, label = label, xlim = xLim, ylim = axis, legend=:bottomleft, lc=:black, title = "t="*string(round(t*time[end];digits=4),"s,   step="*string(floor(Nt*t))))



        pltsol
    end

    gif(anim, "./images/counter_example/anim_"*filename*".gif"; fps = fps)
end