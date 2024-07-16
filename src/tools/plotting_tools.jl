abstract type PlottingMode end
struct SaveMode <: PlottingMode end
struct DisplayMode <: PlottingMode end


function createPlot(::SaveMode)
    plt = plot(size=(750, 600), margin=0.5Plots.cm, legend=:bottomright,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
end

function createPlot(::DisplayMode)
    plt = plot(margin=0.5Plots.cm, legend=:bottomright)
end

function assemblePlot(::SaveMode, pltA)
    n_plots = length(pltA)
    n, m = 2, Int(n_plots/2)
    plot(pltA..., layout=(n, m), size=(750*n_plots, 600*n_plots))
end

function assemblePlot(::DisplayMode, pltA)
    n_plots = length(pltA)
    n, m = 2, Int(n_plots/2)
    plot(pltA..., layout=(n, m))
end