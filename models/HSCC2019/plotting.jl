using Plots                # plotting library
using LaTeXStrings         # to write math labels using GR
using Plots.PlotMeasures   # to specify margins units such as bottom_margin=8mm

# select a plotting backend among: pyplot(), plotly(), gr()
# By default, this script requires that you have installed the GR Plots backend.
# See create_figure_5.jl for recommended setups of other plotting backends.
gr()

if Plots.backend() == Plots.GRBackend()
    # disable graphics output in GR - https://github.com/JuliaPlots/Plots.jl/issues/1182
    ENV["GKSwstype"] = "100"
end
