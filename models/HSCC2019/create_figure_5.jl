# =============================================================================
# To recreate Figure 3, run the following code from the REPL:
#
# julia> include("create_figure_5.jl")
#
# By default, this script requires that you have installed the GR Plots backend.
# See the bottom of this file for recommended setups for other plotting backends.
# =============================================================================
using Plots
using LaTeXStrings         # to write math labels using GR
using Plots.PlotMeasures   # to specify margins units such as bottom_margin=8mm

# select a plotting backend among: pyplot(), plotly(), gr()
gr()

if Plots.backend() == Plots.GRBackend()
    # disable graphics output in GR - https://github.com/JuliaPlots/Plots.jl/issues/1182
    ENV["GKSwstype"] = "100"
end

include("FilteredOscillator.jl")

sol4 = filtered_oscillator(4, ApproximatingDiscretePost(), 20.)
sol_proj4 = get_projection(sol4, [1, 6])

sol16 = filtered_oscillator(16, ApproximatingDiscretePost(), 99.)
sol_proj16 = get_projection(sol16, [1, 18])

plot(sol_proj4, tickfont=font(20, "Times"), guidefontsize=30, xlab=L"x", ylab=L"x_4",
                xtick=[-0.5, 0.0, 0.5], ytick=[-0.5, 0.0, 0.5],
                bottom_margin=8mm, left_margin=8mm)
savefig("Figure 5 Left.pdf")

plot(sol_proj16, tickfont=font(20, "Times"), guidefontsize=30, xlab=L"x", ylab=L"x_{16}",
                xtick=[-0.5, 0.0, 0.5], ytick=[0.0, 0.2, 0.4],
                bottom_margin=8mm, left_margin=8mm)
savefig("Figure 5 Right.pdf")

#
# Commands useful for plotly()
# plot(sol_proj4, tickfont=font(18, "Times"), guidefontsize=18, xlab="x", ylab="x₄")
# x16 = "x"*join(Char.(0x2080 .+ convert.(UInt16, [1, 6])))
# plot(sol_proj16, tickfont=font(18, "Times"), guidefontsize=18, xlab="x", ylab=x16)
# The figures are opened in a browser tab; then click on the "Save as png" button.
#

#
# Commands useful for pyplot()
# Plots.scalefontsizes() # restart to default values
# Plots.scalefontsizes(4) # set font sizes x4 bigger
# x16 = "x"*join(Char.(0x2080 .+ convert.(UInt16, [1, 6])))
# plot(sol_proj4, xlab="x", ylab="x₄")
# savefig("Figure 5 Left.pdf")
# plot(sol_proj4, xlab="x", ylab=x16)
# savefig("Figure 5 Right.pdf")
#
