# =============================================================================
# To recreate Figure 5, run the following code from the REPL:
#
# julia> include("create_figure_5.jl")
#
# By default, this script requires that you have installed the GR Plots backend.
# See the bottom of this file for recommended setups of other plotting backends.
# =============================================================================

include("plotting.jl")

include("FilteredOscillator.jl")

sol4 = filtered_oscillator(4, ApproximatingDiscretePost(), 20.)
sol_proj4 = get_projection(sol4, [1, 6])

sol16 = filtered_oscillator(16, ApproximatingDiscretePost(), 99.)
sol_proj16 = get_projection(sol16, [1, 18])

color = RGB(0.1, 0.2, 1.0)
linecolor = RGBA(1.0, 1.0, 1.0, 0.1)

plot(sol_proj4, color=color, linecolor=linecolor,
                tickfont=font(20, "Times"), guidefontsize=30,
                xlab=L"x\raisebox{-1mm}{\textcolor{white}{.}}",
                ylab=L"x_4\raisebox{2mm}{\textcolor{white}{.}}",
                xtick=[-0.5, 0.0, 0.5], ytick=[-0.5, 0.0, 0.5],
                bottom_margin=8mm, left_margin=6mm)
savefig("Figure 5 Left.pdf")

plot(sol_proj16, color=color, linecolor=linecolor,
                 tickfont=font(20, "Times"), guidefontsize=30,
                 xlab=L"x\raisebox{-1mm}{\textcolor{white}{.}}",
                 ylab=L"x_{16}\raisebox{2mm}{\textcolor{white}{.}}",
                 xtick=[-0.5, 0.0, 0.5], ytick=[0.0, 0.2, 0.4],
                 bottom_margin=8mm, left_margin=6mm)
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
