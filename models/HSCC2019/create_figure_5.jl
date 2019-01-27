# =============================================================================
# To recreate Figure 3, run the following code from the REPL:
#
# julia> include("create_figure_5.jl")
#
# =============================================================================
using Plots, LaTeXStrings

# set the plotting backend among: pyplot(), plotly(), gr()
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

plot(sol_proj4, tickfont=font(18, "Times"), guidefontsize=18, xlab=L"x_1", ylab=L"x_2")
savefig("Figure 5 Left.pdf")

plot(sol_proj16, tickfont=font(18, "Times"), guidefontsize=18, xlab=L"x_1", ylab=L"x_2")
savefig("Figure 5 Right.pdf")

#
# Commands useful for plotly()
# plot(sol_proj4, tickfont=font(18, "Times"), guidefontsize=18)
# plot(sol_proj16, tickfont=font(18, "Times"), guidefontsize=18)
# The figures are opened in a browser tab; then click on the "Save as png" button.
#

#
# Commands useful for pyplot()
# Plots.scalefontsizes() # restart to default values
# Plots.scalefontsizes(4) # set font sizes x4 bigger
# plot(sol_proj4)
# savefig("Figure 5 Left.pdf")
# plot(sol_proj16)
# savefig("Figure 5 Right.pdf")
#
