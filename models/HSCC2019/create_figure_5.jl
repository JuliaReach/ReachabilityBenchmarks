# =============================================================================
# To recreate Figure 3, run the following code from the REPL:
#
# julia> include("create_figure_5.jl")
#
# =============================================================================
using Plots

# disable graphics output in GR - https://github.com/JuliaPlots/Plots.jl/issues/1182
# ENV["GKSwstype"] = "100"

include("FilteredOscillator.jl")

sol4 = filtered_oscillator(4, ApproximatingDiscretePost(), 20.)
sol_proj4 = get_projection(sol4, 7, [1, 6])
plot(sol_proj4, tickfont=font(18, "Times"), guidefontsize=18)
#savefig("Figure 5 Left.pdf")

sol16 = filtered_oscillator(16, ApproximatingDiscretePost(), 99.)
sol_proj16 = get_projection(sol16, 19, [1, 18])
plot(sol_proj16, tickfont=font(18, "Times"), guidefontsize=18)
#savefig("Figure 5 Right.pdf")
