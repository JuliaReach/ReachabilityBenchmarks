# =============================================================================
# To recreate Figure 1, run the following code from the REPL.
#
# julia> include("create_figure_1.jl")
#
# By default, this script requires that you have installed the GR Plots backend.
# See create_figure_5.jl for recommended setups of other plotting backends.
# =============================================================================

using LazySets, LazySets.Approximations

include("plotting.jl")

b = Ball2(zeros(2), 1.)

plot(b, 1e-3, aspectratio=1, alpha=0.3,
              tickfont=font(15, "Times"),
              xtick=[-1.0, -0.5, 0.0, 0.5, 1.0], ytick=[-1.0, -0.5, 0.0, 0.5, 1.0])

plot!(Singleton([1.0, 0.0]), annotations=(1.1, 0.1, text(L"p_1")), color="green")
plot!(Singleton([0.0, 1.0]), annotations=(0.1, 1.1, text(L"p_2")), color="green")
plot!(Singleton([1.0, 1.0]), annotations=(1.09, 1.1, text(L"q")))
plot!(Singleton([0.0, 0.0]), annotations=(0.1, 0.0, text(L"0")), color="green")
plot!(annotations=(1.4, 0.1, text(L"d_1")))
plot!(annotations=(0.1, 1.4, text(L"d_2")))
plot!(annotations=(0.75, 0.8, text(L"ndir")))

plot!(x->x, x->1., -0.8, 1.3, line=1, color="black", linestyle=:dash)
plot!(x->1., x->x, -0.8, 1.3, line=1, color="black", linestyle=:dash)
plot!(x->x+1, x->0., 0.0, 0.4, line=1, color="red", linestyle=:solid, arrow=true)
plot!(x->0., x->x+1, 0.0, 0.4, line=1, color="red", linestyle=:solid, arrow=true)
plot!(x->-x, x->x+1, -1.2, .2, line=1., color="black", linestyle=:dashdot)
plot!(x->x+.6, x->x+.6, -.1, .08, line=1, color="red", linestyle=:solid, arrow=true)

savefig("Figure 1.pdf")
