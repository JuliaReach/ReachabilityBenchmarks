import Reachability
using Plots, Plots.PlotMeasures, LaTeXStrings

include("../ARCH/AFF/Platooning/Platooning.jl")
include("../LinearSwitcher/ls.jl")
include("../ARCH/AFF/Rendezvous/SpacecraftRendezvous.jl")
include("../FilteredOscillator/FilteredOscillator.jl")

function run_platoons()
    println("Platoons")
    opd_decomposed = DecomposedDiscretePost(:out_vars=>[1, 4, 7])
    opd_high = ApproximatingDiscretePost()

    PLAD01_BND42, options_PLAD01_BND42 =
        platooning(; deterministic_switching=true, time_horizon=5.,
                     allowed_distance=42.)
    PLAN01_UNB50, options_PLAN01_UNB50 =
        platooning(; deterministic_switching=true, time_horizon=40.,
                     allowed_distance=50.);

    𝑂_common = Options(:partition=>[(i:i) for i in 1:div(10, 1)])
    𝑂_dense_options_PLAD01_BND42 = merge(options_PLAD01_BND42, Options(:δ => 0.01))
    𝑂_dense_options_PLAN01_UNB50 = merge(options_PLAN01_UNB50, Options(:δ => 0.03))


    t_sol1_d = @elapsed sol1_d = run_platooning(PLAD01_BND42, options_PLAD01_BND42, opd_decomposed, BFFPS19())
    t_sol2_d = @elapsed sol2_d = run_platooning(PLAN01_UNB50, options_PLAN01_UNB50, opd_decomposed, BFFPS19())

    t_sol1_h = @elapsed sol1_h = run_platooning(PLAD01_BND42, options_PLAD01_BND42, opd_high, BFFPSV18())
    t_sol2_h = @elapsed sol2_h = run_platooning(PLAN01_UNB50, options_PLAN01_UNB50, opd_high, BFFPSV18())
    return [("platoon_42", [(t_sol1_d, sol1_d.satisfied), (t_sol1_h, sol1_h.satisfied)]),
                ("platoon_inf", [(t_sol2_d, sol2_d.satisfied), (t_sol2_h, sol2_h.satisfied)])]
end

function run_lm()
    problem, options = linear_switching()

    opC = BFFPS19(:δ=>0.0001, :partition=>[1:2, 3:3, 4:4, 5:5])
    opD = DecomposedDiscretePost(:out_vars=>[1,2], :clustering=>:none)
    t_sol1 = @elapsed sol1 = solve(problem, options, opC, opD)

    opC = BFFPSV18(:δ=>0.0001, :partition=>[1:2, 3:3, 4:4, 5:5])
    opD = ApproximatingDiscretePost()
    t_sol2 = @elapsed sol2 = solve(problem, options, opC, opD)

    return [("lm", [(t_sol1, sol1.satisfied), (t_sol2, sol2.satisfied)])]
end

function run_spacecraft()

    system, options = spacecraft(; abort_time=120.)
    opC = BFFPS19(:partition =>[1:1, 2:2, 3:3, 4:4, 5:5], :δ => 0.04)
    opD = DecomposedDiscretePost(:out_vars=>[1,2,3,4], :clustering=>:none)
    t_sol1_d = @elapsed sol1_d = run_spacecraft(system, options, opC, opD)

    system, options = spacecraft(; abort_time=-1.)
    opC = BFFPS19(:partition =>[1:1, 2:2, 3:3, 4:4, 5:5], :δ => 0.04)
    opD = DecomposedDiscretePost(:out_vars=>[1,2,3,4])
    t_sol2_d = @elapsed sol2_d = run_spacecraft(system, options, opC, opD)

    system, options = spacecraft(; abort_time=120.)
    opC = BFFPSV18(:partition =>[1:1, 2:2, 3:3, 4:4, 5:5], :δ => 0.04)
    opD = ApproximatingDiscretePost()
    t_sol1_h = @elapsed sol1_h = run_spacecraft(system, options, opC, opD)

    system, options = spacecraft(; abort_time=-1.)
    opC = BFFPSV18(:partition =>[1:1, 2:2, 3:3, 4:4, 5:5], :δ => 0.04)
    opD = ApproximatingDiscretePost()
    t_sol2_h = @elapsed sol2_h = run_spacecraft(system, options, opC, opD)

    return [("SpaceCraft_120", [(t_sol1_d, sol1_d.satisfied), (t_sol1_h, sol1_h.satisfied)]),
                ("SpaceCraft_inf", [(t_sol2_d, sol2_d.satisfied), (t_sol2_h, sol2_h.satisfied)])]
end


function run_fo(dimensions, time_steps)
    results = Vector{Tuple{String, Vector{Tuple{Float64, Bool}}}}()
    sizehint!(results, length(dimensions)*length(time_steps))
    for i in dimensions
        for time_step in time_steps
            result_i = Vector{Tuple{Float64, Bool}}()
            problem, options, solver_options = filtered_oscillator(i, 99., true, time_step)

            t_result_d = @elapsed result_d = solve(problem, options, BFFPS19(solver_options), DecomposedDiscretePost(:out_vars=>[1,2]));

            t_result_hl = @elapsed result_hl = solve(problem, options, BFFPSV18(solver_options), LazyDiscretePost());

            push!(results, ("fo($time_step,$i)", [(t_result_d, result_d.satisfied), (t_result_hl, result_hl.satisfied)]))

        end
    end

    return results
end

function print_results(results)
    println("Benchmark name; Decomposed; Lazy")
    for result in results
        name, clocks = result
        d_t, d_s = clocks[1]
        h_t, h_s = clocks[2]
        println("$name, ", round(d_t, digits=2), "($d_s), ",round(h_t, digits=2), "($h_s)")
    end
end

function plot_raw()
    plot(tickfont=font(20, "Times"), guidefontsize=30, aspectratio=1, leg=false)
end

function produce_figure()
	c1 = "orange"
	c2 = "lightblue"
	lower = -0.5
	upper = 0.5
	left = -0.7
	right = 0.7

	problem, options, solver_options = filtered_oscillator(4, 99., false, 0.01)

	plotly()

    plot_raw()
    plot!(xlims=(left, right), ylims=(lower, upper),
		xticks=[-0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6], yticks=[-0.4, -0.2, 0., 0.2, 0.4],
        xlab=L"\raisebox{-1mm}{\textcolor{white}{.}}x",
        ylab=L"\raisebox{1.5mm}{\textcolor{white}{.}}y",
        left_margin=0mm, right_margin=0mm, top_margin=5mm, bottom_margin=10mm,
        size=(1000, 750))

	result_box = solve(problem, options, BFFPS19(solver_options), DecomposedDiscretePost(:out_vars=>[1,2]));
	Xkproj1 = [Approximations.project(set(result_box.Xk[i]), [1, 2], Approximations.BoxDirections) for i in eachindex(result_box.Xk)];

    fig = plot(Xkproj1[1], color=c1)
    for i in 1:length(Xkproj1)
        plot!(fig, Xkproj1[i], color=c1)
    end

	solver_options[:block_options] = Approximations.OctDirections
	solver_options[:ε_proj] = 1e-3
	solver_options[:set_type_proj] = HPolygon
	solver_options[:partition] = [(2*i-1:2*i) for i in 1:div(6, 2)]

	result_oct = solve(problem, options, BFFPS19(solver_options), DecomposedDiscretePost(:overapproximation=>Approximations.OctDirections,:out_vars=>[1,2]));
    Xkproj2 = [Approximations.project(set(result_oct.Xk[i]), [1, 2], Approximations.OctDirections) for i in eachindex(result_oct.Xk)];

    for i in 1:length(Xkproj2)
        plot!(fig, Xkproj2[i], color=c2)
    end

    savefig(fig, "hscc20")
end
produce_figure()