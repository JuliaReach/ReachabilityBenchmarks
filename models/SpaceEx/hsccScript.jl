import Reachability

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
    PLAD01_BND30, options_PLAD01_BND30 =
        platooning(; deterministic_switching=true, time_horizon=20.,
                     allowed_distance=30.)
    PLAN01_UNB50, options_PLAN01_UNB50 =
        platooning(; deterministic_switching=true, time_horizon=40.,
                     allowed_distance=50.);

    𝑂_common = Options(:partition=>[(i:i) for i in 1:div(10, 1)])
    𝑂_dense_options_PLAD01_BND42 = merge(options_PLAD01_BND42, Options(:δ => 0.01))
    𝑂_dense_options_PLAD01_BND30 = merge(options_PLAD01_BND30, Options(:δ => 0.00001))
    𝑂_dense_options_PLAN01_UNB50 = merge(options_PLAN01_UNB50, Options(:δ => 0.03))


    # println("decomposed")
    println("PLAD01_BND42")
    @time sol1_d = run_platooning(PLAD01_BND42, options_PLAD01_BND42, opd_decomposed, BFFPS19())
    println("PLAD01_BND30")
    println("Are safety properties satisfied:", sol1_d.satisfied)
    @time sol2_d = run_platooning(PLAD01_BND30, options_PLAD01_BND30, opd_decomposed, BFFPS19())
    println("Are safety properties satisfied:", sol2_d.satisfied)
    println("PLAN01_UNB50")
    @time sol3_d = run_platooning(PLAN01_UNB50, options_PLAN01_UNB50, opd_decomposed, BFFPS19())
    println("Are safety properties satisfied:", sol3_d.satisfied)

    println("high")
    println("PLAD01_BND42")
    @time sol1_h = run_platooning(PLAD01_BND42, options_PLAD01_BND42, opd_high, BFFPSV18())
    println("Are safety properties satisfied:", sol1_h.satisfied)
    println("PLAD01_BND30")
    @time sol2_h = run_platooning(PLAD01_BND30, options_PLAD01_BND30, opd_high, BFFPSV18())
    println("Are safety properties satisfied:", sol2_h.satisfied)
    println("PLAN01_UNB50")
    @time sol3_h = run_platooning(PLAN01_UNB50, options_PLAN01_UNB50, opd_high, BFFPSV18())
    println("Are safety properties satisfied:", sol3_h.satisfied)
end

function run_lm()
    problem, options = linear_switching()
    println("linear_switchingtoons")

    println("decomposed")
    opC = BFFPS19(:δ=>0.0001, :partition=>[1:2, 3:3, 4:4, 5:5])
    opD = DecomposedDiscretePost(:out_vars=>[1,2])
    @time sol1 = solve(problem, options, opC, opD)
    println("Are safety properties satisfied:", sol1.satisfied)

    println("high")
    opC = BFFPSV18(:δ=>0.0001, :partition=>[1:2, 3:3, 4:4, 5:5])
    opD = ApproximatingDiscretePost()
    @time sol2 = solve(problem, options, opC, opD)
    println("Are safety properties satisfied:", sol2.satisfied)

end

function run_spacecraft()
    println("run_spacecraft")

    println("decomposed")
    system, options = spacecraft(; abort_time=120.)
    opC = BFFPS19(:partition =>[1:1, 2:2, 3:3, 4:4, 5:5], :δ => 0.04)
    opD = DecomposedDiscretePost(:out_vars=>[1,2,3,4])
    @time sol1 = run_spacecraft(system, options, opC, opD)
    println("Deco 120: ", sol1.satisfied)
    # println("Deco 120: ", length(sol1.Xk))

    system, options = spacecraft(; abort_time=-1.)
    opC = BFFPS19(:partition =>[1:1, 2:2, 3:3, 4:4, 5:5], :δ => 0.04)
    opD = DecomposedDiscretePost(:out_vars=>[1,2,3,4])
    @time sol2 = run_spacecraft(system, options, opC, opD)
    println("Deco -1: ", sol2.satisfied)
    # println("Deco 120: ", length(sol2.Xk))

    println("high")
    system, options = spacecraft(; abort_time=120.)
    opC = BFFPSV18(:partition =>[1:1, 2:2, 3:3, 4:4, 5:5], :δ => 0.04)
    opD = ApproximatingDiscretePost()
    @time sol3 = run_spacecraft(system, options, opC, opD)
    println("High 120: ", sol3.satisfied)
    # println("High 120: ", length(sol3.Xk))

    system, options = spacecraft(; abort_time=-1.)
    opC = BFFPSV18(:partition =>[1:1, 2:2, 3:3, 4:4, 5:5], :δ => 0.04)
    opD = ApproximatingDiscretePost()
     @time sol4 = run_spacecraft(system, options, opC, opD)
    println("High -1: ", sol4.satisfied)
    # println("Hugh 120: ", length(sol4.Xk))
end


function run_fo(dimensions, time_steps)
    println("filtered_oscillator")
    for i in dimensions
        for time_step in time_steps
            println("Time step is: ", time_step)
            problem, options, solver_options = filtered_oscillator(i, 99., true, time_step)
            @time result_d = solve(problem, options, BFFPS19(solver_options), DecomposedDiscretePost(:out_vars=>[1,2]));
            println("decomposed i: ", i, " and is satisfied:", result_d.satisfied)

            @time result_h = solve(problem, options, BFFPSV18(solver_options), LazyDiscretePost());
            println("high i: ", i, " and is satisfied:", result_h.satisfied)
            @time result_h = solve(problem, options, BFFPSV18(solver_options), ConcreteDiscretePost());
            println("concrete high i: ", i, " and is satisfied:", result_h.satisfied)
        end
    end
end

run_platoons()
run_lm()
run_spacecraft()
run_fo([16,32,64,128, 256,512,1024], [0.01])
run_fo([64],[0.01, 0.005, 0.001, 0.0005])
