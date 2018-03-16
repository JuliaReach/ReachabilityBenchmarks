#=
Model: building.jl

See also:
bak2017cav_repeatability or
[Building example in Hylaa](https://github.com/stanleybak/hylaa/blob/master/examples/building/building.py)
=#
using Reachability, MAT, Plots

compute(o::Pair{Symbol,<:Any}...) = compute(Options(Dict{Symbol,Any}(o)))

function compute(input_options::Options)
    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "building.mat")
    A = read(file, "A")

    # initial set
    # - x1-x10 are in [0.0002, 0.00025],
    # - x25 is in [-0.0001, 0.0001],
    # - the rest is 0.0
    # 1) Hyperrectangle
    X0 = Hyperrectangle([fill(0.000225, 10); zeros(38)], [fill(0.000025, 10); zeros(14); 0.0001; zeros(23)])
    # 2) Cartesian product of intervals
    # X0 = CartesianProductArray([BallInf(fill(0.000225, 10), 0.000025), Singleton(zeros(14)), BallInf([0.], 0.0001), Singleton(zeros(23))])

    # input set
    B = read(file, "B")
    U = B * BallInf([0.9], .1)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # property: x25 < 6e-3
    p = LinearConstraintProperty(sparsevec([25], [1.0], 48), 6e-3)

    # ===============
    # Problem solving
    # ===============
    if input_options[:mode] == "reach"
        problem_options = Options(:vars => [25],
                                  :partition => [(2*i-1:2*i) for i in 1:24], # 2D blocks
                                  :plot_vars => [0, 25])
    elseif input_options[:mode] == "check"
        problem_options = Options(:vars => [25],
                                  :partition => [(2*i-1:2*i) for i in 1:24], # 2D blocks
                                  :property => p)
    end

    result = solve(S, merge(problem_options, input_options))

    # ========
    # Plotting
    # ========
    if input_options[:mode] == "reach"
        println("Plotting...")
        tic()
        plot(result)
        @eval(savefig(@relpath "building.png"))
        toc()
    end
end # function

# ===================================
# Reach tube computation, dense time
# ===================================

info("warm-up run"; prefix=" ")
compute(:δ => 1e-3, :N => 3, :mode=>"reach", :verbosity => "warn");

info("dense time, 2D blocks Hyperrectangle"; prefix="BENCHMARK SETTINGS: ")
compute(:δ => 1e-3, :T => 20.0, :mode=>"reach");

info("dense time, 2D blocks HPolygon, cf. Table 1 HSCC"; prefix="BENCHMARK SETTINGS: ")
compute(:δ => 1e-3, :T => 20.0, :mode=>"reach",
        :set_type=>HPolygon, :lazy_sih=>false, :ε=>Inf);

info("dense time, 1D blocks Interval"; prefix="BENCHMARK SETTINGS: ")
compute(:δ => 1e-3, :T => 20.0, :mode=>"reach",
        :set_type=>Interval, :partition => [[i] for i in 1:48]);
