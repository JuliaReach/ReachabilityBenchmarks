#=
Model: iss-x182.jl

This is a 270-variable model of component 1r (Russian service module) of the
International Space Station (ISS).

The corresponding SpaceEx model and configuration file are iss.xml and iss.cfg.

NOTE:
This is a two-dimensional model where we have extracted the dynamics of a 2d
block of interest.
=#
include("../../src/LazySets.jl")
include("../../src/Systems.jl")
include("../../src/Approximations.jl")
include("../../src/Reachability.jl")
include("../../src/Transformations.jl")

using LazySets, Systems, Reachability, Approximations, Transformations, MAT, Expokit
using PyPlot
#using Gadfly

function compute(N=3)

# ==========
# 0. Options
# ==========

# algorithm 
algorithm = "no_wrapping_effect"

# plot backend
plot_backend = "pyplot_inline"

# filename of the plot
plot_name = @filename_to_png

# perform intermediate overapproximations?
intermediate_overapproximations = true

# include bloating factor
ct_bloating_factor = true

# time discretization
δ = 0.1

# tolerance
ɛ = 1e-3

# choose one of: "arch_bench_no_input", "arch_bench_with_input"
config = "arch_bench_with_input"

# choose one of: "" (= no transformation), "schur" (= Schur transformation)
transformationMethod = ""

# variables to plot (0 == time)
plot_vars = [0, 2]

# ===========================
# Step 1: System Construction
# ===========================
println("System Construction...")
tic()

A = sparse([0 1; -336.0 -0.1833])
B = sparse([0 0 0; 0.004123 -0.00062444 0.031733])
n = size(A, 1)


# input set
if config == "arch_bench_no_input"
    X0 = BallInf(zeros(n), .0001) # -0.0001 <= xi <= 0.0001 for all i
    U = VoidSet(n)   # [VoidSet(n) for i in 1:N]

elseif config == "arch_bench_with_input"
    X0 = BallInf(zeros(n), .0001) # -0.0001 <= xi <= 0.0001 for all i
    Uint = CartesianProductArray([BallInf([0.05], .05), BallInf([0.9], .1), BallInf([0.95], .05)])
    U = B * Uint

else
    error("undefined configuration of initial states")
end

# instantiate continuous LTI system
D = ContinuousSystem(A, X0, U)

toc()

# ===========================
# Step 2: Time Discretization
# ===========================
println("Time Discretization...")
tic()
Δ = discretize(D, δ, ct_bloating_factor=ct_bloating_factor)
toc()

# ========================
# Step 2.5: Transformation
# ========================
if transformationMethod != ""
    println("Transformation...")
    tic()
    (Δ, transformationMatrix) = transform(transformationMethod, Δ, plot_vars)
    toc()
else
    transformationMatrix = nothing
end

# ====================================
# Step 3: Reachable States Computation
# ====================================
println("Reachable States Computation...")
tic()
Rsets = reach(Δ, ɛ, N, algorithm, intermediate_overapproximations=intermediate_overapproximations)
toc()

# ==================
# Step 4: Projection
# ==================
println("Projection...")
tic()
RsetsProj = project_reach(plot_vars, n, N, δ, ɛ, Rsets, algorithm, transformationMatrix=transformationMatrix)
toc()

# ================
# Step 5: Plotting
# ================
println("Plotting...")
tic()
plot_polygon(RsetsProj, backend=plot_backend, gridlines=true, name=plot_name)
toc()

nothing
end # function
