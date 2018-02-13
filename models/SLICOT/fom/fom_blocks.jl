#=
Feature of the upper part of the FOM model (2x2 blocks with [-1, k; -k, -1] for k = 100, 200, 400).
This model is parametric in the number of blocks K.
=#
include("../../src/LazySets.jl")
include("../../src/Systems.jl")
include("../../src/Approximations.jl")
include("../../src/Reachability.jl")
include("../../src/Transformations.jl")
include("../../src/Benchmarks.jl")
include("../../src/Properties.jl")

using LazySets, Systems, Reachability, Approximations, Transformations, Benchmarks, Properties, MAT, Expokit, PyPlot
#using Gadfly

compute(N::Int64=3, δ::Float64=-1., options::Pair{Symbol,<:Any}...) = compute(N, δ, Options(Dict{Symbol,Any}(options)))

function compute(N::Int64=3, δ::Float64=-1., userProvidedOptions::Options=Options())

# ==========
# 0. Options
# ==========

# model parameter (number of blocks)
K = 3

# use default for unspecified options
options = merge(defaultGlobalOptions, defaultModelOptions(N, K), userProvidedOptions)

if δ < 0.
    δ = options[:δ]
end

# ===========================
# Step 1: System Construction
# ===========================
println("System Construction...")
tic()

# ambient dimension
n = K*2

A = -speye(n)
for i in 1:K
    hundreds = 2^(i-1)*100
    A[2*i-1,2*i] = hundreds
    A[2*i,2*i-1] = -hundreds
end

X0 = BallInf(zeros(n), 0.0001)

B = sparse(1:n, ones(n), fill(10., n))
U = B * BallInf([0.0], 1.0)

# instantiate continuous LTI system
D = ContinuousSystem(A, X0, U)

toc()

# ===========================
# Step 2: Time Discretization
# ===========================
println("Time Discretization...")
tic()
Δ = discretize(D, δ,
    ct_bloating_factor=options[:ct_bloating_factor], pade_expm=options[:pade_expm])
toc()

# ========================
# Step 2.5: Transformation
# ========================
if options[:transformationMethod] != ""
    println("Transformation...")
    tic()
    (Δ, transformationMatrix) = transform(options[:transformationMethod], Δ,
                                options[:plot_vars])
    toc()
else
    transformationMatrix = nothing
end

if options[:mode] == "reach"

# ====================================
# Step 3: Reachable States Computation
# ====================================
println("Reachable States Computation...")
tic()
Rsets = reach(Δ, N; algorithm=options[:algorithm], ɛ=options[:ɛ],
        block=options[:block], blocks=options[:blocks],
        assume_sparse=options[:assume_sparse],
        iterative_refinement=options[:iterative_refinement],
        assume_homogeneous=options[:assume_homogeneous])
toc()

# ==================
# Step 4: Projection
# ==================
println("Projection...")
tic()
if options[:project_output]
    # load the output that we want to observe
    # Example (assuming y = ... in the file named issOutputy3.sage)
    # $ sage outSXmatrixtoJulia.sage 270 issOutputy3.sage MAT 
    M = sparse(ones(n), 1:n, fill(10, n))
else
    M = nothing
end
RsetsProj = project_reach(options[:plot_vars], n, δ, Rsets, options[:algorithm],
    options[:plot_labels], ɛ=options[:ɛ], transformationMatrix=transformationMatrix, M=M)
toc()

# ================
# Step 5: Plotting
# ================
println("Plotting...")
tic()
plot_polygon(view(RsetsProj, options[:plot_indices]),
    backend=options[:plot_backend], gridlines=options[:gridlines],
    name=options[:plot_name], plot_labels=options[:plot_labels])
toc()

elseif options[:mode] == "check"

# =========================
# Step 3: Property checking
# =========================
println("Property Checking...")
tic()
answer = check_property(Δ, N; algorithm=options[:algorithm], ɛ=options[:ɛ],
	 block=options[:block], blocks=options[:blocks],
         assume_sparse=options[:assume_sparse],
         iterative_refinement=options[:iterative_refinement],
         assume_homogeneous=options[:assume_homogeneous], property=options[:property])
toc()
if answer == 0
    println("The property is satisfied!")
    return true
else
    println("The property is violated at index ", answer, " (time point ", answer * δ, ")!")
    return false
end

else
    error("Unsupported mode.")
end

nothing
end # function


#=
Default options for this model (can be overridden by user-defined options)
=#
function defaultModelOptions(N::Int64, K::Int64)::Options
    return Options(
        :assume_sparse=>false,
        :mode=>"reach",
        :algorithm=>"explicit_blocks",
        :δ=>0.01,
        :block=>@block_id(1),
        :blocks=>1:K,
        :property=>LinearConstraintProperty(ones(2*K), 20.), # sum_i x_i < 20
        :project_output=>true,
        :plot_vars=>[0, 1],
        :plot_name=>@filename_to_png,
        :plot_indices=>1 : N #range_last_x_percent(N, 10, 3)
        )
end
