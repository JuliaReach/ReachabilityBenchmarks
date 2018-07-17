# script to run all SLICOT models in several benchmark settings

# Plots is required here, unless never used, to avoid an issue with "world age"
using Reachability, Plots

# model list
models = String[
    "motor",
    "building",
    "pde",
    "heat",
    "iss",
    "beam",
    "mna1",
    "fom",
    "mna5"
   ]

# load models and benchmark scripts
model_library = include(@relpath "model_library.jl")
for model in models
    path = (@relpath "../") * model_library[model]
    include(path)
end
include(@relpath "reach_1D_single.jl")
include(@relpath "reach_1D_all.jl")
include(@relpath "reach_2D_box_two.jl")
include(@relpath "reach_2D_eps_two.jl")
include(@relpath "reach_kD_all.jl")
include(@relpath "check_1D_discrete.jl")
include(@relpath "check_kD_discrete.jl")
include(@relpath "check_1D_dense.jl")
include(@relpath "check_kD_dense.jl")

# create plots (option)?
create_plots = true
if create_plots
    include(@relpath "plot_reach.jl")
end

# run benchmarks


# -- reach_1D_single --
# no modification, optionally create plots
reach_1D_single(models, create_plots)


# -- reach_1D_all --
# remove mna5
models_filtered = filter(e -> e != "mna5", models)
reach_1D_all(models_filtered)


# -- reach_2D_box_two --
# no modification, optionally create plots
reach_2D_box_two(models, create_plots)


# -- reach_2D_eps_two --
# no modification, optionally create plots
reach_2D_eps_two(models, create_plots)


# -- reach_kD_all --
# only use pde, optionally create plots
models_filtered = filter(e -> e ∈ ["pde"], models)
reach_kD_all(models_filtered, create_plots)


# -- check_1D_discrete --
# no modification
check_1D_discrete(models)


# -- check_kD_discrete --
# only use iss and fom; use default partition
models_filtered = filter(e -> e ∈ ["iss", "fom"], models)
check_kD_discrete(models_filtered)


# -- check_1D_dense --
# remove iss and fom; define deltas
models_filtered = filter(e -> e ∉ ["iss", "fom"], models)
delta_map = Dict(
    "motor" => 1e-3,
    "building" => 2e-3,
    "pde" => 3e-4,
    "heat" => 1e-3,
    "beam" => 5e-5,
    "mna1" => 4e-4,
    "mna5" => 3e-1
   )
deltas = Vector{Float64}(length(models_filtered))
for (i, model) in enumerate(models_filtered)
    deltas[i] = delta_map[model]
end
check_1D_dense(models_filtered, deltas)


# -- check_kD_dense --
# only use iss; use default partition
models_filtered = filter(e -> e ∈ ["iss"], models)
delta_map = Dict(
    "iss" => 6e-4,
   )
deltas = Vector{Float64}(length(models_filtered))
for (i, model) in enumerate(models_filtered)
    deltas[i] = delta_map[model]
end
check_kD_dense(models_filtered, deltas)


return nothing
