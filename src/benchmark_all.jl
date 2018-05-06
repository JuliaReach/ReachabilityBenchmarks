# script to run all SLICOT models in several benchmark settings

using Reachability

# model list
models = [
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
include(@relpath "check_1D_discrete.jl")
include(@relpath "check_1D_dense.jl")

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


# -- check_1D_discrete --
# no modification
check_1D_discrete(models)


# -- check_1D_dense --
# remove iss and fom; define deltas
models_filtered = filter(e -> e ∉ ["iss", "fom"], models)
deltas = [
    1e-3,
    2e-3,
    3e-4,
    1e-3,
    5e-5,
    4e-4,
    3e-1
   ]
check_1D_dense(models_filtered, deltas)


return nothing
