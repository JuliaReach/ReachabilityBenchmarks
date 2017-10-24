# ReachabilityBenchmarks

Benchmark suite for reach set computations

## Installation

Clone this repository, and install the Julia package [Reachability.jl](https://github.com/JuliaReach/Reachability.jl),
by following the instruction in the [installation section](https://github.com/JuliaReach/Reachability.jl#installing).

The examples are given as Julia scripts, which you can run by including them in Julia's REPL, e.g.

```julia
julia> include("src/SLICOT/iss")
```

## Usage

The scripts are accompanied by the `compute` function, which calculates reachable states, or a safety property,
depending on the model's options. To use the default options, do

```julia
julia> compute()
```

## Running all the benchmarks

This repository providees a bash script that will execute Julia on all existent SLICIT models. The results are stored
in the folder of each model.
