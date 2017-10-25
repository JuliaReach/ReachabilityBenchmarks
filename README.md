# ReachabilityBenchmarks

Benchmark suite for reach set computations

## Installation

Clone this repository, and install the Julia package [Reachability.jl](https://github.com/JuliaReach/Reachability.jl),
by following the instruction in the [installation section](https://github.com/JuliaReach/Reachability.jl#installing).

The models are given as Julia scripts, which you can run by including them in Julia's REPL, e.g.

```julia
julia> include("models/SLICOT/iss/iss.jl")
```

The models stored in MAT files are loaded using the [MAT.jl](https://github.com/JuliaIO/MAT.jl) Julia package.

There are also bash scripts to run [SpaceEx](http://spaceex.imag.fr/) for most of the models.
These scripts assume that `spaceex` (and some other tools for creating plots) are added to the PATH.

## Usage

The scripts are accompanied by the `compute` function, which calculates reachable states, or checks a safety property,
depending on the model's options. To use the default options, do

```julia
julia> compute()
```

## Running all the benchmarks

This repository provides a bash script that will execute Julia on all existent
[SLICOT](http://slicot.org/20-site/126-benchmark-examples-for-model-reduction) models.
The results are stored in the folder of each model.
The script can also be modified (by just uncommenting the respective parts) to run SpaceEx on all the models.
