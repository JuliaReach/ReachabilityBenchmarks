# This script reads a file in gen format containing polygons and stores them
# in a Julia array. The sequence of polygons can be plotted using a plotting
# backend, such as GR, PyPlot, or Gadfly.
#
# EXAMPLES:
#
# We load the results for the BEAM model, which contain 20000 polygons:
#
#     julia> include("plot_gen.jl")
#     julia> using Reachability
#     julia> import GR
#     julia> Xk = read_gen("SpaceEx_plots/beam_t_x89_20_spaceex.txt");
#     julia> R = ReachSolution(Xk, Options(:plot_vars => [1,2]));
#     julia> plot(R)
#     alternative, more spohisticated plotting command
#     julia> plot(R, color="red", xlabel="t", ylabel = "x89", ylims=(-55000, 55000))

# choose a plotting backend, eg. PyPlot or Gadfly
using Plots, LazySets

"""
    read_gen(filename)

Read a sequence of polygons stored in vertex representation (gen format).

### Input

- `filename` -- path of the file containing the polygons

The x and y coordinates of each vertex are separated by an empty space, and
polygons are separated by empty lines. For example:

```
1.01 1.01
0.99 1.01
0.99 0.99
1.01 0.99

0.908463 1.31047
0.873089 1.31047
0.873089 1.28452
0.908463 1.28452
```

### Output

A sequence of matrices, where the first row is for the x coordinates, and the
second row for the y coordinates.

In the example of above:

```julia
    julia> P = read_gen("test.gen")
    2-element Array{Array{Float64,2},1}:
     [1.01 0.99 0.99 1.01; 1.01 1.01 0.99 0.99]
     [0.908463 0.873089 0.873089 0.908463; 1.31047 1.31047 1.28452 1.28452]
```

### Warning

The input file should end with at least one empty line before the end of files
this is to detect the last polygon.
"""
function read_gen(filename::String)
    Mi = Vector{Vector{Float64}}()
    P = Vector{VPolygon{Float64}}()
    # detects when we finished reading a new polygon, needed because polygons
    # may be separated by more than one end-of-line
    new_polygon = true
    open(filename) do f
        for line in eachline(f)
          if !isempty(line)
              push!(Mi, map(x -> parse(Float64, x), split(line)))
              new_polygon = true
          elseif isempty(line) && new_polygon
              push!(P, VPolygon(Mi))
              Mi = Vector{Vector{Float64}}()
              new_polygon = false
          end
        end
    end
    return P
end

function read_gen_old(filename::String)::Array{Matrix{Float64}, 1}
    Mi = Array{Vector{Float64}, 1}()
    P = Array{Matrix{Float64}, 1}()
    # detects when we finished reading a new polygon, needed because polygons
    # may be separated by more than one end-of-line
    new_polygon = true
    open(filename) do f
        for line in eachline(f)
          if !isempty(line)
              push!(Mi, map(x -> parse(Float64, x), split(line)))
              new_polygon = true
          elseif isempty(line) && new_polygon
              push!(P, hcat(Mi...))
              Mi = Array{Vector{Float64}, 1}()
              new_polygon = false
          end
        end
    end
    return P
end

"""
    plot_polygon(P; backend, [name], [gridlines], [plot_polygon])

Plot one or more polygons given as an array of matrices.

INPUT:

- ``P`` -- a polygon, given as a list of matrices, where the first row is for
           the x coordinates and the second row for the y coordinates
- ``backend``     -- (optional, default: ``'pyplot'``): select the plot backend; valid
  options are:
                 -  ``'pyplot_savefig'`` -- use PyPlot package, save to a file
                 -  ``'pyplot_inline'`` -- use PyPlot package, showing in external program
                 - ``'gadfly'`` -- use Gadfly package, showing in browser
                 - ``''`` -- (empty string), return nothing, without plotting
- ``name``        -- (optional, default: ``'plot.png'``) the filename of the plot
  (if it is saved to disk)
- ``gridlines``   -- (optional, default: false) to display or not gridlines in
                   the output plot
- ``color``       -- (optional, default: ``'red'``) the color of the plotting lines
- ``plot_labels`` -- (optional, default: ``['', '']``) the labels of the axes
"""
function plot_polygon(P::Array{Matrix{Float64}, 1};
                      backend="pyplot_savefig",
                      name="plot.png",
                      gridlines=false,
                      color="red",
                      plot_labels::Vector{String}=["", ""])
    if backend in ["pyplot_inline", "pyplot_savefig"]
        if !isdefined(:PyPlot)
            error("this backend requires that your script loads the PyPlot module")
        end
        eval(Expr(:using, :PyPlot))
        if backend == "pyplot_savefig"
            PyPlot.ioff()  # turn off interactive plotting
            fig = PyPlot.figure()
        end
        gridlines ? grid("on") : grid("off")
        # check if P is iterable
        applicable(start, P) ? nothing : P = [P]
        PyPlot.xlabel(plot_labels[1])
        PyPlot.ylabel(plot_labels[2])
        for pi in P
            xcoords = pi[1, :]
            ycoords = pi[2, :]
            PyPlot.plot(xcoords, ycoords, color=color, linewidth=0.6)
        end
        if backend == "pyplot_savefig"
            PyPlot.savefig(name, bbox_inches="tight")
            PyPlot.close(fig)
        end

    elseif backend in ["gadfly"]
        if !isdefined(:Gadfly)
            error("this backend requires that your script loads the Gadfly module")
        end
        eval(Expr(:using, :Gadfly))
        layers_array = []
        # check if P is iterable
        applicable(start, P) ? nothing : P = [P]
        for pi in P
            xcoords = pi[1, :]
            ycoords = pi[2, :]
            push!(layers_array, Gadfly.layer(x=xcoords, y=ycoords, Geom.polygon(preserve_order=false, fill=true)))
        end
        Gadfly.plot(layers_array...)

    elseif backend == ""
        return nothing
    else
        error("plot backend not valid")
    end
end

"""
    plot_all(path)

Plot all txt files (in gen format) of a specified folder.

INPUT:

- ``path`` -- (optional, default: ``'.'``, the current directory) this is the
              directory where the data files are located

OUTUPT:

One or more files with the extension ``.txt.png``.

To give an idea, for 6 files where each file contains 20000 polygons,
it takes around 30 minutes.
"""
function plot_all(path::String=".")
    files = filter(x -> contains(x, ".txt"), readdir(path))
    for f in files
        println("Plotting ", f)
        plot_polygon(read_gen(f), name=f * ".png", gridlines=true)
    end
end

return nothing
