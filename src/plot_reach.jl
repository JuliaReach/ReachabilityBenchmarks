# plots a reach set

# The next line is needed for a world age issue in Plots:
# https://github.com/JuliaPlots/Plots.jl/issues/1047
# This line must also come before 'using Plots' to prevent an error:
# ERROR: MethodError: no method matching convert(::Type{AssertionError}, ::String)
import GR

using Reachability, Plots

import Reachability.@timing

function plot_reach(result::AbstractSolution, name::String="plot")
    if result.options[:mode] != "reach"
        # ignore plotting command
        return
    end
    info("Plotting...")
    # GR plotting backend does not support subindices like x₁ in savefig
    is_gr = Plots.backend() == Plots.GRBackend()
    @timing begin
        plot(result; use_subindices=!is_gr)
        savefig("$name.png")
    end
end
