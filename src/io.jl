"""
   @relpath(name)

Returns the absolute path to `name` relative to the executing script.

### Input

- `name` -- file name

### Output

A string.

### Notes

This macro is equivalent to `joinpath(@__DIR__, name)`.

The `@relpath` macro is used in model scripts to load data files relative to the
location of the model, wihout having to change directory of your Julia session.

For instance, suppose that the folder `/home/projects/model` contains the script
`my_model.jl`, and suppose that the data file `my_data.dat` located in the same
directory is required to be loaded by `my_model.jl`.
Then,

```julia
# suppose my path working directory is /home/julia/  so we run the script as
# julia -e "include("../projects/models/my_model.jl")"

# in the model file /home/projects/models/my_model.jl we write:
d = open(@relpath "data.dat")
# do stuff with d
```
In this example, the macro `@relpath "data.dat"` evaluates to the string
`/home/projects/models/data.dat`. If the script `my_model.jl` only had
`d = open("dat.dat")`, without `@relpath`, this command would fail as Julia would
have looked for `data.dat` in the *running* directory, giving an error that the file
`/home/julia/data.dat` is not found.

As mentioned above, the interest in using `@relpath` is that the script `my_data.jl`
can now be included from the command line starting Julia from any folder, it isn't
needed to cd to the the folder containing the data file.
"""
macro relpath(name::String)
    __source__.file === nothing && return nothing
    _dirname = dirname(String(__source__.file))
    dir = isempty(_dirname) ? pwd() : abspath(_dirname)
    return joinpath(dir, name)
end