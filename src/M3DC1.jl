__precompile__()
module M3DC1

using HDF5
using StaticArrays
using RecipesBase

const RZsize = (800,600)

include("mesh.jl")
export BBtree, read_mesh, find_element

include("field.jl")
export eval_field

include("line.jl")
export Line, line_integral

end
