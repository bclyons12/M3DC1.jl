__precompile__()
module M3DC1

using HDF5
using StaticArrays
using RecipesBase

const RZsize = (800,600)

include("mesh.jl")
export BBtree, BBnode, find_sibling, insert, read_mesh

include("field.jl")
export eval_field

end
