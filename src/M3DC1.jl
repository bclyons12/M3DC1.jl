__precompile__()
module M3DC1

using HDF5
using StaticArrays
using LazyGrids
using Plots

include("mesh.jl")
export read_mesh

include("field.jl")
export eval_field

end
