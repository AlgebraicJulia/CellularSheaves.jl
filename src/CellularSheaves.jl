""" CellularSheaves.jl is a Julia package for working with cellular sheaves and sheaf Laplacians.
"""
module CellularSheaves

using Reexport

include("BlockSparseArrays/src/BlockSparseArrays.jl")
@reexport using .BlockSparseArrays

# IPM depends only on BlockSparseArrays, and ControlSheaves depends on IPM for the
# local conic safety filters, so it is loaded ahead of the sheaf modules.
include("IPM/IPM.jl")
using .IPM
export IPM

include("network_sheaves/NetworkSheaves.jl")
@reexport using .NetworkSheaves

include("ControlSheaves/ControlSheaves.jl")
using .ControlSheaves
export ControlSheaves

end