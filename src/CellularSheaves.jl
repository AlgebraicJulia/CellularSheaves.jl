""" CellularSheaves.jl is a Julia package for working with cellular sheaves and sheaf Laplacians.
"""
module CellularSheaves

using Reexport

include("BlockSparseArrays/src/BlockSparseArrays.jl")
@reexport using .BlockSparseArrays

include("network_sheaves/NetworkSheaves.jl")
@reexport using .NetworkSheaves

include("ControlSheaves/ControlSheaves.jl")
using .ControlSheaves
export ControlSheaves

# The interior-point method lives in Mumblebee.jl; `CellularSheaves.IPM` is kept
# as an alias so existing `using CellularSheaves.IPM` code keeps working.
using Mumblebee: Mumblebee
const IPM = Mumblebee.IPM
export Mumblebee, IPM

end