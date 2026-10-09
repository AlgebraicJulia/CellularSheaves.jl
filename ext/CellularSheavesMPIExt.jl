# MPI communicator for the boxes of a distributed grid (GridSchwarz): one box
# per rank, halo exchanges between face neighbours, reductions and gathers.
module CellularSheavesMPIExt

using CellularSheaves
using MPI
using CellularSheaves.NetworkSheaves: GridSchwarz
using CellularSheaves.NetworkSheaves.GridSchwarz: BoxCommunicator, BoxLayout, _slabs

struct MPIBoxes <: BoxCommunicator
    comm::MPI.Comm
    buffers::Dict{Any,Any}
end

GridSchwarz.mpi_boxes(comm::MPI.Comm) = MPIBoxes(comm, Dict{Any,Any}())

GridSchwarz.box_count(c::MPIBoxes) = MPI.Comm_size(c.comm)
GridSchwarz.box_rank(c::MPIBoxes) = MPI.Comm_rank(c.comm)
GridSchwarz.box_allreduce(c::MPIBoxes, x, op) = MPI.Allreduce(x, op, c.comm)

function GridSchwarz.box_allgather(c::MPIBoxes, v::AbstractVector)
    local_v = Vector(v)
    counts = MPI.Allgather(length(local_v), c.comm)
    out = similar(local_v, sum(counts))
    MPI.Allgatherv!(local_v, MPI.VBuffer(out, counts), c.comm)
    offsets = cumsum([0; counts])
    return [out[(offsets[r] + 1):offsets[r + 1]] for r in eachindex(counts)]
end

# Device arrays go to MPI directly when the library can read them (CUDA-aware
# MPI); otherwise they are staged through host buffers.
_mpi_ready(x::Array) = true
_mpi_ready(x) = MPI.has_cuda()

function _buffers(c::MPIBoxes, x::AbstractArray, n::Int)
    get!(c.buffers, (typeof(x), n)) do
        device = (similar(x, n), similar(x, n))
        host = _mpi_ready(x) ? device : (Vector{eltype(x)}(undef, n), Vector{eltype(x)}(undef, n))
        (device..., host...)
    end
end

# Send the slab `send` of x to `dest` and receive the slab `recv` from `source`
# (either may be -1: no neighbour on that side).
function _shift!(c::MPIBoxes, x::AbstractArray, send, recv, dest::Int, source::Int, tag::Int)
    dest < 0 && source < 0 && return x
    n = prod(length.(send))
    dsend, drecv, hsend, hrecv = _buffers(c, x, n)
    dest >= 0 && copyto!(dsend, view(x, send...))
    hsend === dsend || (dest >= 0 && copyto!(hsend, dsend))
    MPI.Sendrecv!(hsend, hrecv, c.comm; dest=dest < 0 ? MPI.PROC_NULL : dest, sendtag=tag,
        source=source < 0 ? MPI.PROC_NULL : source, recvtag=tag)
    if source >= 0
        hrecv === drecv || copyto!(drecv, hrecv)
        copyto!(view(x, recv...), reshape(drecv, length.(recv)))
    end
    return x
end

function GridSchwarz.exchange!(c::MPIBoxes, layout::BoxLayout{D}, x::AbstractArray) where {D}
    for d in 1:D
        s = _slabs(layout, d)
        lower, upper = layout.neighbors[d]
        _shift!(c, x, s.send_lower, s.recv_upper, lower, upper, 2d)          # downwards
        _shift!(c, x, s.send_upper, s.recv_lower, upper, lower, 2d + 1)      # upwards
    end
    return x
end

end # module
