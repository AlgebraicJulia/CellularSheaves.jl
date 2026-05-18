module PDESheaves

"""
    PDESheaves

A lightweight module that builds a rectilinear grid, creates a graph whose
vertices correspond to subdomains, and assembles a discretized PDE operator as a
linear operator (or matrix).  The resulting data are wrapped in a
`CellularSheaves.EuclideanSheaf` so that the space of global sections matches
the solution space of the PDE.

The design follows the conventions of the existing `CellularSheaves` codebase:
* Graph construction uses `Graphs.SimpleGraph`.
* Restriction maps project a subdomain solution onto the shared edge (boundary)
  degrees‑of‑freedom via simple averaging; the API can be swapped for more
  sophisticated trace operators.
* The PDE operator is supplied by the user as a function `local_operator(dims,
  Δx, Δy) -> AbstractMatrix` that builds a local stiffness/mass matrix on a
  rectangular cell.
"""
export build_grid_graph, discretize_operator, edge_restriction,
       pde_sheaf

using ..CellularSheaves: EuclideanSheaf, add_vertex_stalk!, add_sheaf_edge!
using Graphs
using LinearAlgebra
using BlockArrays

"""
    build_grid_graph(nx::Int, ny::Int) -> (g::SimpleGraph, vertex_ids::Vector{Int})

Create a rectilinear grid graph with `nx` cells in the x‑direction and `ny`
cells in the y‑direction.  Each cell becomes a vertex; edges connect
4‑connected neighbours.  The function returns the graph and a vector mapping the
vertex index to a unique integer identifier (simply `1:nv(g)` for now).
"""
function build_grid_graph(nx::Int, ny::Int)
    @assert nx > 0 && ny > 0 "grid dimensions must be positive"
    g = SimpleGraph(nx * ny)
    # Helper to convert (i,j) → vertex index (1‑based)
    idx(i, j) = (j - 1) * nx + i
    for j in 1:ny, i in 1:nx
        v = idx(i, j)
        if i < nx
            add_edge!(g, v, idx(i + 1, j))
        end
        if j < ny
            add_edge!(g, v, idx(i, j + 1))
        end
    end
    return g, collect(1:nv(g))
end

"""
    discretize_operator(g::SimpleGraph, subdim::Int,
                       Δx::Float64, Δy::Float64,
                       local_op::Function) -> Vector{Matrix}

Construct a local operator matrix for each vertex (cell).  `local_op` receives
the cell dimension (e.g. number of degrees‑of‑freedom per cell) and the cell
spacing `Δx, Δy`, and returns a square matrix of size `subdim`.
The function returns a vector of matrices, one per vertex, ready to be used as
vertex stalk dimensions in `EuclideanSheaf`.
"""
function discretize_operator(g::SimpleGraph, subdim::Int,
                            Δx::Float64, Δy::Float64,
                            local_op::Function)
    # For a regular grid the same local operator applies everywhere.
    op = local_op(subdim, Δx, Δy)
    return fill(op, nv(g))
end

"""
    edge_restriction(g::SimpleGraph, subdim::Int) -> Vector{Matrix}

Create restriction maps for each edge.  For a simple trace operator we take
the average of the two neighbouring cell values, which for a scalar field is the
identity of size `subdim`.  For vector‑valued fields you may replace this with a
more elaborate projection.  The return value is a vector of matrices ordered as
`(src→edge, dst→edge)` for each edge in the order returned by `edges(g)`.
"""
function edge_restriction(g::SimpleGraph, subdim::Int)
    # Identity restriction (no reduction) – edge stalk dimension equals vertex
    # stalk dimension.  Users can later replace with a true trace if needed.
    I = Matrix{Float64}(I, subdim, subdim)
    maps = Matrix{Float64}[]
    for e in edges(g)
        push!(maps, I)   # src → edge
        push!(maps, I)   # dst → edge
    end
    return maps
end

"""
    pde_sheaf(nx::Int, ny::Int, subdim::Int,
             Δx::Float64, Δy::Float64,
             local_op::Function) -> EuclideanSheaf{Float64}

High‑level constructor that builds a rectilinear grid, assembles the associated
graph, creates vertex stalks of dimension `subdim`, and installs the
restriction maps using `edge_restriction`.  The returned `EuclideanSheaf`
can be used directly with `nullspace_ldlt`, `nearest_global_section`, etc.
"""
function pde_sheaf(nx::Int, ny::Int, subdim::Int,
                  Δx::Float64, Δy::Float64,
                  local_op::Function)
    g, _ = build_grid_graph(nx, ny)
    # Vertex stalk dimensions – all identical for a regular grid
    stalks = repeat([subdim], nv(g))
    sheaf = EuclideanSheaf{Float64}(stalks)
    # Add vertices to the underlying graph (done in EuclideanSheaf ctor)
    # Add edges with restriction maps
    for e in edges(g)
        i, j = src(e), dst(e)
        rm = edge_restriction(g, subdim)  # returns a flat vector; we use first two entries
        # rm vector order matches pushes: src→edge, dst→edge
        add_sheaf_edge!(sheaf, i, j, rm[1], rm[2])
    end
    return sheaf
end

end # module PDESheaves
