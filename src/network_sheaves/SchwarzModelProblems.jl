# Standard model problems for testing and benchmarking the Schwarz methods:
# finite-difference Poisson problems on grid domains, including the notched
# rectangle with its two re-entrant corners.
module SchwarzModelProblems

export GridDomain, unit_square, notched_rectangle, poisson_matrix, box_partition, grid_values

using ArgCheck: @argcheck
using SparseArrays

"""
    GridDomain(predicate, h, width, height)

A planar domain in ``[0, \\mathrm{width}] \\times [0, \\mathrm{height}]``,
sampled on the uniform grid ``(a h, b h)`` of spacing `h`. A grid point is an
unknown (a dof) when it lies strictly inside the box and `predicate(x, y)`
holds. All other grid points are on the boundary, where homogeneous Dirichlet
conditions are imposed.

Fields: `h`, `width`, `height`, `inside` (a `BitMatrix` over the interior grid,
`true` for dofs), `index` (dof number of each interior grid point, `0` if it is
not a dof), and `points` (the ``(x, y)`` coordinates of each dof).
"""
struct GridDomain
    h::Float64
    width::Float64
    height::Float64
    inside::BitMatrix
    index::Matrix{Int}
    points::Vector{NTuple{2,Float64}}
end

function GridDomain(predicate, h::Real, width::Real, height::Real)
    @argcheck h > 0 && width > h && height > h
    nx = round(Int, width / h) - 1
    ny = round(Int, height / h) - 1
    @argcheck (nx + 1) * h ≈ width && (ny + 1) * h ≈ height "width and height must be multiples of h"
    inside = BitMatrix(undef, nx, ny)
    index = zeros(Int, nx, ny)
    points = NTuple{2,Float64}[]
    for b in 1:ny, a in 1:nx
        x, y = a * h, b * h
        inside[a, b] = predicate(x, y)
        if inside[a, b]
            push!(points, (x, y))
            index[a, b] = length(points)
        end
    end
    return GridDomain(Float64(h), Float64(width), Float64(height), inside, index, points)
end

"""
    unit_square(m) -> GridDomain

The unit square with ``m \\times m`` interior grid points, ``h = 1/(m+1)``.
"""
unit_square(m::Integer) = GridDomain((x, y) -> true, 1 / (m + 1), 1.0, 1.0)

"""
    notched_rectangle(m; width=2, height=1, notch_width=height/4, notch_depth=height/2) -> GridDomain

The notched rectangle: ``[0, \\mathrm{width}] \\times [0, \\mathrm{height}]``
with a slot of width `notch_width` and depth `notch_depth` cut from the middle
of the top edge, sampled with ``h = \\mathrm{height}/(m+1)``. The bottom of the
slot has two *re-entrant* corners of angle ``3\\pi/2``. Near them the solution
of the Poisson problem behaves like ``r^{2/3}`` and its gradient is unbounded.
This is the standard test of how a method copes with corner singularities.

When the slot cuts through a subdomain boundary, the subdomains next to it get
an irregular shape, with corners of their own.
"""
function notched_rectangle(m::Integer; width::Real=2.0, height::Real=1.0,
                           notch_width::Real=height / 4, notch_depth::Real=height / 2)
    @argcheck 0 < notch_width < width && 0 < notch_depth < height
    h = height / (m + 1)
    tol = h / 4
    in_notch(x, y) = abs(x - width / 2) <= notch_width / 2 + tol && y >= height - notch_depth - tol
    return GridDomain((x, y) -> !in_notch(x, y), h, width, height)
end

"""
    poisson_matrix(dom::GridDomain) -> SparseMatrixCSC

The 5-point finite-difference discretization of ``-\\Delta u`` on `dom` with
homogeneous Dirichlet conditions on its boundary (including the boundary of any
notch). The matrix is symmetric positive definite and an M-matrix.
"""
function poisson_matrix(dom::GridDomain)
    nx, ny = size(dom.inside)
    I, J, V = Int[], Int[], Float64[]
    scale = 1 / dom.h^2
    for b in 1:ny, a in 1:nx
        k = dom.index[a, b]
        k == 0 && continue
        push!(I, k); push!(J, k); push!(V, 4scale)
        for (da, db) in ((1, 0), (-1, 0), (0, 1), (0, -1))
            a2, b2 = a + da, b + db
            (1 <= a2 <= nx && 1 <= b2 <= ny) || continue
            k2 = dom.index[a2, b2]
            k2 == 0 && continue
            push!(I, k); push!(J, k2); push!(V, -scale)
        end
    end
    n = length(dom.points)
    return sparse(I, J, V, n, n)
end

"""
    box_partition(dom::GridDomain, px, py) -> Vector{Int}

Cut the bounding box of `dom` into ``p_x \\times p_y`` equal boxes and label
every dof with its box. Boxes that contain no dof (for example inside a notch)
are dropped and the rest are numbered consecutively, row by row.
"""
function box_partition(dom::GridDomain, px::Integer, py::Integer)
    @argcheck px >= 1 && py >= 1
    box(x, y) = (clamp(ceil(Int, y / dom.height * py), 1, py) - 1) * px +
                clamp(ceil(Int, x / dom.width * px), 1, px)
    raw = [box(x, y) for (x, y) in dom.points]
    used = sort!(unique(raw))
    relabel = Dict(b => i for (i, b) in enumerate(used))
    return [relabel[b] for b in raw]
end

"""
    grid_values(dom::GridDomain, u) -> Matrix{Float64}

Place the dof values `u` on the interior grid for plotting, with `NaN` outside
the domain. Rows run along ``x`` and columns along ``y``.
"""
function grid_values(dom::GridDomain, u::AbstractVector)
    @argcheck length(u) == length(dom.points)
    V = fill(NaN, size(dom.inside))
    for (ab, k) in pairs(dom.index)
        k == 0 || (V[ab] = u[k])
    end
    return V
end

end
