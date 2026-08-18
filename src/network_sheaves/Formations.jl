module Formations

using LinearAlgebra
using ..EuclideanSheaves: EuclideanSheaf, add_sheaf_edge!
using ArgCheck: @argcheck

export se3_translation_matrix, se3_rotation_matrix, se3_affine_matrix, affine_translation_matrix, build_escort_topology, build_escort_ring, build_escort_clique

"""
    se3_translation_matrix(d::AbstractVector)

Returns the 4x4 homogeneous affine translation matrix.
"""
function se3_translation_matrix(d::AbstractVector)
    @assert length(d) == 3 "Translation vector must be 3-dimensional"
    [I(3) -d; 0 0 0 1]
end

"""
    affine_translation_matrix(d::AbstractVector)

Returns the `(n+1)x(n+1)` homogeneous affine translation matrix for a translation
vector `d` of arbitrary dimension `n`. Generalizes `se3_translation_matrix` (the
`n=3` case) to any stalk dimension, for non-SE(3) formations (e.g. planar agents).
"""
function affine_translation_matrix(d::AbstractVector)
    n = length(d)
    [I(n) -d; zeros(1, n) 1]
end

"""
    se3_rotation_matrix(R::AbstractMatrix)

Returns the 4x4 homogeneous affine rotation matrix.
"""
function se3_rotation_matrix(R::AbstractMatrix)
    @assert size(R) == (3, 3) "Rotation matrix must be 3x3"
    [R zeros(3); 0 0 0 1]
end

"""
    se3_rotation_matrix(; θx::Real=0.0, θy::Real=0.0, θz::Real=0.0)

Returns the 4x4 homogeneous affine rotation matrix corresponding to Euler angles (XYZ convention).
"""
function se3_rotation_matrix(; θx::Real=0.0, θy::Real=0.0, θz::Real=0.0)
    Rx = [1.0 0.0 0.0; 0.0 cos(θx) -sin(θx); 0.0 sin(θx) cos(θx)]
    Ry = [cos(θy) 0.0 sin(θy); 0.0 1.0 0.0; -sin(θy) 0.0 cos(θy)]
    Rz = [cos(θz) -sin(θz) 0.0; sin(θz) cos(θz) 0.0; 0.0 0.0 1.0]
    se3_rotation_matrix(Rz * Ry * Rx)
end

"""
    se3_affine_matrix(R::AbstractMatrix, d::AbstractVector)

Returns the 4x4 homogeneous affine matrix with rotation `R` and translation `d`.
"""
function se3_affine_matrix(R::AbstractMatrix, d::AbstractVector)
    @assert length(d) == 3 "Translation vector must be 3-dimensional"
    @assert size(R) == (3, 3) "Rotation matrix must be 3x3"
    [R -d; 0 0 0 1]
end

"""
    se3_affine_matrix(d::AbstractVector; θx::Real=0.0, θy::Real=0.0, θz::Real=0.0)

Returns the 4x4 homogeneous affine matrix using Euler angles and translation `d`.
"""
function se3_affine_matrix(d::AbstractVector; θx::Real=0.0, θy::Real=0.0, θz::Real=0.0)
    R_hom = se3_rotation_matrix(; θx=θx, θy=θy, θz=θz)
    R = R_hom[1:3, 1:3]
    se3_affine_matrix(R, d)
end

"""
    build_escort_topology(kind::Symbol, n_agents::Int, target_node::Int, radius::Float64;
                          observers=1:n_agents, D::Int=4, affine::Bool=true) -> EuclideanSheaf{Float64}

Constructs and returns a `EuclideanSheaf` (with stalk dimension `D`, default 4) for an `n_agents`
escort formation around `target_node`.

An escort formation bundles two concerns that are mathematically independent, and this function
keeps them so:

- **Geometry** — where each agent sits. Regardless of `kind`, agent `i` is placed at angle
  `2π(i-1)/n_agents` and distance `radius` in the plane spanned by the first two translation
  coordinates (this is exactly what `build_escort_ring` has always done).
- **Consensus topology** — which pairs of agents are directly wired together by a sheaf edge.
  `kind` selects this graph:

  - `:ring`   — cycle: agent `i` shares an edge with agent `i % n_agents + 1`. A `2`-agent ring
    is a degenerate 2-cycle (it produces two parallel edges between the same pair of agents), but
    the formation is still well-defined and rigid.
  - `:path`   — open chain: agent `i` shares an edge with `i + 1`, for `i in 1:n_agents-1`.
  - `:star`   — hub-and-spoke: agent `1` shares an edge with every agent `i in 2:n_agents`.
  - `:clique` — all-to-all: every pair `i < j` shares an edge.

Every edge constraint has the algebraic form "agent `i` = centre + `d_i`" — a shared centre plus a
per-agent translation offset — so the constraint system is globally realizable for *any connected*
choice of `kind`. Consequently, for all four topologies above, the resulting sheaf's space of
exact global sections is exactly `D`-dimensional, parameterized by the formation centre: the
formation is rigid no matter which agents happen to be directly wired together. This is the
property the hierarchical layered-control architecture depends on (see
`docs/issues/007-nested-layered-systems-design.md`, §3.2).

When `affine=true` (default), stalks use `D-1` homogeneous-affine translation coordinates plus one
homogeneous row (e.g. `D=4` for SE(3): 3D translation), and restriction maps are
`affine_translation_matrix` offsets — this recovers the original SE(3) escort ring for `D=4`. When
`affine=false`, stalks are `D` plain (non-homogeneous) Euclidean coordinates; a purely linear
restriction map cannot represent a translation at all, so every restriction map is the identity (a
pure consensus topology) and `radius` must be `0.0`.

`observers` names which agents (local indices `1:n_agents`) are pinned to `target_node`.
"""
function build_escort_topology(kind::Symbol, n_agents::Int, target_node::Int, radius::Float64;
                               observers=1:n_agents, D::Int=4, affine::Bool=true)
    @argcheck kind in (:ring, :path, :star, :clique) "kind must be one of :ring, :path, :star, :clique (got :$kind)"
    @argcheck affine || radius == 0.0 "Non-affine (linear) stalks cannot represent a nonzero radius offset; set affine=true or radius=0.0"
    @argcheck all(1 .<= o .<= n_agents for o in observers) "observers must be within 1:n_agents"
    min_agents = kind == :ring ? 2 : 1
    @argcheck n_agents >= min_agents "n_agents must be >= $min_agents for kind=:$kind"

    total_nodes = max(n_agents, target_node)
    sheaf = EuclideanSheaf{Float64}(fill(D, total_nodes))

    restriction_matrix = if affine
        trans_dim = D - 1
        # Translation offset for agent i in the world frame, placed in the plane
        # spanned by the first two translation coordinates (zero elsewhere)
        function angular_offset(i)
            angle = (i - 1) * 2π / n_agents
            d = zeros(trans_dim)
            trans_dim >= 1 && (d[1] = cos(angle) * radius)
            trans_dim >= 2 && (d[2] = sin(angle) * radius)
            return d
        end
        i -> affine_translation_matrix(angular_offset(i))
    else
        # Non-affine (plain linear) stalks cannot represent a translation at all,
        # so every restriction map is simply the identity (pure consensus).
        i -> Matrix{Float64}(I, D, D)
    end

    # Consensus edges: which agents are directly wired together, per `kind`.
    consensus_edges = if kind == :ring
        [(i, i % n_agents + 1) for i in 1:n_agents]
    elseif kind == :path
        [(i, i + 1) for i in 1:(n_agents - 1)]
    elseif kind == :star
        [(1, i) for i in 2:n_agents]
    else # :clique
        [(i, j) for i in 1:n_agents for j in (i + 1):n_agents]
    end

    for (i, j) in consensus_edges
        add_sheaf_edge!(sheaf, i, j, restriction_matrix(i), restriction_matrix(j))
    end

    # Pin observers to the target
    F_target = Matrix{Float64}(I, D, D)
    for i in observers
        add_sheaf_edge!(sheaf, i, target_node, restriction_matrix(i), F_target)
    end

    return sheaf
end

"""
    build_escort_ring(n_agents::Int, target_node::Int, radius::Float64; observers=1:n_agents, D::Int=4, affine::Bool=true)

Constructs and returns a `EuclideanSheaf` for an `n_agents` escort *ring* around `target_node` —
a thin wrapper around [`build_escort_topology`](@ref)`(:ring, ...)`. See that docstring for the
full explanation of the geometry/topology split, `D`, and `affine`.
"""
build_escort_ring(n_agents::Int, target_node::Int, radius::Float64; observers=1:n_agents, D::Int=4, affine::Bool=true) =
    build_escort_topology(:ring, n_agents, target_node, radius; observers, D, affine)

"""
    build_escort_clique(n_agents::Int, target_node::Int, radius::Float64; observers=1:n_agents, D::Int=4, affine::Bool=true)

Constructs and returns a `EuclideanSheaf` for an `n_agents` escort *clique* (all-to-all consensus)
around `target_node` — a thin wrapper around [`build_escort_topology`](@ref)`(:clique, ...)`. See
that docstring for the full explanation of the geometry/topology split, `D`, and `affine`.
"""
build_escort_clique(n_agents::Int, target_node::Int, radius::Float64; observers=1:n_agents, D::Int=4, affine::Bool=true) =
    build_escort_topology(:clique, n_agents, target_node, radius; observers, D, affine)

end # module
