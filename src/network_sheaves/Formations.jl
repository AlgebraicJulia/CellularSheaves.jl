module Formations

using LinearAlgebra
using ..EuclideanSheaves: EuclideanSheaf, add_sheaf_edge!
using ArgCheck: @argcheck

export se3_translation_matrix, se3_rotation_matrix, se3_affine_matrix, affine_translation_matrix, build_escort_ring, build_escort_clique

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
    build_escort_ring(n_agents::Int, target_node::Int, radius::Float64; observers=1:n_agents, D::Int=4, affine::Bool=true)

Constructs and returns a `EuclideanSheaf` (with stalk dimension `D`, default 4) for an `n_agents` escort ring around `target_node`.
The consensus edges form a topological cycle graph (ring) with translation offsets placed in the plane spanned by the
first two translation coordinates. The target pinning edges are added for the agents specified in `observers`.

When `affine=true` (default), stalks use `D-1` homogeneous-affine translation coordinates plus one homogeneous row
(e.g. `D=4` for SE(3): 3D translation), and restriction maps are `affine_translation_matrix` offsets — this recovers
the original SE(3) escort ring for `D=4`. When `affine=false`, stalks are `D` plain (non-homogeneous) Euclidean
coordinates; a purely linear restriction map cannot represent a translation at all, so every restriction map is the
identity (a pure consensus ring) and `radius` must be `0.0`.
"""
function build_escort_ring(n_agents::Int, target_node::Int, radius::Float64; observers=1:n_agents, D::Int=4, affine::Bool=true)
    @argcheck affine || radius == 0.0 "Non-affine (linear) stalks cannot represent a nonzero ring radius offset; set affine=true or radius=0.0"
    @argcheck all(1 .<= o .<= n_agents for o in observers) "observers must be within 1:n_agents"
    total_nodes = max(n_agents, target_node)
    sheaf = EuclideanSheaf{Float64}(fill(D, total_nodes))

    restriction_matrix = if affine
        trans_dim = D - 1
        # Translation offset for agent i in the world frame, placed in the plane
        # spanned by the first two translation coordinates (zero elsewhere)
        function ring_offset(i)
            angle = (i - 1) * 2π / n_agents
            d = zeros(trans_dim)
            trans_dim >= 1 && (d[1] = cos(angle) * radius)
            trans_dim >= 2 && (d[2] = sin(angle) * radius)
            return d
        end
        i -> affine_translation_matrix(ring_offset(i))
    else
        # Non-affine (plain linear) stalks cannot represent a translation at all,
        # so every restriction map is simply the identity (pure consensus).
        i -> Matrix{Float64}(I, D, D)
    end

    # Create the ring consensus edges
    for i in 1:n_agents
        j = i % n_agents + 1
        add_sheaf_edge!(sheaf, i, j, restriction_matrix(i), restriction_matrix(j))
    end

    # Pin observers to the target
    for i in observers
        F_target = Matrix{Float64}(I, D, D)
        add_sheaf_edge!(sheaf, i, target_node, restriction_matrix(i), F_target)
    end

    return sheaf
end

"""
    build_escort_clique(n_agents::Int, target_node::Int, radius::Float64; observers=1:n_agents)

Constructs and returns a `EuclideanSheaf` (with D=4 stalk dimension) for an `n_agents` escort formation around `target_node`.
The consensus edges form a topological clique (all-to-all) with the appropriate SE(3) translation offsets. 
The target pinning edges are added for the agents specified in `observers`.
"""
function build_escort_clique(n_agents::Int, target_node::Int, radius::Float64; observers=1:n_agents)
    total_nodes = max(n_agents, target_node)
    sheaf = EuclideanSheaf{Float64}(fill(4, total_nodes))
    
    # Create the clique consensus edges (all pairs)
    for i in 1:n_agents
        for j in (i+1):n_agents
            angle_i = (i - 1) * 2π / n_agents
            angle_j = (j - 1) * 2π / n_agents
            
            di = [cos(angle_i), sin(angle_i), 0.0] * radius
            dj = [cos(angle_j), sin(angle_j), 0.0] * radius
            
            Fi = se3_translation_matrix(di)
            Fj = se3_translation_matrix(dj)
            
            add_sheaf_edge!(sheaf, i, j, Fi, Fj)
        end
    end
    
    # Pin observers to the target
    for i in observers
        angle_i = (i - 1) * 2π / n_agents
        di = [cos(angle_i), sin(angle_i), 0.0] * radius
        
        Fi = se3_translation_matrix(di)
        F_target = [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0] # Identity
        
        add_sheaf_edge!(sheaf, i, target_node, Fi, F_target)
    end
    
    return sheaf
end

end # module
