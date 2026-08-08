struct HSDScaling{T} <: AbstractScaling{T}
    cscl::FVector{T}   # length n: per-column scaling (constant on each cone block)
    rscl::FVector{T}   # length m: per-row scaling
    bscl::FScalar{T}   # HSD embedding border scale d₂ (Stage 2; 1 = unused)
end

"""
    HSDScaling{T}(n, m)

Construct an identity (trivial) scaling with all entries equal to 1.
"""
function HSDScaling{T}(n::Int, m::Int) where {T}
    cscl = FVector{T}(undef, n)
    rscl = FVector{T}(undef, m)
    bscl = FScalar{T}(undef)

    fill!(cscl, one(T))
    fill!(rscl, one(T))
    fill!(bscl, one(T))

    return HSDScaling{T}(cscl, rscl, bscl)
end

"""
    update!(scaling, c, g)

Stage 2 of the two-stage equilibration of the HSD embedding operator

    [ H   -Bᵀ  -c ]
    [ B    0   -g ]
    [ cᵀ   gᵀ   0 ] .

Stage 1 (Ruiz on the interior `[H, B]`, done by `equilibrate!`) is assumed
complete, so `c`, `g` arrive already interior-scaled. This computes the
closed-form scalar `d₂` that equilibrates the `[c; g]` border row/column to unit
∞-norm, applies it in place (`c .*= d₂`, `g .*= d₂`), and stores it in
`scaling.bscl`. The embedding's corner entry is zero, so the general formula
`d₂ = min(1/‖[c;g]‖∞, 1/√|γ|)` reduces to `d₂ = 1/‖[c;g]‖∞`.

O(n): a couple of ∞-norm passes over `c` and `g`. Reusable across a chain
(the interior scaling never changes; only this border scalar does) — hence
callable from both `init` and `reinit!`.
"""
function update!(scaling::HSDScaling{T}, c::AbstractVector{T}, g::AbstractVector{T}) where {T}
    β = max(norm(c, Inf), norm(g, Inf))

    if iszero(β)
        d = one(T)
    else
        d = inv(β)
    end

    rmul!(c, d)
    rmul!(g, d)

    scaling.bscl[] = d
    return d
end

"""
    unscale!(p, d, y, τ, κ, scaling)

HSD-specific unscaling. Applies the interior inverse scaling to `p`, `d`, `y`
(as `unscale!(p, d, y, scaling)`), and maps the internal homogenizer `τ` and
slack `κ` back through the Stage-2 border scale `d₂`, returning the user-frame
`(τ, κ) = (d₂·τ, κ/d₂)`. The complementarity product `τκ` is preserved.
"""
function unscale!(p::AbstractVector, d::AbstractVector, y::AbstractVector, τ::T, κ::T, scaling::HSDScaling{T}) where {T}
    unscale!(p, d, y, scaling)

    d₂ = scaling.bscl[]
    return d₂ * τ, κ / d₂
end
