"""
    IPMProblem{T, I, C}

A primal-dual pair of conic quadratic programs.

```math
\\begin{aligned}
\\text{(P)} \\quad & \\min_p \\quad \\tfrac{1}{2} p^\\top Q p + c^\\top p \\\\
& \\text{s.t.} \\quad B p = g, \\quad p \\in K
\\end{aligned}
\\qquad
\\begin{aligned}
\\text{(D)} \\quad & \\max_{p, y} \\quad -\\tfrac{1}{2} p^\\top Q p + g^\\top y \\\\
& \\text{s.t.} \\quad Q p + c - B^\\top y = d, \\quad d \\in K^*
\\end{aligned}
```

Solve using [`solve`](@ref).
"""

struct IPMProblem{T, I, C <: AbstractCone}
    Q::BlockSparseMatrix{T, I}
    B::BlockSparseMatrix{T, I}
    c::FVector{T}
    g::FVector{T}
    K::FVector{C}

    function IPMProblem{T, I, C}(Q::BlockSparseMatrix, B::BlockSparseMatrix, c::FVector, g::FVector, K::FVector) where {T, I, C <: AbstractCone}
        @assert nrows(B) == length(g)
        @assert ncols(B) == ncols(Q) == length(c)
        @assert nvtxs(B) == nvtxs(Q) == length(K)

        for v in vtxs(B)
            @assert ncols(B, v) == ncols(Q, v)
        end

        return new{T, I, C}(Q, B, c, g, K)
    end
end

"""
    IPMProblem(Q, B, c, g, K)

Construct an [`IPMProblem`](@ref).
"""
function IPMProblem(Q::BlockSparseMatrix{T, I}, B::BlockSparseMatrix{T, I}, c::AbstractVector{T}, g::AbstractVector{T}, K::AbstractVector{C}) where {T, I, C <: AbstractCone}
    return IPMProblem{T, I, C}(Q, B, c, g, K)
end

function IPMProblem{T, I, C}(Q::BlockSparseMatrix, B::BlockSparseMatrix, c::AbstractVector, g::AbstractVector, K::AbstractVector) where {T, I, C <: AbstractCone}
    if !(c isa FVector{T})
        c = FVector{T}(c)
    end

    if !(g isa FVector{T})
        g = FVector{T}(g)
    end

    if !(K isa FVector{C})
        K = FVector{C}(K)
    end

    return IPMProblem{T, I, C}(Q, B, c, g, K)
end

function showproblem(io::IO, prob::IPMProblem; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "(P)  min   ½ pᵀ Q p + cᵀ p")
    println(io, pad, "     s.t.  B p = g,  p ∈ K")
    println(io)
    println(io, pad, "(D)  max  -½ pᵀ Q p + gᵀ y")
    println(io, pad, "     s.t.  Q p + c - Bᵀ y ∈ K*")
    println(io)
    @printf(io, "%sQ: %6d × %-6d  c: %d\n", pad, size(prob.Q, 1), size(prob.Q, 2), length(prob.c))
    @printf(io, "%sB: %6d × %-6d  g: %d\n", pad, size(prob.B, 1), size(prob.B, 2), length(prob.g))
    println(io)

    nnoc = npos = nsoc = nsdp = nexp = npow = 0

    for k in prob.K
        if k isa CofreeCone
            nnoc += 1
        elseif k isa PositiveCone
            npos += 1
        elseif k isa SecondOrderCone
            nsoc += 1
        elseif k isa SemidefiniteCone
            nsdp += 1
        elseif k isa ExponentialCone
            nexp += 1
        elseif k isa PowerCone
            npow += 1
        end
    end

    println(io, pad, "K:")
    nnoc > 0 && @printf(io, "%s    CofreeCone:       %d\n", pad, nnoc)
    npos > 0 && @printf(io, "%s    PositiveCone:     %d\n", pad, npos)
    nsoc > 0 && @printf(io, "%s    SecondOrderCone:  %d\n", pad, nsoc)
    nsdp > 0 && @printf(io, "%s    SemidefiniteCone: %d\n", pad, nsdp)
    nexp > 0 && @printf(io, "%s    ExponentialCone:  %d\n", pad, nexp)
    npow > 0 && @printf(io, "%s    PowerCone:        %d\n", pad, npow)
    return
end

function Base.show(io::IO, ::MIME"text/plain", prob::T) where {T <: IPMProblem}
    println(io, T, ":")
    return showproblem(io, prob; indent=2)
end

@kwdef struct IPMSettings{T}
    verbose::Bool = false
    step_frac::T = 0.99
    feas_tol::T = 1e-8
    gap_tol::T = 1e-8
    itmax::Int = 100
    near_factor::T = 1000.0
    stall_tol::T = 1e-6
    forcing_ceil::T = 0.3
    forcing_frac::T = 1.0
    refine_itmax::Int = 10
    refine_stall::T = 0.5
    scale_itmax::Int = 10
    kkt::KKTSettings{T} = UzawaSettings{T}()
end

function showsettings(io::IO, set::IPMSettings; indent::Integer=0)
    pad = " "^indent
    @printf(io, "%sverbose:      %8s  feas_tol:     %8.2e\n", pad, set.verbose, set.feas_tol)
    @printf(io, "%sstep_frac:    %8.2e  gap_tol:      %8.2e\n", pad, set.step_frac, set.gap_tol)
    @printf(io, "%sitmax:        %8d  stall_tol:    %8.2e\n", pad, set.itmax, set.stall_tol)
    @printf(io, "%snear_factor:  %8.2e  forcing_ceil: %8.2e\n", pad, set.near_factor, set.forcing_ceil)
    @printf(io, "%sforcing_frac: %8.2e  refine_stall: %8.2e\n", pad, set.forcing_frac, set.refine_stall)
    @printf(io, "%srefine_itmax: %8d  scale_itmax:  %8d\n", pad, set.refine_itmax, set.scale_itmax)
    println(io)
    println(io, pad, "kkt: ", typeof(set.kkt))
    showsettings(io, set.kkt; indent=indent+4)
    return
end

function Base.show(io::IO, ::MIME"text/plain", set::T) where {T <: IPMSettings}
    println(io, T, ":")
    return showsettings(io, set; indent=2)
end
