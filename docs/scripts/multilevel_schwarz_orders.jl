# Why the two V-cycle orders of MultilevelSchwarz differ.
#
#   julia --project=docs docs/scripts/multilevel_schwarz_orders.jl
#
# Every correction of either order is a lifted box correction
#     C_{i,ℓ} = P_ℓ R_{i,ℓ}ᵀ S(R_{i,ℓ} A_ℓ R_{i,ℓ}ᵀ) R_{i,ℓ} Q_ℓ
# (identical in both orders by base change). A V-cycle over corrections C₀, C₁, …
# is the noncommutative polynomial p(C₀, C₁, …; A) of B ← B + Cₖ (I − A B). Since
# Rᵢᵀ Xᵢ Aᵢ Yᵢ Rᵢ = (Rᵢᵀ Xᵢ Rᵢ) A (Rᵢᵀ Yᵢ Rᵢ) for Aᵢ = Rᵢ A Rᵢᵀ:
#     multigrid in each box, then Schwarz  = Σᵢ p(C_{i,0}, C_{i,1}, …; A)
#     Schwarz on each level of the hierarchy = p(Σᵢ C_{i,0}, Σᵢ C_{i,1}, …; A)
# The difference is the cross terms C_{i,ℓ} A C_{j,ℓ'} (i ≠ j) of p's products,
# each containing the box coupling Rᵢ A Rⱼᵀ. This script checks the identities,
# shows the difference vanishing linearly with the couplings between boxes, and
# compares the two preconditioners.

using CellularSheaves
using LinearAlgebra
using SparseArrays
using Printf
using Random

function random_stencil_matrix(rng, n)
    L = LinearIndices(n)
    rows, cols, vals = Int[], Int[], Float64[]
    for I in CartesianIndices(n)
        d = 0.1
        for j in 1:length(n), s in (-1, 1)
            c = rand(rng)
            d += c
            k = I[j] + s
            1 <= k <= n[j] || continue
            push!(rows, L[I]); push!(cols, L[Base.setindex(Tuple(I), k, j)...]); push!(vals, -c)
        end
        push!(rows, L[I]); push!(cols, L[I]); push!(vals, d)
    end
    return sparse(rows, cols, vals, prod(n), prod(n))
end

function compose(A, corrections)
    B = zeros(size(A))
    for C in corrections
        B = B + C * (I - A * B)
    end
    return B
end

# The lifted corrections C[i][ℓ + 1] of every box on every level.
function lifted_corrections(A, n, boxes, factors, smoother, coarsest)
    L = length(factors)
    sizes, covers = [n], [boxes]
    P, Q, As = [sparse(1.0I, prod(n), prod(n))], [sparse(1.0I, prod(n), prod(n))], Any[A]
    for r in factors
        E, T = prolongation_matrix(sizes[end], r), transfer_matrix(sizes[end], r)
        push!(P, P[end] * E); push!(Q, T * Q[end]); push!(As, T * As[end] * E)
        push!(covers, [coarse_box(U, r) for U in covers[end]])
        push!(sizes, coarse_points(sizes[end], r))
    end
    return [[begin
                 R = box_restriction(sizes[ℓ + 1], covers[ℓ + 1][i])
                 S = ℓ < L ? smoother : coarsest
                 Matrix(P[ℓ + 1] * R' * S(R * As[ℓ + 1] * R') * R * Q[ℓ + 1])
             end for ℓ in 0:L] for i in eachindex(boxes)]
end

vorder(L) = [0:L; (L - 1):-1:0]
reldiff(X, Y) = norm(X - Y) / max(norm(X), norm(Y))
contraction(B, A) = maximum(abs, eigvals(Matrix(I - B * A)))

rng = Random.MersenneTwister(11)
n, factors = (16, 12), [(2, 2), (2, 2)]
L = length(factors)
A = random_stencil_matrix(rng, n)
overlapping = [(x, y) for x in (1:8, 5:12, 9:16) for y in (1:8, 5:12)]
disjoint = [(x, y) for x in (1:8, 9:16) for y in (1:4, 5:12)]

println("1. Both V-cycles are the same polynomial, of the box corrections or of their sums")
for (name, ω) in (("Jacobi", 1.0), ("Jacobi, ω = 1/4", 0.25))
    sm = M -> ω * Diagonal(1 ./ diag(M))
    co = M -> ω * inv(Matrix(M))
    C = lifted_corrections(A, n, overlapping, factors, sm, co)
    V1 = schwarz_of_multigrid(A, n, overlapping, factors; composition = :multiplicative, smoother = sm, coarsest_solve = co)
    V2 = multigrid_of_schwarz(A, n, overlapping, factors; composition = :multiplicative, smoother = sm, coarsest_solve = co)
    p_each = sum(compose(A, [C[i][ℓ + 1] for ℓ in vorder(L)]) for i in eachindex(C))
    p_sum = compose(A, [sum(C[i][ℓ + 1] for i in eachindex(C)) for ℓ in vorder(L)])
    @printf("   %-16s Σᵢ p(Cᵢ) vs schwarz_of_multigrid %.1e, p(Σᵢ Cᵢ) vs multigrid_of_schwarz %.1e, orders differ by %.3f\n",
        name, reldiff(p_each, V1), reldiff(p_sum, V2), reldiff(V1, V2))
end

println("\n2. Disjoint boxes, couplings between boxes scaled by ε: the difference is O(ε)")
box_of = zeros(Int, prod(n))
for (i, U) in enumerate(disjoint)
    box_of[findnz(box_restriction(n, U))[2]] .= i
end
rows, cols, vals = findnz(A)
inside = box_of[rows] .== box_of[cols]
for ε in (1.0, 0.1, 0.01, 0.0)
    Aε = sparse(rows, cols, ifelse.(inside, vals, ε .* vals), size(A)...)
    V1 = schwarz_of_multigrid(Aε, n, disjoint, factors; composition = :multiplicative)
    V2 = multigrid_of_schwarz(Aε, n, disjoint, factors; composition = :multiplicative)
    @printf("   ε = %-5g difference %.2e\n", ε, reldiff(V1, V2))
end

println("\n3. One level, two sweeps of Schwarz with exact box solves: no multigrid needed")
C = lifted_corrections(A, n, overlapping, Tuple{Int,Int}[], identity, M -> inv(Matrix(M)))
Cs = [C[i][1] for i in eachindex(C)]
each = sum(compose(A, [Ci, Ci]) for Ci in Cs)
@printf("   Σᵢ p(Cᵢ) = Σᵢ Cᵢ (a box solved twice is solved once): %.1e\n", reldiff(each, sum(Cs)))
@printf("   p(Σᵢ Cᵢ) vs Σᵢ p(Cᵢ): %.3f (the second sweep sees the neighbours' corrections)\n",
    reldiff(compose(A, [sum(Cs), sum(Cs)]), each))

println("\n4. As preconditioners: spectral radius of I − B A (overlapping boxes)")
for ω in (1.0, 0.5, 0.25)
    sm = M -> ω * Diagonal(1 ./ diag(M))
    co = M -> ω * inv(Matrix(M))
    ρ = map((:additive, :multiplicative)) do comp
        (contraction(schwarz_of_multigrid(A, n, overlapping, factors; composition = comp, smoother = sm, coarsest_solve = co), A),
         contraction(multigrid_of_schwarz(A, n, overlapping, factors; composition = comp, smoother = sm, coarsest_solve = co), A))
    end
    @printf("   ω = %-4g additive: %.3f (both orders)   V-cycle: in each box %.3f, on the hierarchy %.3f\n",
        ω, ρ[1][1], ρ[2][1], ρ[2][2])
end
