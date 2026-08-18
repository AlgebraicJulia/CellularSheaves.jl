# =============================================================================
#
# # Frictional Contact Dynamics
#
# A collection of rigid bodies move with velocities v = (v₁, …, vₘ), vᵢ ∈ ℝ⁶,
# so the system has 6m degrees of freedom. Bodies that contact each other
# experience Coulomb friction, and there are n points of contact. We discretize
# the time interval into steps [tₖ, tₖ₊₁] ⊆ ℝ. Impacts make v jump, so a
# contact transmits no force in the ordinary sense; over one step it transmits
# an impulse. At each step the n contact points have relative velocities and
# impulses
#
#   u = (u₁, …, uₙ), uⱼ ∈ ℝ³
#   r = (r₁, …, rₙ), rⱼ ∈ ℝ³,
#
# each with a normal component along the contact normal and a tangential
# component in the contact plane:
#
#   uⱼ = (unⱼ, utⱼ) ∈ ℝ × ℝ²
#   rⱼ = (rnⱼ, rtⱼ) ∈ ℝ × ℝ².
#
# These are constrained to satisfy the contact equation
#
#   uⱼ + Φⱼ(u) ∈ Kⱼ*,  uⱼ + Φⱼ(u) ⊥ rⱼ ∈ Kⱼ,  (1)
#
# where Kⱼ is a friction cone with coefficient μⱼ
#
#   Kⱼ  = { rⱼ : μⱼ rnⱼ ≥    ‖rtⱼ‖ }
#   Kⱼ* = { uⱼ :    unⱼ ≥ μⱼ ‖utⱼ‖ }
#
# and Φⱼ(u) = (μⱼ‖utⱼ‖, 0) is a function called De Saxcé's correction.
# The momentum balance and the kinematics close the system:
#
#   M v + f = Hᵀ r,   H v + w = u,            (2)
#
# where M ∈ ℝ⁶ᵐˣ⁶ᵐ is the mass matrix, H ∈ ℝ³ⁿˣ⁶ᵐ maps global to relative
# velocities, and f and w collect what is known at the start of the step.
#
# To simulate the system defined by (1) and (2), we solve at each time step a
# sequence of second-order cone programs. Defining P to be the block
# diagonal matrix
#
#        [ P₁       ]            [ 1        ]
#   P  = [    ⋱     ] ,   Pⱼ  =  [    μⱼ    ]
#        [       Pₙ ]            [       μⱼ ]
#
# and re-parametrizing
#
#   u ← P u,   H ← P H,   w ← P(w + φ(s)),   r ← P⁻¹ r,
#
# we solve
#
#   min   1/2 vᵀ M v + fᵀ v
#   v, u
#
#   s.t.  u - H v = w                         (3)
#               u ∈ Lⁿ
#
# and then update sⱼ ← μⱼ‖utⱼ‖. This is Cadoux's fixed point iteration.
#
# # Data
#
# The benchmarks use data from FCLIB. To download it, run
#
#   git clone --filter=blob:none --no-checkout \
#       https://github.com/FrictionalContactLibrary/fclib-library
#   cd fclib-library && git sparse-checkout set Global/siconos && git checkout master
#
# and point FCLIB_DIR at the checkout.
#
# # Reference
#
#  - Acary, Armand, Nguyen & Shpakovych,
#    Second-order cone programming for frictional contact mechanics using
#    interior point algorithm,
#    Optim. Methods Softw. 39(3), 2024. doi:10.1080/10556788.2023.2296438
#
# =============================================================================

using AppleAccelerate
using LinearAlgebra, Printf, SparseArrays
using HDF5

using CellularSheaves.IPM
using CellularSheaves.IPM: CofreeCone, SecondOrderCone, IPMSettings, HSDSettings
using CellularSheaves.BlockSparseArrays: block, blocksparse, colrange, rowrange

include("../utils.jl")

const FCLIB_DIR = get(ENV, "FCLIB_DIR", joinpath(@__DIR__, "..", "..", "data", "fclib-library"))

# -----------------------------------------------------------------------------
# FCLIB Global reader
# -----------------------------------------------------------------------------

# An fclib matrix is {m, n, nz, p, i, x} with nz ≥ 0 triplet, -1 csc, -2 csr,
# and 0-based indices.
function fclib_matrix(g)
    m = only(read(g, "m")); n = only(read(g, "n")); nz = only(read(g, "nz"))
    p = read(g, "p"); i = read(g, "i"); x = read(g, "x")

    if nz >= 0
        return sparse(i .+ 1, p .+ 1, x, m, n)
    elseif nz == -1
        return SparseMatrixCSC(m, n, p .+ 1, i .+ 1, x)
    elseif nz == -2
        return sparse(transpose(SparseMatrixCSC(n, m, p .+ 1, i .+ 1, x)))
    end

    return error("unknown fclib sparse type nz = $nz")
end

# (M, H, w, f, mu) in the notation of (1). The stored H is m × 3n, so it is the
# transpose of the reference's H.
function read_fclib_global(path::AbstractString)
    return h5open(path, "r") do fid
        g = fid["fclib_global"]
        M = fclib_matrix(g["M"])
        H = sparse(permutedims(fclib_matrix(g["H"])))
        v = g["vectors"]
        (M, H, read(v, "w"), read(v, "f"), read(v, "mu"))
    end
end

# -----------------------------------------------------------------------------
# Stalk partition
# -----------------------------------------------------------------------------

# The bodies are the connected components of M. For the rigid body families this
# recovers the 6 × 6 blocks exactly; a deformable instance would come back as a
# single component, which is correct but useless.
function bodypartition(M::SparseMatrixCSC)
    n = size(M, 2)
    parent = collect(1:n)

    function find(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end

        return i
    end

    for j in 1:n, p in nzrange(M, j)
        a = find(M.rowval[p]); b = find(j)
        a != b && (parent[a] = b)
    end

    label = Dict{Int, Int}()
    cols = Vector{Int}[]

    for j in 1:n
        r = find(j)
        v = get!(label, r) do
            push!(cols, Int[])
            length(cols)
        end
        push!(cols[v], j)
    end

    return cols
end

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------

# Stalks 1..nb are the body velocities (cofree), nb+c is the contact velocity
# u_c (second order cone). One 3-row block per contact,
#
#   u_c - Pμ H_c v = Pμ w_c
#
# touching its own cone stalk and the one or two bodies the contact joins.
function build_contact(M, H, w, f, mu)
    nc = length(mu)
    bodycols = bodypartition(M)
    nb = length(bodycols)

    stalk = zeros(Int, size(M, 2))

    for (v, cols) in pairs(bodycols), j in cols
        stalk[j] = v
    end

    row_ids = Int[]
    col_ids = Int[]
    blks = Matrix{Float64}[]

    Ht = sparse(transpose(H))   # column access to H's rows
    I3 = Matrix{Float64}(I, 3, 3)

    for c in 1:nc
        rows = 3c - 2:3c
        P = Diagonal([1.0, mu[c], mu[c]])

        # bodies touched by this contact
        touched = Int[]

        for i in rows, p in nzrange(Ht, i)
            v = stalk[Ht.rowval[p]]
            v in touched || push!(touched, v)
        end

        for v in touched
            push!(row_ids, c); push!(col_ids, v)
            push!(blks, -(P * Matrix(H[rows, bodycols[v]])))
        end

        push!(row_ids, c); push!(col_ids, nb + c); push!(blks, copy(I3))
    end

    # dimensions given explicitly: a body in no contact would otherwise be
    # assigned a zero-dimensional stalk
    B = blocksparse(row_ids, col_ids, blks, fill(3, nc),
                    vcat(length.(bodycols), fill(3, nc)))

    # Q: the mass matrix on each body stalk, zero on the cone stalks
    Q = IPM.allocblockdiag(B)
    fill!(Q, 0)

    for v in 1:nb
        copyto!(block(Q, v, v, v), Matrix(M[bodycols[v], bodycols[v]]))
    end

    # (1) is ½ vᵀMv + fᵀv and IPMProblem is ½ pᵀQp - cᵀp, so c = -f.
    c = zeros(size(B, 2))

    for v in 1:nb
        c[colrange(B, v)] .= -f[bodycols[v]]
    end

    g = zeros(3nc)

    for c_ in 1:nc
        g[rowrange(B, c_)] .= Diagonal([1.0, mu[c_], mu[c_]]) * w[3c_ - 2:3c_]
    end

    cones = vcat(IPM.AbstractCone[CofreeCone() for _ in 1:nb],
                 IPM.AbstractCone[SecondOrderCone() for _ in 1:nc])

    return IPMProblem(Q, B, c, g, cones)
end

function build_contact_jump(M, H, w, f, mu; optimizer)
    nv = size(M, 2)
    nc = length(mu)
    P = Diagonal(vec([ones(nc) mu mu]'))

    Hs = P * H
    ws = P * w

    model = Model(optimizer)
    @variable(model, v[1:nv])
    @variable(model, u[1:3nc])

    @constraint(model, u .- Hs * v .== ws)

    for i in 1:nc
        @constraint(model, u[3i - 2:3i] in JuMP.SecondOrderCone())
    end

    @objective(model, Min, 0.5 * dot(v, M * v) + dot(f, v))
    return model
end

# -----------------------------------------------------------------------------
# Reduced (Delassus) form — Mosek's preferred formulation
# -----------------------------------------------------------------------------

# Eliminating v from (1) leaves a QP in the impulses, min ½ r̂ᵀ Ŵ r̂ + q̂ᵀ r̂ over
# r̂ ∈ L³ⁿ, with the Delassus operator Ŵ = Ĥ M⁻¹ Ĥᵀ. Ŵ is only PSD in exact
# arithmetic; formed as the product it comes back numerically indefinite (665
# negative eigenvalues, λmin ≈ -0.018, on Chute-1885), which MOI cannot Cholesky
# into a cone. So we never form the product: with M_b = L_b L_bᵀ per body,
# G = Ĥ L⁻ᵀ gives Ŵ = G Gᵀ, PSD by construction, and the quadratic goes in
# factored form s ≥ ½ ‖Gᵀ r̂‖². Mosek is faster here than on (1); Clarabel is not.

pmu(mu) = Diagonal(vec([ones(length(mu)) mu mu]'))

# Ŵ = G Gᵀ with G = Ĥ L⁻ᵀ, assembled body by body so nothing densifies and the
# product Ĥ M⁻¹ Ĥᵀ is never formed. Returns (G, q̂).
function delassus_gram(M, H, w, f, mu)
    nv = size(M, 2)
    nc = length(mu)
    P = pmu(mu)
    Hh = P * H
    wh = P * w

    Is = Int[]; Js = Int[]; Vs = Float64[]
    z = zeros(nv)

    for cols in bodypartition(M)
        L = cholesky(Symmetric(Matrix(M[cols, cols]))).L

        Hb = Hh[:, cols]
        rows = unique(Hb.rowval)                    # only the contacts on this body
        isempty(rows) || let
            X = Matrix(Hb[rows, :]) / UpperTriangular(Matrix(L)')
            for (jj, j) in pairs(cols), (ii, i) in pairs(rows)
                iszero(X[ii, jj]) && continue
                push!(Is, i); push!(Js, j); push!(Vs, X[ii, jj])
            end
        end

        z[cols] .= LowerTriangular(Matrix(L)) \ f[cols]
    end

    G = sparse(Is, Js, Vs, 3nc, nv)
    return G, wh - G * z
end

# min ½ ‖Gᵀ r‖² + q̂ᵀ r  s.t.  r_c ∈ L³, quadratic in factored (rotated cone)
# form so the frontend never needs a Cholesky of the singular Ŵ.
function build_reduced_jump(G::SparseMatrixCSC, q::AbstractVector, nc::Integer; optimizer)
    nv = size(G, 2)

    model = Model(optimizer)
    @variable(model, r[1:3nc])
    @variable(model, t[1:nv])
    @variable(model, s)

    @constraint(model, t .== transpose(G) * r)
    @constraint(model, [s; 1.0; t] in JuMP.RotatedSecondOrderCone())

    for c in 1:nc
        @constraint(model, r[3c - 2:3c] in JuMP.SecondOrderCone())
    end

    @objective(model, Min, s + dot(q, r))
    return model
end

# -----------------------------------------------------------------------------
# Benchmark
# -----------------------------------------------------------------------------

# Three instances from the two Global families the reference tabulates, spanning
# its size range: BoxStacks tops out at 557 contacts, Chute at ~3200.
const CASES = [
    ("BoxStacks-552", "Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-552-8.hdf5"),
    ("Chute-1885",    "Global/siconos/Chute/Chute-ndof-9216-nc-1885-2782.hdf5"),
    ("Chute-3248",    "Global/siconos/Chute/Chute-ndof-12672-nc-3248-3892.hdf5"),
]

function benchmark(; tol = 1.0e-8, clarabel = false, mosek = false, runhsd = true, cases = CASES)
    @printf("%-11s %6s %6s %6s %9s %9s %9s %9s\n",
            "problem", "bodies", "cts", "stalks", "IPM", "HSD", "Clarabel", "Mosek")

    for (name, rel) in cases
        path = joinpath(FCLIB_DIR, rel)

        if !isfile(path)
            @printf("%-11s  missing: %s\n", name, rel)
            continue
        end

        M, H, w, f, mu = read_fclib_global(path)
        prob = build_contact(M, H, w, f, mu)

        iset = IPMSettings{Float64}(; feas_tol = tol, gap_tol = tol)
        hset = HSDSettings{Float64}(; feas_tol = tol, gap_tol = tol)

        ipm = measure_ipm(prob, iset)
        hsd = runhsd ?
            measure_ipm(prob, hset) :
            (t = NaN, status = "-", obj = NaN)

        # Clarabel on the velocity form (1); Mosek on the reduced (Gram) form.
        cla = clarabel ?
            measure_jump(() -> build_contact_jump(M, H, w, f, mu;
                                                   optimizer = clarabel_opt(; tol))) :
            (t = NaN, status = "-", obj = NaN)

        mos = if mosek
            G, q = delassus_gram(M, H, w, f, mu)
            measure_jump(() -> build_reduced_jump(G, q, length(mu);
                                                   optimizer = mosek_opt(; tol)))
        else
            (t = NaN, status = "-", obj = NaN)
        end

        @printf("%-11s %6d %6d %6d %s %s %s %s\n",
                name, length(prob.K) - length(mu), length(mu), length(prob.K),
                fmt_time(ipm; width = 9), fmt_time(hsd; width = 9),
                fmt_time(cla; width = 9), fmt_time(mos; width = 9))
    end

    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = parse_args()
    benchmark(; tol = args.tol, clarabel = args.clarabel, mosek = args.mosek, runhsd = args.hsd)
end

# =============================================================================
# Sample run: 2026-08-15 (--clarabel --mosek --no-hsd), tol = 1e-8.
# IPM / Clarabel solve the velocity form (1); Mosek the reduced (Gram) form.
#
#   problem       bodies    cts stalks       IPM   HSD  Clarabel     Mosek
#   BoxStacks-552    450    552   1002    26.1ms     —    28.1ms    17.6ms
#   Chute-1885      6144   1885   8029   107.2ms     —   159.4ms   112.5ms
#   Chute-3248      8448   3248  11696   198.8ms     —   400.6ms   235.7ms
# =============================================================================
