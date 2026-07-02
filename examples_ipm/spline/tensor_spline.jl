######################################################################
# tensor_spline.jl
#
# Shape-constrained TENSOR-PRODUCT Bézier surfaces on a (possibly
# lossy) grid of patches — the corpus's first 2D instance, and the
# first whose sheaf H¹ (= ker δᵀ, the dual-coset dimension) is large,
# structured, and predicted by a formula the solve can confirm.
#
# Base graph: one vertex per ACTIVE cell of an Mx×My grid, stalk
# ℝ^{(n+1)²} (the bidegree-(n,n) Bernstein control net, x-index
# fastest). Continuity edges between adjacent active cells match
# transversal boundary jets:
#     vertical shared line:   kron(I, jetₓ)  (x-jets, orders 0..k_e)
#     horizontal shared line: kron(jet_y, I)
# Per-edge smoothness k_e makes the grid LOSSY with real semantics:
#     k_e = k   uniform C^k surface           (CAD/spline standard)
#     k_e < k   crease / ridge line           (terrain, feature lines)
#     k_e = -1  fault: edge dropped           (geological faults)
#     cell off  polyomino domain / holes      (irregular domains)
#
# H¹ bookkeeping (verified in tensor_spline_oracle.py): rank deficiency
# of B decomposes per cycle but NOT topologically —
#     corner 4-cycle (4 active cells, 4 live segments):
#         (min(kx_lo,kx_hi)+1) · (min(ky_lo,ky_hi)+1)     [(k+1)² uniform]
#     unit hole (8-ring), uniform order k:
#         max(2k−n+1, 0)²      — transport around the ring must pass
#         THROUGH patches, and left/right jet spans of a degree-n patch
#         overlap in only (2k−n+1)₊ indices per coordinate.
# So in the practical regime n > 2k holes are free; trees/faulted grids
# have H¹ = 0; the full grid is the maximal (stress-test) case. The
# redundant rows are CONSISTENT (g exact); per duality.md §7 the dual y
# is a coset mod ker δᵀ and the Uzawa CG pins its basepoint (min-norm).
# The 2×2 default demo IS the probe for that solver path: run it first.
#
# Objective: least-squares data fit + λr·thin-plate roughness
# ∫(f_xx² + 2f_xy² + f_yy²), all Kronecker products of the univariate
# Grams. mode = :intensity swaps the data term for a spatial-Poisson
# log-likelihood with ExponentialCone leaves, reusing the mle_spline.jl
# construction verbatim (φ = kron of Bernstein evals) — the
# Alizadeh–Papp multi-dimensional arrival-rate application.
#
# Shape: :nonneg (orthant on the stalk — sufficient, tightens under
# degree elevation; NO exact 2D variant exists, see tensor_spline.md),
# :monotone_x / :monotone_y (difference leaves kron(I,Δ) / kron(Δ,I)),
# :free. Convexity (PSD-Hessian SOC leaves) documented, deferred.
#
# NOTE: written against the CellularSheaves.IPM PR-67 API; not executed
# here (no Julia toolchain). tensor_spline_oracle.py is the numerical
# ground truth for every value quoted in tensor_spline.md.
######################################################################

include("nonnegative_spline_lp.jl")   # univariate Bernstein/jet numerics

if !@isdefined(rand_poisson)
    """Poisson(λ) via cumulative exponentials (log-space Knuth)."""
    function rand_poisson(rng, λ::Float64)
        k = -1; t = 0.0
        while t < λ
            t += -log(rand(rng)); k += 1
        end
        return k
    end
end

# ---- tensor numerics -------------------------------------------------

"""Un-normalized ∫₀¹ (f⁽ᵃ⁾)² Bernstein Gram, (n+1)×(n+1)."""
deriv_gram(n::Int, a::Int) = a == 0 ? bernstein_gram(n) :
    (falling(n, a)^2) .* (diffop(n, a)' * bernstein_gram(n - a) * diffop(n, a))

"""Trace-normalized thin-plate Gram ∫∫ f_xx² + 2f_xy² + f_yy² on an
hx×hy cell, in vec-layout with the x-index fastest: kron(G_y, G_x)."""
function thinplate_gram(n::Int, hx::Float64, hy::Float64)
    G0, G1, G2 = deriv_gram(n, 0), deriv_gram(n, 1), deriv_gram(n, 2)
    G = hx * hy .* (kron(G0, G2) ./ hx^4 .+ 2 .* kron(G1, G1) ./ (hx^2 * hy^2) .+
                    kron(G2, G0) ./ hy^4)
    G ./= tr(G)
    return G
end

"""Ground truth: two bumps (:nonneg/:free/intensity) or a monotone ramp."""
truth2(shape::Symbol, x, y) =
    shape === :monotone_x ? (0.5 + 0.5 * tanh(4 * (x - 0.5))) * (0.8 + 0.2 * y) :
    shape === :monotone_y ? (0.5 + 0.5 * tanh(4 * (y - 0.5))) * (0.8 + 0.2 * x) :
    exp(-((x - 0.3)^2 + (y - 0.35)^2) / (2 * 0.12^2)) +
        0.6 * exp(-((x - 0.7)^2 + (y - 0.75)^2) / (2 * 0.10^2))

"""Cell index and local coordinate of t ∈ [0,1] on a uniform M-partition."""
locate01(t, M) = (i = clamp(floor(Int, t * M) + 1, 1, M);
                  (i, clamp(t * M - (i - 1), 0.0, 1.0)))

# ---- instance --------------------------------------------------------

struct TensorSplineInstance
    Mx::Int; My::Int; n::Int; k::Int
    mode::Symbol                              # :regress | :intensity
    shape::Symbol                             # :nonneg | :monotone_x | :monotone_y | :free
    active::BitMatrix                         # Mx×My cell mask
    kx::Matrix{Int}                           # (Mx−1)×My vertical lines; −1 = fault
    ky::Matrix{Int}                           # Mx×(My−1) horizontal lines
    obs::Vector{NTuple{2, Float64}}           # :intensity events
    w::Vector{Float64}
    Qp::Dict{Tuple{Int, Int}, Matrix{Float64}}  # :regress data Grams
    cp::Dict{Tuple{Int, Int}, Vector{Float64}}
    λr::Float64; ε::Float64
end

function sample_active(rng, active, Mx, My)
    while true
        x, y = rand(rng), rand(rng)
        a, _ = locate01(x, Mx); b, _ = locate01(y, My)
        active[a, b] && return (x, y)
    end
end

"""
    generate_tensor_instance(; Mx, My, n, k, mode, shape,
                               holes, faults, creases, N, σ, λr, ε, seed)

holes   :: Vector{(a,b)}           — deactivate cells (polyomino domain)
faults  :: Vector{(a,b,dir)}       — drop the edge from cell (a,b), dir ∈ (:E,:N)
creases :: Vector{(a,b,dir,ke)}    — lower that edge to C^{ke}
:regress accumulates per-cell data-fit Grams (as the 1D generator does);
:intensity samples an inhomogeneous Poisson process on the active region
with E[#events] ≈ N by thinning.
"""
function generate_tensor_instance(;
        Mx::Int = 2, My::Int = 2, n::Int = 4, k::Int = 2,
        mode::Symbol = :regress, shape::Symbol = :nonneg,
        holes = Tuple{Int, Int}[], faults = Tuple{Int, Int, Symbol}[],
        creases = Tuple{Int, Int, Symbol, Int}[],
        N::Int = 100 * Mx * My, σ::Float64 = 0.05,
        λr::Float64 = 1e-4, ε::Float64 = 1e-8, seed::Int = 1)
    @assert n ≥ max(k + 1, 2) && mode in (:regress, :intensity)
    @assert shape in (:nonneg, :monotone_x, :monotone_y, :free)
    rng = MersenneTwister(seed)
    ns = n + 1; d = ns * ns

    active = trues(Mx, My)
    for (a, b) in holes; active[a, b] = false; end
    kx = fill(k, Mx - 1, My); ky = fill(k, Mx, My - 1)
    setedge!(a, b, dir, v) = dir === :E ? (kx[a, b] = v) : (ky[a, b] = v)
    for (a, b, dir) in faults;       setedge!(a, b, dir, -1); end
    for (a, b, dir, ke) in creases;  (@assert 0 ≤ ke ≤ k); setedge!(a, b, dir, ke); end

    obs = NTuple{2, Float64}[]; w = Float64[]
    Qp = Dict{Tuple{Int, Int}, Matrix{Float64}}()
    cp = Dict{Tuple{Int, Int}, Vector{Float64}}()

    if mode === :regress
        for b in 1:My, a in 1:Mx
            active[a, b] || continue
            Qp[(a, b)] = zeros(d, d); cp[(a, b)] = zeros(d)
        end
        for _ in 1:N
            (x, y) = sample_active(rng, active, Mx, My)
            z = truth2(shape, x, y) + σ * randn(rng)
            a, u = locate01(x, Mx); b, v = locate01(y, My)
            φ = kron(bernstein_eval(n, v), bernstein_eval(n, u))
            Qp[(a, b)] .+= 2.0 .* (φ * φ')
            cp[(a, b)] .-= 2.0 .* z .* φ
        end
    else
        grid = 0.0025:0.005:0.9975                       # active-region stats
        vals = Float64[]
        for x in grid, y in grid
            a, _ = locate01(x, Mx); b, _ = locate01(y, My)
            active[a, b] && push!(vals, truth2(:free, x, y))
        end
        area = count(active) / (Mx * My)
        smax = maximum(vals); sint = (sum(vals) / length(vals)) * area
        λmax = N * smax / sint
        for _ in 1:rand_poisson(rng, λmax * area)        # thinning
            (x, y) = sample_active(rng, active, Mx, My)
            rand(rng) * λmax ≤ N * truth2(:free, x, y) / sint && push!(obs, (x, y))
        end
        w = ones(length(obs))
    end
    return TensorSplineInstance(Mx, My, n, k, mode, shape, active, kx, ky,
                                obs, w, Qp, cp, λr, ε)
end

# ---- per-patch objective (shared by the native and JuMP builders) ----

function patch_QC(inst::TensorSplineInstance)
    ns = inst.n + 1; d = ns * ns
    hx, hy = 1.0 / inst.Mx, 1.0 / inst.My
    # normalized intensity: λ = s·f, s = Σw (see mle_spline.md §3.1)
    sc = inst.mode === :intensity ? sum(inst.w) : 1.0
    Qbase = sc^2 .* (inst.λr .* thinplate_gram(inst.n, hx, hy) .+
                     inst.ε .* Matrix{Float64}(I, d, d))
    Q = Dict{Tuple{Int, Int}, Matrix{Float64}}()
    C = Dict{Tuple{Int, Int}, Vector{Float64}}()
    for b in 1:inst.My, a in 1:inst.Mx
        inst.active[a, b] || continue
        Qi = copy(Qbase); ci = zeros(d)
        if inst.mode === :regress
            Qi .+= inst.Qp[(a, b)]; ci .+= inst.cp[(a, b)]
        else
            ci .+= sc * hx * hy / d               # the compensator s·∫f (linear)
        end
        Q[(a, b)] = Qi; C[(a, b)] = ci
    end
    return Q, C
end

# ---- H¹ predictor ----------------------------------------------------

"""
Predicted dim ker δᵀ (verified in tensor_spline_oracle.py):
Σ over interior corner points with 4 active cells and 4 live segments of
(min(kx_lo, kx_hi)+1)·(min(ky_lo, ky_hi)+1), plus max(2k−n+1,0)² per unit
hole whose 8-ring is active with uniform-order edges.
"""
function h1_predict(inst::TensorSplineInstance)
    tot = 0
    for b in 1:(inst.My - 1), a in 1:(inst.Mx - 1)
        (inst.active[a, b] && inst.active[a + 1, b] &&
         inst.active[a, b + 1] && inst.active[a + 1, b + 1]) || continue
        kxlo, kxhi = inst.kx[a, b], inst.kx[a, b + 1]
        kylo, kyhi = inst.ky[a, b], inst.ky[a + 1, b]
        min(kxlo, kxhi, kylo, kyhi) ≥ 0 || continue
        tot += (min(kxlo, kxhi) + 1) * (min(kylo, kyhi) + 1)
    end
    for b in 2:(inst.My - 1), a in 2:(inst.Mx - 1)     # unit holes
        inst.active[a, b] && continue
        ring = all(inst.active[a + da, b + db]
                   for da in -1:1, db in -1:1 if !(da == 0 && db == 0))
        ring && (tot += max(2 * inst.k - inst.n + 1, 0)^2)
    end
    return tot
end

# ---- native model ----------------------------------------------------

"""
    build_tensor_problem(inst) :: (IPMProblem, info)

Vertex order: active patches (a fastest, then b); monotone difference
leaves (same order); intensity leaves (observation order). Edge order:
shape edges, observation edges, then continuity (per cell: E then N).
"""
function build_tensor_problem(inst::TensorSplineInstance)
    Mx, My, n, k = inst.Mx, inst.My, inst.n, inst.k
    ns = n + 1; d = ns * ns
    Ins = Matrix{Float64}(I, ns, ns)
    Rj(ke) = jet_operator(n, 1, ke, 1.0, :right)
    Lj(ke) = jet_operator(n, 1, ke, 1.0, :left)

    row_ids = Int[]; col_ids = Int[]; blocks = Matrix{Float64}[]
    col_cone = AbstractCone[]; col_dim = Int[]
    rhs = Dict{Int, Vector{Float64}}()
    new_vertex!(dim, cone) = (push!(col_dim, dim); push!(col_cone, cone); length(col_dim))
    ec = Ref(0); new_edge!() = (ec[] += 1; ec[])
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v); push!(blocks, Matrix{Float64}(mat)))

    # patch vertices
    patch = Dict{Tuple{Int, Int}, Int}()
    for b in 1:My, a in 1:Mx
        inst.active[a, b] || continue
        cone = inst.shape === :nonneg ? PositiveCone() : CofreeCone()
        patch[(a, b)] = new_vertex!(d, cone)
    end

    # monotone difference-leaves: kron(I,Δ)P − leaf = 0, leaf ≥ 0
    if inst.shape === :monotone_x || inst.shape === :monotone_y
        D1 = diffop(n, 1)
        D = inst.shape === :monotone_x ? kron(Ins, D1) : kron(D1, Ins)
        m = size(D, 1)
        for b in 1:My, a in 1:Mx
            haskey(patch, (a, b)) || continue
            lv = new_vertex!(m, PositiveCone()); e = new_edge!()
            place!(e, patch[(a, b)], D)
            place!(e, lv, -Matrix{Float64}(I, m, m))
        end
    end

    # intensity observation leaves (u,y,t) ∈ ExponentialCone (mle_spline.jl)
    leaf = Int[]; obsmap = Tuple{Int, Vector{Float64}}[]
    if inst.mode === :intensity
        for (x, y) in inst.obs
            a, u = locate01(x, Mx); b, v = locate01(y, My)
            φ = kron(bernstein_eval(n, v), bernstein_eval(n, u))
            pv = patch[(a, b)]
            push!(obsmap, (pv, φ))
            lv = new_vertex!(3, ExponentialCone()); push!(leaf, lv)
            e = new_edge!()
            place!(e, pv, vcat(φ', zeros(1, d)))       # φᵀP − u = 0
            place!(e, lv, [-1.0 0.0 0.0; 0.0 1.0 0.0]) # …and y = 1
            rhs[e] = [0.0, 1.0]
        end
    end

    # continuity edges (per cell: E then N)
    for b in 1:My, a in 1:Mx
        haskey(patch, (a, b)) || continue
        if a < Mx && haskey(patch, (a + 1, b)) && inst.kx[a, b] ≥ 0
            ke = inst.kx[a, b]; e = new_edge!()
            place!(e, patch[(a, b)], kron(Ins, Rj(ke)))
            place!(e, patch[(a + 1, b)], -kron(Ins, Lj(ke)))
        end
        if b < My && haskey(patch, (a, b + 1)) && inst.ky[a, b] ≥ 0
            ke = inst.ky[a, b]; e = new_edge!()
            place!(e, patch[(a, b)], kron(Rj(ke), Ins))
            place!(e, patch[(a, b + 1)], -kron(Lj(ke), Ins))
        end
    end

    # assemble
    B = blocksparse(row_ids, col_ids, blocks)
    nvtx = length(col_dim)
    @assert nvtxs(B) == nvtx "column blocks not contiguous 1..nvtx"
    g = zeros(size(B, 1)); for (e, v) in rhs; g[rowrange(B, e)] .= v; end

    Qd, Cd = patch_QC(inst)
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for ((a, b), v) in patch
        block(Q, v, v, v) .= Qd[(a, b)]
        c[colrange(B, v)] .= Cd[(a, b)]
    end
    for (j, lv) in enumerate(leaf)
        c[colrange(B, lv)[3]] = -inst.w[j]
    end

    cones = AbstractCone[col_cone[v] for v in 1:nvtx]
    prob = IPMProblem(Q, B, c, g, cones)
    info = (patchcols = Dict(cell => colrange(B, v) for (cell, v) in patch),
            leafcols = [colrange(B, lv) for lv in leaf],
            obsmap = obsmap, nvtx = nvtx, h1 = h1_predict(inst))
    return prob, info
end

tensor_settings() = IPMSettings{Float64}(
    kkt = UzawaSettings{Float64}(raug = 1e5), feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)

"""Mean-square error of the fitted surface against the truth on the active region."""
function fitted_mse(inst, res, info; ngrid = 60)
    n = inst.n; ns = n + 1
    tot = 0.0; cnt = 0
    for gx in 1:ngrid, gy in 1:ngrid
        x, y = (gx - 0.5) / ngrid, (gy - 0.5) / ngrid
        a, u = locate01(x, inst.Mx); b, v = locate01(y, inst.My)
        inst.active[a, b] || continue
        P = reshape(res.p[info.patchcols[(a, b)]], ns, ns)
        f = bernstein_eval(n, u)' * P * bernstein_eval(n, v)
        tot += (f - truth2(inst.shape, x, y))^2; cnt += 1
    end
    return tot / cnt
end

"""
Build+solve demo. THE DEFAULT (2×2, n=4, k=2) IS THE RANK-DEFICIENCY
PROBE: B has h1 = 9 consistent redundant rows (one corner). Run this
before anything larger; watch the Uzawa CG for stagnation.
"""
function tensor_demo(; kwargs...)
    inst = generate_tensor_instance(; kwargs...)
    prob, info = build_tensor_problem(inst)
    solve(prob, tensor_settings())                     # JIT warmup
    t = @elapsed (res = solve(prob, tensor_settings()))
    pres = isempty(res.history) ? NaN : res.history.pres[end]
    dres = isempty(res.history) ? NaN : res.history.dres[end]
    mincoef = minimum(minimum(res.p[r]) for r in values(info.patchcols))
    @printf("2D %-9s %-10s  %dx%d n=%d k=%d  vtx=%d  H1(pred)=%d  %.1f ms  it=%d  %s\n",
            string(inst.mode), string(inst.shape), inst.Mx, inst.My, inst.n, inst.k,
            info.nvtx, info.h1, 1e3t, res.niter, res.status)
    mse = inst.mode === :regress ? fitted_mse(inst, res, info) : NaN
    @printf("    pres=%.2e  dres=%.2e  min-coeff=%.2e  fit-mse=%.2e\n",
            pres, dres, mincoef, mse)
    return inst, res, info
end

# =====================================================================
# JuMP reference model (Mosek comparison)
# =====================================================================

function build_jump_tensor(inst::TensorSplineInstance, optimizer)
    Mx, My, n, k = inst.Mx, inst.My, inst.n, inst.k
    ns = n + 1; d = ns * ns
    Ins = Matrix{Float64}(I, ns, ns)
    Rj(ke) = jet_operator(n, 1, ke, 1.0, :right)
    Lj(ke) = jet_operator(n, 1, ke, 1.0, :left)
    Qd, Cd = patch_QC(inst)

    model = Model(optimizer); set_silent(model)
    P = Dict{Tuple{Int, Int}, Vector{VariableRef}}()
    for b in 1:My, a in 1:Mx
        inst.active[a, b] || continue
        P[(a, b)] = inst.shape === :nonneg ?
            @variable(model, [1:d], lower_bound = 0.0) : @variable(model, [1:d])
    end
    if inst.shape === :monotone_x || inst.shape === :monotone_y
        D1 = diffop(n, 1)
        D = inst.shape === :monotone_x ? kron(Ins, D1) : kron(D1, Ins)
        for Pv in values(P); @constraint(model, D * Pv .>= 0); end
    end
    for b in 1:My, a in 1:Mx
        haskey(P, (a, b)) || continue
        if a < Mx && haskey(P, (a + 1, b)) && inst.kx[a, b] ≥ 0
            ke = inst.kx[a, b]
            @constraint(model, kron(Ins, Rj(ke)) * P[(a, b)] .== kron(Ins, Lj(ke)) * P[(a + 1, b)])
        end
        if b < My && haskey(P, (a, b + 1)) && inst.ky[a, b] ≥ 0
            ke = inst.ky[a, b]
            @constraint(model, kron(Rj(ke), Ins) * P[(a, b)] .== kron(Lj(ke), Ins) * P[(a, b + 1)])
        end
    end
    obj = sum(0.5 * P[cell]' * Qd[cell] * P[cell] + Cd[cell]' * P[cell] for cell in keys(P))
    if inst.mode === :intensity
        N = length(inst.obs)
        @variable(model, t[1:N])
        for (j, (x, y)) in enumerate(inst.obs)
            a, u = locate01(x, Mx); b, v = locate01(y, My)
            φ = kron(bernstein_eval(n, v), bernstein_eval(n, u))
            @constraint(model, [t[j], 1.0, φ' * P[(a, b)]] in MOI.ExponentialCone())
        end
        obj -= sum(inst.w[j] * t[j] for j in 1:N)
    end
    @objective(model, Min, obj)
    return model, P
end

# =====================================================================
# Benchmark: Sheaf IPM vs Mosek (2D factorization path)
# =====================================================================

function run_tensor_benchmark(; optimizer = nothing, dual_optimizer = nothing,
                               solver_name::String = "Mosek", nwarmup::Int = 2, nruns::Int = 5)
    optimizer === nothing && error("Pass optimizer, e.g. run_tensor_benchmark(optimizer=Mosek.Optimizer)")
    # (mode, shape, Mx, My, holes, faults, raug) — raug tuned per case from sweep
    cases = [
        (:regress, :nonneg,      4,  4, Tuple{Int, Int}[], Tuple{Int, Int, Symbol}[], 1e5),
        (:regress, :nonneg,      8,  8, Tuple{Int, Int}[], Tuple{Int, Int, Symbol}[], 1e5),
        (:regress, :nonneg,     16, 16, Tuple{Int, Int}[], Tuple{Int, Int, Symbol}[], 1e5),
        (:regress, :monotone_x,  8,  8, Tuple{Int, Int}[], Tuple{Int, Int, Symbol}[], 1e5),
        (:regress, :nonneg,      8,  8, [(4, 4), (5, 5)],  Tuple{Int, Int, Symbol}[], 1e5),  # holes
        (:intensity, :nonneg,    4,  4, Tuple{Int, Int}[], Tuple{Int, Int, Symbol}[], 1e1),
    ]
    sname = rpad(solver_name, 9)
    sname_d = rpad(solver_name * "-D", 9)
    println("="^100)
    println("Tensor-Spline Benchmark (2D grid, lossy variants): Sheaf IPM vs $solver_name")
    println("="^100)
    if dual_optimizer !== nothing
        @printf("%-10s %-11s %6s %6s %5s %5s %5s %9s %9s %9s %7s %7s\n",
                "Mode", "Shape", "raug", "grid", "|V|", "H1", "n", "IPM(ms)", sname, sname_d, "P/IPM", "D/IPM")
    else
        @printf("%-10s %-11s %6s %6s %5s %5s %5s %9s %9s %8s\n",
                "Mode", "Shape", "raug", "grid", "|V|", "H1", "n", "IPM(ms)", sname, "speedup")
    end
    println("-"^100)
    for (mode, shape, Mx, My, holes, faults, raug) in cases
        inst = generate_tensor_instance(; mode, shape, Mx, My, holes, faults, N = 50 * Mx * My)
        prob, info = build_tensor_problem(inst)
        settings = IPMSettings{Float64}(
            kkt = UzawaSettings{Float64}(raug = raug),
            feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)
        for _ in 1:nwarmup; solve(prob, settings); end
        t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])
        for _ in 1:nwarmup
            m, _ = build_jump_tensor(inst, optimizer); optimize!(m)
        end
        t_mosek = minimum([begin
            m, _ = build_jump_tensor(inst, optimizer)
            @elapsed optimize!(m)
        end for _ in 1:nruns])

        if dual_optimizer !== nothing
            for _ in 1:nwarmup
                m, _ = build_jump_tensor(inst, dual_optimizer); optimize!(m)
            end
            t_dual = minimum([begin
                m, _ = build_jump_tensor(inst, dual_optimizer)
                @elapsed optimize!(m)
            end for _ in 1:nruns])
            @printf("%-10s %-11s %6.0e %3dx%-3d %5d %5d %5d %9.1f %9.1f %9.1f %6.2fx %6.2fx\n",
                    string(mode), string(shape), raug, Mx, My, info.nvtx, info.h1, inst.n,
                    t_ipm * 1000, t_mosek * 1000, t_dual * 1000,
                    t_mosek / t_ipm, t_dual / t_ipm)
        else
            @printf("%-10s %-11s %6.0e %3dx%-3d %5d %5d %5d %9.1f %9.1f %7.2fx\n",
                    string(mode), string(shape), raug, Mx, My, info.nvtx, info.h1, inst.n,
                    t_ipm * 1000, t_mosek * 1000, t_mosek / t_ipm)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmark_utils.jl"))
    opts = parse_benchmark_args(ARGS)

    if opts.mosek
        using MosekTools
    else
        using Clarabel
    end
    optimizer, dual_optimizer = get_optimizers(opts; lp_only = false)
    if !opts.mosek
        dual_optimizer = nothing  # skip dualization (ExponentialCone not supported)
    end
    solver_name = opts.mosek ? "Mosek" : "Clarabel"
    print_benchmark_config(opts; lp_only = false)

    run_tensor_benchmark(; optimizer, dual_optimizer, solver_name,
                         nwarmup = opts.nwarmup, nruns = opts.nruns)
end
