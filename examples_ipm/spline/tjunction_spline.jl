######################################################################
# tjunction_spline.jl
#
# TWO-LEVEL ADAPTIVE REFINEMENT of the tensor-product surface: selected
# coarse cells are replaced by their four half-size children, producing
# T-junctions — and the corpus's first genuinely NON-IDENTITY,
# NON-SQUARE-SYMMETRIC restriction maps with classical semantics.
#
# At a T-junction a coarse patch abuts two fine children along one
# line. The edge says: the coarse boundary jet, RESTRICTED to the
# child's sub-segment and reparametrized, equals the child's boundary
# jet. Restriction-to-half-segment in Bernstein coordinates is the
# de Casteljau subdivision matrix S (splitting at t = ½), so a vertical
# T-edge reads
#
#     kron(S_half, jetₓ(T = h))·P_coarse = kron(I, jetₓ(T = h/2))·P_fine
#
# with PHYSICAL jet scaling (jet_operator's T argument finally earns
# its keep: coarse and fine cells have different widths, and C^k must
# hold in physical coordinates). Coarse–coarse lines keep one full
# edge; splitting them into two subdivided half-edges would impose the
# same constraint set with (k+1)(n+1) redundant rows per line — a
# bookkeeping trap the enumerator avoids.
#
# H¹ (measured FIRST in tjunction_spline_oracle.py, then stated). The
# JUNCTION RULE
#
#     dim ker δᵀ = (k+1)² × #{interior fine-mesh vertices where
#                             ≥ 3 distinct patches meet}
#
# — fine corners, T-points (coarse + two children), and mixed coarse
# corners all contribute the same (k+1)² block — is EXACT in every
# tested configuration with n > 2k, and subsumes the uniform-grid
# corner rule (4 patches at a corner). For n ≤ 2k a refined region
# enclosed by a coarse ring behaves like a PARTIALLY FILLED HOLE:
# measured extras appear, bounded by (2k−n+1)₊² per coarse ring not
# spanned by pure-coarse corner cycles (+1 per single-cell ring at
# n = 2k for k = 2 and 3; +3 of the bound 4 at (n,k) = (3,2), the
# classical cubic-C² family; +0 for an enclosed 2×2 block). The exact
# n ≤ 2k law is open; h1_predict_tj returns the junction rule — a
# lower bound there — and the oracle measures the truth per instance.
#
# Scope: mode = :regress, shape ∈ (:nonneg, :free), full rectangle,
# uniform k. Monotone leaves, the :intensity likelihood, and the lossy
# dials of tensor_spline.jl compose per-patch and per-edge without
# interaction — kept out of v1 to isolate the T-junction machinery.
# T-spline linear-independence / analysis-suitability issues do NOT
# arise here: no reduced global basis is ever constructed; we constrain
# raw per-patch polynomials by exact jet agreement.
#
# NOTE: written against the CellularSheaves.IPM PR-67 API; not executed
# here. tjunction_spline_oracle.py is the numerical ground truth.
######################################################################

include("tensor_spline.jl")   # tensor numerics + univariate closure + IPM

# ---- de Casteljau subdivision at t = 1/2 -----------------------------

"""(S_L, S_R): Bernstein coefficients of the left/right half, reparametrized.
S_L[i,j] = C(i,j)c^j(1−c)^{i−j} (j ≤ i);  S_R[i,j] = C(n−i,j−i)c^{j−i}(1−c)^{n−j}."""
function subdiv_half(n::Int; c::Float64 = 0.5)
    SL = zeros(n + 1, n + 1); SR = zeros(n + 1, n + 1)
    for i in 0:n, j in 0:n
        j ≤ i && (SL[i + 1, j + 1] = binomial(i, j) * c^j * (1 - c)^(i - j))
        j ≥ i && (SR[i + 1, j + 1] = binomial(n - i, j - i) * c^(j - i) * (1 - c)^(n - j))
    end
    return SL, SR
end

# ---- instance --------------------------------------------------------

"""Patch keys: (:c, a, b) coarse on the Mx×My grid; (:f, p, q) fine on the
doubled 2Mx×2My grid (children of refined cells only)."""
struct TJunctionInstance
    Mx::Int; My::Int; n::Int; k::Int
    shape::Symbol                                   # :nonneg | :free
    refined::BitMatrix                              # Mx×My
    Qp::Dict{Tuple{Symbol, Int, Int}, Matrix{Float64}}
    cp::Dict{Tuple{Symbol, Int, Int}, Vector{Float64}}
    λr::Float64; ε::Float64
end

"""Owning patch key and its local coordinates for a physical point."""
function tj_owner(refined, Mx, My, x, y)
    p, uf = locate01(x, 2 * Mx); q, vf = locate01(y, 2 * My)
    a, b = cld(p, 2), cld(q, 2)
    if refined[a, b]
        return (:f, p, q), uf, vf
    else
        _, uc = locate01(x, Mx); _, vc = locate01(y, My)
        return (:c, a, b), uc, vc
    end
end

"""Mark the `nref` cells with the largest mean truth for refinement."""
function refine_by_truth(Mx, My; nref::Int = max(1, (Mx * My) ÷ 4))
    score = [(sum(truth2(:free, (a - 1 + s) / Mx, (b - 1 + t) / My)
                  for s in 0.1:0.2:0.9, t in 0.1:0.2:0.9), a, b)
             for a in 1:Mx, b in 1:My]
    refined = falses(Mx, My)
    for (_, a, b) in sort(vec(score); rev = true)[1:nref]
        refined[a, b] = true
    end
    return refined
end

function generate_tjunction_instance(;
        Mx::Int = 3, My::Int = 3, n::Int = 4, k::Int = 2,
        shape::Symbol = :nonneg, refine = Tuple{Int, Int}[],
        refined::Union{Nothing, BitMatrix} = nothing,
        N::Int = 100 * Mx * My, σ::Float64 = 0.05,
        λr::Float64 = 1e-4, ε::Float64 = 1e-8, seed::Int = 1)
    @assert n ≥ max(k + 1, 2) && shape in (:nonneg, :free)
    rng = MersenneTwister(seed)
    ns = n + 1; d = ns * ns
    R = refined === nothing ? falses(Mx, My) : copy(refined)
    for (a, b) in refine; R[a, b] = true; end

    Qp = Dict{Tuple{Symbol, Int, Int}, Matrix{Float64}}()
    cp = Dict{Tuple{Symbol, Int, Int}, Vector{Float64}}()
    for b in 1:My, a in 1:Mx
        if R[a, b]
            for q in (2b - 1):(2b), p in (2a - 1):(2a)
                Qp[(:f, p, q)] = zeros(d, d); cp[(:f, p, q)] = zeros(d)
            end
        else
            Qp[(:c, a, b)] = zeros(d, d); cp[(:c, a, b)] = zeros(d)
        end
    end
    for _ in 1:N
        x, y = rand(rng), rand(rng)
        z = truth2(shape, x, y) + σ * randn(rng)
        key, u, v = tj_owner(R, Mx, My, x, y)
        φ = kron(bernstein_eval(n, v), bernstein_eval(n, u))
        Qp[key] .+= 2.0 .* (φ * φ'); cp[key] .-= 2.0 .* z .* φ
    end
    return TJunctionInstance(Mx, My, n, k, shape, R, Qp, cp, λr, ε)
end

# ---- shared per-patch objective and edge enumeration -----------------

function tj_patch_QC(inst::TJunctionInstance)
    n = inst.n; ns = n + 1; d = ns * ns
    hx, hy = 1.0 / inst.Mx, 1.0 / inst.My
    QC = inst.λr .* thinplate_gram(n, hx, hy) .+ inst.ε .* Matrix{Float64}(I, d, d)
    QF = inst.λr .* thinplate_gram(n, hx / 2, hy / 2) .+ inst.ε .* Matrix{Float64}(I, d, d)
    Q = Dict{Tuple{Symbol, Int, Int}, Matrix{Float64}}()
    C = Dict{Tuple{Symbol, Int, Int}, Vector{Float64}}()
    for (key, Qdata) in inst.Qp
        Q[key] = (key[1] === :c ? QC : QF) .+ Qdata
        C[key] = inst.cp[key]
    end
    return Q, C
end

"""
    tj_edges(inst) :: Vector{(keyL, AL, keyR, AR)}   meaning AL·P_L = AR·P_R

Enumeration order (mirrored by the oracle): vertical coarse lines
(b outer, a inner; one full edge if both sides coarse, else two
sub-segment edges, lower first); interior vertical half-lines of
refined cells; horizontal coarse lines; interior horizontal half-lines.
"""
function tj_edges(inst::TJunctionInstance)
    Mx, My, n, k = inst.Mx, inst.My, inst.n, inst.k
    ns = n + 1
    hx, hy = 1.0 / Mx, 1.0 / My
    Ins = Matrix{Float64}(I, ns, ns)
    SL, SR = subdiv_half(n); S = (SL, SR)
    RC(T) = jet_operator(n, 1, k, T, :right)
    LC(T) = jet_operator(n, 1, k, T, :left)
    E = Tuple{Tuple{Symbol, Int, Int}, Matrix{Float64},
              Tuple{Symbol, Int, Int}, Matrix{Float64}}[]
    Rf = inst.refined

    for b in 1:My, a in 1:(Mx - 1)                         # vertical lines
        Lr, Rr = Rf[a, b], Rf[a + 1, b]
        if !Lr && !Rr
            push!(E, ((:c, a, b), kron(Ins, RC(hx)), (:c, a + 1, b), kron(Ins, LC(hx))))
        else
            for jh in 1:2
                q = 2b - 2 + jh
                kL, AL = Lr ? ((:f, 2a, q), kron(Ins, RC(hx / 2))) :
                              ((:c, a, b), kron(S[jh], RC(hx)))
                kR, AR = Rr ? ((:f, 2a + 1, q), kron(Ins, LC(hx / 2))) :
                              ((:c, a + 1, b), kron(S[jh], LC(hx)))
                push!(E, (kL, AL, kR, AR))
            end
        end
    end
    for b in 1:My, a in 1:Mx                               # interior vertical
        Rf[a, b] || continue
        for jh in 1:2
            q = 2b - 2 + jh
            push!(E, ((:f, 2a - 1, q), kron(Ins, RC(hx / 2)),
                      (:f, 2a, q), kron(Ins, LC(hx / 2))))
        end
    end
    for b in 1:(My - 1), a in 1:Mx                         # horizontal lines
        Dr, Ur = Rf[a, b], Rf[a, b + 1]
        if !Dr && !Ur
            push!(E, ((:c, a, b), kron(RC(hy), Ins), (:c, a, b + 1), kron(LC(hy), Ins)))
        else
            for ih in 1:2
                p = 2a - 2 + ih
                kD, AD = Dr ? ((:f, p, 2b), kron(RC(hy / 2), Ins)) :
                              ((:c, a, b), kron(RC(hy), S[ih]))
                kU, AU = Ur ? ((:f, p, 2b + 1), kron(LC(hy / 2), Ins)) :
                              ((:c, a, b + 1), kron(LC(hy), S[ih]))
                push!(E, (kD, AD, kU, AU))
            end
        end
    end
    for b in 1:My, a in 1:Mx                               # interior horizontal
        Rf[a, b] || continue
        for ih in 1:2
            p = 2a - 2 + ih
            push!(E, ((:f, p, 2b - 1), kron(RC(hy / 2), Ins),
                      (:f, p, 2b), kron(LC(hy / 2), Ins)))
        end
    end
    return E
end

"""H¹ predictor — the junction rule: (k+1)² per interior fine-mesh vertex
where ≥ 3 distinct patches meet. EXACT for n > 2k (oracle-verified); a
lower bound for n ≤ 2k when coarse rings enclose refined regions."""
function h1_predict_tj(inst::TJunctionInstance)
    Mx, My = inst.Mx, inst.My
    own(p, q) = (a = cld(p, 2); b = cld(q, 2);
                 inst.refined[a, b] ? (:f, p, q) : (:c, a, b))
    cnt = 0
    for q in 1:(2My - 1), p in 1:(2Mx - 1)
        owners = Set([own(p, q), own(p + 1, q), own(p, q + 1), own(p + 1, q + 1)])
        length(owners) ≥ 3 && (cnt += 1)
    end
    return cnt * (inst.k + 1)^2
end

# ---- native model ----------------------------------------------------

"""Vertex order: coarse patches (b outer, a inner, unrefined cells),
then fine patches (q outer, p inner, refined parents); edges per tj_edges."""
function build_tjunction_problem(inst::TJunctionInstance)
    ns = inst.n + 1; d = ns * ns

    row_ids = Int[]; col_ids = Int[]; blocks = Matrix{Float64}[]
    col_cone = AbstractCone[]; col_dim = Int[]
    new_vertex!(dim, cone) = (push!(col_dim, dim); push!(col_cone, cone); length(col_dim))
    ec = Ref(0); new_edge!() = (ec[] += 1; ec[])
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v); push!(blocks, Matrix{Float64}(mat)))

    cone() = inst.shape === :nonneg ? PositiveCone() : CofreeCone()
    vtx = Dict{Tuple{Symbol, Int, Int}, Int}()
    for b in 1:inst.My, a in 1:inst.Mx
        inst.refined[a, b] && continue
        vtx[(:c, a, b)] = new_vertex!(d, cone())
    end
    for q in 1:(2 * inst.My), p in 1:(2 * inst.Mx)
        inst.refined[cld(p, 2), cld(q, 2)] || continue
        vtx[(:f, p, q)] = new_vertex!(d, cone())
    end

    for (kL, AL, kR, AR) in tj_edges(inst)
        e = new_edge!(); place!(e, vtx[kL], AL); place!(e, vtx[kR], -AR)
    end

    B = blocksparse(row_ids, col_ids, blocks)
    nvtx = length(col_dim)
    @assert nvtxs(B) == nvtx "column blocks not contiguous 1..nvtx"
    g = zeros(size(B, 1))                                  # exact: all-zero rhs

    Qd, Cd = tj_patch_QC(inst)
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for (key, v) in vtx
        block(Q, v, v, v) .= Qd[key]
        c[colrange(B, v)] .= Cd[key]
    end
    prob = IPMProblem(Q, B, c, g, AbstractCone[col_cone[v] for v in 1:nvtx])
    info = (patchcols = Dict(key => colrange(B, v) for (key, v) in vtx),
            nvtx = nvtx, ncols = size(B, 2), h1 = h1_predict_tj(inst))
    return prob, info
end

"""Fitted MSE against the truth (owner-aware evaluation)."""
function tj_mse(inst, res, info; ngrid = 80)
    n = inst.n; ns = n + 1
    tot = 0.0
    for gx in 1:ngrid, gy in 1:ngrid
        x, y = (gx - 0.5) / ngrid, (gy - 0.5) / ngrid
        key, u, v = tj_owner(inst.refined, inst.Mx, inst.My, x, y)
        P = reshape(res.p[info.patchcols[key]], ns, ns)
        f = bernstein_eval(n, u)' * P * bernstein_eval(n, v)
        tot += (f - truth2(inst.shape, x, y))^2
    end
    return tot / ngrid^2
end

function tjunction_demo(; kwargs...)
    inst = generate_tjunction_instance(; refine = [(2, 2)], kwargs...)
    prob, info = build_tjunction_problem(inst)
    solve(prob, tensor_settings())
    t = @elapsed (res = solve(prob, tensor_settings()))
    pres = isempty(res.history) ? NaN : res.history.pres[end]
    @printf("TJ %-8s  %dx%d(+%d ref) n=%d k=%d  vtx=%d dof=%d  H1(pred)=%d  %.1f ms it=%d %s\n",
            string(inst.shape), inst.Mx, inst.My, count(inst.refined), inst.n, inst.k,
            info.nvtx, info.ncols, info.h1, 1e3t, res.niter, res.status)
    @printf("    pres=%.2e  fit-mse=%.2e\n", pres, tj_mse(inst, res, info))
    return inst, res, info
end

"""The adaptive payoff, at a FIXED data budget: the coarse mesh underfits,
the uniform-fine mesh overfits (25 coefficients per patch vs ~N/(2M)² points),
and adaptive refinement beats both by spending coefficients only where the
truth has structure."""
function adaptive_demo(; M::Int = 4, n::Int = 4, k::Int = 2, N::Int = 2400, seed::Int = 1)
    for (label, Mx, My, refined) in
        (("coarse  $(M)x$(M)", M, M, falses(M, M)),
         ("adaptive$(M)x$(M)", M, M, refine_by_truth(M, M; nref = M * M ÷ 4)),
         ("fine    $(2M)x$(2M)", 2M, 2M, falses(2M, 2M)))
        inst = generate_tjunction_instance(; Mx, My, n, k, refined = BitMatrix(refined),
                                           N, seed)
        prob, info = build_tjunction_problem(inst)
        res = solve(prob, tensor_settings())
        @printf("  %-14s dof=%5d  H1=%3d  fit-mse=%.3e  it=%d %s\n",
                label, info.ncols, info.h1, tj_mse(inst, res, info), res.niter, res.status)
    end
end

# =====================================================================
# JuMP reference (shares tj_patch_QC and tj_edges)
# =====================================================================

function build_jump_tjunction(inst::TJunctionInstance, optimizer)
    d = (inst.n + 1)^2
    Qd, Cd = tj_patch_QC(inst)
    model = Model(optimizer); set_silent(model)
    P = Dict(key => (inst.shape === :nonneg ?
                     @variable(model, [1:d], lower_bound = 0.0) :
                     @variable(model, [1:d])) for key in keys(Qd))
    for (kL, AL, kR, AR) in tj_edges(inst)
        @constraint(model, AL * P[kL] .== AR * P[kR])
    end
    @objective(model, Min,
               sum(0.5 * P[key]' * Qd[key] * P[key] + Cd[key]' * P[key] for key in keys(P)))
    return model, P
end

# =====================================================================
# Benchmark: Sheaf IPM vs Mosek (T-junction LP path)
# =====================================================================

function run_tjunction_benchmark(; optimizer = nothing, dual_optimizer = nothing, solver_name::String = "Mosek", nwarmup::Int = 2, nruns::Int = 5)
    optimizer === nothing && error("Pass optimizer, e.g. run_tjunction_benchmark(optimizer=Mosek.Optimizer)")
    # (Mx, My, refine_list, label, raug)
    cases = [
        (3, 3, [(2, 2)], "3x3 center",  1e5),
        (4, 4, [(2, 2), (3, 3)], "4x4 diag",  1e5),
        (4, 4, [(2, 2), (2, 3), (3, 2), (3, 3)], "4x4 block",  1e5),
        (6, 6, [(3, 3), (3, 4), (4, 3), (4, 4)], "6x6 block",  1e5),
        (8, 8, [(i, j) for i in 3:6, j in 3:6], "8x8 4x4blk", 1e5),
    ]
    sname = rpad(solver_name, 9)
    sname_d = rpad(solver_name * "-D", 9)
    println("="^95)
    println("T-Junction Benchmark (adaptive LP): Sheaf IPM vs $solver_name")
    println("="^95)
    if dual_optimizer !== nothing
        @printf("%-12s %6s %5s %5s %5s %9s %9s %9s %7s %7s\n",
                "Config", "raug", "DOF", "|V|", "H1", "IPM(ms)", sname, sname_d, "P/IPM", "D/IPM")
    else
        @printf("%-12s %6s %5s %5s %5s %9s %9s %8s\n",
                "Config", "raug", "DOF", "|V|", "H1", "IPM(ms)", sname, "speedup")
    end
    println("-"^95)
    for (Mx, My, refine_list, label, raug) in cases
        inst = generate_tjunction_instance(; Mx, My, refine = refine_list, N = 100 * Mx * My)
        prob, info = build_tjunction_problem(inst)
        settings = IPMSettings{Float64}(
            kkt = UzawaSettings{Float64}(raug = raug),
            feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)
        for _ in 1:nwarmup; solve(prob, settings); end
        t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])
        for _ in 1:nwarmup
            m, _ = build_jump_tjunction(inst, optimizer); optimize!(m)
        end
        t_mosek = minimum([@elapsed begin
            m, _ = build_jump_tjunction(inst, optimizer); optimize!(m)
        end for _ in 1:nruns])

        if dual_optimizer !== nothing
            for _ in 1:nwarmup
                m, _ = build_jump_tjunction(inst, dual_optimizer); optimize!(m)
            end
            t_dual = minimum([@elapsed begin
                m, _ = build_jump_tjunction(inst, dual_optimizer); optimize!(m)
            end for _ in 1:nruns])
            @printf("%-12s %6.0e %5d %5d %5d %9.1f %9.1f %9.1f %6.2fx %6.2fx\n",
                    label, raug, info.ncols, info.nvtx, info.h1,
                    t_ipm * 1000, t_mosek * 1000, t_dual * 1000,
                    t_mosek / t_ipm, t_dual / t_ipm)
        else
            @printf("%-12s %6.0e %5d %5d %5d %9.1f %9.1f %7.2fx\n",
                    label, raug, info.ncols, info.nvtx, info.h1,
                    t_ipm * 1000, t_mosek * 1000, t_mosek / t_ipm)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmark_utils.jl"))
    opts = parse_benchmark_args(ARGS)

    if opts.mosek
        using MosekTools
        optimizer = Mosek.Optimizer
        dual_optimizer = Dualization.dual_optimizer(Mosek.Optimizer)
        solver_name = "Mosek"
    else
        using OSQP
        optimizer = OSQP.Optimizer
        dual_optimizer = Dualization.dual_optimizer(OSQP.Optimizer)
        solver_name = "OSQP"
    end
    println("Solver: $solver_name")
    println("Runs: $(opts.nruns), Warmup: $(opts.nwarmup)\n")

    run_tjunction_benchmark(; optimizer, dual_optimizer, solver_name,
                            nwarmup = opts.nwarmup, nruns = opts.nruns)
end
