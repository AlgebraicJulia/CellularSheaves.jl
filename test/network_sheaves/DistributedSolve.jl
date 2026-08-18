using Test
using CellularSheaves
using CliqueTrees.Multifrontal
using LinearAlgebra
using Graphs
using SparseArrays
using Distributed
using Random

@testset "DistributedSolve" begin

    # A small diagonal shift keeps the sheaf Laplacian strictly positive
    # definite (it's only PSD, with a 1-D nullspace, for a connected graph
    # with identity restriction maps) so it can be factored with a plain
    # ChordalCholesky. The shift doesn't change the clique tree, which is
    # the only thing this test is actually exercising.
    function spd_sheaf_laplacian(g)
        s = sheaf_from_graph(g, 1, d -> Matrix{Float64}(I, d, d); symmetric_edges=true)
        return sparse(Matrix(sheaf_laplacian_matrix(s)) + 1e-6I)
    end

    function check_topology(g)
        M = spd_sheaf_laplacian(g)
        n = size(M, 1)
        F = cholesky!(ChordalCholesky(M), NoPivot())
        L = F.L

        # one workspace, reused across every solve below — the intended usage
        # (build once per factorization, re-solve every step)
        ws = TreeWorkspace(L, Float64)

        for _ in 1:5
            b = rand(n)

            expected_fwd = copy(b)
            ldiv!(L, expected_fwd)
            actual_fwd = copy(b)
            tree_forward_ldiv!(L, actual_fwd)
            @test actual_fwd ≈ expected_fwd
            # workspace forward matches the recursive forward (hence ldiv!)
            ws_fwd = copy(b)
            tree_forward_ldiv!(L, ws_fwd, ws)
            @test ws_fwd ≈ expected_fwd

            expected_full = copy(expected_fwd)
            ldiv!(L', expected_full)
            actual_full = copy(actual_fwd)
            tree_backward_ldiv!(L, actual_full)
            @test actual_full ≈ expected_full
            # workspace backward matches too, completing the full solve
            ws_full = copy(ws_fwd)
            tree_backward_ldiv!(L, ws_full, ws)
            @test ws_full ≈ expected_full
        end

        # the workspace path hoists the per-node Dict/collect/findfirst work
        # out of the solve, so it allocates far less than the recursive path
        # (the reason to prefer it when re-solving every step)
        b = rand(n)
        r1 = copy(b); tree_forward_ldiv!(L, r1)          # warm up / compile
        w1 = copy(b); tree_forward_ldiv!(L, w1, ws)
        r2 = copy(b); w2 = copy(b)
        rec_alloc = @allocated tree_forward_ldiv!(L, r2)
        ws_alloc = @allocated tree_forward_ldiv!(L, w2, ws)
        @test ws_alloc < rec_alloc
    end

    @testset "wide/star topology" begin
        check_topology(star_graph(10))
    end

    @testset "narrow/path topology" begin
        check_topology(path_graph(12))
    end

    @testset "single-vertex edge case" begin
        check_topology(path_graph(1))
    end

    # star and path are deliberately the two extremes (widest vs. narrowest
    # possible tree) — good for exercising edge cases, but nothing in
    # tree_forward_ldiv!/tree_backward_ldiv! is specific to either shape.
    # These exercise a general irregular graph and a 2D grid (the latter
    # relevant to this repo's multi-agent tracking examples) to confirm that.
    function random_connected_graph(n, seed)
        rng = Random.MersenneTwister(seed)
        g = Graphs.erdos_renyi(n, 0.2, rng = rng)
        for i in 1:(n - 1)
            add_edge!(g, i, i + 1)
        end
        return g
    end

    @testset "general topology: random irregular graph" begin
        check_topology(random_connected_graph(16, 42))
    end

    @testset "general topology: 2D grid graph" begin
        check_topology(Graphs.grid([4, 5]))
    end

    function check_partition(g, target_workers)
        M = spd_sheaf_laplacian(g)
        F = cholesky!(ChordalCholesky(M), NoPivot())
        L = F.L
        S = L.S
        n = nv(S.res)

        p = partition_tree(L, target_workers)

        # every node assigned to exactly one worker
        @test sort(vcat(p.chunks...)) == collect(1:n)
        @test all(!=(0), p.owner)

        # cut edges are provably minimal for a single-rooted tree: workers - 1
        nroots = count(j -> iszero(S.pnt[j]), 1:n)
        @test length(p.cut_edges) == length(p.chunks) - nroots

        # every chunk induces one connected rooted subtree
        for members in p.chunks
            mset = Set(members)
            local_roots = [m for m in members if iszero(S.pnt[m]) || !(S.pnt[m] in mset)]
            @test length(local_roots) == 1
            seen = Set{Int}()
            stack = [local_roots[1]]
            while !isempty(stack)
                x = pop!(stack)
                push!(seen, x)
                for c in neighbors(S.chd, x)
                    c in mset && push!(stack, c)
                end
            end
            @test seen == mset
        end

        # each worker's extracted slice matches the original factor data
        # bit-for-bit, and every entry is covered by exactly one worker
        totalD = 0
        totalL = 0
        for w in 1:length(p.chunks)
            wd = worker_factorization(L, p, w)
            totalD += length(wd.Dval)
            totalL += length(wd.Lval)
            for (k, j) in enumerate(wd.ids)
                Dp, Dp2 = S.Dptr[j], S.Dptr[j + 1] - 1
                @test L.Dval[Dp:Dp2] == wd.Dval[wd.Dptr[k]:wd.Dptr[k + 1] - 1]
                Lp, Lp2 = S.Lptr[j], S.Lptr[j + 1] - 1
                @test L.Lval[Lp:Lp2] == wd.Lval[wd.Lptr[k]:wd.Lptr[k + 1] - 1]
            end
        end
        @test totalD == S.Dptr[n + 1] - 1
        @test totalL == S.Lptr[n + 1] - 1
    end

    @testset "partitioning: wide/star collapses toward one chunk" begin
        for tw in (1, 2, 3)
            check_partition(star_graph(10), tw)
        end
        # only fragments once the threshold drops below a single node's weight
        check_partition(star_graph(10), 9)
    end

    @testset "partitioning: narrow/path splits cleanly" begin
        for tw in (1, 2, 3, 4)
            check_partition(path_graph(12), tw)
        end
    end

    @testset "partitioning: single-vertex edge case" begin
        check_partition(path_graph(1), 1)
    end

    @testset "partitioning: general topology" begin
        for tw in (1, 2, 3, 4)
            check_partition(random_connected_graph(16, 42), tw)
            check_partition(Graphs.grid([4, 5]), tw)
        end
    end

    @testset "distributed_tree_solve: real multi-process run" begin
        active = Int[]
        for p in workers()
            try
                remotecall_fetch(myid, p)
                push!(active, p)
            catch
                try rmprocs(p; waitfor=1.0) catch; end
            end
        end
        needed = 5 - length(active)
        if needed > 0
            added = addprocs(needed; exeflags = ["--project=$(Base.active_project())"])
            append!(active, added)
        end
        pids = active[1:5]
        @everywhere pids using CellularSheaves

        function check_distributed(g, target_workers)
            M = spd_sheaf_laplacian(g)
            n = size(M, 1)
            F = cholesky!(ChordalCholesky(M), NoPivot())
            L = F.L

            for _ in 1:3
                b = rand(n)
                expected = copy(b)
                ldiv!(L, expected)
                ldiv!(L', expected)

                actual = distributed_tree_solve(L, b, target_workers; pids = pids)
                @test actual ≈ expected
            end
        end

        # path splits cleanly into exactly 2 chunks -> exercises real
        # cross-process RemoteChannel messages at the one cut edge
        check_distributed(path_graph(12), 2)

        # star collapses to 1 chunk at this target -> exercises the
        # single-worker, zero-cut-edge path (no messages needed at all)
        check_distributed(star_graph(10), 2)

        # a general graph -> exercises multiple cut edges / multiple
        # simultaneous cross-process channels at once, not just one
        check_distributed(random_connected_graph(16, 42), 3)
        check_distributed(Graphs.grid([4, 5]), 3)
    end
end
