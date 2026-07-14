# Benchmark script for symnorm_impl
#
# Usage:
#   julia --project=. src/BlockSparseArrays/bench/symnorm_bench.jl

using LinearAlgebra
using Printf
using Random

using CellularSheaves.BlockSparseArrays
using CellularSheaves.BlockSparseArrays: blocksparse, vtxs, ncols

# Build a block-diagonal matrix (like H in the IPM)
function make_block_diagonal(nblocks::Int, blocksize::Int; T=Float64)
    II = Int[]
    JJ = Int[]
    VV = Matrix{T}[]

    for v in 1:nblocks
        push!(II, v)
        push!(JJ, v)
        # Random SPD block
        A = randn(T, blocksize, blocksize)
        push!(VV, A' * A + LinearAlgebra.I * blocksize)
    end

    return blocksparse(II, JJ, VV)
end

# Build a block-diagonal with varying block sizes (like E11's H with pair + SOC stalks)
function make_mixed_block_diagonal(; T=Float64, npair=256, nsoc=257, pair_size=24, soc_size=6)
    II = Int[]
    JJ = Int[]
    VV = Matrix{T}[]

    # Pair stalks (larger blocks)
    for v in 1:npair
        push!(II, v)
        push!(JJ, v)
        A = randn(T, pair_size, pair_size)
        push!(VV, A' * A + LinearAlgebra.I * pair_size)
    end

    # SOC stalks (smaller blocks)
    for v in 1:nsoc
        push!(II, npair + v)
        push!(JJ, npair + v)
        A = randn(T, soc_size, soc_size)
        push!(VV, A' * A + LinearAlgebra.I * soc_size)
    end

    return blocksparse(II, JJ, VV)
end

# Build a path graph with off-diagonal blocks (like B'B structure)
function make_path_graph(nvertices::Int, blocksize::Int; T=Float64)
    II = Int[]
    JJ = Int[]
    VV = Matrix{T}[]

    for v in 1:nvertices
        # Diagonal block
        push!(II, v)
        push!(JJ, v)
        A = randn(T, blocksize, blocksize)
        push!(VV, A' * A + LinearAlgebra.I * blocksize)

        # Off-diagonal block to next vertex (lower triangle)
        if v < nvertices
            push!(II, v + 1)
            push!(JJ, v)
            push!(VV, randn(T, blocksize, blocksize))
        end
    end

    return blocksparse(II, JJ, VV)
end

function run_bench()
    Random.seed!(42)
    nruns = 100
    nwarmup = 5

    println("=" ^ 70)
    println("  symnorm_impl benchmark")
    println("=" ^ 70)
    println()

    cases = [
        # (name, matrix)
        ("E11-like (256 pair + 257 soc)", make_mixed_block_diagonal(npair=256, nsoc=257, pair_size=24, soc_size=6)),
        ("uniform 24x24, 256 blocks", make_block_diagonal(256, 24)),
        ("uniform 24x24, 512 blocks", make_block_diagonal(512, 24)),
        ("uniform 6x6, 512 blocks", make_block_diagonal(512, 6)),
        ("uniform 48x48, 128 blocks", make_block_diagonal(128, 48)),
        ("path graph 24x24, 256 vtx", make_path_graph(256, 24)),
        ("path graph 12x12, 512 vtx", make_path_graph(512, 12)),
    ]

    println("Matrix info:")
    println("-" ^ 70)
    for (name, A) in cases
        nv = length(collect(vtxs(A)))
        nc = ncols(A)
        nval = length(A.val)
        @printf("  %-35s  vtx=%4d  cols=%5d  vals=%7d\n", name, nv, nc, nval)
    end
    println()

    println("Timings ($nruns runs each, excluding $nwarmup warmup):")
    println("-" ^ 70)
    @printf("  %-35s  %10s  %10s  %10s\n", "case", "sym norm", "raw norm", "ratio")
    println("-" ^ 70)

    for (name, A) in cases
        SA = Symmetric(A, :L)

        # Warmup
        for _ in 1:nwarmup
            norm(SA)
            norm(A)
        end

        # Timed runs
        t_sym = @elapsed for _ in 1:nruns
            norm(SA)
        end

        t_raw = @elapsed for _ in 1:nruns
            norm(A)
        end

        t_sym_us = t_sym / nruns * 1e6
        t_raw_us = t_raw / nruns * 1e6
        ratio = t_sym_us / t_raw_us

        @printf("  %-35s  %8.1f μs  %8.1f μs  %8.2fx\n", name, t_sym_us, t_raw_us, ratio)
    end

    println()
    println("Breakdown for E11-like case:")
    println("-" ^ 70)

    A = cases[1][2]
    SA = Symmetric(A, :L)

    # Time components
    t_normInf = @elapsed for _ in 1:nruns
        maximum(abs, A.val)
    end

    t_sum_sq = @elapsed for _ in 1:nruns
        sum(abs2, A.val)
    end

    t_sqrt_sum = @elapsed for _ in 1:nruns
        sqrt(sum(abs2, A.val))
    end

    # Re-time the norms for this specific case
    for _ in 1:nwarmup
        norm(SA)
        norm(A)
    end
    t_sym_e11 = @elapsed for _ in 1:nruns
        norm(SA)
    end
    t_raw_e11 = @elapsed for _ in 1:nruns
        norm(A)
    end

    @printf("  maximum(abs, vals):      %8.1f μs\n", t_normInf / nruns * 1e6)
    @printf("  sum(abs2, vals):         %8.1f μs\n", t_sum_sq / nruns * 1e6)
    @printf("  sqrt(sum(abs2, vals)):   %8.1f μs\n", t_sqrt_sum / nruns * 1e6)
    @printf("  norm(A) [current]:       %8.1f μs\n", t_raw_e11 / nruns * 1e6)
    @printf("  norm(Symmetric(A,:L)):   %8.1f μs\n", t_sym_e11 / nruns * 1e6)
end

run_bench()
