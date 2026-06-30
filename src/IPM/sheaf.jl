"""
    sheaf(I, J, V)

Construct a coboundary matrix from restriction map data.

Arguments:
- I: source vertex indices (one per restriction map)
- J: target edge indices (one per restriction map)
- V: restriction maps (each V[k] maps from vertex stalk I[k] to edge stalk J[k])

Returns B, the coboundary matrix as BlockSparseMatrix.
Sign convention: for edge e = (u, v) with u < v, B[e, u] = -F_u, B[e, v] = +F_v.
"""
function sheaf(I::Vector{Int}, J::Vector{Int}, V::Vector{<:AbstractMatrix})
    map = Dict{Tuple{Int, Int}, Int}(); e = 0

    CI = Int[]
    CJ = Int[]
    CV = Matrix{Float64}[]

    for (i, j, M) in zip(I, J, V)
        if !haskey(map, (i, j))
            if haskey(map, (j, i))
                map[i, j] = map[j, i]
            else
                map[i, j] = e += 1                
            end
        end

        if i < j
            sign = -1
        else
            sign =  1
        end

        push!(CI, map[i, j])
        push!(CJ,        i)
        push!(CV, sign * M)
    end

    return blocksparse(CI, CJ, CV)
end
