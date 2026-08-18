function conedegree(cones::AbstractVector, B::BlockSparseMatrix)
    ν = 0

    for v in vtxs(B)
        ν += degree(cones[v], ncols(B, v))
    end

    return ν
end

function scalenorm(v::AbstractVector, s::AbstractVector)
    out = zero(promote_eltype(v, s))

    for i in eachindex(v, s)
        out += abs2(v[i] / s[i])
    end

    return sqrt(out)
end
