function conedegree(cones::AbstractVector, B::BlockSparseMatrix)
    ν = 0

    for v in vtxs(B)
        ν += degree(cones[v], ncols(B, v))
    end

    return ν
end
