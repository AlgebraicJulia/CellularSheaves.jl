############################################################################################
# Edge Iterator
############################################################################################

struct SheafEdgeIter{T, I} <: AbstractEdgeIter
    S::Sheaf{T, I}
end

function Base.show(io::IO, iter::SheafEdgeIter)
    print(io, "SheafEdgeIter $(length(iter))")
    return
end

############################
# Abstract Graph Interface #
############################

function Graphs.edgetype(::Sheaf{T, I}) where {T, I}
    return SimpleEdge{I}
end

function Graphs.edges(S::Sheaf)
    return SheafEdgeIter(S)
end

#######################
# Iteration Interface #
#######################

function Base.iterate(iter::SheafEdgeIter{T, I}, (v, e)::Tuple{I, I} = (one(I), one(I))) where {T, I}
    S = iter.S

    vp1 = v + one(I)
    ep1 = e + one(I)

    while e ≤ narcs(S)
        @inbounds u = S.tgt[e]
        @inbounds enext = S.xsrc[vp1]

        @inbounds while enext ≤ e
            v = vp1; vp1 = v + one(I); enext = S.xsrc[vp1]
        end

        if u > v
            edge = SimpleEdge{I}(v, u)

            if enext ≤ ep1
                v = vp1
            end

            return (edge, (v, ep1))
        end

        e = ep1
        ep1 = e + one(I)

        if enext ≤ e
            v = vp1; vp1 = v + one(I)
        end
    end

    return
end

function Base.length(iter::SheafEdgeIter)
    return ne(iter.S)
end

function Base.eltype(::Type{SheafEdgeIter{T, I}}) where {T, I}
    return SimpleEdge{I}
end

function Base.in(edge, iter::SheafEdgeIter)
    return has_edge(iter.S, edge)
end
