module ConicSheaves

export ConicSheaf, vertex_cone, vertex_cones, edge_vector, set_vertex_cone!,
    set_edge_vector!, optimization_problem

using Graphs
using AutoHashEquals: @auto_hash_equals
using LinearAlgebra
using JuMP
using MathOptInterface: AbstractVectorSet, Reals

using ..SheafInterface
import ..SheafInterface: vertex_stalks, edge_stalks, edge_stalk_dimensions,
    add_vertex_stalk!, add_sheaf_edge!, underlying_graph,
    get_vertex_stalk, get_edge_stalk, get_restriction_map
using ..EuclideanSheaves: UnorderedPair

@auto_hash_equals struct ConicSheaf{T} <: AbstractNetworkSheaf
    vertex_stalks::Vector{Int}
    edge_stalks::Dict{UnorderedPair{Int}, Int}
    underlying_graph::Graph
    restriction_maps::Dict{Pair{Int}, Matrix{T}}
    vertex_cones::Vector{AbstractVectorSet}
    edge_vectors::Dict{UnorderedPair{Int}, Vector{T}}
end

function ConicSheaf{T}(vertex_stalks::Vector{Int}) where T
    edge_stalks = Dict{UnorderedPair{Int}, Int}()
    underlying_graph = Graph(length(vertex_stalks))
    restriction_maps = Dict{Pair{Int}, Matrix{T}}()
    vertex_cones = AbstractVectorSet[Reals(d) for d in vertex_stalks]
    edge_vectors = Dict{UnorderedPair{Int}, Vector{T}}()
    return ConicSheaf{T}(vertex_stalks, edge_stalks, underlying_graph, restriction_maps, vertex_cones, edge_vectors)
end

function vertex_stalks(s::ConicSheaf)
    return s.vertex_stalks
end

function edge_stalks(s::ConicSheaf)
    return s.edge_stalks
end

function underlying_graph(s::ConicSheaf)
    return s.underlying_graph
end

function get_vertex_stalk(s::ConicSheaf, v::Int)
    return s.vertex_stalks[v]
end

function get_edge_stalk(s::ConicSheaf, v1::Int, v2::Int)
    return s.edge_stalks[UnorderedPair(v1, v2)]
end

function edge_stalk_dimensions(s::ConicSheaf)
    return [get_edge_stalk(s, src(e), dst(e)) for e in edges(underlying_graph(s))]
end

function get_restriction_map(s::ConicSheaf, v1::Int, v2::Int)
    return s.restriction_maps[v1 => v2]
end

function vertex_cones(s::ConicSheaf)
    return s.vertex_cones
end

function vertex_cone(s::ConicSheaf, v::Int)
    return s.vertex_cones[v]
end

function edge_vector(s::ConicSheaf, v1::Int, v2::Int)
    stored = s.edge_vectors[UnorderedPair(v1, v2)]

    if v1 < v2
        return stored
    else
        return -stored
    end
end

function add_vertex_stalk!(s::ConicSheaf, stalk_size::Int; cone::AbstractVectorSet = Reals(stalk_size))
    push!(s.vertex_stalks, stalk_size)
    push!(s.vertex_cones, cone)
    add_vertex!(s.underlying_graph)
    return length(s.vertex_stalks)
end

function add_sheaf_edge!(s::ConicSheaf{T}, v1::Int, v2::Int, rm1::AbstractMatrix, rm2::AbstractMatrix;
                         b::AbstractVector = zeros(T, size(rm1, 1))) where T
    stalk_size = size(rm1, 1)
    add_edge!(s.underlying_graph, v1, v2)
    key = UnorderedPair(v1, v2)
    s.edge_stalks[key] = stalk_size
    s.restriction_maps[v1 => v2] = Matrix{T}(rm1)
    s.restriction_maps[v2 => v1] = Matrix{T}(rm2)

    if v1 < v2
        s.edge_vectors[key] = Vector{T}(b)
    else
        s.edge_vectors[key] = -Vector{T}(b)
    end

    return ne(s.underlying_graph)
end

function set_vertex_cone!(s::ConicSheaf, v::Int, cone::AbstractVectorSet)
    s.vertex_cones[v] = cone
    return cone
end

function set_edge_vector!(s::ConicSheaf{T}, v1::Int, v2::Int, b::AbstractVector) where T
    key = UnorderedPair(v1, v2)

    if v1 < v2
        s.edge_vectors[key] = Vector{T}(b)
    else
        s.edge_vectors[key] = -Vector{T}(b)
    end

    return b
end

function optimization_problem(s::ConicSheaf{T}, Q, f, optimizer) where T
    n = sum(s.vertex_stalks)

    colptr = [0; cumsum(s.vertex_stalks)] .+ 1

    model = Model(optimizer)
    @variable(model, x[1:n])

    for j in eachindex(s.vertex_stalks)
        K = s.vertex_cones[j]
        jstrt = colptr[j]
        jstop = colptr[j + 1] - 1

        if !(K isa Reals)
            @constraint(model, x[jstrt:jstop] in K)
        end
    end

    for e in edges(s.underlying_graph)
        i = src(e)
        j = dst(e)
        istrt = colptr[i]
        istop = colptr[i + 1] - 1
        jstrt = colptr[j]
        jstop = colptr[j + 1] - 1
        be = s.edge_vectors[UnorderedPair(i, j)]

        @constraint(model, get_restriction_map(s, i, j) * x[istrt:istop] -
                           get_restriction_map(s, j, i) * x[jstrt:jstop] .== be)
    end

    if isnothing(Q) && isnothing(f)
        @objective(model, Min, 0)
    elseif isnothing(Q)
        @objective(model, Min, -dot(f, x))
    elseif isnothing(f)
        @objective(model, Min, 0.5 * x' * Q * x)
    else
        @objective(model, Min, 0.5 * x' * Q * x - dot(f, x))
    end

    return model, x
end

end
