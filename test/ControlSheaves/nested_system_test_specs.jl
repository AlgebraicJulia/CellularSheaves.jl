# Shared `NestedSystemSpec` builders for the NestedSystems tower-compiler tests (Issue 009)
# and its downstream solver/DSL tests (Issues 010, 012). Defined once here rather than
# duplicated per test file.
#
# Assumes `CellularSheaves.ControlSheaves.NestedSystems` is already `using`-ed by the includer.

"""
    two_level_spec() -> NestedSystemSpec

Two flat teams, each observing its own target. `depth == 2`: `H_N` (raw agents + targets) and
`H_0` (each team collapsed to its `D`-dimensional centre) directly, no intermediate refinement.
"""
function two_level_spec()
    team1 = LeafTeam(:team1, :ring, 3, 1.0)
    team2 = LeafTeam(:team2, :ring, 3, 1.0)
    root = RefinedSystem(:root, AbstractSystemNode[team1, team2])
    targets = [TargetSpec(:t1), TargetSpec(:t2)]
    observations = [Observation([1], 1), Observation([2], 2)]
    return NestedSystemSpec(root, targets, observations, 3, true)
end

"""
    three_level_spec() -> NestedSystemSpec

A single subsystem, itself refined into two leaf teams: `root -> mid -> {teamA, teamB}`.
`depth == 3`.
"""
function three_level_spec()
    teamA = LeafTeam(:teamA, :ring, 3, 1.0)
    teamB = LeafTeam(:teamB, :star, 2, 1.0)
    mid = RefinedSystem(:mid, AbstractSystemNode[teamA, teamB], [(1, 2)])
    root = RefinedSystem(:root, AbstractSystemNode[mid])
    targets = [TargetSpec(:t1)]
    observations = [Observation([1, 1], 1)]
    return NestedSystemSpec(root, targets, observations, 3, true)
end

"""
    irregular_spec() -> NestedSystemSpec

An irregular tree: the left child is refined into two leaf teams, the right child is a bare
leaf team. `depth == 3` (the left branch needs one more pushforward step than the right).
"""
function irregular_spec()
    leftA = LeafTeam(:leftA, :ring, 3, 1.0)
    leftB = LeafTeam(:leftB, :path, 2, 1.0)
    leftchild = RefinedSystem(:leftchild, AbstractSystemNode[leftA, leftB], [(1, 2)])
    rightchild = LeafTeam(:rightchild, :clique, 3, 1.0)
    root = RefinedSystem(:root, AbstractSystemNode[leftchild, rightchild])
    targets = [TargetSpec(:t1), TargetSpec(:t2)]
    observations = [Observation([1], 1), Observation([2], 2)]
    return NestedSystemSpec(root, targets, observations, 3, true)
end

"""
    shared_target_spec() -> NestedSystemSpec

Two flat teams that both observe the *same* single target.
"""
function shared_target_spec()
    team1 = LeafTeam(:team1, :ring, 3, 1.0)
    team2 = LeafTeam(:team2, :ring, 3, 1.0)
    root = RefinedSystem(:root, AbstractSystemNode[team1, team2])
    targets = [TargetSpec(:shared)]
    observations = [Observation([1], 1), Observation([2], 1)]
    return NestedSystemSpec(root, targets, observations, 3, true)
end

"""
    two_target_team_spec() -> NestedSystemSpec

One refined system split into two rigid sub-teams wired by an internal consensus edge, with each
sub-team observing a *different* target (`D == 4`). The whole refined system collapses to a single
rigid vertex by `H_0`, so pulling the two targets apart is a demand the tower cannot meet: the
direct solve stretches the internal edge, the hierarchical solve cannot. This is the strict case
of the energy-gap inequality (Issue 010).
"""
function two_target_team_spec()
    teamA = LeafTeam(:teamA, :ring, 3, 1.0)
    teamB = LeafTeam(:teamB, :ring, 3, 1.0)
    mid = RefinedSystem(:mid, AbstractSystemNode[teamA, teamB], [(1, 2)])
    root = RefinedSystem(:root, AbstractSystemNode[mid])
    targets = [TargetSpec(:t1), TargetSpec(:t2)]
    observations = [Observation([1, 1], 1), Observation([1, 2], 2)]
    return NestedSystemSpec(root, targets, observations, 4, true)
end

"""
    flat_equivalent_spec() -> NestedSystemSpec

The degenerate case with no refinement to speak of: a single rigid team observing a single
target, `D == 4` affine. `depth == 2`, so the tower is just `H_N -> H_0` — structurally the same
two-level pipeline the flat `Layered` code implements. The answer here is known in closed form
(see the golden test), which makes it an oracle rather than a self-comparison.
"""
function flat_equivalent_spec()
    team = LeafTeam(:team, :ring, 5, 0.3)
    root = RefinedSystem(:root, AbstractSystemNode[team])
    targets = [TargetSpec(:t1)]
    observations = [Observation([1], 1)]
    return NestedSystemSpec(root, targets, observations, 4, true)
end

"""
    default_targets(spec::NestedSystemSpec) -> Vector{Vector{Float64}}

Distinct, well-separated boundary values, one per target of `spec`. When `spec.affine`, the last
coordinate is the homogeneous `1.0` row that affine restriction maps expect.
"""
function default_targets(spec::NestedSystemSpec)
    D = spec.D
    n_free = spec.affine ? D - 1 : D
    return [begin
                v = zeros(Float64, D)
                for i in 1:n_free
                    v[i] = t * (i == 1 ? 2.0 : -1.0 / i)
                end
                spec.affine && (v[D] = 1.0)
                v
            end
            for t in 1:length(spec.targets)]
end

"""
    degenerate_spec() -> NestedSystemSpec

An over-constrained/malformed team: a `:ring` formation with only one agent (`build_escort_topology`
requires at least two agents for a ring). The tower cannot represent this spec.
"""
function degenerate_spec()
    bad_team = LeafTeam(:bad, :ring, 1, 1.0)
    root = RefinedSystem(:root, AbstractSystemNode[bad_team])
    targets = [TargetSpec(:t1)]
    observations = [Observation([1], 1)]
    return NestedSystemSpec(root, targets, observations, 3, true)
end
