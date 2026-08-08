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
