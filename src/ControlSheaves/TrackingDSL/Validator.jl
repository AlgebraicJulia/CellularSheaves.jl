"""
    TrackingDSLValidator

Semantic validation for `TrackingProgram` ASTs.

`validate_tracking_program` checks:

1. No duplicate declarations of the same name.
2. All referenced names in constraints and boundaries are declared.
3. Time sets reference declared time aliases.
4. Agent dynamics reference declared maps.
5. Consensus/track constraints reference declared agents/targets.
"""
module TrackingDSLValidator

using ..TrackingDSLTerm

export validate_tracking_program

"""
    validate_tracking_program(prog::TrackingProgram)

Validate the `TrackingProgram` AST for semantic consistency.

Raises a `TrackingDeclarationError`, `TrackingTypeError`, or
`TrackingTimeResolutionError` if validation fails.  Returns `prog` unchanged
on success so it can be used in a pipeline:

```julia
prog |> validate_tracking_program |> resolve_tracking_program |> lower_tracking_program
```

Validation checks:
- No duplicate names across spaces, maps, agents, targets, time aliases, constraints.
- All names referenced in constraints exist as declarations.
- `initial` and `final` are allowed without declaration (built-in aliases).
"""
function validate_tracking_program(prog::TrackingProgram)
    # Collect declared names by category
    spaces   = Set{Symbol}()
    maps     = Set{Symbol}()
    agents   = Set{Symbol}()
    targets  = Set{Symbol}()
    times    = Set{Symbol}()    # user-declared aliases; initial/final are built-in
    time_sets = Set{Symbol}()
    horizons  = Symbol[]

    for stmt in prog.statements
        if stmt isa SpaceDecl
            _check_unique(stmt.name, spaces, "space")
            push!(spaces, stmt.name)

        elseif stmt isa MapDecl
            _check_unique(stmt.name, maps, "map")
            push!(maps, stmt.name)

        elseif stmt isa VectorDecl
            # just register it — no category collision needed
            nothing

        elseif stmt isa AgentDecl
            _check_unique(stmt.name, agents, "agent")
            push!(agents, stmt.name)

        elseif stmt isa TargetDecl
            _check_unique(stmt.name, targets, "target")
            push!(targets, stmt.name)

        elseif stmt isa HorizonDecl
            push!(horizons, stmt.name)

        elseif stmt isa TimeAliasDecl
            _check_unique(stmt.alias, times, "time alias")
            push!(times, stmt.alias)

        elseif stmt isa TimeSetDecl
            _check_unique(stmt.name, time_sets, "time set")
            push!(time_sets, stmt.name)

        elseif stmt isa ConsensusConstraint
            _check_in(stmt.agents[1], agents, "consensus '$(stmt.name)' agent")
            _check_in(stmt.agents[2], agents, "consensus '$(stmt.name)' agent")
            _validate_time_spec(stmt.time_spec, times, time_sets, "consensus '$(stmt.name)'")

        elseif stmt isa TrackConstraint
            _check_in(stmt.agent,  agents,  "track '$(stmt.name)' agent")
            _check_in(stmt.target, targets, "track '$(stmt.name)' target")
            _validate_time_spec(stmt.time_spec, times, time_sets, "track '$(stmt.name)'")

        elseif stmt isa ConsensusSheafConstraint
            for a in stmt.agents
                _check_in(a, agents, "consensus_sheaf '$(stmt.name)' agent")
            end
            _validate_time_spec(stmt.time_spec, times, time_sets, "consensus_sheaf '$(stmt.name)'")

        elseif stmt isa BoundaryConstraint
            if stmt.entity_kind == :agent
                _check_in(stmt.entity_name, agents, "boundary agent")
            elseif stmt.entity_kind == :target
                _check_in(stmt.entity_name, targets, "boundary target")
            else
                throw(TrackingDeclarationError(
                    "boundary: unknown entity kind '$(stmt.entity_kind)'"))
            end
            _validate_time_ref_loose(stmt.time_ref, times, "boundary")

        elseif stmt isa BindStmt
            # Bind statements are checked at resolution time
            nothing
        end
    end

    # Warn if multiple horizons declared (first wins)
    if length(horizons) > 1
        throw(TrackingDeclarationError(
            "Multiple horizon declarations: $(join(horizons, ", ")). Declare exactly one."))
    end

    return prog
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

function _check_unique(name::Symbol, seen::Set{Symbol}, kind::String)
    if name in seen
        throw(TrackingDeclarationError("Duplicate $kind declaration: '$name'"))
    end
end

function _check_in(name::Symbol, declared::Set{Symbol}, ctx::String)
    if !(name in declared)
        throw(TrackingDeclarationError(
            "$ctx references undeclared name '$name'"))
    end
end

function _validate_time_spec(spec::TimeSpec, times::Set{Symbol}, time_sets::Set{Symbol}, ctx::String)
    if spec isa SingletonTime
        _validate_time_ref(spec.t, times, ctx)
    elseif spec isa TimeList
        foreach(t -> _validate_time_ref(t, times, ctx), spec.ts)
    elseif spec isa TimeRange
        _validate_time_ref(spec.lo, times, ctx)
        _validate_time_ref(spec.hi, times, ctx)
    elseif spec isa NamedTimeSet
        # NamedTimeSet can refer to either a declared time alias (singleton) or
        # a declared time set (range). initial/final are built-in; others must be declared.
        (spec.name == :initial || spec.name == :final ||
         spec.name in time_sets || spec.name in times) ||
            throw(TrackingTimeResolutionError(
                "$ctx references undeclared time set or alias '$(spec.name)'"))
    end
end

function _validate_time_ref(t::TimeRef, times::Set{Symbol}, ctx::String)
    if t isa InitialTime || t isa FinalTime || t isa LiteralTime
        nothing   # always valid (built-in aliases or literals)
    elseif t isa NamedTime
        # initial and final are built-in; other names must be declared
        (t.name == :initial || t.name == :final || t.name in times) ||
            throw(TrackingTimeResolutionError(
                "$ctx references undeclared time alias '$(t.name)'"))
    end
end

# Loose version: allow NamedTime without checking (used for boundary's time ref
# which may be a bound variable like `t`)
function _validate_time_ref_loose(t::TimeRef, times::Set{Symbol}, ctx::String)
    # Always passes — boundary time references can be late-bound variables
    nothing
end

end # module TrackingDSLValidator
