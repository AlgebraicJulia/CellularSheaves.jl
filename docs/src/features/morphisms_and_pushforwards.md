# Morphisms and Pushforwards

This page highlights the APIs for mapping sheaves across graphs and comparing sheaf data under those maps.

## When to use this page

Use this guide when you want to:

- Relate two graphs with an explicit map.
- Verify sheaf-level compatibility under graph maps.
- Transport sections and operators along a map.

## Graph Homomorphisms

- Main types and combinators: `GraphHomomorphism`, `compose`, `fiber_vertices`, `fiber_edges`, `cross_edges`, `graph_pushout`.
- See: [Graph Homomorphisms](../api/graph_homomorphisms.md).

## Sheaf Morphisms

- Core APIs: `SheafMorphism`, `ComplexMorphism`, `is_morphism`, `make_sheaf_morphism_spec`, `id`.
- See: [Sheaf Morphisms](../api/sheaf_morphisms.md).

## Pushforward and Transfer Maps

- Main APIs: `all_fiber_bases`, `fiber_section_basis`, `pushforward_sheaf`, `pushforward_transfer_map`.
- See: [Pushforwards](../api/pushforwards.md).

## Pushout Construction

- Main APIs: `SheafSpan`, `pushout_sheaf`.
- See: [Pushouts](../api/pushouts.md).

```julia
# Sketch: build a graph map f : G -> H and pushforward a sheaf.
# (Use your concrete graph and sheaf constructors from examples/tests.)
f = GraphHomomorphism(g_dom, g_cod, vmap, emap)

Sf = pushforward_sheaf(s, f)
T = pushforward_transfer_map(s, f)
```

Typical workflow:

1. Define a `GraphHomomorphism`.
2. Build or validate the induced sheaf map with `SheafMorphism` / `is_morphism`.
3. Compute `pushforward_sheaf` and `pushforward_transfer_map` for downstream computations.
