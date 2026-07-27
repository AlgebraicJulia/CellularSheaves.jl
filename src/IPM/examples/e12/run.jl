include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-25 (--quick --mosek)  [Mosek benchmarked DUAL — ~1.1x faster than primal]
# -----------------------------------------------------------------------------
# Z=4    dof=3040   n1=2364   blk=65    IPM   65.6ms  HSD   72.6ms (1.11x)  Cla  200.3ms (3.05x)  Msk  188.9ms (2.88x)
# Z=8    dof=7004   n1=5512   blk=65    IPM  167.5ms  HSD  179.9ms (1.07x)  Cla  528.3ms (3.15x)  Msk  501.6ms (2.99x)
# Z=16   dof=14932  n1=11808  blk=65    IPM  381.8ms  HSD  444.0ms (1.16x)  Cla 1151.6ms (3.02x)  Msk 1282.6ms (3.36x)
# IPM: DOF^1.11  HSD: DOF^1.14  Clarabel: DOF^1.10  Mosek: DOF^1.20  (dual shaves ~12% off primal at Z=16; sheaf still leads)
# =============================================================================
