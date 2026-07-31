# The alpha window: two laws and a controller

Status: candidate of record. Every law, constant, and controller component
below is either derived, measured, or replay-ablated against its alternative;
nothing is deployed, and nothing ships without the owner's decision.

Theory document. Supersedes the measurement-heavy floor_law.md / ceiling_law.md
as the statement of record; those remain as evidence appendices. All quantities
are decades (log10). All sensor inputs are computed by the solver during its
normal work; nothing here requires a probe, a sweep, or an extra solve.

## 1. Setting

Each IPM iteration solves KKT systems through the augmented operator
F = (1/alpha) H + BtB: factor F, then CRAIG to the inner tolerance
eps = max(force_tol, floor_tol). The usable alpha window is bounded by cost on
one side and arithmetic on the other:

  - FLOOR (low alpha): the solve becomes expensive -- the target deepens and
    iterations grind. Failure is gradual, announced by counts.
  - CEILING (high alpha): kappa(F) ~ alpha grows until finite precision cannot
    deliver the accuracy eps demands. Failure is abrupt, announced late.

The window is where one cheap solve suffices. Both edges move every iteration
(the floor rides the trajectory, the ceiling rides the tolerance), so a fixed
alpha is wrong eventually; the laws below turn each step's own byproducts into
calibrated distances to both edges.

## 2. The start-gap identity

The master fact, exact up to roundoff, for every role's base solve:

    r0 = (1/alpha) * || S (y0 - lambda*) ||        S = B F^-1 Bt

The initial CRAIG residual is the distance from the solve's WARM START to its
own true multiplier, deflated by alpha. Consequences:

  - every role's r0 falls one decade per decade of alpha (the 1/alpha
    prefactor is universal), so residual-decades convert to alpha-decades
    at slope 1: a residual gauge IS an alpha altimeter;
  - the numerator N = r0 * alpha names itself per role by the start:

        predictor  (cold start, y0 = 0)      N = ||S Dy||, the step's size
        corrector  (same-step start)         N = size of the centering
                                             correction (~1/3 of the step;
                                             coupled to the predictor)
        woodbury   (cross-step start, hsd)   N = inter-iteration drift of its
                                             solution (decoupled, volatile)

  - the identity is start-policy-contingent: warm-starting the predictor
    lowers the floor itself and breaks any sensor built on ||Dy||, but never
    breaks a sensor built on r0, which measures the operative gap whatever
    its source. The residual form is the invariant one.

The hsd bordered (3-row, tau) system is solved by elimination into exactly
these 2-row solves plus a scalar border lift, with r0 harvested pre-border and
the corrector warm-started from the predictor's pre-border solution. At the
level where the residuals live, hsd IS ipm; that is why one law serves both.

## 3. The floor law

The floor sits where the pending workload exceeds a fixed budget:
alpha_floor ~ N / eps. Measuring workload directly:

    mg   = max(log r0_p, log r0_c) - log eps      (decades still to grind)

    (mg deliberately includes the corrector: its start-gap is the centering
    correction's workload -- coupled to the predictor by the same-step warm
    start, binding the max on ~10% of rows. Corrector-skipped mode: mg <-
    log r0_p - log eps alone, constants recalibrated -- pre-registered
    prediction c -> +0.3, the universal budget -- plus a supplementary
    breakdown check, since the abstention gate's sensitivity lives in the
    corrector's refinement passes.)

    d_lo = -mg + c - 0.05 * relu(pbase - 1)       if state == 0 and npass <= 1
    ABSTAIN                                        otherwise

    c_ipm = 0.3      c_hsd = 0      (eta = 0.05, baseline 1 structural)

Readings:

  - slope 1: one decade of pending work costs one decade of floor margin
    (confirmed by tolerance-dial intervention; a mechanism, not a schedule
    artifact);
  - c is the workload budget the min-cost edge tolerates. It is
    tolerance-invariant and, measured with full role accounting, solver-
    universal at ~0.3: c_hsd = 0 because hsd's third (woodbury) solve spends
    ~0.3 of the budget invisibly to the two-role mg. The shipped-in-doc
    constant is a shrinkage estimator: the hidden workload's mean is absorbed
    rather than measured, because its per-row reading is drift-noise;
  - the count term is the exponent-drift meter: below the floor, spectral
    tail mass flattens r0's slope (r0 = sqrt((N/alpha)^2 + W_tail)), and each
    predictor iteration past the first counts one crossed-over direction --
    0.05 decades of descent per iteration, twenty per decade;
  - the gate: refinement passes firing twice is the ceiling's symptom; the
    floor gauge saturates there and abstains rather than guesses.

## 4. The ceiling law

The predictor's logged dual residual p_d is the breakdown quantity itself,
measured. Its structure over alpha:

    p_d = max( ramp, peg )
    ramp ~ u * alpha * sigma2max / h     (attainability: slope 1 in alpha)
    peg  ~ machine floor                 (the measurement's own roundoff)

The ceiling is where the ramp crosses eps -- refinement cannot beat its own
attainability floor (each pass must halve the residual and stalls when it
can't). Slope 1 again converts residual-decades to alpha-decades:

    d_hi >= log eps - log p_d + c2       if state == 0 and pbase <= 3
    ABSTAIN                               otherwise

    c2_ipm = -0.3     c2_hsd = -0.4

Semantics: ONE number, one-sided by construction. Because measured p_d is a
max, it never under-reads the ramp, so the emitted d_hi is exact when the
sensor is on its ramp and a valid (slack) lower bound when it has bottomed
out on the peg -- and no consumer needs to know which. There is no pinned
detector in the control path.

Two edges: the labeled ceiling is the ECONOMIC edge (counts leave the
plateau); p_d senses the HARD edge (arithmetic breakdown), which sits above
it by a variable gap whose mean c2 absorbs. Consequent bias profile: mild
optimism at the economic edge (costs iterations, not correctness); ipm turns
pessimistic past the edge (safe); hsd stays optimistic -- consumers treating
the edge as binding subtract a MARGIN of 0.5 for hsd. The margin is a margin,
not a sensor: above the ceiling, counts and p_d are redundant gauges of one
breakdown process, so no count-correction analog exists.

The gate: base-iteration grinding is the floor's symptom; it lifts p_d off
the peg for the wrong reason (the struggling solve's own amplified residual),
so the gauge abstains below the floor.

## 5. The mirror

    FLOOR:    primal residual ABOVE tolerance  ->  decades above the floor
    CEILING:  dual residual  BELOW tolerance   ->  decades below the ceiling

Same conversion (slope 1), same shape (margin + small budget constant), and
cross-gated: each law's ABSTAIN is triggered by the other law's home-territory
symptom (refinement passes ~ ceiling; base grinding ~ floor). An abstention is
therefore simultaneously a refusal and a referral to the other gauge. Coverage
is complementary: the floor gauge answers everywhere below its abstention (its
sensor rides 1/alpha over the whole range); the ceiling gauge answers exactly
(as equality) on its ramp and one-sidedly elsewhere.

Both laws err pessimistic in their hard regimes (endgame floors; far-above
ceilings on ipm), and the one anti-safe bias (hsd near-ceiling) carries an
explicit margin.

## 6. The controller

Because the laws are calibrated altimeters, control is deadbeat: read, place,
no probing. Changing alpha is free (F is refactored every iteration anyway).
Every component below was chosen by replay ablation against its alternative;
the record follows the loop.

    floor(step) -> d_lo | ABSTAIN
    ceil(step)  -> d_hi | ABSTAIN        (d_hi one-sided; may sit on its peg)

    CAP  = 1.0 (ipm) / 1.5 (hsd)         = 1 decade wanted margin + measured
                                           near-edge bias (0 / 0.5)
    controller(step, alpha):             # after every IPM step
        if state != 0 or floor == ABSTAIN:   return alpha / 10^2   # truck on
        F = alpha_log - d_lo
        if ceil == ABSTAIN:                  return 10^(F + 1.5)
              # ceiling abstains only below the floor: climb out
        H = alpha_log + d_hi
        return 10^( min( (F + H)/2,  H - CAP ) )

    There is NO peg detection anywhere: the ceiling's one-sided number feeds
    the midpoint directly. When the sensor sits on its machine-precision peg
    the number is a slack lower bound and the midpoint lands conservatively
    below true center; when on the ramp it is exact. Replay-verified: deleting
    the former pegged-fallback branch changes nothing fleet-wide except
    eliminating the one-gauge climbing phase (immediate vault instead).

    cold start: the library's own anchored default,
    alpha = aaug + raug*(||H||+||Q||)/||B||^2 (raug = 1e7); the first reading
    corrects even a six-decade-wrong anchor in one step.

Ablation record (replay on oracle fleets; each component isolated):

  - PLACEMENT = midpoint of [F, H - CAP]: beats fixed floor-offset by ~10% on
    narrow/tilted windows (the floor edge carries a universal one-iteration
    lip by the window construction; midpoint clears it), ties elsewhere;
    simpler (no target constant). Floor-offset survives only as the fallback
    when the ceiling gauge is pegged or absent (F + 1.5).
  - CAP = bias + 1: realized margin = nominal cap minus the solver's
    near-edge bias, decade for decade (measured in closed loop); cap size is
    otherwise free in steady state. 0.5 realizes ZERO margin on hsd; 2+
    buys nothing more. No separate economic setback: the calibrated H is
    already economic-edge-referenced.
  - ENDGAME LATCH DEMOTED TO DIAGNOSTIC (second vestigial mechanism caught
    by re-audit after the midpoint adoption; the first was the pegged
    fallback). Its rationale -- stop chasing collapsing floors -- belonged
    to the offset policy; midpoint never chases the floor, and the gauges
    track a collapsing window on their own (replay: identical descent to
    within one decade, latch cost +4 for no reduction in violations or
    breaches). As an actuator it also carried unbounded-retreat risk on any
    false fire, being latched. The detector survives as a report only: the
    count-rise signal (one true fire, zero false on the fleet) announces
    approaching precision exhaustion. The drift-innovation alternative had
    already been rejected (false-fired on every clean problem).
  - RECOVERY = truck on: a breached step is degraded, not repeated; the
    correction applies to the next iteration. Mid-step re-solves (REDO) were
    struck for implementation cost before testing.
  - DROPPED FOR CAUSE (36-cell ablation, all cells identical): informed
    descent on floor-abstain (the band where it could beat a blind drop is
    1-2 decades wide -- one drop crosses it; above it state=1 blinds both
    gauges anyway); per-iteration jump cap and deadband (big jumps happen
    only off the first reading into wide early windows, where the laws are
    accurate; settled placements are naturally small). Carry-and-drift died
    with its consumers; stall-perch died with carry; both-abstain needs no
    special case (same drop).
  - PEGGED-FALLBACK DELETED (caught as an inconsistency with the law's
    bound semantics): bound-as-location is safe by construction and replay-
    identical; the peg concept survives nowhere in the control path.
  - NO positional policy survives a one-iteration total window collapse
    (terminal cliff); the endgame retreat softens entry, truck-on recovers.

Diagnostics beyond alpha: both ceiling edges (economic from counts, hard from
the ramp) and the window width d_lo + d_hi -- width shrinking toward ~2
decades is an early precision-exhaustion warning.

Deliberately absent from v1: the Ritz/quadrature cost model (r0 IS the exact
spectral quadrature at the current alpha; the ten-node caricature earns its
keep only for counterfactual alphas -- placement optimization and foresight,
registered as v2), probing, gradient estimation, any memory at all
(the endgame detector reports but does not steer), and the ritz_beta
tripwire (registered add-on).

## 7. Constants of record

    quantity                          value      status
    floor budget c (ipm)               0.3       calibrated; = universal budget
    floor budget c (hsd)               0.0       = universal minus hidden
                                                 woodbury workload (~0.3)
    descent rate eta                   0.05      shared; 20 iters per decade
    count baseline                     1         structural
    floor gate                         npass<=1  structural
    ceiling offset c2 (ipm/hsd)     -0.3/-0.4    calibrated (mean two-edge gap)
    ceiling gate                    pbase<=3     chosen; alternative: floor
                                                 law's own d_lo >= 1
    controller CAP (ipm/hsd)        1.0/1.5      replay-derived: bias + 1
    pegged fallback offset             1.5       chosen (floor-referenced)
    truck-on drop                      2         chosen
    endgame detector (diagnostic only) rising, min+2   1 true fire, 0 false;
                                                 reports, never steers

All calibrated constants are round to within the label noise; none is
load-bearing to the second decimal. Everything above is a candidate pending
fresh-seed certification and the registered alternatives; nothing is decided
here beyond the replay evidence stated.

## 8. Assumptions and domain

  - Predictor solves start cold (the identity's N = step size, and the floor
    labels themselves, are start-policy-dependent).
  - Corrector present (the floor gate's sensitivity lives in the corrector's
    refinement passes; corrector-skipped operation needs a supplementary
    breakdown check and recalibration -- the ceiling sensor itself is
    predictor-side and survives).
  - Spectrum compressible; X04-class (near-dependent constraint rows)
    problems are outside both laws' domain and belong to presolve. v1 does
    NOT police this assumption -- on such problems both laws are confidently
    wrong and the controller will follow them. A ritz_beta-based tripwire is
    the registered add-on for detecting the regime in-flight.
  - Production tolerance schedule for the constants; the slope-1 forms are
    schedule-independent (intervention-tested).
  - Endgame iterations: floors collapse benignly; the laws go pessimistic
    and the controller stops descending.

## 9. Open register

Pending, in rough order of leverage:

  - FRESH-SEED CERTIFICATION: both laws' forms were selected iterating on the
    oracle fleets; constants are stable but the forms carry selection debt.
    One scoring pass on fresh seeds converts every number above from "lightly
    flattered" to final. Registered challengers ride along: the norm_dy floor
    sensor for ipm (beat the residual form by 0.035 in CV; cold-start-
    contingent by the start-gap identity), the sum-combiner for multi-role
    margins, the ceiling below-zone resolution (domain semantics vs cross-
    gate vs instrumentation), and hsd ceiling gate npass<=2.
  - CORRECTOR-SKIP RUNS: sharp pre-registered prediction -- single-role mg on
    the corrector-skipped ipm recalibrates to c ~ +0.3, the universal budget,
    in a third role configuration. Also required before that mode ships: a
    supplementary breakdown check (the floor gate's sensitivity lives in the
    corrector's refinement passes).
  - INSTRUMENTATION (one regen): diag(H)/diag(BtB) ratio columns (deterministic
    pinned-regime ceiling intercept) and the pre-border ||Dy|| for hsd (makes
    the hsd norm-law well-posed for the first time; validity derived from
    source before any data existed to fit it on).
  - REMAINING FORCING LEGS: the 0.3 dial and the endgame-onset prediction
    (moves earlier with larger forcing_frac iff the crossover is
    accuracy-channel).
  - OPEN QUESTION: per-problem spread of the floor budget c (0.0 .. 3.0;
    not w-bookkeeping, not tolerance). Last candidate standing: c tracks the
    decades one marginal iteration buys at the edge (contraction-rate
    property).
  - V2 REGISTER: the Ritz/quadrature cost model for placement optimization
    and foresight (counterfactual alphas; the exact integral is already
    consumed at the current alpha); the ritz_beta dense-spectrum tripwire;
    a persistence-hardened drift-innovation endgame detector, should the
    count detector ever prove insufficient.

## 10. Fleet replay summary (the controller's evidence)

On the replay fleet, the five-line controller: matches or halves the cost of
the anchored fixed default (26 vs 53 on the rising-floor exemplar), holds
96-100% window residency wherever windows exist, recovers from a six-decade
bad start in two reflexive drops and one vault (total penalty ~18 CRAIG
iterations), retreats ahead of endgame collapse when its latch fires (one
true fire, zero false, entry overshoot halved), and absorbs the one
unavoidable event -- a one-iteration total window collapse -- by trucking on.
Its three verbs: cruise, retreat, recover. Everything it does is the two
gauges plus five lines; everything it refuses to do, a gauge refused first.
