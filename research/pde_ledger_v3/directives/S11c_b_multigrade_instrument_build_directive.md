# S11c-b — build directive: the background-order MULTIGRADE INSTRUMENT for the 20 genuine cross-engine cases

## Role — a MEASUREMENT, not an adjudication
The committed adjudication layer (`scripts/S11c_b_adjudicated_comparison.py`) routed exactly twenty core
cross-engine cases as genuine independent differences (twelve `FLAG`, eight `RESIDUAL_BULK`). This instrument
computes, for each of those twenty aligned operand pairs, the `(eta_bg, sigma_W)` background-bookkeeper
multigrade of the SymPy operand, the Wolfram operand, and their residual. It **states no conclusion** about
which engine is correct, which grades "should" be present, or whether either engine over- or under-retains.
It emits computed decomposition objects; the interpretation is written elsewhere, against the spec, by hand.

The two bookkeepers are `eta_bg` (the zero-jet contrast bookkeeper) and `sigma_W` (the first-background-jet
bookkeeper), both dimensionless SymPy symbols already present in the aligned operands after the layer's
Bridge D. The instrument treats every other symbol — `W_0`, `L_W`, `rho_br`, the profile jets, the trial
fields, coupling coefficients — as an atom.

## The twenty cases — exact keys (copy verbatim; do not add, drop, or re-derive the set)
Twelve `FLAG` core cases:
- `SLAB_OPERATOR_TERM_ORIGINS  (OBJECT=ADVECTIVE, BRANCH∈{LAB_HELD,MATERIAL_ADVECTED}, DENSITY∈{RHO4_CONSTANT,RHOBR_CONSTANT})`  — 4
- `SLAB_OPERATOR_TERM_ORIGINS  (OBJECT=KINETIC,   BRANCH∈{LAB_HELD,MATERIAL_ADVECTED}, DENSITY∈{RHO4_CONSTANT,RHOBR_CONSTANT})`  — 4
- `ADMISSIBILITY_OPERATOR_OPERAND (OBJECT=BODY_FORCE, BRANCH∈{LAB_HELD,MATERIAL_ADVECTED}, DENSITY∈{RHO4_CONSTANT,RHOBR_CONSTANT}, DOF=THETA)`  — 4

Eight `RESIDUAL_BULK` cases:
- `COUPLING_KERNEL (BRANCH∈{LAB_HELD,MATERIAL_ADVECTED}, DENSITY∈{RHO4_CONSTANT,RHOBR_CONSTANT}, SECTOR=TRANSVERSE_TO_THICKNESS)` — 4
- `COUPLING_KERNEL (OBJECT=ADJOINTNESS_OPERAND_FORWARD, BRANCH∈{…}, DENSITY∈{…}, SECTOR=TRANSVERSE_TO_THICKNESS)` — 4

⛔ The `FLAG SLAB_OPERATOR (OBJECT=DIVERGENCE_ROUTE_FIXTURE)` line the layer prints is an instrument
self-test fixture, **not** one of these twenty. Do not include it.

## Reuse the reviewed alignment — do NOT re-derive any physics (call the layer's own helpers)
The aligned operands already exist inside the committed, twice-reviewed layer. Consume them; do not rebuild
the energy, the operator, the coupling kernel, or the admissibility balance. Bind the layer as `A`, its
comparator as `C`, its engine module as `P` — the same three the layer imports.

- Parse the committed transcripts the layer parses by default: `C.load_py(C.DEFAULT_PY)` and
  `C.load_wl(C.DEFAULT_WL)`.
- Rename map: `active_names = dict(A.WL_TO_PY_RENAME)`. ⚠ This object lives on the adjudication layer
  (`A.WL_TO_PY_RENAME`, which is `H.WL_TO_PY_RENAME` from the handcoded module) — it does **not** exist on the
  comparator `C`. Do not invent one and do not import it from `C`.
- Build the case set with the layer's own case iterator (`A._family_cases(family, py_tags, wl_tags, state)`)
  for the families that contain the twenty keys, and select exactly the twenty keys above. Construct the
  `state` the same way `A`'s `main` does.
- For each case, build the aligned pair **exactly as the layer's `main` does — bridge_d is applied ONCE, by
  `_bridge_d`, not inside `transform`**:
  `pre = A.transform(case.operand_a, active_names, bridge_a=True, bridge_d=False, collapse=None)` and the same
  for `operand_b`; then `left = A._bridge_d(pre_a)`, `right = A._bridge_d(pre_b)`. ⛔ Do **not** pass
  `bridge_d=True` to `transform` (the layer never does; that would double-apply or diverge from the reviewed
  path). None of the twenty families are energy-gated, so this non-gated path is the one the layer routed.
- ⛔ Do **not** call `A._classify_case`, the Euler operator, or the divergence classifier — this instrument
  needs only the aligned operands, not the routing.
- The residual must be the SAME object the layer routed as `A_minus_B`:
  - For the three scalar families (`ADVECTIVE`, `ADMISSIBILITY_OPERATOR_OPERAND`, `COUPLING_KERNEL`) the
    aligned operands are scalar SymPy expressions; use `A._arithmetic_residual(left, right)`.
  - For `KINETIC` the operands are heteromorphic containers (SymPy `left` is a nested tuple; Wolfram `right`
    is a labelled `Association`), so they do **not** share container paths. Reproduce the routed residual with
    the layer's adapter: `pairs = A._kinetic_pairs(case.family, case.key, left, right)` (fail loudly if it
    returns `None`), then `A._kinetic_residual(pairs)`. The pairs carry the SEMANTIC component labels
    (`U_MOMENTUM_ROWS[0..2]`, `THICKNESS_ROW`) that align the SymPy and Wolfram components.

## The object to compute — the exact `(eta_bg, sigma_W)` multigrade on a fixed grade window
The operands are heterogeneous: some are scalar SymPy expressions, and the kinetic case is a container of
component expressions. **Every scalar leaf is graded against the same two symbols the operands actually
contain**: bind `eta_bg = P.eta_bg` and `sigma_W = P.sigma_W` (both created `real=True`). ⚠ A plain
`sp.Symbol("eta_bg")` is a *different* atom and would silently mis-grade every bookkeeper occurrence into the
`(0,0)` cell — assert `eta_bg is P.eta_bg` and `sigma_W is P.sigma_W` at import.

For a scalar leaf `expr`, define its multigrade over the fixed inclusive rectangular window
`a ∈ [0, N], b ∈ [0, N]` (`N` a module constant; set `N = 4`, comfortably above the observed operand degrees)
by **exact coefficient extraction — not `sp.series`** (whose truncation-order semantics are ambiguous and
whose expansion of a bookkeeper-in-denominator leaf is unbounded):

    c[a,b] = ( expr.diff(eta_bg, a).diff(sigma_W, b).subs({eta_bg: 0, sigma_W: 0}) ) / (factorial(a)*factorial(b))

This is exact and deterministic for every leaf, polynomial or rational (e.g. the coupling leaves that carry
`1/(1 + eta_bg*w1_profile)`). Every `c[a,b]` MUST be free of `eta_bg` and `sigma_W` (verify; a nonfree
coefficient means the wrong symbol identity or a bad extraction).

Then emit the **exact remainder** that the window does not capture — nothing is silently dropped:

    R = sp.cancel( sp.together( expr - Σ_{a,b∈[0,N]²} c[a,b] * eta_bg**a * sigma_W**b ) )

and emit `R`'s computed leading-grade information: scan `d[a,b] = R.diff(eta_bg,a).diff(sigma_W,b).subs(0)`
for `a,b ∈ [0, 2N]` and emit the set of lowest-total-order `(a,b)` at which `d[a,b] ≠ 0` (empty if `R == 0`).
For a correct extraction, `d[a,b] == 0` for every `(a,b)` inside the window `[0,N]²` — emit that fact as a
computed object (`R` has no window content), so a builder cannot hide coefficients by dumping the leaf into
`R`.

Emit, per case, three multigrade objects: `MULTIGRADE_A` (SymPy operand), `MULTIGRADE_B` (Wolfram operand),
`MULTIGRADE_RESIDUAL` (residual), each with its per-leaf remainder `R` and leading-grade set.

### Pinned output format (so two faithful builders emit differenceable stdout)
Serialise every payload with `C.serialise`. Use ONE nesting convention and ONE grade-key spelling:
- A multigrade is a `C.Association` keyed **first** by the leaf path, **then** by the grade pair.
- The leaf path key is: for a scalar-family operand, the single key `TextAtom("ROOT")`; for the `KINETIC`
  container, the `_kinetic_pairs` semantic label (e.g. `TextAtom("U_MOMENTUM_ROWS[0]")`, `TextAtom("THICKNESS_ROW")`).
  A/B/residual multigrades for a kinetic case MUST use these SAME pair labels as keys — never the raw SymPy
  tuple index or the raw Wolfram Association label — so their grades line up leaf-for-leaf.
- The grade-pair key is the string `f"{a},{b}"` (e.g. `TextAtom("1,0")`). The value is the coefficient object.
- Emit the remainder and leading-grade set under sibling keys `TextAtom("REMAINDER")` and
  `TextAtom("REMAINDER_LEADING_GRADES")` at the leaf-path level.

## Guards — PRINT the residual object, never assert it (the three clauses)
1. **A script PRINTS computed objects; it never states a conclusion.** Every emitted payload is a SymPy
   object or a plain path/grade key — never a sentence describing a result, never a verdict on an engine.
2. **PRINT the residual; do not assert it.** Emit each guard object; do not `assert` it zero and do not exit
   nonzero on it. Compute → emit.
3. **Interpretation belongs to the step record**, not to this script.

Emit these guard objects per leaf (computed, printed, not asserted) — each normalised so a correct
decomposition yields the exact zero object:
- `RECONSTRUCTION = sp.cancel(sp.together( leaf - Σ_{a,b∈[0,N]²} c[a,b]·eta_bg**a·sigma_W**b - R ))`.
- `WINDOW_CLEAN(a,b) = R.diff(eta_bg,a).diff(sigma_W,b).subs({eta_bg:0, sigma_W:0})` for every `(a,b)∈[0,N]²`
  (a correct remainder carries no window content, so each is the zero object).
- `GRADE_DIFFERENCE(path, a, b) = MULTIGRADE_A[path][a,b] − MULTIGRADE_B[path][a,b] − MULTIGRADE_RESIDUAL[path][a,b]`
  for every leaf path and `(a,b)` — normalise each; a correct residual makes each the zero object. Because
  the kinetic multigrades share the `_kinetic_pairs` labels, this differences matched leaves, not mismatched
  containers.
- `REMAINDER_DIFFERENCE(path) = sp.cancel(sp.together( R_A[path] − R_B[path] − R_residual[path] ))`.

These guards check the arithmetic of the decomposition itself (that coefficient extraction round-trips, that
the residual is the elementwise operand difference, and that the kinetic leaves are aligned). They reference
**no** expected physics value and **no** grade population — they are the analogue of a parser round-trip, not
an acceptance on the physics, and strengthening them to exact-zero introduces no rule-5 leak.

## The ONLY place physical symbols are combined by hand is upstream (already reviewed)
This instrument constructs no action, no ansatz, no operator. Every operand it touches is reached by
importing and calling the reviewed layer. It introduces no hand-typed CAS object standing in for a computed
one; its only new computation is coefficient extraction. A tag name here names the case, the leaf path, and
the grade key only — never a coefficient value, sign, or order.

## ⛔ Rule 5 — withhold every expected value and every classification
Do not encode, in the script or its tags, which grades `(a,b)` any operand or residual populates, how many
grades are populated, that any residual is "first order" or "higher order", that any operand is zero, that
either engine is missing/wrong/right, or any coefficient expression. The instrument discovers the grades by
computation and prints them. There is no physics acceptance criterion and none may be written: a genuine
difference in the emitted grades is the informative outcome and must not be tuned toward.

## Acceptance — value-free, mechanical, exact
- Runs to completion over all twenty cases (no case silently skipped; if a key is not found, fail loudly).
- Asserts `eta_bg is P.eta_bg` and `sigma_W is P.sigma_W` at import; every emitted coefficient is free of
  those two symbols.
- For each case emits `MULTIGRADE_A`, `MULTIGRADE_B`, `MULTIGRADE_RESIDUAL` (with per-leaf `REMAINDER` and
  `REMAINDER_LEADING_GRADES`), and the four guard objects.
- For the KINETIC cases, the A, B, and residual multigrades are keyed by the `_kinetic_pairs` semantic labels
  (the same leaf-path set across all three), so `GRADE_DIFFERENCE` differences matched leaves; an instrument
  that emits path-heterogeneous kinetic multigrades fails this clause.
- Every `RECONSTRUCTION`, `WINDOW_CLEAN(a,b)`, `GRADE_DIFFERENCE`, and `REMAINDER_DIFFERENCE` object
  normalises to exactly the zero object (checked mechanically on our side). These are algebraic
  self-consistency checks, not physics acceptances — an instrument that grades in the wrong variables, merges
  container components, or dumps a leaf into its remainder fails at least one of them.
- Uses the single pinned output format above (one nesting convention, `"a,b"` grade keys), so the stdout is
  deterministically differenceable.
- Confirms it consumed the layer's aligned operands (print the `A`/`C`/`P` module file paths and the
  transcript paths it parsed, so the reviewer can see it re-derived nothing).
- Pure SymPy; no Wolfram kernel. Parsing the two large committed transcripts dominates runtime (a few
  minutes); foreground is fine. ⛔ Never wrap in a shell `timeout`.

## Deliverable
`research/pde_ledger_v3/scripts/S11c_b_background_multigrade.py`. Its stdout is the deliverable; it is
captured and differenced on our side. The script does not decide anything.
