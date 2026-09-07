# Independent physics review — the blind Wolfram N6 engine (a SCRIPT)

## Artifact
`research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl` (current working tree; a Codex-written blind
Mathematica engine, 1004 lines). It reproduces, independently and blind, the per-engine c2 N6 checks the SymPy side
established: (1) the **carrier reconcile** `C_E − C_M`, and (2) the **source-naturality** residual
`R_cov = ms − ms_pred` (the commuting square), at fixed anchoring `α∈{LAB_HELD, MATERIAL_ADVECTED}` and density
`ρ∈{RHO4_CONSTANT, RHOBR_CONSTANT}`, with able-to-fail FORM controls.

## What to check
Does this engine **faithfully COMPUTE** the two checks with genuine, able-to-fail controls — or does any load-bearing
result sit in a hand-typed payload, a tautology, or a check whose corruption does not move it? Specifically:
- **Carrier**: `C_E` (Eulerian pressure-slot coefficients) and `C_M` (material face-fold coefficients) are each reached
  by computation; `CARRIER_BRIDGE_RESIDUAL = C_E − C_M` is a genuine coefficient-level difference of two INDEPENDENT
  routes (Eulerian slab rows vs the material face builder with its internal covector map), ⛔ not one derived from the
  other. The three-way affine split `R_N6 = I(ΔC,ms) + B(C_M,ΔS) + B(ΔC,ΔS)` is exact (`SPLIT_CHECK` samplewise-zero),
  CARRIER uses `ms` (not `es`), SOURCE/CROSS are closed-response-only (no bare `−C·p`).
- **Source-naturality**: `ms = source(μ_M, V_M)` (material pull-back-then-vary) vs `ms_pred = source(μ_E∘Φ, V_E)`
  (the re-derived Eulerian `μ_E` with the prolonged field map Φ substituted — vary-then-substitute-Φ). NON-circular:
  `ms_pred` must NOT be built from the material pullback, from `ms`, or by subtracting Eulerian−material. Φ is prolonged
  to the jet order PRESENT in the re-derived `μ_E` (rank-2), via a live derivative chain, with a pre-substitution domain
  census. `R_COV_INCREMENT` is a closed-response contraction (⛔ not the affine `I`, which re-adds the bare term).
- **`μ_M` vs `μ_E∘Φ`** are two DISTINCT constructions (pull-back-then-vary vs the field-map substitution), kept
  separately parameterizable — conflating them makes `R_cov` tautological.
- **PIT**: an EXACT finite-field polynomial-identity probe (⛔ not a floating-point residual), both momentum legs
  sampled independently (`k_out ≠ k_in`), joint singular rejection, real-branch cell before modular reduction,
  per-cell coverage, a WL-derived δ = min(1, D/(N−E)) with bad-prime handling. ⛔ No residual-zero exit / assert.
- **No physics disposition / no VERDICT**; operands then residual; emission not conditional on a payload's value.
- **The three caveats** (Φ physical correctness; face-velocity transform; extracted/omitted-block leakage) are carried
  as open, not silently closed.

## What you are handed
- The `.wl` above.
- The physics you derive from independently: `directives/S11c_c2_SHARED_PHYSICS.md` §5c/§§1–2, the route-2 construction
  spec `_measurements/S11c_c2_N6_route2_spec_astra.md`, and the sibling specs (`S11c_a`/`S11c_b`/`S11c_c1`/`S11b`
  `_SHARED_PHYSICS.md`).
- The cleared SymPy instruments `scripts/S11c_c2_N6_reconcile_sympy.py` and `scripts/S11c_c2_N6_covariance_sympy.py`
  — the OBJECTS the WL engine reproduces. ⚠ The WL engine must be INDEPENDENTLY correct; ⛔ do NOT treat "matches the
  SymPy representation" as the test (they may legitimately differ in form). The test is whether the WL engine's own
  construction is correct and its controls bite.
- The build directive `directives/S11c_c2_N6_wl_build_directive.md` is the spec the engine implements — use it to see
  what was commissioned, ⛔ but derive the physics yourself; a code-to-directive match is not a physics check.

## Required method (this is a SCRIPT — derive independently, then ablate)
Form your own view of what `C_E − C_M` and `R_cov` should be from the physics BEFORE reading the engine. Then:

⛔⛔ **A FORM ablation is MANDATORY, not optional — it is the only thing that catches the worst defect.** Change the
STRUCTURE of a load-bearing object (flip a sign AND an off-diagonal, collapse two independent symbols into one), re-run,
and report the LITERAL diff. A COEFFICIENT rescale tests arithmetic; only a FORM change tests physics.

Run these one-sided corruptions (each in a /tmp COPY; report the literal stdout diff of the affected tag AND that the
others are unchanged):
1. **Carrier FORM knife** — corrupt the material covector/normal map feeding `C_M` (a FORM change), holding `ms`
   uncorrupted ⇒ `CARRIER_BRIDGE_RESIDUAL`/`CARRIER_CHANNEL` MUST move; `C_E`, the Eulerian operand, and
   `SOURCE_BRIDGE_RESIDUAL` MUST stay fixed. (If it does not move, the carrier reconcile is decorative.)
2. **Φ-coefficient knife** — make the ACTUAL material path use `2·a_ρ` (or drop `h_α`) while the prediction keeps the
   declared Φ ⇒ `R_COV` MUST move even though the `a_ρ+h_α` truth table is unchanged. Use `RHOBR_CONSTANT`/`LAB_HELD`.
3. **θ-independent junk knife** — add `κ_j·J_μ·e_W` to `μ_M` only (prediction unchanged) at
   `MATERIAL_ADVECTED.RHO4_CONSTANT` ⇒ `R_COV` MUST move.
4. **Non-circularity probe** — force `ms_pred` to be built from the material pullback (make the prediction circular) ⇒
   `R_COV` should collapse to structural zero (proving the shipped `ms_pred` is genuinely independent). If breaking the
   Eulerian route also moves the material route, they were never independent.
5. **PIT probe** — set `k_in = k_out` on an off-diagonal object and report whether an off-diagonal discrepancy is
   hidden; confirm the shipped δ is WL-derived (degree/exclusion from the circuit), not a transferred constant.

⭐ Ask of every emitted result: **WHICH LINE COMPUTED THIS?** Give the line number or report it as an uncomputed
hand-typed payload. Report any `assert`/exit that precedes the value it guards, any conclusion emitted as an
unconditional literal, any tautological residual (`A := B/C` then `A·C − B`), and any place a corruption to a
load-bearing map leaves the output byte-identical.

## Mathematica run discipline (both legs identical — these bind YOU)
⛔ Wrap EVERY kernel run in `timeout 600`. A 600 s hit is a FAILED ablation — report it and move on. ⛔ NEVER raise the
timeout, and ⛔ never run more than ONE kernel at a time (the licence has TWO seats and another leg may be running).
⛔ Copy the `.wl` to /tmp and ablate the COPY; ⛔ never modify the working tree. ⭐ Save every ablation script AND its
literal stdout to named absolute paths under /tmp, and report those paths. If a full 4-case run is too slow, ablate a
SINGLE case (per the knife's stated `(α,ρ)`) — a per-case run is ~500 s.

## Physics filter
Report a finding only if it catches a way the physics could be wrong (a control that cannot bite, a hand-typed result,
a circular or tautological residual, a frozen varying field, a wrong PIT that hides a discrepancy). ⛔ Do not report
"the engine would be wrong on a different input."

## Output
For each finding: the `.wl` line, what is wrong, the ablation script + its literal stdout showing it, and the minimal
fix. If a control genuinely bites (corruption moves its target, others fixed), say so with the literal diff. End with:
BUILD CLEAR, or DEFECTS with the blocking items.
