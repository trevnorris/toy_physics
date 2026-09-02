# Consultation: the S11c-b bulk-row jet-depth cross-engine disagreement — is it real, which engine is right, what next?

You are one of two independent advisors (the other is a different engine) asked to assess where an orchestrator's
cross-engine reconciliation stands and to give input on how to proceed. This is a dual-engine CAS ledger: a SymPy
engine (PY) and a blind Wolfram engine (WL) each build the same brane physics independently; a disagreement
between them is the measurement, not a bug to paper over. ⛔ Do NOT try to make the disagreement go away with
careful prose — if you think one engine is wrong, say which and prove it. ⛔ Do NOT pre-judge which engine is
correct because of its role.

## The object

The "brane slab operator" is the Euler-Lagrange (EL) operator of a brane energy density that is quadratic in the
dynamical fields (a displacement field `u` (3 comps), a width/thickness field `e_W`, a phase field `theta`) with
VARIABLE background coefficients that depend on a background width profile `W_bg(x)` and its spatial derivatives.
The "bulk momentum EL rows" (a.k.a. strong U-momentum rows / `U_INTERNAL` / `BULK_ENERGY` slot) are the momentum
sector of that operator from the bulk energy (excluding the separate kinetic and face/boundary contributions):

    row_a  =  dL/du_a  -  d_i ( dL/d(d_i u_a) )        (sum over i; d_i = spatial divergence)

`W_bg(x)` is kept LIVE (a genuine function of position) and differentiated; only AFTER the full EL + divergence is
the result reduced to the RETAINED GRADE. The retained grade keeps terms LINEAR in the shape parameter
`sigma_W ≡ dW/W_0` and linear in `eta_bg` (both truncated to power <= 1); higher powers are dropped. KEY FACT
(verify it): a background jet `d^n W_bg` of ANY order n carries exactly ONE power of `sigma_W` (order n differs
from order 1 only by powers of `1/L_W`), so the retained grade does NOT distinguish an order-2 (Hessian) jet from
an order-3 jet — a single jet of any order survives; a product of two jets is `sigma_W^2` and is dropped.

## The disagreement (measured; commands + literal outputs below — re-run and check them)

The two engines DISAGREE on the maximum background-jet order (order = number of spatial derivatives on `W_bg`,
i.e. `Total` of the jet multi-index) that appears in the retained-grade bulk momentum EL rows:

- **WL: order 3.**  10 order-3 atoms, e.g. `widthProfileJet[2,1,0]`, `[3,0,0]`, `[1,1,1]`.
- **PY: order 2.**  Order-3 is NEVER generated; the rows are identical whether PY is allowed depth 2, 3, or 4.

### PY measurement — `reconciliation_probes/tower_depth_probe.py` (+ `.out`)
Runs PY's own depth control (`S11c_b_brane_operator_sympy_audit.py:task_tower_depth_control`, `:3731`) in
isolation (it is skipped in the committed `.out` under PRIMARIES_ONLY) for case (LAB_HELD, RHO4_CONSTANT,
EULERIAN), building `retained_grade(live_strong_rows(construct_energy(branch, background_depth=d).density, …, d))`
for d in {1,2,3,4}. Literal stdout:
```
RESIDUAL_21_HESSIAN_SANITY_NONZERO_LEAF_COUNT: 5 of 5      (depth2-depth1: instrument is LIVE)
RESIDUAL_21_HESSIAN_SANITY_BG_JET_ORDERS_PRESENT: [1, 2]
RESIDUAL_32_DECISIVE_NONZERO_LEAF_COUNT: 0 of 5           (depth3-depth2: IDENTICAL)
RESIDUAL_43_TERMINATION_NONZERO_LEAF_COUNT: 0 of 5       (depth4-depth3: IDENTICAL)
UNREDUCED_DEPTH3_ORDER>=3_ATOMS (0): []                   (order-3 not even generated pre-reduction)
```
PY caps at `STRONG_ROW_JET_DEPTH = 2` (`sympy_audit.py:59`); the tower is built by `total_derivative` (`:696`),
the momentum row by `operator_from_density`/`euler_derivative` (`:1968`,`:752`). Note the EL applies a single
`total_derivative` (the `-d_i` term); a single `total_derivative` with depth>=2 maps a first jet -> Hessian but
does NOT map a Hessian -> order-3 (the `range(2, background_depth)` loop, `:725`).

### WL measurement — `reconciliation_probes/wl_provenance_probe.wl` (+ `.out`)
Loads the committed WL engine defs (reduced scale, `basisRepresentativeIndices={16}`), evaluates
`evaluatedModel["EULERIAN","LAB_HELD","RHO4_CONSTANT"]`, reports the engine's own `backgroundJetOrder`
(`mathematica…:1857`) per provenance slot. Literal stdout:
```
FULL_U_MOMENTUM_ROWS_MAX_JET_ORDER: 3
ORIGIN_KINETIC_U_MAX_JET_ORDER: 0        (0 jet atoms)
ORIGIN_BULK_ENERGY_U_MAX_JET_ORDER: 3    (all 10 order-3 atoms are HERE — the bulk EL slot)
ORIGIN_FACE_U_MAX_JET_ORDER: 0           (0 jet atoms)
FULL_MINUS_SLOTSUM_U (expect 0 vector): {0, 0, 0}
```
The order-3 sits entirely in the BULK slot (WL builds `U_MOMENTUM_ROWS = kineticU + rows["U_INTERNAL"] +
faceRows["U_FACE"]`, `mathematica…:1170`; the live EL is `constrainedRowsWithLiveEnergyEL`, reduced by
`finalBackgroundReduction`, `:1310`,`:1405`). So it is NOT an object-boundary artifact between the two engines'
"strong row" definitions.

### #89b context (cuts both ways — do not let it decide the answer for you)
#89b just FIXED a WL bug where WL wrongly FROZE the background (it reduced the operator while an outer
`Inactive[Div]` was still held, dropping mixed/higher jets); the fix keeps `W_bg` live and differentiates, which
RESTORED the order-3 jets. So either (a) WL's order-3 is the correct consequence of finally keeping `W_bg` live,
and PY's hard depth-2 cap is the SAME freeze un-fixed (a jet/derivative order held fixed is itself the finding);
or (b) the un-freeze over-restored an order-3 the physical retained grade does not carry, and PY's depth-2 is
right. The freeze-check `_measurements/S11c_b_89b_py_sibling_freeze_check.md` established PY reduces in the right
ORDER (reduce-last) but did NOT settle the DEPTH cap — that is this question.

## What we need from you

1. **Verify the measurements.** Re-run either probe, or derive the numbers your own way; report literal output.
   Do you agree it is a genuine disagreement located in the BULK momentum EL rows (not a definitional/boundary
   mismatch, not a reduction-order artifact)?
2. **Which engine's bulk-row jet depth is physically CORRECT — and prove it.** Derive, from the structure of a
   quadratic brane energy density `L(u, d u, W_bg(x))` with the retained grade above, the MAXIMUM order of a
   single background jet `d^n W_bg` that the bulk momentum EL row `dL/du - d_i(dL/d(d_i u))` must carry. The crux:
   what is the highest background-jet order carried by the `d_i u`-COEFFICIENTS in the energy density (call it p)?
   The row then reaches order p+1 via the EL divergence. Is p=1 (grad only -> row order 2, PY) or p=2 (Hessian ->
   row order 3, WL)? Ground this in what the brane energy density actually contains (you may read either engine's
   energy-basis construction: PY `construct_energy`/`retained_energy_basis`; WL `constructEnergyData` and the
   selected records), and say whether one engine drops or invents a background differentiation.
   ⛔ A prose "I derived it and got X" is worth nothing. If you compute, save your script AND its literal stdout to
   named absolute paths and cite them. If you argue by hand, show the term-level algebra.
3. **Recommendation.** Given where this sits: is PY's `STRONG_ROW_JET_DEPTH=2` a rule-17 freeze that must be
   raised (WL correct), or is WL over-generating (PY correct), or is it genuinely undecidable from here and needs
   a specific further computation (name it)? What is the single best next step?

## Constraints
- ⛔ If you run a CAS kernel (SymPy is free; Mathematica has a 2-seat licence — assume the other advisor may hold
  one, so wrap any `wolframscript` call in `timeout 600` and never run more than one kernel at a time). Copy any
  file you ablate to /tmp and ablate the copy; never modify the working tree.
- Save scripts + literal stdout to named absolute paths and cite them. Report findings as: agreement/disagreement
  on (1), your derivation + verdict on (2), and a concrete recommendation on (3).
- If you think the orchestrator's framing above is wrong, say so and why — that is the most useful thing you can do.
