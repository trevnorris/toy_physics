# Independent physics review — S11c-b #88 blast-radius INSTRUMENT (Codex-written script)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_88_blast_radius.py`

A diagnostic SymPy instrument. For each strong DOF operator row (U_MOMENTUM_1/2/3, MU_THETA=θ,
THICKNESS_EW=e_W) it computes, at retained grade (η,σ_W ≤ 1), the difference between the
physically-correct variable-coefficient brane operator (the completed §3a spurion basis, 15
invariants/source, with the Euler–Lagrange total-divergence step differentiating the background
spurion — the "Hessian-retaining dx") and the operator the current engine emits (which freezes
the background at first order). It prints residuals, a three-driver decomposition, an
absorbability span, and four structural controls. It must NOT state any conclusion.

## What to check (the physics that must be right)
The instrument's job is to answer: *does correcting the frozen-spurion §3a basis (engine count
26 → correct 40) change each strong operator row, and is that change genuinely new operator
structure (non-absorbable into the frozen coefficients) rather than a relabeling?* A wrong
instrument would (a) mis-build the corrected or frozen operator so the residual is an artifact,
(b) let a re-parametrization masquerade as new structure (or vice versa), (c) freeze the Hessian
where it should retain it (or the reverse), (d) apply the retained-grade truncation wrongly, or
(e) silently compare the instrument against itself instead of the committed engine.

## What you are handed
- The artifact above.
- The engine it reuses: `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`
  (key: `dx`:616, `DERIVATIVE_MAP`:611-613, `operator_from_density`:1459-1511, `operator_dx`:1850-1874,
  `background_dx`:2122-2137, `first_shape_series`:713-725, `construct_energy`:1228, `enumerate_new_candidates`:1073,
  `live_basis_substitution`:1178, `basis_euler_signatures`:936, `quotient_independent_indices`:972).
- The build directive the artifact was built from:
  `research/pde_ledger_v3/directives/S11c_b_88_blast_radius_build_directive.md`
  (read it to know the intended object — but the artifact can satisfy the directive and still be
  wrong; your job is the physics, not directive-conformance).
- The spec `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` §1d, §3a.
- Settled context: `directives/_measurements/S11c_b_86_reference_result.md` (corrected basis = 40).
  ⛔ You are NOT given, and must NOT try to infer, the expected per-row answer — the settled
  cross-engine family verdicts are withheld on purpose.

## Required method — SCRIPT branch (derive independently; ablate; show literal stdout)
Copy the artifact to `/tmp` and work on the copy; never modify the working tree.

1. **Derive one row from first principles yourself, before trusting the artifact.** Pick one row
   (e.g. MU_THETA). Independently: build the completed density (frozen engine density + the
   omitted spurion invariants), compute the corrected EL with a Hessian-retaining dx and the
   frozen EL with the global dx, truncate with `first_shape_series`, and form the residual. Write
   your own script; save it + its literal stdout to named absolute paths and report them. Compare
   your residual's structure to the artifact's `RESIDUAL` for that row. A prose derivation is
   discarded — script + stdout or it did not happen.
2. **FORM ablations (mandatory — a coefficient rescale is not enough).** On the /tmp copy, make
   each of these structural changes and report the literal diff in the emitted tags:
   - Zero the second-background-jet map entries (make `hessian_dx` ≡ `frozen_dx`). Expected: the
     residual and all three drivers collapse toward 0, and `CONTROL_ENGINE` still passes (it uses
     the frozen dx). If the residual does NOT collapse, the Hessian retention is not what drives
     the result — report it.
   - Swap the two sources' live jets (`grad_W`↔`grad_mu`) in the completion, or zero one source.
     Expected: `CONTROL_JACOBIAN` bites (a source template loses its correct jet). If it passes
     silently, the Jacobian control is decorative.
   - Corrupt the driver reconstruction (drop one driver). Expected: `CONTROL_RECON` assert fires.
   - Perturb the frozen-density extraction so it no longer matches the engine. Expected:
     `CONTROL_ENGINE` assert fires. If it does not, the control is comparing the instrument to
     itself (tautological) — this is the highest-value thing to check.
3. **Is `CONTROL_ENGINE` genuinely independent?** Confirm it extracts the frozen row from the
   committed `operator_from_density`/`build_operator` output (ε-stripped, inertia-free) and
   asserts equality with the instrument's own `EL_FROZEN` — not `A − A`. Quote the lines.
4. **Absorbability span correctness.** The `S11CB_88_ABSORB_*` object must test constant-coefficient
   absorbability over the JOINT monomial basis (DOF jets AND background-profile jets/bookkeepers),
   with the scalar field being background-INDEPENDENT constants only. Check that background factors
   (`grad_W`, Hessian atoms, profile values, free coefficients) are NOT left in the matrix entries
   as the field — if they are, `RANK_GAIN` is meaningless (a first jet `g·q` and a Hessian `H·q`
   would show gain 0). Verify with a one-line CAS check on the artifact's own construction.
5. **No conclusion, no answer-target.** Confirm the only asserts are the four structural controls
   (`CONTROL_ENGINE`, `CONTROL_RECON`, `CONTROL_JACOBIAN` termcount) — none asserting a residual is
   nonzero, a rank gain is positive, or any per-row physics value. Confirm no tag PRINTS a verdict
   word or a hand-typed conclusion. Report any `assert` that precedes the value it guards.
6. **Retained-grade + Hessian-order.** Confirm the second background derivative enters at σ_W¹ (so
   it survives η,σ_W ≤ 1) and matches the committed `operator_dx`/`background_dx` normalization
   (`sigma_W·w1_profile_dij/L_W`, `mu_R·sigma_W·m1_profile_dij/(W0·L_W)`). Confirm the truncation is
   applied AFTER the dx/EL, not before (a truncate-then-differentiate order would drop real terms).

## Ablation-sandbox constraints (both legs identical)
- Copy the artifact to /tmp and ablate the COPY. Never modify the working tree.
- Save every ablation/derivation script AND its literal stdout to named absolute paths; report them.
- The artifact is pure SymPy (no Mathematica kernel); no timeout wrapper is needed, but if any run
  exceeds a few minutes, report it rather than waiting indefinitely.

## Physics filter
Report a finding only if it catches a way the instrument would answer the #88 question wrongly —
a mis-built operator, a control that cannot bite, an absorbability test that mis-classifies, a
frozen/corrected dx error, a truncation-order error, or a smuggled conclusion. Do not report
"it would be wrong on a different physics problem" or style.
