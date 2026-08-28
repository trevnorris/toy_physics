# Independent review — S11c-b cross-engine comparator (a build leg; the instrument, not the physics)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_cross_engine_comparator.py`
(with its synthetic fixtures `/var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S11c_b_cross_engine_comparator.py`).

This is a **comparator**, not a physics engine. It does not derive physics. It reads two blind engines'
already-committed tag streams and, for every `S11CB_*` family, keys each case on a typed axis key, applies a
fixed set of name/CAS folds, and prints `operand_A` (SymPy), `operand_B` (Wolfram), and `A − B` per joined
case, plus per-family accounting. **It must decide nothing** (no PASS/FAIL/VERDICT/target).

## What to check — the ONE failure mode that matters
A comparator earns its keep only if it **cannot manufacture false agreement and cannot silently hide a real
disagreement.** Every finding you report must catch one of these:
1. **A fold that MASSAGES two genuinely different operands into agreement** — e.g. a global
   `mu_theta_L/mu_theta_M → mu_theta` collapse, or *naming/folding a specific energy-basis representative
   identity*, so a real per-branch or first-background-jet difference is reconciled away and prints as `0`.
2. **A silent 0-join or silent drop** — a family (or a case within it) that extracts 0 joined cases, or an
   unpaired/duplicate/parse-failed/axis-mismatch case, without emitting an accounting line naming it.
3. **A mis-key** — two distinct physical cases colliding onto one key (a duplicate silently overwriting), or
   one case's operand compared against the wrong sibling because an axis was positional-guessed instead of
   typed.
4. **A residual that is zero by construction rather than by measurement** — `A − B` short-circuited, or an
   operand compared against something derived from itself.
5. **The comparator deciding anything** — any PASS/FAIL/VERDICT/FINAL_STATUS, any hard-coded zero/nonzero
   "expected" target, any emission gated on a payload's *value*.

The governing rule (do not drift from it): **basic mechanical pass (exact zeros + hand-verified renames) is
fine; everything else must be compared AS-IS and FLAGGED — never blanket-collapsed.** A `X(args) → X`
rename that strips real argument-dependence, or a canonicalizer that reduces an operand modulo a structure
that carries physics, MANUFACTURES agreement. That is the exact defect this leg exists to catch.

## The reconciliation folds the comparator is REQUIRED to get right (from the build brief)
- **μ_θ per-branch registry, argument-preserving**: `{LAB_HELD: mu_theta_L, MATERIAL_ADVECTED: mu_theta_M}`,
  consulted ONLY where the key pins the branch. ⛔ A global `mu_theta_L/M → mu_theta` collapse is the defect.
- **Energy basis is a non-unique QUOTIENT — compared RAW and FLAGGED.** The two engines may pick different
  basis representatives modulo total in-plane divergences. ⛔ The comparator must NOT name, assume, or fold
  any representative identity — a variable-coefficient IBP generates first-background-jet terms that are
  PHYSICS, so a pre-named collapse could mask them. A representative difference must SURVIVE as a documented
  `axis_set_mismatch`/nonzero residual, not be reconciled to 0.
- **Coupling-kernel ADJOINTNESS residual is a density MODULO compact-support in-plane IBP.** The comparator
  must surface BOTH the raw adjoint density AND its divergence-reduced form (or reduce before a zero-test).
  ⛔ It must NOT read a nonzero adjoint density as "non-adjoint" and it must NOT let the reducer zero a
  genuinely non-adjoint density: the real operator is non-self-adjoint via the dissipative `Λ_X` face
  response, and that must survive the reducer as nonzero.
- **CONTROL residuals compared raw** (UNIFORM_LIMIT, CONTROL_FORM, CONTROL_INDEPENDENCE, HOMOGENEITY,
  REP_INVARIANCE): each carries its own RESIDUAL operand; compare A/B/A−B, no target.
- **Axis-typed key** `(OBJECT, BRANCH, DENSITY, SECTOR, SOURCE, DIRECTION, DOF)`, each token TYPED (reject two
  values for one axis; never positional-guess). SECTOR must reconcile the PY/WL spellings
  `TRANSVERSE_TO_THICKNESS ≡ TRANSVERSE_TO_THETA_EW_UL` (and the mirror) as an axis rename, arg-preserving.

## Definition of done (check these EMPIRICALLY, by running — not by reading)
- Every emitted S11CB family prints an accounting line `{join, py_only, wl_only, duplicate_key, parse_failed,
  axis_set_mismatch}`, each with `join>0` OR a documented unpaired/mismatch + reason. **No family silently
  extracts 0.**
- It prints operand A, B, A−B **before** any guard; asserts nothing on the measured payloads (synthetic
  asserts live only in the separate test file); exits 0 on disagreement; nonzero only on operational failure.
- A `RUN_ACCOUNTING` summary line reports `families`, `families_with_join`, `families_with_unpaired`,
  `parse_failed`, `duplicate_key`.

## What you are handed
- The comparator + its test file (paths above).
- Its two read-only committed inputs:
  - PY: `/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_b_brane_operator_sympy_audit.out`
    (one-line `PY_S11CB_<Q>: <srepr>`).
  - WL: `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11c_b_brane_operator_mathematica_audit.out`
    (multi-line `WL_S11CB_<Q>` `<| … |>` associations with `Inactive[…]`).
- The build brief it was built to:
  `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_b_comparator_build_directive.md`.
- ⛔ You are NOT handed a "correct" comparison result. There is none. Do not invent an expected residual for
  any family; a genuine `A − B ≠ 0` is a finding for later adjudication, NOT a comparator bug.

## Required method — this is a SCRIPT; ABLATE, do not merely read
Code-reading alone has repeatedly missed real defects here. For every load-bearing fold, **ablate it and
report the literal diff in output.**

1. **Fast structural ablations on the synthetic fixtures** (the test file runs in <1s). For each fold below,
   corrupt it on a `/tmp` COPY of the comparator and re-run the fixtures; report which tests move and the
   literal diff:
   - **FORM ablation (MANDATORY, the only thing that has ever caught the worst defect):** change the
     *structure* of a load-bearing object, not a coefficient. Concretely, at minimum:
     (a) Force the μ_θ registry to the global collapse `mu_theta_L/M → mu_theta`. Does a case that SHOULD
         show a per-branch difference now residual to 0? If nothing moves, the registry is not actually
         branch-consulted — report it.
     (b) Insert a representative-identity fold on the ENERGY_BASIS operands (reduce them modulo a total
         in-plane divergence before comparison). Does a first-background-jet difference vanish? The comparator
         must have NO such reducer on the energy basis; confirm by grepping the real code path and by the
         fixture behaviour.
     (c) Break the SECTOR axis rename (map `TRANSVERSE_TO_THICKNESS` to a wrong sibling). Do coupling-kernel
         joins collapse to py_only/wl_only? That proves the rename is load-bearing and correct.
     (d) Corrupt ONE operand in one fixture (one-sided): confirm `A − B` moves and the other operand does
         NOT — i.e. the residual is an independent typed recursion, not zero-by-construction.
   - Report the literal before/after for each.
2. **One real end-to-end run** on the two committed `.out` files (the builder measured ~18 min wall; budget
   for it). Confirm empirically: every family prints its accounting line; the `RUN_ACCOUNTING` summary is
   present; no family shows `join=0` without a documented unpaired/mismatch reason; no `PASS/FAIL/VERDICT`
   token appears; exit code is 0. Save the run's stdout to a named absolute path and report the path.
3. **Adjointness reducer probe:** hand the reducer (i) a true total-divergence density (must reduce to 0) and
   (ii) a genuinely non-adjoint density with a `Λ_X`-like face term (must stay nonzero). Report both literal
   results. A reducer that zeros case (ii) is masking real physics — a finding.
4. **Independent re-key of a sample:** pick ~5 cases per family from the raw `.out` payloads and key them BY
   HAND from the axis tokens; confirm the comparator joined the same PY/WL pair. ⛔ Do NOT use the
   comparator's own `make_key` as ground truth — that is circular. A collision you find by hand that the
   comparator did not flag as `duplicate_key` is a finding.
5. For any claim you make, name the LINE that computes it (or report it as uncomputed). Report any `assert`
   that precedes the value it guards — an assert before the emit hides a defect by turning it into a crash.

## Physics filter
Report a finding only if it catches a way the COMPARISON could manufacture false agreement, hide a real
disagreement, silently drop/mis-key a case, or emit a verdict. Do NOT report "the comparator would be wrong
on a different input" or "engine A and engine B disagree here" (a real residual is expected and is not a bug).

## Ablation sandbox & operational constraints (identical for both legs)
- Copy the comparator to `/tmp` and ablate the COPY. ⛔ NEVER modify the working tree.
- This is pure Python — no Mathematica kernel, no licence seat. Run ablations **FOREGROUND**; do NOT use any
  background job or monitor loop; report directly when done. (A background-monitor leg has stalled here twice.)
- The one real run is ~18 min; run it once, foreground, and move on. Do not re-run it in a loop.
- Save every ablation script AND its literal stdout to named absolute paths, and report those paths. A prose
  "I checked and it's fine" with no script + stdout is discarded.

## Report
Separately, per fold: the ablation you ran, its literal diff, and your verdict (sound / manufactures
agreement / hides disagreement / silently drops / decides something). Then the real-run accounting summary
and its saved path. Then any unaddressed definition-of-done item. Keep it to findings, not narration.
