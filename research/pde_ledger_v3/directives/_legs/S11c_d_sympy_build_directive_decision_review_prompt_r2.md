# Independent physics review — S11c-d SymPy BUILD DIRECTIVE (decision review, ROUND 2)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_sympy_build_directive.md`

## What to check
Orchestrator-written **THIN SymPy BUILD DIRECTIVE** for step **S11c-d** (the G2 TRIGGER gate). Physics authority is
the CLEARED spec `directives/S11c_d_SHARED_PHYSICS.md` v10 (`399a8516`). Governs the **SymPy engine only**. The
directive must POINT at the spec (⛔ never restate/drift the physics) and fix only the build-mechanical layer + the
deferred §1c Fourier-reduction **element census against the real rows** + leak discipline. It is physics-bearing, so
**any finding that changes what is computed or what may be claimed is a must-fix**, not a nit.

**Round 1 found four must-fix defects, now folded.** Verify each fold is **correct, complete, and introduced no new
defect** — AND independently re-scan the whole directive (a fold can breed a fresh defect):

1. **F1 — the census was `VALUE`-only; each `(α,ρ)` case of both closed rows is a FIVE-slot payload.** The fold adds
   (§2 case-structure + §3 census) the slots `VALUE, MULTIGRADE, DIMENSION_L_T_M, COMPUTED_BRANCH_BINDINGS,
   FOURIER_PROFILE_BINDINGS`, and classifies the two extra slots: `FOURIER_PROFILE_BINDINGS` = c2's own 3-D hat
   definitions (transfer-only) → **not-to-bind** (same hazard class as `dtn_kernel`/a typed `(2π)` map; the engine
   computes its own reduction of the applied hats, incl. middle-leg); `COMPUTED_BRANCH_BINDINGS` = ONM 3-D roots at
   `k_out/k_in/k_mid` → reduce/carry at all three sites. **Verify against the real rows:** (a) is the five-slot
   payload real (each slot 8× = 4 cases × 2 rows)? (b) is the **classification physically right** — is
   `FOURIER_PROFILE_BINDINGS` genuinely a ready-made 3-D convention map that must NOT be bound, and does it really
   cover transfer only (so relying on it would leave middle-leg hats unreduced)? (c) is `COMPUTED_BRANCH_BINDINGS`
   correctly a reduce-at-all-three-sites operand? (d) **is any OTHER convention-bearing 3-D element, in ANY of the
   five slots, still MISSED** by the census? (e) does the directive avoid typing c2's numeric 3-D factor as the 1-D
   answer?

2. **Q3 — the candidate reduction map leaked through prohibitions.** The fold removed the typed factors
   (`[L_W/(2π)]`, `(2π)²L_W`) and the assembled `A_3D ≡ …` relation, genericizing to "a typed `(2π)`/normalization
   map." **Verify:** no typed reduction factor/relation remains anywhere (the flux-step's tangential `δ²(Q_∥)` is a
   legitimate method object to strip, per spec §1c — not a reduction answer; confirm that distinction holds). Does
   any remaining prohibition still name a real expected shape?

3. **Q4 — the withheld criterion named "the `O(1)` grating reductio".** The fold genericized §6 to "the falsification
   acceptance criterion (numeric magnitude bound), orchestrator-side," dropping the order token. **Verify:** the
   withheld criterion's value/order is not disclosed anywhere; leak discipline otherwise clean (c2 debt qualitative,
   no status counts, no expected value/sign/order/parity/grade/baseline as a builder target).

4. **Q7 — bound Riesz data was marked EMIT-only, contradicting the decision-list `resonances / local spectrum`
   handoff.** The fold moves `S11CD_BOUND_POLE_SET_AND_RIESZ_DATA` + `S11CD_BOUND_SPECTRAL_OVERLAP` into the EXPORT
   list (with the recursive closure + no-S11c-e-manifest caveat) and removes them from emit-only. **Verify:** the
   export membership now matches the declared boundary; nothing else that S11c-e binds is left emit-only, and nothing
   is exported that e cannot bind.

## What you are handed (no do-not-read list — what a leg must not use, it is not given)
- The artifact above.
- The physics authority: `directives/S11c_d_SHARED_PHYSICS.md` (source of truth).
- The real export files + loader (verify computationally, from `/var/projects/toy_physics/research/pde_ledger_v3`;
  Python has sympy): `scripts/S11c_b_exports.py`, `scripts/S11c_c1_exports.py`, `scripts/S11c_c2_exports.py`,
  `scripts/ledger_fold.py`.
- Precedents: `directives/S11c_c2_sympy_build_directive.md`, `directives/S11c_c1_sympy_build_directive.md`,
  `directives/S11c_decisions.md`.
- Skills + governance: `.claude/skills/build/SKILL.md`, `.claude/skills/review-legs/SKILL.md`, `CLAUDE.md`.
- The orchestrator's grounding (NOT a source of truth — re-derive yourself): the round-1 review record
  `directives/_measurements/S11c_d_sympy_build_directive_decision_review.md` and the census
  `directives/_measurements/S11c_d_sympy_build_directive_census.md`.

## Required method
**DOCUMENT that PINS FACTS — use both.** (1) Read the cleared spec first, form your own view of what it establishes /
defers to the build directive / withholds, then read the directive; quote spec-vs-directive for every finding.
(2) **Computationally verify** the pinned facts and the folds (⛔ a prose claim is discarded — run it; save the script
+ its literal stdout to a named absolute path; report those paths). Cover the same universe as round 1: the 3-parent
fold + `IMPORT_KEYS` rule + open-vs-closed provenance hazard; the whole-row Fourier/measure/DiracDelta census
**across all five payload slots**; the recipe-creep guard (WHAT not HOW); leak discipline; the HELD-PHYSICS
pointers; the script clauses + builder-lane; EMIT-vs-EXPORT.

## Physics filter
Report a finding only if it catches a way the **physics or the build could go wrong** — an incorrect/incomplete
census (any slot), a wrong/silent import binding, a leaked expected value, a misrepresentation/drop of a cleared-spec
control, recipe-creep, a broken builder-lane/script obligation, or a **new defect introduced by a round-1 fold**. For
each, state whether it changes what is computed or may be claimed (must-fix) or is a one-pass nit, and give
spec-vs-directive or command+stdout. If nothing survives, say the directive is CLEAR (sound).

## Ablation sandbox
Copy any file you introspect to `/tmp` and work on the copy; ⛔ never modify the working tree. Save every script + its
literal stdout to named absolute paths and report those paths.

## Bounds
Write your report and exit. ⛔ Do not spawn agents or supervisor orchestration. Finish IN-TURN, foreground blocking
only.
