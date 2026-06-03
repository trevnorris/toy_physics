---
unit_id: 236
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T22:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 236

## Per-finding outcomes

### F1 — missing_verification_script (subtype missing_mathematica)

**Classification:** resolved

**What changed:**
Codex added the new second-engine script
`mathematica/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_mathematica_audit.wl`
(330 lines, untracked) plus its committed output. The SymPy `.py` is unchanged
(no diff, no git status entry); only the SymPy `.txt` was refreshed by the
orchestrator re-run. The captured `stage_236_diff.patch` is empty because both
the `.wl` and its output are untracked (`git diff` shows tracked changes only),
not because nothing was written — the file is present on disk and ran to exit 0.

**Assessment — independence (genuinely different route, not a transliteration):**
- **M4 (the central independence requirement):** the `.py` hardcodes
  `L_rm_dep = [[0,-1,1],[0,-1,0]]` (line 89) and builds the projectors as
  `S_rm_dep · diag · L_rm_dep`. The `.wl` does NOT posit `L_rm_dep`. It takes the
  plane block `sRmDep[[{2,3},All]]` and *solves* the recovery with
  `LinearSolve[planeBlock, {dK,dMu}]`, cross-checks it against an independent
  `Solve[sRmDep·{uNt,uEta} == {0,dK,dMu}]`, then reconstructs the left-inverse
  matrix via `Coefficient` extraction (`matrixFromLinearAction`) and only afterward
  compares to the expected `L`. This is exactly the directive-prescribed alternative.
- **M3:** the `.py` forms `C_rm_dep = S_rm_dep · M_rm` (matrix product) vs a
  hardcoded matrix. The `.wl` builds the compiler by *acting on basis vectors*
  (`sRmDep·(mRm·#) & /@ IdentityMatrix[2]`, transposed back to a matrix) AND adds a
  product-vs-basis-action agreement check on a generic input. Different construction.
- **M5:** the `.py` checks idempotency/complementarity by repeating `P^2 - P = 0`
  matrix products. The `.wl` derives the projectors from the *solved* recovery map
  (`sRmDep·(selector·(lRecovered·deltaVec))`), checks idempotency/annihilation by
  *acting on a generic vector*, and additionally verifies them *spectrally*
  (`Eigenvalues` sorted to `{0,0,1}`, `MatrixRank == 1`) — the eigenvalue route the
  directive named as acceptable. Strictly broader, demonstrably different.
- Variable names (camelCase `epsEtaStar/qNt` vs snake_case), choreography
  (basis-action + LinearSolve + spectral), and assertion mechanism (coefficient
  extraction from symbolic linear actions vs matrix subtraction) all differ. Not a port.

The edit fully addresses F1; no collateral changes; directive `## Applied: F1`
reports `deviation: none` with a single `files_changed` entry, which matches disk.

## Exec log assessment

**SymPy:** exit=0. Notable lines: `Delta_mu - Delta_Keta = q_nt`;
`||y_eta||^2 = 2*q_eta**2`; `All Stage 236 symbolic checks passed.` (unchanged
script; refreshed output only).

**Mathematica:** exit=0. 37 `PASS:`, 0 `FAIL:`, no Solve/LinearSolve/FullSimplify
warnings or `$Failed`/`Indeterminate`. Notable lines:
`M4 LinearSolve and Solve recovery agree = {0, 0}` → PASS (independent recovery);
`Recovered {q_nt, q_eta} from {0, Delta_Keta, Delta_mu} = {-dK + dMu, -dK}`
(matches `L_rm_dep` solved, not asserted); `M5 P_nt eigenvalues` → PASS,
`M5 P_nt rank = 0` residual (rank 1, correct); `M8 Delta_mu - Delta_Keta residual = 0`
and `M8 Delta_mu equals Delta_Keta iff q_nt is zero = True`;
`||y_eta||^2 = 2*qEta^2`; `All Stage 236 Mathematica checks passed.`

PASS distribution maps cleanly to the manifest: M1:1, M2:2, M3:3, M4:4, M5:10,
M6:2, M7:3, M8:8, M9:4 = 37. Every M1–M9 item — including the static-strip
`Delta_mu − Delta_Keta = q_nt`, the iff-equivalence, the equal-drift ray
`{0,-E1,-E1} = (R1/c_eta){0,1,1}` (checked both via R1- and E1-substitution), the
three norm identities `2 q_eta^2 = 2 E1^2 = 2 R1^2/c_eta^2`, and the static-only
restoration gap (`y_rm + Delta_y_static = y_eta` ≠ 0, vs `y_rm + Delta_y_orbit = 0`)
— is covered by an explicit non-vacuous check. No silent parser-skip.

**Non-vacuous FAIL path:** confirmed. `fail[]` calls `Exit[1]`; `expectZero`
routes to `fail` whenever `zeroResidualQ` is false, and `zeroResidualQ` requires
every flattened entry to be `=== 0` after `cleanResidual` (FullSimplify∘Together∘
Expand, with ConditionalExpression stripping). Printed residuals are the
post-simplification values (`{0,0,0}`, `= 0`), so they reflect real reductions, not
`True` placeholders. `expectTrue` (M8) uses `Equivalent[...]` and would print a
non-True symbolic form on failure. A sign flip in `S_rm_dep` or the
`-1/(1-eps)` compiler entry would surface as a nonzero residual and exit nonzero.

**Output freshness:** confirmed. `.wl` mtime 1780460063 < its output `.txt`
1780460233 (output regenerated after the script). SymPy output `.txt` 1780460233 >
`.py` 1778522211 (refreshed). Log body byte-matches the committed
Mathematica output `.txt`.

## Material-change assessment

`material_change`: false. The fix only adds an independent verification engine;
it derives no new result, changes no constant, and the `.py`-derived quantities are
unchanged. Downstream units (>236) do not depend on anything the `.wl` produced.
No narrow or broad re-audit warranted on account of this stage.

## Side observations (non-blocking)

- The repo shows the full batch VII.2 (stages 231–242) in flight; for stage 236
  specifically the only touched files are the new `.wl`, its output, and the
  refreshed SymPy `.txt`. The stage-236 `.py`, `paper/`, and the stage-236 `notes/`
  file are untouched. (Modified notes are 231/232 — out of this unit's scope.)
- The `.wl` notes/prose self-labeling is N/A here (verifier is scripts-only); the
  known "Stage 253" renumber drift was confirmed by the auditor as a stale label,
  not a finding. Canonical 236 is ground truth.

## Verdict justification

The single finding F1 (missing Mathematica engine) is fully resolved by a new
`.wl` that reaches every M1–M9 conclusion via a genuinely independent route —
LinearSolve/Solve-derived recovery instead of a hand-asserted `L_rm_dep`,
basis-action compiler construction, and spectral (eigenvalue/rank) projector
checks — rather than transliterating the `.py`'s `S_rm_dep·diag·L_rm_dep`
choreography. Both engines exit 0, the FAIL path is non-vacuous, 37/37 PASS map to
the manifest with no silent skips, outputs are fresh, and there are no out-of-scope
or regressive edits. Verdict: verified; material_change: false.
