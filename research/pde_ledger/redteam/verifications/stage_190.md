---
unit_id: 190
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-01T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 190

## Per-finding outcomes

### F1 — missing_mathematica (dual-engine gap)

**Classification:** resolved

**What changed:**
A new independent-route Mathematica audit was created at
`mathematica/moving_throat_pde_stage190_direct_defect_vs_dressing_split_mathematica_audit.wl`
(283 lines), covering all five load-bearing sections of the SymPy reference
(`scripts/moving_throat_pde_stage190_direct_defect_vs_dressing_split_sympy_audit.py`).
The SymPy script was NOT modified (mtime 1778522333, unchanged; only its `.txt`
was re-run for the freshness gate). The committed Mathematica output
(`mathematica/output/...mathematica_audit.txt`) shows 41 PASS / 0 FAIL and the
exec log ends `# exit_code: 0`.

**Assessment:** Correct, substantive, and a genuinely independent route. The
`.wl` re-derives — rather than re-types — every conclusion via Mathematica-native
machinery on a different decomposition:

- **Sec II (slippage ledger):** strongest independence signal. SymPy uses a
  one-parameter multiplicative drift (substitute each primitive
  `x -> x*exp(s*x1)`, then `d/ds log(expr)|_{s=0}`). The `.wl` instead applies a
  logarithmic Euler operator,
  `Total[MapThread[#2*#1*D[Log[monomial], #1]&, {vars, drifts}]]`
  (`logEuler`, lines 127–134) — Euler's homogeneous-function operator
  `Σ_i drift_i · x_i · ∂log f/∂x_i`. Different mechanism, same five closed forms.
- **Sec III (direct-defect split):** SymPy postulates the closed forms `Atr`,
  `Sigma_nt` and checks the identity. The `.wl` *derives* the projection
  coefficient via `Coefficient` extraction —
  `atrProjected = Coefficient[xiDirect, sigmaChi]/Coefficient[trackingCoordinate, sigmaChi]`
  (line 167), then forms `sigmaNtProjected = xiDirect - atrProjected*trackingCoordinate`
  (line 168) and independently confirms `Coefficient[sigmaNtProjected, sigmaChi] == 0`
  (line 182, a check absent from SymPy). Genuine reconstruction, not transliteration.
- **Sec IV (compiler):** SymPy hand-writes the three inverse-reconstruction closed
  forms. The `.wl` inverts the 3×3 triangular compiler natively with
  `LinearSolve[compiler, {thetaVar, xiVar, dressVar}]` (line 207) and compares
  each component to the target; adds a `Det`-based determinant-formula check
  (line 202) not present in SymPy.
- **Sec V (no-go filter):** SymPy uses explicit scalar combinations. The `.wl`
  uses matrix projectors (`projectors.xLane`, lines 241–249) and a
  `DiagonalMatrix[{4, 4/5}]` invariant kernel with vector·matrix·vector products
  (lines 257–259), and checks the no-linear-term claims via
  `Coefficient[Normal[Series[...]], epsAx, 1]` (lines 267–274) rather than
  `diff(...).subs(eps,0)`.
- **Sec I (support-blindness):** the structural rebuild
  (`loadedTransfer = lambda0*(1-epsEta)/loadedTarget`) mirrors SymPy's
  decomposition — this is the one section closest to a port — but it is NOT a bare
  re-type: it adds a new mass-factorization check `M_tr - M_mix*S(zeta;eps)`
  (line 94) and replaces SymPy's weak "spoiled residual ≠ 0" guard with a stronger
  **exact numerical witness**, `(spoiledResidual /. {eps->1/3, zeta->1/2}) + 46/133 == 0`
  (lines 114–117). I independently confirmed the spoiled residual evaluates to
  exactly `-46/133` at that point, so this is a real falsifiable assertion, not X−X.

No collateral edits: no other script files were touched (`stage_190_diff.patch`
is empty because the `.wl` is a brand-new file).

## Exec log assessment

**SymPy:** exit=0. Reference engine, unchanged. Output reproduces all section
residuals as `0` (e.g. `Xi_direct - (A_tr Sigma_tr + Sigma_nt) = 0`,
`I[x,x] - 7/10 eps^2 x1^2 = 0`).

**Mathematica:** exit=0 (`# exit_code: 0` in `stage_190_mathematica.log`). 41 PASS,
0 FAIL. Notable lines:
- `PASS: Sigma_chi - (gamma1 + c1 - kappaU)` (Euler-operator route)
- `PASS: projected A_tr` and `PASS: nontracking residual has no Sigma_chi`
  (Coefficient-projection route)
- `PASS: Sigma_tr inverse by LinearSolve` (native matrix inversion)
- `PASS: spoiled exact witness at eps=1/3,zeta=1/2` (exact-value control)
- `PASS: I[x,x] - 7/10 eps^2 x1^2` (projector/kernel route)

**Output freshness:** confirmed. `.wl` mtime 1780335170; `.wl` output mtime
1780335323 (newer → regenerated post-creation). SymPy output mtime 1780335323
(re-run for the gate; the `.py` itself at 1778522333 is untouched).

## Material-change assessment

`material_change`: false. The SymPy reference engine was not modified, so no
derived result changed. The `.wl` only adds an independent confirmation of
already-established conclusions; downstream units depending on stage-190 outputs
are unaffected.

## Cross-engine agreement

Yes. Every conclusion derived by the `.wl` matches the SymPy result: the five
slippage closed forms, the central decomposition `Xi_direct = A_tr*Sigma_tr + Sigma_nt`,
the explicit `A_tr`/`C_tr`, the triangular compiler determinant and its three
inverse reconstructions plus three rigidity theorems, the support-blindness of
`T^2`/`R_target`/`N_*` (with the spoiled-packet control firing), and the quadratic
grouped invariant `I = (7/10)eps^2` with no linear feed-down. The `.wl` derives
`A_tr` and `Sigma_nt` by projection and reaches the same closed forms SymPy
postulates, and inverts the compiler by `LinearSolve` to the same recovered
formulas SymPy hand-writes.

## Side observations (non-blocking)

- **Banner provenance label.** The SymPy script still prints "STAGE 173" banners
  (a pre-renumbering artifact already recorded as cosmetic in the original audit's
  observation (a)). The new `.wl` correctly prints "STAGE 190", so the two outputs
  carry different headers. Cosmetic only; both sets of math are unambiguously
  stage-190's. Not a finding.
- **Off-diagonal derivative checks** (`dXi/dSigma_eta` etc., `.wl` lines 203–205)
  are zero by construction in both engines — the same guaranteed-by-construction
  status the original audit assigned to SymPy's A15. The `.wl` does at least build
  `xiExpr`/`dressExpr` from the actual `compiler.normalVector` product, so the
  matrix layout is genuinely exercised; the load-bearing triangularity content is
  carried by the substantive Coefficient-projection (Sec III) and LinearSolve
  inversion (Sec IV) checks. No new tautology was introduced by Codex.

## Verdict justification

The sole finding (the dual-engine gap) is fully resolved. Codex produced a
genuinely independent Mathematica route — logarithmic Euler operators for the
ledger, `Coefficient`-based projection for the direct-defect split, native
`LinearSolve`/`Det` matrix algebra for the compiler inversion, and
projector/kernel matrix products plus `Series`/`Coefficient` for the no-go filter
— not a line-by-line port of the SymPy closed forms. It even strengthens one
guard into an exact numerical witness (`-46/133`, independently confirmed). Both
engines exit 0 with all checks passing, outputs are fresh, the SymPy reference is
untouched, and every cross-engine conclusion agrees. Verdict: verified.
