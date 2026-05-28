---
unit_id: 157
batch: IV.6
verifier_model: claude-opus-4-7
verify_date: 2026-05-28T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 157

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`.wl` lines 102-105 replace the old `deltaCFromTangent = -16 sigmaStar (dR /. dg -> gp dr)` projector with an independent Solve on a fresh symbol set (`dCsym, dKsym`): `solDeltaC = Solve[{dCsym - 9 sigmaStar dKsym == 0, 5 dCsym - 72 sigmaStar dKsym == 0}, {dCsym, dKsym}]; deltaCIndep = FullSimplify[dCsym /. First[solDeltaC]]; expectZero["delta C from canonical-even Solve", deltaCIndep]`. The literal `-16 sigmaStar` multiplier no longer appears as an `expectZero` target.

**Assessment:**
Directive option A applied. The new check routes the `delta C = 0` claim through a fresh 2×2 homogeneous Solve (not through the closed-form `-16 sigmaStar` projector shared with Python), so it genuinely exercises whether the system's determinant is nonzero rather than re-evaluating an already-proven zero. Non-tautological. The dE2/dE4 coefficient literals (27, 243, 9, 72, 5) remain shared between engines, but the directive explicitly accepted option A over option B for that reason.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
`.py` lines 112-114 now read `sol_deltaC = sp.solve([sp.Eq(dE2,0), sp.Eq(dE4,0)], [deltaC, dkappa], dict=True)[0]; deltaC_from_pair = sp.simplify(sol_deltaC[deltaC]); expect_zero("delta C from canonical-even Solve", deltaC_from_pair)`. The trivial `-16 sigma_star * dR.subs(dg, gp*dr)` assertion is removed. The `.wl` mirrors with the F1 Solve-based check carrying the same banner string.

**Assessment:**
Directive option B applied in both engines, with mirrored wording. The new assertion is no longer a literal multiplication of a known-zero quantity; it now exercises the rank of the canonical-even pair. Output transcripts confirm the new banner `delta C from canonical-even Solve = 0` in both `.txt` files. Note: sympy invokes `sp.solve` on the same system twice (line 107 for `even_preservation`, line 112 for `sol_deltaC`) — mildly redundant but harmless and consistent with the directive.

### F3 — stale_output

**Classification:** resolved

**What changed:**
No source edit required; refreshed outputs supersede the old mtimes.

**Assessment:**
sympy output mtime 11:30 > script mtime 10:03; mathematica output mtime 11:31 > script mtime 10:03. Both fresh.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `R(g_*) - 1/4 = 0`
- `tangent motion keeps delta R = 0 = 0`
- `canonical-even preservation solutions = [{deltaC: 0, delta_kappa: 0}]`
- `delta C from canonical-even Solve = 0`
- `Stage 158 tangent expansion packet = 0`

**Mathematica:** exit=0. Notable lines:
- `PASS: R(g_*) - 1/4`
- `PASS: tangent motion keeps delta R = 0`
- `canonical-even preservation solutions = {{deltaC -> 0, dKappa -> 0}}`
- `delta C from canonical-even Solve = 0` followed by `PASS: delta C from canonical-even Solve`
- `Stage 157 Mathematica audit passed.`

**Output freshness:** Confirmed. Both `.txt` mtimes (11:30, 11:31 2026-05-28) are newer than their respective script mtimes (10:03 2026-05-28).

## Material-change assessment

`material_change`: false.

The edits replace a trivial assertion with a substantive equivalent that asserts the same outcome (delta C = 0 under canonical-even preservation). No constants, no derived numerical results, no symbolic identities downstream of this stage are altered. The carry-forward tuple (`r_F1`, `g_*`, `Sigma_0^can`, `T_hat_can`, `Pi_can`, `S_can`) is unchanged. Downstream units do not depend on the wording of the assertion.

## Side observations (non-blocking)

- The Stage 2 banner in the sympy output still reads "Carry-forward numerical basepoint from Stages 138-139" while the .wl banner reads "Stages 155-156". The directive didn't ask for either; the auditor noted the 138-139 phrasing as pre-existing. Not a verification blocker.
- The new sympy block re-Solves the same 2×2 system already Solved one line earlier for `even_preservation`. Functionally fine; cosmetically could reuse `even_preservation[0][deltaC]`. Not blocking.

## Verdict justification

All three findings are resolved. Both engines exit 0 with the new substantive `delta C from canonical-even Solve` check appearing in both transcripts. The literal `-16 sigmaStar` projector is gone from the `.wl`, satisfying F1's verification criterion. The replacement assertion is non-tautological (a homogeneous 2×2 Solve outcome rather than a multiplication-by-zero), satisfying F2. Output mtimes are post-fix, satisfying F3. No regressions in either log; numerical residuals unchanged.
