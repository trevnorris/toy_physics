---
unit_id: 171
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-08T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 171

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Codex made an additive edit to the `.wl` only (diff touches a single file, `mathematica/...stage171...audit.wl`):
1. `wl:59-63` — added `Clear[eps2]; bCombSeries = Coefficient[Normal[Series[(b2 + b0/9) /. {c -> c + eps2*dc, w -> w + eps2*dw}, {eps2,0,1}]], eps2];` then `expectZero["BdG obstruction bundle (series route)", bCombSeries - bCombFormula];` immediately after the existing D[]-route `BdG obstruction bundle` check (wl:58).
2. `wl:153-154` — defined `kScalar = (k1 - b01 - z01)/9 + (-m1 - b21 - z21);` and `gScalar = n01 - p0*(k1 - b01 - z01);` outside the Do-loop.
3. `wl:170-176` — inside the loop, added `weak-axisymmetric K scalar route` (`eps*lamVal*kScalar - kRebuilt`) and `weak-axisymmetric G scalar route` (`eps*lamVal*gScalar - gRebuilt`) checks.

No existing checks were removed; no `.py` change (not in git status); no prose touched.

**Assessment:**
The edit matches the directive exactly with no collateral changes.

- **bCombSeries** is genuinely distinct-mechanism: it recovers the first variation of `(b2 + b0/9)` by Taylor expansion in a fresh `eps2` perturbation (`Series` + `Coefficient`), as opposed to the existing `D[]`-summed total differential (`bCombExact`). Both `c` and `w` perturbations contribute (neither `δc` nor `δw` coefficient is identically zero), and it is compared against the hand-typed `bCombFormula`, so a wrong typed coefficient leaves a nonzero residual — non-tautological. `Clear[eps2]` guards against a leaked prior binding.
- **kScalar/gScalar** rebuild the `(𝔎₁,𝔊₁)` amplitudes from the Section-1 split grouping and compare against the loop's `kRebuilt`/`gRebuilt`, which use a different grouping of the same `(k1,m1,b01,b21,z01,z21,n01,p0)` symbols. The residual vanishes only if the split algebra is correct — a distinct route from the literally-typed `kMicro`/`gMicro`, hence additive independent value.

All 12 new checks (BdG series route + K/G scalar route × {1, 1/2, −1}) print `= 0` / `PASS` in the exec log. This closes the F1 transliteration finding: the `.wl` now carries independent-mechanism cross-routes on the BdG and weak-axisymmetric sections, joining the pre-existing `Series` routes on Z/N.

## Exec log assessment

**SymPy:** exit=n/a. No sympy log was captured for this stage (the F1 fix was `.wl`-only and the `.py` was untouched, so no sympy re-run was required).

**Mathematica:** exit=0. Notable lines:
- `BdG obstruction bundle (series route) = 0` / `PASS: BdG obstruction bundle (series route)`
- `weak-axisymmetric K scalar route lambda=1 = 0` / `PASS` (and `=1/2`, `=-1`)
- `weak-axisymmetric G scalar route lambda=1 = 0` / `PASS` (and `=1/2`, `=-1`)
- `Stage 171 Mathematica audit passed.` / `# exit_code: 0`

All pre-existing checks (split, BdG D-route, Maxwell primitives, δZ/δN, Z/N bundles + their Series routes, weak-axisym obstruction) still print `= 0` / `PASS` — no regressions.

**Output freshness:** confirmed. The saved `mathematica/output/...stage171...audit.txt` mtime is 2026-06-08 17:13:08, newer than the `.wl` mtime 2026-06-08 16:48:00, and it appears `M` (modified) in git status — regenerated post-fix.

## Material-change assessment

`material_change`: false.

The edit is purely additive verification scaffolding inside the `.wl` audit script — new `expectZero` cross-route checks against the *same* paper-typed targets. No derived constant, formula, or carry-forward result changed (the three carry-forward summary lines are identical). No downstream unit depends on the audit script's internal check count, so nothing > 171 is materially affected.

## Side observations (non-blocking)

- The added Series-route comment lines and the `eps2` reuse for the BdG route mirror the established pattern already used for the Z/N bundles (wl:129-135, 144-148), so the file is internally consistent.

## Verdict justification

The single low-severity `mathematica_transliteration` finding (F1) is resolved: Codex added exactly the two requested distinct-mechanism routes (a `Series`/`Coefficient` BdG-bundle check in a fresh `eps2` perturbation vs. the existing `D[]`-summation, and `kScalar`/`gScalar` split-derived weak-axisymmetric checks vs. the literal `kMicro`/`gMicro`), all additive with no existing checks removed, the `.py` untouched, all 12 new checks plus all prior checks passing, and the Mathematica script exits 0 with the output `.txt` regenerated. The new routes are genuinely distinct-mechanism and non-tautological. No regressions and no material change. Verdict: verified.
