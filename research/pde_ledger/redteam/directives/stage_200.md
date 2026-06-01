---
unit_id: 200
batch: V.3
created_at: 2026-06-01T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-01T18:30:20Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 200

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT touch the SymPy script (it is correct and substantive), and do NOT touch paper.tex, notes/, or any prose document. Do NOT edit the cosmetic stage-number banner strings ("STAGE 183") — they are out of scope.

After editing, RUN the Mathematica script (`math -script /var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl`) and iterate until it exits 0 with all in-file checks passing.

F1 and F2 overlap: applying F2's Section III rewrite discharges the Section III portion of F1. Apply F2 first, then the remaining Section I portion of F1.

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl:236-238`

**Issue:** The Section III actual/target ratios are written in already-collapsed form built directly from the `m·orbit` definitions (`TActual=mT TOrbit`, `KetaActual=mK KetaOrbit`, `muActual=mMu muOrbit`). As a result the three Section III `expectZero` checks reduce to algebraic `Log[a^b]=b Log[a]` identities and never exercise the Ctr/Cnt/epsEta monomials' exponent structure. The SymPy script (sympy:304-311) instead rebuilds the full monomial with the actual values substituted and divides by an independent target symbol (`Ctr_target`, etc.), which genuinely tests the chart-conversion law. Make the Mathematica route match the SymPy substantive route.

**Required change:**
Replace the three lines

```
CtrActualRatio = normalizeExpr[(TActual/TOrbit)^(1 + chi0s)];
CntActualRatio = normalizeExpr[(muActual/muOrbit) (KetaOrbit/KetaActual) (TActual/TOrbit)^(-fStar)];
epsEtaActualRatio = normalizeExpr[KetaOrbit/KetaActual];
```

with

```
CtrActualRatio = normalizeExpr[ctrMonomial[gf, cf, KUf, TActual, chi0s, deltaUs, L]/CtrTarget];
CntActualRatio = normalizeExpr[cntMonomial[lamf, gf, KUf, KetaActual, KWf, muActual, TActual, eStar, fStar, L, sigma]/CntTarget];
epsEtaActualRatio = normalizeExpr[epsEtaMonomial[cf, KUf, KetaActual]/epsEtaTarget];
```

Notes for safe application:
- Keep `TOrbit`, `KetaOrbit`, `muOrbit` (math:222-230) and `TActual`, `KetaActual`, `muActual` (math:232-234) exactly as they are. `TOrbit` is the orbit solve that makes `ctrMonomial[gf, cf, KUf, TOrbit, ...]` equal `CtrTarget`, so dividing the actual-value monomial by `CtrTarget` gives the genuine ratio.
- The symbols `CtrTarget`, `CntTarget`, `epsEtaTarget`, `mT`, `mK`, `mMu`, `gf`, `cf`, `KUf`, `KWf`, `lamf`, `sigma`, `eStar`, `fStar`, `chi0s`, `deltaUs`, `L` are all already declared positive/real in `$Assumptions` (math:84-114). Do NOT add new symbol declarations.
- Leave the three `expectZero` calls (math:240-245) unchanged: their targets `Log[CtrActualRatio] - (1 + chi0s) Log[mT]`, `Log[CntActualRatio] - (Log[mMu] - Log[mK] - fStar Log[mT])`, and `Log[epsEtaActualRatio] + Log[mK]` will now test the real monomial conversions.
- Leave `qMismatch`/`qMismatchExpected`/`Dmis` (math:247-254) unchanged.

After the rewrite, the three Section III ratios must reduce (under `normalizeExpr` / `PowerExpand` with the positivity assumptions) to `mT^(1+chi0s)`, `mMu/(mK mT^fStar)`, and `1/mK` respectively, so each `expectZero` still evaluates to 0 / PASS — but now because the monomial exponents are correct, not by construction.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 200` and confirm: (a) math:236 now reads `CtrActualRatio = normalizeExpr[ctrMonomial[gf, cf, KUf, TActual, chi0s, deltaUs, L]/CtrTarget]`; (b) the three Section III `expectZero` lines print `= 0` / `PASS`; (c) the printed `q(mismatch)` matrix is unchanged from the committed output (`{(1 + chi0s) Log[mT], -Log[mK] + Log[mMu] - fStar Log[mT], -Log[mK]}`); (d) the script exits 0.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl`
- summary: Section III actual/target ratios now divide the full helper monomials by the independent targets before checking the mismatch chart.
- deviation: Also parenthesized the existing `cntMonomial` helper body so Mathematica evaluates all factors used by the required helper-based quotient.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl:118-130` (Section I ratio construction)

**Issue:** The `.wl` mirrors the `.py` derivation route rather than deriving independently. The most concentrated symptom (Section III) is fixed by F2. The remaining transliteration symptom this directive targets is Section I, where Mathematica hardcodes the pre-reduced monomial ratios instead of building them from the helper functions the way SymPy does (sympy:118-119 divides two `ctr_monomial`/`cnt_monomial`/`eps_eta_monomial` calls after applying the ratio substitution). Making Section I derive the ratios from the helper functions gives the second engine an independent route into the load-bearing `Mderived` Jacobian.

**Required change:**
Replace the Section I ratio block

```
CtrRatio = FullSimplify[
  PowerExpand[(rg rc/rU)^(1 + deltaUs) (rT/rU)^(1 + chi0s)],
  Assumptions -> $Assumptions
];
CntRatio = FullSimplify[
  PowerExpand[
    (rla^2 rmu/(rK rW^2))
      (rg^2 rla^2/(rU rW))^eStar
      (rT/rU)^(-fStar)
  ],
  Assumptions -> $Assumptions
];
epsRatio = FullSimplify[PowerExpand[rc^2/(rK rU)], Assumptions -> $Assumptions];
```

with a construction that divides the monomials evaluated at scaled-vs-unscaled arguments, mirroring SymPy's `(Ctr2/Ctr1).subs(ratio_subs)` route. Introduce a `ratioSubs` rule list mapping the "2" symbols to `r * "1"` symbols, build `Ctr1/Ctr2/Cnt1/Cnt2/eps1/eps2` from the helper functions, and form the ratios:

```
ratioSubs = {
  lam2 -> rla lam1, c2 -> rc c1, gam2 -> rg gam1, KU2 -> rU KU1,
  Keta2 -> rK Keta1, KW2 -> rW KW1, mu2 -> rmu mu1, T2 -> rT T1
};
Ctr1 = ctrMonomial[gam1, c1, KU1, T1, chi0s, deltaUs, L];
Ctr2 = ctrMonomial[gam2, c2, KU2, T2, chi0s, deltaUs, L];
Cnt1 = cntMonomial[lam1, gam1, KU1, Keta1, KW1, mu1, T1, eStar, fStar, L, sigma];
Cnt2 = cntMonomial[lam2, gam2, KU2, Keta2, KW2, mu2, T2, eStar, fStar, L, sigma];
eps1 = epsEtaMonomial[c1, KU1, Keta1];
eps2 = epsEtaMonomial[c2, KU2, Keta2];
CtrRatio = FullSimplify[PowerExpand[(Ctr2/Ctr1) /. ratioSubs], Assumptions -> $Assumptions];
CntRatio = FullSimplify[PowerExpand[(Cnt2/Cnt1) /. ratioSubs], Assumptions -> $Assumptions];
epsRatio = FullSimplify[PowerExpand[(eps2/eps1) /. ratioSubs], Assumptions -> $Assumptions];
```

Notes for safe application:
- The "1" and "2" symbols (`lam1,c1,gam1,KU1,Keta1,KW1,mu1,T1,lam2,...,T2`) and the ratio symbols (`rla,rc,rg,rU,rK,rW,rmu,rT`) are all already declared positive in `Clear[...]` (math:69-82) and `$Assumptions` (math:84-114). Do NOT add new declarations.
- Everything downstream of this block (`ratioToLogs`, `qPair`, `Dvec`, `Mderived`, `Mexpected`, the two Section I `expectZero` calls at math:167-168) stays unchanged; the rewritten `CtrRatio`/`CntRatio`/`epsRatio` must simplify to the same reduced forms (`(rc rg rT (...)^deltaUs (rT/rU)^chi0s)/rU^2`, etc.) so the Jacobian and comparisons are unaffected.
- Do NOT alter Sections II, IV, V. Section III is handled by F2.

**Claim manifest** (the rewritten Section I must independently verify):
- M1: the log-ratio Jacobian of the three primitive monomial ratios equals the carried Stage 192 matrix `M_*` — i.e. row 0 `[0, 1+deltaUs, 1+deltaUs, -(2+chi0s+deltaUs), 0, 0, 0, 1+chi0s]`, row 1 `[2(1+eStar), 0, 2 eStar, fStar-eStar, -1, -(2+eStar), 1, -fStar]`, row 2 `[0, 2, 0, -1, -1, 0, 0, 0]` (`Mderived - Mexpected == 0`).
- M2: `qPair - Mderived . Dvec == 0`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 200` and confirm: (a) math Section I builds `CtrRatio` etc. from `ctrMonomial`/`cntMonomial`/`epsEtaMonomial` quotients with a `ratioSubs` rule (no hand-collapsed literal); (b) `derived M_* - carried Stage 192 matrix` prints the all-zero matrix and `q^(2<-1) - M_* Delta x^(2<-1)` prints `{0,0,0}`; (c) the script exits 0 with all PASS lines intact.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl`
- summary: Section I now builds the primitive monomial ratios from scaled and unscaled helper-function calls via `ratioSubs`.
- deviation: Same helper-body parenthesization as F2 was required so `cntMonomial` carries the `eStar` and `fStar` factors into the quotient.
