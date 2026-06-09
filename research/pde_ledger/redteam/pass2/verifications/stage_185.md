---
unit_id: 185
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T01:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 185

## Per-finding outcomes

### F1 — mathematica_transliteration (USER-AUTHORIZED full re-author; CHECKPOINT)

**Classification:** resolved

**What changed:**
The `.wl` monomial-exponent compilation was re-authored to derive the primitive
exponent vectors instead of hand-coding them.
- New machinery (`wl:117-148`): `primitiveVars`/`primitiveRatios`/`primitiveDrifts`/`logVars`,
  plus `monomialExponentVector[m]` (`wl:125-132`) — substitutes each primitive var → `Exp[logVar]`,
  takes `Log`, `PowerExpand`s, and reads off `Coefficient[logForm, logVar_i]`; `ratioFromExponentVector`
  (`wl:134-137`) and `driftFromExponentVector` (`wl:139-142`). Only the FIVE base kernel monomials are
  declared in genuine primitive form (`wl:144-148`: `chiMonomial=gammaVar*cetaUVar/kuVar`,
  `deltaUMonomial=tuVar/kuVar`, `epsWMonomial=gammaVar^2*lamWVar^2/(kuVar*kweffVar)`, `zMonomial`,
  `epsEtaMonomial`) — these match the reference kernel defs (`wl:83-101`), NOT the M_* rows.
- Composites built FROM those bases: `ctrMonomial = chiMonomial^(1+deltaUs)*deltaUMonomial^(1+chi0s)`
  (`wl:183`), `cntMonomial = zMonomial*epsWMonomial^eStar*deltaUMonomial^(-fStar)` (`wl:194`).
  Exponent vectors then EXTRACTED automatically (`wl:184,195,204`).
- Det-minor and Solve re-routed through the compiled drifts via
  `ledgerSigmaTr/Nt/Eta = sigmaTrCompiled/...` (`wl:243-251`, `252-255`).

**Assessment:**
Correct and addresses the finding at the checkpoint bar. Contrast with the `.py` (py:157-179),
which hand-writes the literal primitive powers `KU_ratio^(-(2+chi0s+deltaUs))`,
`lamW_ratio^(2+2*E_star)`, `KU_ratio^(F_star-E_star)`, etc. In the new `.wl` those exponents are
NEVER typed; they emerge mechanically from `monomialExponentVector`. A wrong base-monomial
structure or wrong composite power propagates to a wrong extracted exponent vector and breaks
`ctrRatio - ctrRatioPrimitive` (wl:188) / `sigmaTrCompiled - sigmaTr` (wl:189) /
`cntRatio - cntRatioPrimitive` (wl:199). This is a genuine independent derivation, not a relabeled
copy. Keeping `firstRatioDrift` (wl:48) is acceptable: it is only the slope operator on the
(now independently constructed) ratios, while the compiled drift uses the separate
`driftFromExponentVector` dot-product — the load-bearing piece (exponent provenance) is derived.
No collateral edit beyond the authorized re-author; `.py` untouched.

## Exec log assessment

**SymPy:** exit=0. No re-run needed (`.py` unchanged); output unchanged in substance, all residuals `= 0`.
- `det M_*^(tau,keta,mu) - (1+chi0s) = 0` (L76)
- `tracking/dressing/nontracking substitution = 0` (L81-83)

**Mathematica:** exit=0. Every check PASS:
- `det M_*^(tau,keta,mu) - (1+chi0s) = 0` → `PASS` (L108-109) — re-derived through `ledgerSigma*` = compiled drifts.
- `C_tr,* / C_nt,* ratio from primitive coordinates = 0` → PASS (L61-62, L71-72) — exercises the derived exponent vectors.
- `d ln C_tr,*/C_nt,*/epsilon_eta (primitive compiler) - Sigma_* = 0` → PASS (L63-64, L73-74, L83-84).
- `tracking/dressing/nontracking substitution` → PASS (L114-119).

**Output freshness:** confirmed. `.wl` mtime 2026-06-08 22:44 < `.txt` 2026-06-09 00:25; sympy `.txt`
2026-06-09 00:25. Both regenerated post-fix; exec-log dates (00:24/00:25) match.

## Material-change assessment

`material_change`: false. The re-author is method-only; every one of the 17 deliverable values is
preserved. Σ_tr, Σ_nt, Σ_eta, E_*, F_*, the three monomials and their primitive exponents, the
triangular law (Θ1, Ξ1, R1+Ξ1), `det = 1+χ0s`, and the solved triple (τ1, κη, μ1) all print
identically (algebraically) to the prior run — verified against the auditor's reconciliation table
(report L116-138, 17/17 MATCH). No downstream unit's consumed value changes. (Stage 186 consumes the
M_* setup; det = 1+χ0* is unchanged.)

## Side observations (non-blocking)

- `monomialExponentVector` relies on `PowerExpand` over real-positive primitives; sound here because
  `$Assumptions` declares all `*Ref > 0` and the `logVar`s are introduced as free coefficients — the
  `Coefficient` extraction is exact for these pure power monomials. Not a finding.
- The base monomial defs (wl:144-148) duplicate the kernel structure already in `chi00`/`deltaU0`/
  `epsW0`/`zratio0`/`epseta0` (wl:83-101); the `ratio from primitive coordinates` checks (wl:188,199,209)
  cross-tie the two, so a divergence would be caught. Healthy redundancy, not duplication risk.

## Verdict justification

`verified`. The single authorized finding (full re-author of the second engine) is resolved at the
checkpoint bar: the monomial-exponent compilation is no longer a hand-coded mirror differentiated by
the same operator — `monomialExponentVector` genuinely DERIVES the primitive exponent vectors from the
base-kernel monomial definitions via log-coordinate `Coefficient` extraction, so a wrong exponent would
now fail (the `.py` instead hand-types those exponents). Keeping `firstRatioDrift` is acceptable since
the exponent provenance — the load-bearing part — is now independent. All seven deliverable groups
still verify; `det ∂(Σ_tr,Σ_nt,Σ_eta)/∂(τ1,κη,μ1) = 1+χ0*` is re-derived through the compiled drifts
(wl:246-251, PASS); the zero-defect solve and three back-substitutions PASS; `math -script` exits 0;
`.py` unchanged; outputs fresh; all 17 deliverable values preserved.
