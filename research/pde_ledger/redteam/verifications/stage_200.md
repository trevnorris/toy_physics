---
unit_id: 200
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-01T19:05:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 200

This is a CHECKPOINT (trust anchor MTDC-T9.6). Reviewed at the higher bar: both
engines must establish the load-bearing results by genuinely substantive checks,
and the second engine must be an independent route, not a transliteration.

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/...stage200..._mathematica_audit.wl:119-131` (diff hunk 2). The Section I
ratio block previously hardcoded the pre-collapsed literals
`CtrRatio = (rg rc/rU)^(1+deltaUs)(rT/rU)^(1+chi0s)`, `CntRatio = (rla^2 rmu/(rK rW^2))(...)^eStar(rT/rU)^(-fStar)`,
`epsRatio = rc^2/(rK rU)`. It now introduces a `ratioSubs` rule
(`lam2->rla lam1, ..., T2->rT T1`), builds `Ctr1/Ctr2/Cnt1/Cnt2/eps1/eps2` from the
`ctrMonomial`/`cntMonomial`/`epsEtaMonomial` helper functions at unscaled-vs-scaled
arguments, and forms `CtrRatio = FullSimplify[PowerExpand[(Ctr2/Ctr1) /. ratioSubs]]`
(likewise Cnt/eps). This mirrors SymPy's `(Ctr2/Ctr1).subs(ratio_subs)` route
(sympy:118-120). No new symbols declared (all "1"/"2"/r-symbols pre-declared at
math:70-114).

**Assessment:**
Genuine independent route, not a re-typed literal. `Ctr2/Ctr1` is a ratio of two
distinct monomials in independent symbols; applying `ratioSubs` cancels the "1"
symbols and yields `(rc rg rT ((rc rg)/rU)^deltaUs (rT/rU)^chi0s)/rU^2` (math log
L15), which is exactly the SymPy reduced form (sympy output L10-16) — same answer,
now DERIVED from the helper monomials rather than asserted. The resulting `qPair`
feeds the `D[]` Jacobian and the load-bearing comparisons close: `Mderived - Mexpected`
is the all-zero 3x8 matrix (M1, math log L24-25) and `qPair - Mderived.Dvec = {0,0,0}`
(M2, math log L26-27). Crucially, if the `ctrMonomial`/`cntMonomial` exponents were
wrong, the derived ratio would change and `Mderived` would no longer equal the carried
Stage-192 literal `Mexpected` — so the second engine now independently constrains the
monomial exponent structure feeding the Jacobian. Deviation (cntMonomial body
parenthesized, diff hunk 1) is correct and necessary: without the wrapping parens the
`:=` RHS would bind only the first factor and drop the `eStar`/`fStar` factors; the log
confirms `Cnt_2/Cnt_1` carries `(rg rla)^(2 eStar)` and `(rU/rT)^fStar` (math log L17),
so the full monomial is present. Benign — does not change any derived result.

### F2 — tautological_check (HIGH)

**Classification:** resolved

**What changed:**
`mathematica/...stage200..._mathematica_audit.wl:237-239` (diff hunk 3). The three
Section III ratios previously built from the collapsed `m·orbit` pieces
(`(TActual/TOrbit)^(1+chi0s)`, `(muActual/muOrbit)(KetaOrbit/KetaActual)(TActual/TOrbit)^(-fStar)`,
`KetaOrbit/KetaActual`) are replaced by full-monomial-over-independent-target quotients:
`CtrActualRatio = normalizeExpr[ctrMonomial[gf,cf,KUf,TActual,chi0s,deltaUs,L]/CtrTarget]`,
`CntActualRatio = .../CntTarget`, `epsEtaActualRatio = epsEtaMonomial[cf,KUf,KetaActual]/epsEtaTarget`.
The three `expectZero` targets (math:241-246) are unchanged in form, matching SymPy
sympy:304-315 exactly. `TOrbit/KetaOrbit/muOrbit` (the orbit solve, math:223-231) and
`TActual/KetaActual/muActual` (math:233-235) are unchanged.

**Assessment:**
Genuinely de-tautologized — the checks now exercise the monomial exponents, NOT a
`Log[a^b]=b Log[a]` collapse. Concretely: `TOrbit` is the orbit solve that inverts
`ctrMonomial`, so `ctrMonomial[gf,cf,KUf,TOrbit,...] = CtrTarget`. Because `tU` enters
`ctrMonomial` ONLY through the factor `(Pi^2 tU/(L^2 KUf))^(1+chi0s)`, substituting
`TActual = mT·TOrbit` gives `ctrMonomial[...,mT·TOrbit,...] = mT^(1+chi0s)·CtrTarget`,
and dividing by the INDEPENDENT symbol `CtrTarget` yields `mT^(1+chi0s)`. The check
`Log[CtrActualRatio] - (1+chi0s)Log[mT] = 0` therefore passes only because the monomial's
thermal-factor exponent is `(1+chi0s)`. Falsification test: had the exponent been mistyped
(e.g. `chi0s` instead of `1+chi0s`), the ratio would be `mT^chi0s` and the residual
`-Log[mT] != 0` -> FAIL. The old hand-built `(TActual/TOrbit)^(1+chi0s)` could never fail
regardless of the monomial; the new form can. Same logic holds for epsEta (exercises that
Keta enters with exponent -1: `epsEtaMonomial[cf,KUf,mK·KetaOrbit]/epsEtaTarget = 1/mK`,
genuinely testing the sign) and for Cnt (exercises the `mMu/(mK mT^fStar)` structure
through the orbit-solved `muOrbit`/`KetaOrbit`). All three print `= 0` / PASS (math log
L43-48), `q(mismatch)` is unchanged at `{(1+chi0s)Log[mT], -Log[mK]+Log[mMu]-fStar Log[mT], -Log[mK]}`
(math log L50), and `M_* Dmis - q(mismatch) = {0,0,0}`. The Section III agreement is now
independent confirmation of the chart-conversion law, no longer vacuous.

## Exec log assessment

**SymPy:** exit=0. The substantive reference engine, untouched. Section I
`derived M_* - carried Stage 192 matrix` all-zero (output L55-60); Section III
`Ctr(actual)/Ctr_* - m_T^(1+chi0_*) = 0`, `Cnt(actual)/Cnt_* - m_mu/(m_K m_T^F_*) = 0`,
`epsEta(actual)/epsEta_* - 1/m_K = 0` (output L182-184). All sections close.

**Mathematica:** exit=0. 7 PASS lines, 0 FAIL. Section I `Mderived - Mexpected`
all-zero 3x8 (log L24-25), `qPair - Mderived.Dvec = {0,0,0}` (L26-27); Section III
three chart checks PASS (L43-48). Reduced forms match SymPy byte-for-byte where both
compute the same quantity.

**Cross-engine:** the load-bearing results agree across engines AND the agreement is
now non-vacuous in Section III (the prior defect). Section I reduced ratios identical;
`Mderived = Mexpected` in both; mismatch chart `q` identical; Sections II/IV/V residuals
all zero in both.

**Output freshness:** confirmed. `.wl` mtime 2026-06-01 12:29:43; `.txt` outputs
(both math and sympy) 2026-06-01 12:42:20 — both regenerated post-fix (newer than the
edited script).

## Material-change assessment

`material_change`: false. Every edit was to the Mathematica engine's derivation ROUTE,
not to any value. The reduced Section I ratios, `Mderived`, the Section III mismatch
chart `q`, and all printed residuals are identical to the committed pre-fix output (and
to SymPy). No derived constant changed; no downstream unit depends on a changed result.

## Side observations (non-blocking)

- The in-script banner still reads "STAGE 183" and the ledger header "STAGE 183 LEDGER";
  these are the pre-renumbering cosmetic artifacts the auditor and directive explicitly
  placed out of scope. Not a math issue; not blocking.
- Section IV's cocycle checks remain near-definitional (log-additivity of composed
  ratios), as the auditor noted ("partial" / near-definitional). This was not a directed
  finding and the check is not vacuous (it uses the real `Mderived`); no action needed.

## Verdict justification

`verified`. Both directed findings are resolved at the checkpoint bar. F2 (HIGH): the
Section III checks are now genuinely substantive — `CtrActualRatio` is the full
`ctrMonomial[...,TActual,...]` divided by the independent `CtrTarget`, so the residual
`Log[...] - (1+chi0s)Log[mT]` vanishes only because the monomial's exponent structure is
correct, and a wrong exponent would make it FAIL (verified by reasoning through the orbit
solve, not just trusting the PASS). F1: Section I now derives `CtrRatio`/`CntRatio`/`epsRatio`
from helper-monomial quotients via `ratioSubs`, an independent route into the `Mderived`
Jacobian; `Mderived - Mexpected` and `qPair - Mderived.Dvec` are all-zero. No remaining
transliteration or vacuous check in Sections I-V. Both engines exit 0, outputs are fresh,
and the cross-engine agreement is now independent (not echo). The SymPy reference engine
was untouched and remains substantive throughout.
