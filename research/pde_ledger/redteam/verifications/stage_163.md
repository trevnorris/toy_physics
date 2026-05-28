---
unit_id: 163
batch: IV.6
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-28T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 163

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Orchestrator applied edits directly (no `## Applied: F1` block). In
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.wl`:
- Lines 38-40: implicit-function route added — `gPrimeImplicit = -D[fComp,r]/D[fComp,g] /. g -> gMinus` with `expectZero["gPrime matches implicit-function route", gPrime - gPrimeImplicit]`.
- Lines 61-81: chain-rule Series route — `rExpr = lam/Sqrt[Ks*Kq]`, `gExpr = gq*Sqrt[Ks]/(gs*Sqrt[Kq])`, `pertRule` log-perturbations, `Series[..., {eps,0,1}]` then `Coefficient[..., eps]`, with three new `expectZero` lines for `delta r`, `delta g`, and `delta_perp` via series.
- Line 83: `Clear[Ks, Kq, lam, gs, gq, eps, pertRule, rExpr, ...]` restores clean scope for downstream sections.
- Original assertions (delta F, delta R, microscopic identity, delta C, delta Pi split) preserved verbatim. SymPy file untouched.

**Assessment:**
Edits match the directive exactly. The two new routes are structurally distinct: (a) implicit-function derivative uses `D[fComp,g]` and `D[fComp,r]` partials never invoked in the SymPy script; (b) `Series[...,{eps,0,1}]`/`Coefficient` has no SymPy counterpart in the audit. Neither is tautological — both could fail under sign or coefficient error in the original `gPrime` or hand-written `deltaG`/`deltaR`. The `Clear` block correctly restores scope before `sigmaStar`/`dPi` sections. No collateral edits beyond the directive's scope.

## Exec log assessment

**SymPy:** exit=0. Unchanged; all original PASS lines (delta F, delta R, microscopic, delta C, delta Pi) reproduce.

**Mathematica:** exit=0. New PASS lines all present:
- line 6: `PASS: gPrime matches implicit-function route`
- line 16: `PASS: delta r series matches hand form`
- line 18: `PASS: delta g series matches hand form`
- line 20: `PASS: microscopic delta_perp via series route`

Pre-existing PASS lines (delta F, delta R, microscopic identity, delta C linkage, delta Pi split) all still appear at lines 8, 10, 13, 26, 28.

**Output freshness:** Mathematica output mtime 1779989544 > script mtime 1779984351. SymPy output mtime 1779989435 > script mtime 1779945093. Both outputs regenerated post-fix.

## Material-change assessment

`material_change`: false. The edit adds independent-derivation assertions; no derived results change. Downstream units not affected.

## Side observations (non-blocking)

- Script banner still reads `STAGE 163 — OFF-FAMILY NORMAL COORDINATE` correctly (the auditor flagged the SymPy banner as `STAGE 146` in its verdict justification, but that is a SymPy-only label issue and not part of F1's scope).

## Verdict justification

F1 is fully resolved. The Mathematica script now carries two genuinely independent derivation routes (implicit-function ratio for `g_-'`; chain-rule Log/Series expansion for the microscopic `delta_perp`) alongside the original parallel assertions. Both new routes use Mathematica machinery with no SymPy mirror in the audit, so cross-engine agreement is no longer structurally guaranteed. Both engines exit 0 and outputs are fresh. No regressions, no collateral edits.
