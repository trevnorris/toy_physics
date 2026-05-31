---
unit_id: 184
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-30T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 184

## Per-finding outcomes

### F1 — tautological_check (HIGH)

**Classification:** resolved

**What changed:**
Four drift quantities in the Mathematica audit were re-routed from hardcoded restatements of their comparison targets to genuine series extractions from the branch composites (`mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl`):
- Line 53: `dlnTtr = FullSimplify[SeriesCoefficient[Log[tTr/tTr0], {small, 0, 1}], ...]` (was `FullSimplify[-cStar*theta1, ...]`)
- Line 60: `dlnNtr = FullSimplify[SeriesCoefficient[Log[nTr/nTr0], {small, 0, 1}], ...]` (was `FullSimplify[xi1 + bStar*theta1, ...]`)
- Line 65: `dlnEpsEta = FullSimplify[SeriesCoefficient[Log[epsEtaVar/epsEta], {small, 0, 1}], ...]` (was `FullSimplify[sigmaEta, ...]`)
- Line 71: `dlnEcomp = FullSimplify[SeriesCoefficient[Log[eComp/eComp0], {small, 0, 1}], ...]` (was `FullSimplify[-epsEta*sigmaEta/(1 - epsEta), ...]`)

The diff (`stage_184_diff.patch`) touches exactly these four lines and nothing else. No collateral edits. The `expectZero` lines, composite definitions, zero-map block, and banner were all left unchanged (the optional STAGE 167 → 184 banner relabel was correctly not done, and not gated on).

**Assessment:**
The fix is correct and the new checks are genuinely non-tautological. I traced each derivation against the composite definitions that carry the `small` parameter:
- `rTr = rTr0*Exp[small*theta1]` (line 44), `t2 = t20*Exp[small*xi1]` (line 45), `epsEtaVar = epsEta*(1 + small*sigmaEta)` (line 46), `rTarget = lam0*(1 - epsEtaVar)/t2` (line 47). All four carry `small`.
- `dlnTtr`: `Log[rTr^(-cStar)/rTr0^(-cStar)] = -cStar*small*theta1` → first-order coeff `-cStar*theta1`. Compared at line 55 against `sigmaTr` (independently defined at line 40 as `-cStar*theta1`). If `cStar` or the `-cStar` exponent in `tTr` were wrong, this would now FAIL — the SeriesCoefficient pulls `cStar` from the composite exponent while `sigmaTr` carries it as the closed form, so they are two distinct routes to the same constant.
- `dlnNtr`: `Log[(t2*rTr^bStar)/(t20*rTr0^bStar)] = small*xi1 + bStar*small*theta1` → `xi1 + bStar*theta1`. Compared at line 62 against independently-defined `sigmaNT` (line 41). Genuinely load-bearing in `bStar`.
- `dlnEpsEta`: `Log[(epsEta*(1+small*sigmaEta))/epsEta] = Log[1+small*sigmaEta]` → first-order coeff `sigmaEta`. Compared against `sigmaEta`. Now derived from `epsEtaVar`, not restated.
- `dlnEcomp`: `eComp = rTarget*t2/lam0 = 1 - epsEtaVar = 1 - epsEta*(1+small*sigmaEta)`, `eComp0 = 1 - epsEta`; `Log[eComp/eComp0]` first-order coeff `-epsEta*sigmaEta/(1-epsEta)`. Matches the directive's expected value and the complement assertion at line 73 closes to 0.

The four drift-law assertions (55, 62, 67, 73), the three mirror checks (76-78), and the three zero-map checks (83-85) now exercise the composites rather than `expr - expr`. The printed values are unchanged from the pre-fix transcript and match the directive's expected post-fix values verbatim. SymPy was correctly left untouched (no `.py` in the diff; SymPy script mtime unchanged at 1778522333).

### F2 — insufficient_verification (MEDIUM)

**Classification:** resolved

**What changed:**
No separate edit (per directive, F2 shares F1's root cause and fix). The four composite definitions that were previously dead code are now consumed:
- `tTr`/`tTr0` (51-52) consumed by `dlnTtr` (line 53)
- `nTr`/`nTr0` (58-59) consumed by `dlnNtr` (line 60)
- `epsEtaVar` (46) consumed by `dlnEpsEta` (line 65)
- `eComp`/`eComp0` (69-70) consumed by `dlnEcomp` (line 71)

**Assessment:**
Confirmed by direct read of the post-fix script: every composite now appears inside a `Log[.../...]` SeriesCoefficient call, so none remains write-only. The Mathematica engine now independently exercises D2–D5 from the branch composites. Resolving F1 resolved F2 as the directive intended.

## Exec log assessment

**SymPy:** exit=0. The SymPy script was not part of this directive and was not edited. The log (`stage_184_sympy.log`) re-confirms the genuine series-derived checks already present pre-audit: `delta ln T_* - Sigma_tr = 0`, `delta ln N_* - Sigma_nt = 0`, `delta ln eps_eta - Sigma_eta = 0`, `selected-branch complement identity = 0`. exit_code: 0.

**Mathematica:** exit=0. Notable lines from `stage_184_mathematica.log`:
- `delta ln T_* = -(((1 + chi0)*(1 + deltaU)*(1 + chi0 + deltaU)*theta1)/(chi0*deltaU))` ... `PASS: delta ln T_* - Sigma_tr`
- `delta ln N_* = (2*(1 + chi0 + deltaU)*theta1)/deltaU + xi1` ... `PASS: delta ln N_* - Sigma_nt`
- `delta ln eps_eta = sigmaEta` ... `PASS: delta ln eps_eta - Sigma_eta`
- `delta ln[(R_target T^2)/Lambda0] = (epsEta*sigmaEta)/(-1 + epsEta)` ... `PASS: selected-branch complement identity`
- All three zero-map checks PASS. `Stage 184 Mathematica audit passed.` exit_code: 0.

The printed drift values are produced by the SeriesCoefficient operation yet are unchanged from pre-fix, exactly as required — a passing run that is now also non-vacuous.

**Output freshness:** Confirmed. Mathematica output mtime 1780126869 > script mtime 1780126166 (output regenerated post-fix). SymPy output mtime 1780126862 > script mtime 1778522333 (script untouched, output re-run). The saved `.txt` (`mathematica/output/...stage184...txt`) matches the exec log line-for-line with all 11 PASS lines.

## Material-change assessment

`material_change`: false.

The fix only strengthens the second-engine (Mathematica) corroboration of D2–D5. No forward-carried constant, composite definition, or derived value changed — the printed drift values are identical to the pre-fix transcript. The SymPy engine already established the math pre-audit. No downstream unit depends on a changed result, so no `upstream_stale` propagation is warranted on substantive grounds.

## Side observations (non-blocking)

- The banner on both the `.wl` (line 26) and the SymPy `.py` still reads "STAGE 167 — EXACT BRANCH-INVARIANT COORDINATES" though this is stage 184. The directive marked this an optional, non-gating cosmetic relabel; Codex correctly left it (editing the SymPy banner was prohibited under this directive, and the Mathematica relabel was optional). Cosmetic only — does not affect any assertion. Not a basis to fail verification.

## Verdict justification

Both findings are `resolved`. The diff matches the directive's four prescribed substitutions exactly with no collateral edits; the SymPy script was correctly left untouched. I verified by tracing each derivation that the four drift quantities are now genuinely produced by `SeriesCoefficient[Log[composite/reference], {small, 0, 1}]` of composites that carry the `small` parameter, and are compared against independently-defined closed forms (`sigmaTr`, `sigmaNT`, `sigmaEta`) — so the drift-law assertions are no longer `expr - expr` tautologies and would fail if `cStar`/`bStar` or a composite exponent were wrong (F1). The previously dead composites are now load-bearing (F2). Both exec logs exit 0 with all PASS lines, outputs are fresh, and the printed values are unchanged. Verdict: verified.
