---
unit_id: 178
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

# Verification — unit 178

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage178_outgoing_port_coloading_mathematica_audit.wl:103–110` adds a new independent cross-derivation of `nu_r` inside Section 4, immediately after the existing `expectZero["nu_r - [kappa1 + sigma_r]", ...]` (line 101):

```
nuFromData = FullSimplify[
  Coefficient[Normal[Series[Log[pA^2/dA^2], {eps, 0, 1}]], eps*lam],
  Assumptions -> $Assumptions
];
expectZero["nu via log-data vs slippage", nuFromData - nuExpected];
Print["nuFromData = ", fmt[nuFromData]];
```

The diff matches the directive's "required change" verbatim. The existing §1–§4 checks (`pExpected`/`dExpected`, `nuDirect`, `nuExpected`) were left untouched. No collateral edits beyond this insert and the F2 banner line.

**Assessment:**
The edit correctly addresses the finding and the new check is genuinely independent and non-tautological. Traced the data flow:

- **Independent inputs.** `nuFromData` is built from `pA` (line 62) and `dA` (line 63), the component-level port quantities assembled directly from the five physical drifts (`oU, gW, rr, gU, oW`). It does not reference `pExpected`, `dExpected`, the weights `alpha/beta/chi/zeta`, nor the §4 slippage pieces `mR/iR/hR`.
- **Independent extraction route.** It uses `Coefficient[Normal[Series[Log[pA^2/dA^2], {eps,0,1}]], eps*lam]` — a logarithmic-slope extraction. This differs from §1 (which used a linear ratio-minus-one extraction `(Normal[Series[nAr]]/n0r - 1)/(eps*lam)` on the *abstract* `pAr/dAr` symbols) and from §3 (linear series on `pA/p`, `dA/delta`). The log route never retraces the `2*(pExpected - dExpected)` factoring the auditor flagged as the shared SymPy choreography.
- **Independent target.** The residual `nuFromData - nuExpected` compares the log-data route against `nuExpected`, the Stage-177 slippage assembly (`kappa1 + 2*mR + 2*iRExpr*iR/(1+iRExpr) + 2*hRExpr*hR/(1-hRExpr)`), not against `nuDirect = 2*(pExpected-dExpected)`. So a construction error in the `pExpected`/`dExpected` factoring would now be caught by this path.
- **Non-vacuous.** The exec log (mathematica line 52 / saved txt line 47) shows `nuFromData` resolves to a non-trivial rational form `(-2*ou2*(oU+oW)*ow2 + 4*r^2*rr)/(ou2*ow2 - r^2) + (2*(gw*(gW+oU)*ou2 + gu*r*(gU+rr)))/(gw*ou2 + gu*r)` and the residual against the slippage form prints `0` (PASS). This is a real cancellation between two differently-built expressions, not a `0 − 0` identity.

This is NOT an X−X port. The finding's core concern — that a shared construction error would be reproduced identically rather than caught — is now resolved by a path from premises to `nu_r` that does not reuse the SymPy variable choreography. The `$Assumptions` in scope at line 105 are the §3 set (reals + `ou2,ow2,r,gu,gw > 0`), appropriate for the simplification; the eps·lam coefficient of the log series is well-defined because `eps` and `lam` occur only as the product in `pA`/`dA`, and the exec confirms a finite correct result with exit 0.

### F2 — stale_output (banner label defect)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage178_outgoing_port_coloading_sympy_audit.py:34`: `banner("STAGE 161 …")` → `banner("STAGE 178 — OUTGOING-PORT CO-LOADING THEOREM")`.
- `mathematica/moving_throat_pde_stage178_outgoing_port_coloading_mathematica_audit.wl:26`: `banner["STAGE 161 …"]` → `banner["STAGE 178 — OUTGOING-PORT CO-LOADING THEOREM"]`.

Both edits match the directive exactly; no other change to either file.

**Assessment:**
Correct. Both refreshed saved transcripts now read `STAGE 178 — OUTGOING-PORT CO-LOADING THEOREM` at line 3 (sympy txt line 3, mathematica txt line 3). No assertion depends on the banner string; both scripts still exit 0. Pure provenance-label fix, math unaffected.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 178 — OUTGOING-PORT CO-LOADING THEOREM` (banner now correct).
- `nu_r - 2(p_r-d_r) = 0`, `Xi_1 - (nubar_N - kappa_1) = 0`, `p_r formula = 0`, `d_r formula = 0`, `alpha+beta-1 = 0`, `chi-zeta-1 = 0`, `nu_r - [kappa1 + sigma_r] = 0` — all residuals vanish.

**Mathematica:** exit=0. Notable lines:
- `STAGE 178 — OUTGOING-PORT CO-LOADING THEOREM` (banner now correct).
- `PASS: nu_r - [kappa1 + sigma_r]` (existing §4 check).
- `nu via log-data vs slippage = 0` / `PASS: nu via log-data vs slippage` — the new F1 independent check, residual zero.
- `nuFromData = (-2*ou2*(oU + oW)*ow2 + 4*r^2*rr)/(ou2*ow2 - r^2) + (2*(gw*(gW + oU)*ou2 + gu*r*(gU + rr)))/(gw*ou2 + gu*r)` — non-trivial, confirming non-vacuity.
- `Stage 178 Mathematica audit passed.`

**Output freshness:** confirmed. Script mtime `1780125319` for both; saved outputs newer — sympy txt `1780126749`, mathematica txt `1780126758`. Outputs regenerated post-fix. Both saved transcripts carry the corrected banner and (for `.wl`) the new `nu via log-data vs slippage` PASS line.

## Material-change assessment

`material_change`: false.

F1 adds an additional independent confirmation of the *same* `nu_r` identity (the residual is zero — no derived result changed). F2 is a banner-label string only. No constant, no symbolic result, no paper-stated value is altered. No downstream unit depends on either edit.

## Side observations (non-blocking)

- The git diff is tightly scoped: exactly the two named files, +11/−2 lines, consisting solely of the F1 insert and the two banner edits. No scope creep.
- The new check shares the raw inputs `pA`/`dA` with §3 (as the directive intended — §3 confirms those raw inputs have the claimed slopes `pExpected`/`dExpected`); independence is in the extraction route and the comparison target, which is the correct interpretation of the second-engine policy. Non-blocking.

## Verdict justification

Both findings are resolved. F1's new Mathematica-native `Log[pA^2/dA^2]` Series/Coefficient route is genuinely independent of the SymPy hand-form: it draws from the component-level port data, uses a distinct logarithmic extraction, and compares against the slippage assembly rather than the shared `pExpected`/`dExpected` factoring — so it is not an X−X port and is non-vacuous (it produces a non-trivial rational form whose residual against the slippage target cancels to zero). F2 corrected the misattributed banner in both engines and refreshed transcripts confirm it. Both engines exit 0, all in-file checks pass, outputs are fresh, the diff shows no regressions or scope creep, and nothing material changed downstream. Verdict: verified.
