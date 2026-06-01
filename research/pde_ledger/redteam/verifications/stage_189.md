---
unit_id: 189
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-01T12:10:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 189 (iteration 2 re-verify)

This OVERWRITES the iteration-1 verification, which wrongly passed the Section II
selected-branch product as a genuine check. Iteration 1's residual tautology was:
`Rtarget_oneport = Lambda0*(1-epseta)/T2_direct` back-defined R_target, so
`Rtarget_oneport*T2_direct - Lambda0*(1-epseta)` was identically 0 (definitional,
rank-2 compatibility — it DEFINES R_target and cannot verify it). Iteration 2 (per
Claude+Codex consult 019e843e, CONCUR) demoted that identity to a printed definition
and added a genuine direct-slope bridge in both engines. This re-verify focuses on the
iter-2 delta (F1 re-fix + F4 mirror); F2 and F3 were already resolved in iter-1 and are
re-confirmed unchanged.

## Per-finding outcomes

### F1 — tautological_check (Section II selected-branch identity) — ITER-2 RE-FIX

**Classification:** resolved

**What changed:**
- SymPy `scripts/...sympy_audit.py:89-136`. The tautological assertion is GONE. Lines
  100-107 now only PRINT three things: `T2_direct` (the one-port continuum form), the
  R_target DEFINITION `Rtarget_definition = Lambda_0(1-epseta)/T2_direct` (a `sp.pprint`,
  NOT an `expect_zero`), and the explicit `Lambda0_val = 27 pi^2 G c_s^5/(20 a^5 c^5)`.
  There is no `expect_zero("R_target * T_A^2 - Lambda_0 (1-epseta)", ...)` anywhere.
- The new genuine bridge is lines 109-136. The concrete coherent continuum shape
  `T2_coh = Zw(1+chi0)^2/(OmW2(1-epssplit)^2)` (line 111, with
  `epssplit = epsW*(1 - 2/11 deltaU/(1+deltaU))`) is perturbed along a path parameter `s`
  (`T2_coh_pert`, lines 126-134), and the first-order log-slope
  `direct_slope = diff(log(T2_coh_pert/T2_coh), s)|_{s=0}` (line 135) is asserted equal to
  `epsA*lambdaA*Xi1_closed` (line 136).

**Assessment — is the bridge genuinely non-tautological? (THE KEY QUESTION):**
Yes. I traced exactly how `Xi1_closed` is built (lines 119-125):
- `eps1 = diff(epssplit, epsW)*epsW1 + diff(epssplit, deltaU)*deltaU1` (line 119) — the
  first-order variation of the `epssplit` MAP via chain rule on `epssplit`'s own
  definition. It depends on `epssplit`, NOT on `T2_coh` or any of its derivatives.
- `Xi1_closed = zetaZ - omegaW + 2*chi1/(1+chi0) + 2*eps1/(1-epssplit)` (lines 120-125) —
  a HAND-WRITTEN grouped closed form assembled from the independent input perturbation
  amplitudes (`zetaZ`, `omegaW`, `chi1`) and `eps1`. It is the analytically-written
  log-derivative of the transfer shape's component structure.

  `Xi1_closed` does NOT reference `direct_slope`, `T2_coh_pert`, `s`, or any derivative of
  `T2_coh`. It is built BEFORE the perturbation path is even constructed (lines 126-135).
  The two sides of line 136 are computed by genuinely independent routes: the LHS
  (`direct_slope`) is SymPy machine-differentiating `log(T2_coh_pert/T2_coh)`; the RHS is
  the hand-grouped `Xi1_closed`. The check passes ONLY because the concrete `T2_coh`'s
  actual log-slope equals the independently-grouped expression. It would FAIL if the
  continuum closed form were wrong (e.g. wrong power on `(1-epssplit)`, wrong factor of 2
  on `chi1`, wrong `epssplit` coefficient 2/11) OR if `Xi1_closed`'s grouping were wrong.
  This is the stage-181 pattern (`scripts/...stage181...:76` eps-split + Xi1 def, then
  the Xi1-vs-T^2-derivative check). Non-tautological. Confirmed.

  (The input amplitudes `zetaZ, chi1, omegaW, epsW1, deltaU1` appear in BOTH the
  perturbation path and in `Xi1_closed` — this is correct and necessary: they are the
  shared independent perturbation inputs. The test is whether the slope's functional
  grouping of those shared inputs matches the hand-written grouping, which is a real
  algebraic content check, not a renaming.)

- Mathematica mirror `mathematica/...mathematica_audit.wl:101-130`: SAME de-tautologization
  (lines 105-108 PRINT `t2OnePort`, `rTargetDefinition`, `lambda0Value`; no `Rtarget*T2`
  assertion) and SAME genuine bridge. `xi1Closed` (lines 117-120) is built from
  `eps1Bridge = D[epsSplit,epsW]*epsW1 + D[epsSplit,deltaU]*deltaU1` and the hand grouping,
  independent of `directSlope` (line 126). Crucially the `.wl` builds `t2Coherent` by
  SUBSTITUTION `t2OnePort /. {rho->chi0, epsW->epsSplit}` (line 111) rather than re-typing
  the closed form (SymPy re-types it), and writes `epsSplit` in a different surface form
  (`2*deltaU/(11*(1+deltaU))` vs SymPy's `2/11*deltaU/(1+deltaU)`) — genuinely independent,
  not a transliteration.

### F2 — tautological_check (Section I defect-notation/compatibility) — re-confirmed (iter-1)

**Classification:** resolved

**What changed (unchanged since iter-1):**
SymPy lines 65-87 introduce independent defect symbols `Theta1, Xi1, Sigma_eta` and the
observable substitution `obs_sub = {dr: Theta1, dn: Xi1 + Bstar*Theta1, de: Sigma_eta}`,
then check `dlnT2.subs(obs_sub) - Xi1 == 0`, `dlnOneMinus.subs - (R1 + dlnT2.subs) == 0`,
and the derived compatibility relation. The old `Xi - (dn - Bstar*dr)`, `Rcal - dlnRtarget`,
`(Rcal+Xi) - dlnOneMinus`, and the line-63 `trf[2]+trf[0]-trf[1]` compatibility-row checks
are gone (only the genuine matrix check at line 62 remains). The `.wl` mirror uses a
different route entirely: a finite-log Jacobian `D[transferLogFunctions, obsVars]/.zeroObs`
for the matrix (lines 65-83) and exponential paths `Exp[s*theta1]` for the transfer
identities (lines 85-99). Non-tautological in both engines.

### F3 — paper_misalignment (stale STAGE 172 banner) — re-confirmed (iter-1)

**Classification:** resolved

**What changed (unchanged since iter-1):** Authorized relabel per settled
canonical-stage-number convention (direction a). SymPy banner line 35 reads
`"STAGE 189 — TRANSFER-SHAPE / OUTGOING-PREFACTOR COMPILER"` and line 223 reads
`"STAGE 189 LEDGER"`; output `.txt` prints "STAGE 189" throughout. The "Stage 188"
upstream references are left as-is (188 is canonical branch-observable source). String-only;
no math change. The notes-prose legacy numbers (240/239) are out of red-team scope.

### F4 — missing_mathematica — re-confirmed + ITER-2 mirror

**Classification:** resolved

**What changed:** The `.wl` exists and independently re-verifies every load-bearing
assertion of the corrected SymPy script via different native primitives: finite-log
Jacobian (Section I matrix), exponential perturbation paths (Section I identities),
substitution-built coherent shape + path-derivative slope (Section II bridge), `D[...,{w,n}]`
Taylor-coefficient extraction (Sections III/V `u2,u4,P0,P2,P4`), and `Solve` for the
constant-prefactor branch `N2_const,N4_const` (Section VI). It is NOT a line-by-line port:
different variable names, different derivation routes (derivatives/Coefficient/Solve vs
SymPy's `series`), different surface forms. The iter-2 mirror correction (printed R_target
definition + non-tautological slope bridge) is present. All conclusions agree with the
SymPy engine.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `dln T^2 - Xi_1 = 0`, `dln(1-epseta) - (R_1 + Xi_1) = 0`, `compatibility: ... = 0`
- `direct-slope bridge dln T_A^2 - epsilon lambda_A Xi_1 = 0` (the genuine iter-2 bridge)
- The old `Lambda_0 (1-epseta) / R_target - T_A^2` line and the iter-1
  `R_target * T_A^2 - Lambda_0 (1-epseta)` product line are BOTH absent.
- Reaches "STAGE 189 LEDGER".

**Mathematica:** exit=0. Notable lines:
- `PASS: direct-slope bridge dln T_A^2 - epsilon lambda_A Xi_1` (mirror of the genuine bridge)
- `PASS: coherent local D/N specialization` (confirms substitution-built shape matches)
- All Sections I–VI PASS; ends with `Exit[0]`. No `Rtarget*T2` product assertion present.

**Output freshness:** confirmed. Both `.txt` outputs (mtime 2026-06-01 11:53:34) are newer
than their scripts (sympy `.py` 11:48:06; `.wl` 11:48:34) — regenerated post iter-2 fix.

## Material-change assessment

`material_change`: false. The iter-2 edits demoted a definitional identity to a print and
added a new (passing) verification; no derived numeric/symbolic RESULT that a downstream
unit consumes changed. The verified compiler formulas, constants (`Lambda_0`, `A`, `B`,
`G_5`, `P0_target`), and identities are unchanged. No downstream re-audit needed on account
of this fix.

## Side observations (non-blocking)

- In Section II the input perturbation amplitudes (`zetaZ, chi1, omegaW, epsW1, deltaU1`)
  are shared between the perturbation path and `Xi1_closed`. This is correct and is what
  makes the slope a function of physical inputs; flagging only to note it was scrutinized
  and is not a hidden circularity — `Xi1_closed`'s GROUPING of those inputs is the load.
- A15 (`Gamma5 - a^5 P0/(27 c_s^5)`, line 212) remains a mild self-restatement
  (Gamma5 := G5out*P0 := a^5 P0/(27 c_s^5)). The original auditor flagged this as "mild"
  and did not raise a finding; not in scope to introduce one here.

## Verdict justification

`verdict: verified`. The iteration-1 residual tautology is eliminated in BOTH engines: the
selected-branch product `expect_zero` is gone, replaced by a printed R_target definition.
The new direct-slope bridge is genuinely non-tautological — `Xi1_closed` is assembled by
hand from the independent `epssplit`-map variation and the input perturbation amplitudes,
built before the perturbation path and never touching `direct_slope` or the `T2_coh`
derivative, so the check passes only because the concrete continuum closed form's actual
log-slope equals the independently-grouped expression and would fail on a wrong closed form
or wrong grouping. The `.wl` mirrors the correction via a genuinely independent route
(substitution-built coherent shape, finite-log Jacobian, exp-paths, derivative/Solve
compilers) and agrees with SymPy on every conclusion. Both engines exit 0; outputs are
fresh. F2/F3/F4 remain resolved from iteration 1. All four findings resolved.
