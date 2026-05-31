---
unit_id: 177
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-30T07:25:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 3
material_change: false
---

# Verification — unit 177

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Three edits, all matching the directive exactly (confirmed against `stage_177_diff.patch` and the current file state):

1. New independent factorization check inserted at `.wl:61-63` (between the `d ln H` check and the `Portwise outgoing-defect amplitude` banner):
   ```
   banner["Load-factor factorization (independent check)"];
   expectZero["load-factor factorization lambda0^2/k = M^2 (1+I)^2/(1-H)^2",
     lambda0^2/k - mCal^2*(1 + iCal)^2/(1 - hCal)^2];
   ```
2. Per-port amplitude re-anchored to the factored form at `.wl:66-68`. The old single line `sigmaExact = FullSimplify[D[Log[(lambdaP^2/kP)/(lambda0^2/k)], eps]...]` was replaced with a build from the perturbed `(M,I,H)` factors:
   ```
   lambdaSqOverKP = mCalP^2*(1 + iCalP)^2/(1 - hCalP)^2;
   lambdaSqOverK0 = mCal^2*(1 + iCal)^2/(1 - hCal)^2;
   sigmaExact = FullSimplify[D[Log[lambdaSqOverKP/lambdaSqOverK0], eps] /. eps -> 0, ...];
   ```
3. Banner string corrected in both scripts: `.wl:26` and `.py:32` now read `STAGE 177` (was `STAGE 160`).

**Assessment:**
All three edits are correct and complete; nothing beyond what the directive named was touched (diff is exactly the 3 prescribed hunks). 

- The factorization check is non-tautological. `lambda0 = (ou2*gw + r*gu)/(ou2*ow2 - r^2)` and `mCal/iCal/hCal` are defined independently at `.wl:32-35`; the residual `lambda0^2/k - mCal^2(1+iCal)^2/(1-hCal)^2` is only zero because the algebra closes (verified: `mCal^2(1+iCal)^2/(1-hCal)^2 = (gw^2/(k ow2^2))((ou2 gw + r gu)/(ou2 gw))^2 / ((ou2 ow2 - r^2)/(ou2 ow2))^2 = (ou2 gw + r gu)^2/(k (ou2 ow2 - r^2)^2) = lambda0^2/k`). It would fail if `lambda0` were defined inconsistently with the `(M,I,H)` decomposition — exactly the setup error the second-engine policy is meant to catch. The output log line 20-21 shows `= 0` / `PASS`.
- The re-anchored amplitude is a genuinely different intermediate route, not an echo of the SymPy raw-`lambdaP` path. SymPy differentiates `log(lambdaP^2/kP / (lambda0^2/k))`; Mathematica now differentiates `log(M_P^2(1+I_P)^2/(1-H_P)^2 / [M^2(1+I)^2/(1-H)^2])`. The two routes reach the same `sigma_r` only through the factorization (independently checked in edit 1), so the engines now cross-check rather than transliterate. The target `lam*sigmaR` is unchanged, so the equality is still a real constraint on the slope coefficients — non-tautological. Output line 26-27 confirms `Sigma_{A,r} = lambda_A sigma_r = 0` / `PASS`.
- Banner fix is cosmetic and correct in both engines (sympy log line 8, mathematica log line 8).

### F2 — insufficient_verification (informational)

**Classification:** resolved (informational; no edit expected or applied)

**What changed:** Nothing, as directed. The directive (and report) explicitly instruct NOT to add a `sum_r rho_r sigma_r == Xi_1` collapse assertion because it reduces to linearity of a sum (trivial-by-linearity), i.e. any such check would be tautological. The diff confirms no edit was made to the block 3 / `.py:91-100` / `.wl:66-74` region.

**Assessment:** The "no edit" disposition is sound. I independently confirm the reasoning: with the per-port amplitude `Sigma_{A,r} = eps*lam*sigma_r` already verified (block 2), the grouped collapse `Xi_load = sum_r rho_r Sigma_r = eps*lam*sum_r rho_r sigma_r = eps*lam*Xi_1` is a pure linearity step over a weighted sum and carries no adversarial signal. Adding it would create a passing-but-tautological assertion, which would itself be a `partial` outcome. Leaving it documented is the correct call.

### F3 — insufficient_verification (informational)

**Classification:** resolved (informational; no edit expected or applied)

**What changed:** Nothing, as directed. The prefactor-slope block (`.py:111-122` / `.wl:91-101`) was left untouched (confirmed: that region does not appear in the diff).

**Assessment:** The "no edit" disposition is sound. The block defines `xiLoad = n1 - d1` and `pA = n0 e^{eps lam n1}/(d0 e^{eps lam d1})`, then asserts `pSlope/p0 - xiLoad == 0`. The order-eps log-slope of that quotient is `n1 - d1` by construction, so the assertion cannot fail — it is the definitional "log of a quotient = difference of logs." Strengthening D6 would require re-deriving Stage 241 (out of scope) or the tautological F2 sum. Documenting it as a known non-adversarial definitional identity, rather than dressing it up, is correct; adding nothing avoids manufacturing a false sense of verification.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 177 — WEAK-AXISYMMETRIC OUTGOING-SLIPPAGE COLLAPSE` (banner fixed)
- `weak-axisymmetric d ln M/I/H = 0`
- `Sigma_{A,r} = lambda_A sigma_r = 0`
- `grouped trace vanishes = 0`, `a_Xi - eps Xi1/4 = 0`, `b_Xi - 3 eps Xi1/4 = 0`, `b_Xi - 3 a_Xi = 0`

(SymPy was not given the factorization check by the directive — F1 part 1/2 are Mathematica-only; SymPy received only the banner fix. This is consistent with the directive, which deliberately differentiates the engines.)

**Mathematica:** exit=0. Notable lines:
- `STAGE 177 — WEAK-AXISYMMETRIC OUTGOING-SLIPPAGE COLLAPSE` (banner fixed)
- `load-factor factorization lambda0^2/k = M^2 (1+I)^2/(1-H)^2 = 0` / `PASS` (new independent check)
- `Sigma_{A,r} = lambda_A sigma_r = 0` / `PASS` (now via factored route)
- all trace/anomaly and prefactor checks `= 0` / `PASS`
- `Stage 177 Mathematica audit passed.`

**Output freshness:** confirmed fresh. Both committed `.txt` outputs were regenerated post-fix:
- sympy script mtime `2026-05-30 01:09:47`, sympy output mtime `2026-05-30 01:38:56` (newer)
- mathematica script mtime `2026-05-30 01:09:47`, mathematica output mtime `2026-05-30 01:39:04` (newer)
The committed mathematica `.txt` contains the new `load-factor factorization ... = 0` / `PASS` lines and the `STAGE 177` banner; the committed sympy `.txt` carries the `STAGE 177` banner. Outputs match the captured exec logs.

## Material-change assessment

`material_change`: false.

All edits are additive or cosmetic. The new factorization check verifies an identity (`lambda0^2/k = M^2(1+I)^2/(1-H)^2`) that was already implicitly assumed; it changes no derived result. The re-anchored amplitude reaches the identical `sigma_r` through a different route and the assertion target (`lam*sigmaR`) is unchanged. The banner fix is a label only. No downstream-consumable quantity (slopes, sigma_r, lane pattern, anisotropy values, Xi_1 = P1/P0 tie) was altered. No downstream units need re-auditing on account of unit 177.

## Side observations (non-blocking)

- The factorization identity verified in F1 is the bridge the appendix (eq:app-part05-load-factor-factorization) asserts; with edit 2 the Mathematica `sigma_r` no longer touches the raw `lambdaP` symbol (`.wl:47`) at all, relying instead on the factored form plus the new factorization check to tie back to `lambda0`. This is the intended differentiation and is sound, but note the raw-`lambdaP` definition at `.wl:47` is now used only by the SymPy-mirrored... actually it is now unused in Mathematica except as documentation. Not a defect; flagging for awareness only. Do not block on this.
- The two engines now genuinely differ in the amplitude derivation route (raw quotient in SymPy vs factored-(M,I,H) in Mathematica), which is the strongest practical mitigation of the transliteration concern available within this stage's symbol scope.

## Verdict justification

F1 was applied exactly as the directive prescribed, with no collateral edits (diff = the 3 named hunks). The new Mathematica factorization check is non-tautological and passes, the re-anchored per-port amplitude now derives `sigma_r` through a genuinely independent factored route while still constraining it against the unchanged `lam*sigmaR` target, and both banners read STAGE 177. F2 and F3 were correctly left as informational no-edit findings: I independently confirm that any collapse or prefactor-slope assertion would be tautological-by-linearity / definitional and that documenting them (rather than adding a passing-but-empty check) is the right disposition. Both engines exit 0 with every in-file check passing, and both committed outputs are fresh and consistent with the exec logs. Verdict: verified.
