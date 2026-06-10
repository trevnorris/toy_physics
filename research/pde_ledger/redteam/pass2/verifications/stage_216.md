---
unit_id: 216
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 216

The sole finding in the original report is non-script: F1 = `paper_misalignment` (the card's `\stagefield{Verification}` line reads "Mathematica audit: none yet" while a passing `.wl` exists). The directive directs Codex to apply NOTHING and routes the card-text-lag to the user (paper-side, P4-51 deferral). No source edit was directed and none was needed; consequently there is no `stage_216_diff.patch` (correctly absent). Verification confirms (A) outputs clean/fresh, (B) the audit disposition holds on the refreshed artifacts, (C) the card-text-lag deferral is correctly classed, and (D) `material_change: false`.

## Per-finding outcomes

### F1 — paper_misalignment (paper_missing_script_claim, stale verification metadata)

**Classification:** resolved (correctly routed — non-script, USER-DEFERRED to P4-51)

**What changed:**
Nothing in scripts. `paper/stages/stage_216.tex` is off-limits to the red-team (scripts-only scope). The directive (`directives/stage_216.md:14-17`) explicitly states there are NO non-paper_misalignment findings to apply and that "Codex applies nothing on this unit until the user chooses a direction." The directive carries only the resolve-block (directions (a) name the `.wl` as stage 218 does / (b) leave the card) and no `## Applied:` block — consistent with zero Codex invocation on this unit. The git status shows no modification to any `scripts/`, `scripts/output/`, `mathematica/`, or `mathematica/output/` stage-216 file (only the new untracked directive), confirming no source edit landed.

**Assessment:**
Correctly classified and correctly deferred. The `.wl` exists, is dated 2026-06-02, and passes all six check families ("All Stage 216 Mathematica audit checks passed.", math log L79), so the card understates verification coverage but contradicts no math result. This is a documentation-metadata staleness from the pass-1 dual-engine retrofit, not a math disagreement. Routing it to the user (P4-51 paper-cleanup) is the right disposition — it lies outside the scripts-only red-team scope and any fix is a paper-owner decision between directions (a) and (b). Non-blocking.

## Disposition re-confirmation (independence + non-tautology, on the refreshed artifacts)

- **Genuinely independent `.wl`:** confirmed. The `.wl` derives every load-bearing object by a different route than the `.py`'s posit-and-verify choreography. M1/M2 gradient optimum: derived via the constrained Lagrange stationarity system (`Solve[Join[stationarity, {aVec.aVec==1}], …, Reals]`, positive-μ branch) — math log L14-21 prints the derived `{kLam/Sqrt[S], …}` point and `Sqrt[S]` slope with positive-branch dominance, vs the `.py` positing `a_grad = k/‖k‖` (sympy log L13-20). M4 leverage/slack: encoded as matrix quadratic forms `aᵀ(J−I)a` / `aᵀ(nI−J)a` (math log L50-53) vs the `.py`'s explicit-monomial sum (sympy log L29-30). M5 bound: proved spectrally via `Eigenvalues[J−I]={4,−1,−1,−1,−1}` and `Max[eig]−4=0` (math log L58,63) — a strictly stronger, different argument than the `.py`'s single-point substitution. M6 bracket: solved via `Solve[…==0, tau, Reals]` selecting the smaller positive root (math log L69-77) vs the `.py`'s posit-then-verify-quadratic. Routes differ at the construction level (derive-vs-posit, matrix-vs-monomial, spectral-vs-point); no transliteration.
- **Non-tautological:** confirmed. Each posit-verify residual CAN fail (a wrong `a_grad` would make `k_grad−k_norm≠0`; the `½κτ²−kτ+H0` residual is not `0=0` by construction). The leverage identities expand to genuine polynomial cancellations. The `.wl`'s eigenvalue/Lagrange-Solve checks are independent derivations, not echoes.
- **0 value misalignments:** confirmed. The report reconciled 14 deliverable values (symbolic forms + the 4/3/2/1 ladder + structural integers) against notes/appendix, all MATCH; the stage emits no free numeric constants of its own (the Bézout budget 162/324/750/1500 belongs to 217/218). The refreshed outputs carry the same forms in both engines: gradient point, max slope `√S`, face gaps `k_axis²`, leverage/slack identities `=0`, barycenter leverage `4` (+ top eigenvalue `4` in `.wl`), bracket `2H0/(k+√(k²−2H0κ))`.

## Exec log assessment

**SymPy:** exit=0. Notable lines: `||a_grad||^2 = 1` (L14); `w_sigma - ((sum a_i)^2 - sum a_i^2) = 0` (L29); `5*sum a_i^2 - (sum a_i)^2 - sum_(p<q)(a_p-a_q)^2 = 0` (L30); `w_sigma(a_eq) = 4` (L32); certified bracket `tau = 2*H0/(k + sqrt(-2*H0*kappa + k**2))` (L38); support ceiling `5` (L46-48). Every identity resolves to `0`/expected form; `# exit_code: 0` (L50).

**Mathematica:** exit=0. Every check prints `PASS` and no `FAIL`: `PASS: M1 gradient-optimal unit norm` (L17); `PASS: M2 positive Lagrange value minus Sqrt[S]` (L19); `PASS: M4 off-diagonal quadratic identity` (L51); `off-diagonal leverage eigenvalues = {4, -1, -1, -1, -1}` (L58) with `PASS: M5 top quadratic-form eigenvalue minus 4` (L64); `PASS: M6 solved smaller root minus stated bracket form` (L71). Closes "All Stage 216 Mathematica audit checks passed." (L79); `# exit_code: 0` (L80).

**Output freshness:** confirmed fresh. Both `.txt` outputs carry mtime 2026-06-09 16:51:54, newer than their scripts (`.py` and `.wl` both 2026-06-02 11:23). The orchestrator's independent re-run refreshed them (216's outputs were already current — likely byte-identical), and git shows no modification to the committed outputs. No `stage_216_diff.patch` exists, consistent with zero Codex source edits.

## Material-change assessment

`material_change`: false. No source code changed (no Codex edit was directed; the diff patch is correctly absent). The only filesystem change is the regenerated committed `.txt` outputs (transcript refresh, no derived-result change) and the new untracked red-team directive. No downstream unit (> 216) is affected.

## Side observations (non-blocking)

None beyond the single reported finding. The auditor correctly did not flag M3 (face gap = missing square) as transliteration despite the shared trivial subtraction — it is a corollary, not a load-bearing object. I concur and add nothing.

## Verdict justification

`verified`. The lone finding is a low-severity stale-metadata `paper_misalignment` (card says "Mathematica audit: none yet" while a passing `.wl` exists), correctly classed as a card-text-lag and USER-DEFERRED to P4-51; the directive directs Codex to apply nothing and no source edit landed (no diff patch, git status clean on all script/output files). Both engines pass on the refreshed artifacts (sympy exit 0, mathematica exit 0, every check PASS, no FAIL), and the saved `.txt` outputs are fresh (mtime Jun 9 16:51 > script mtime Jun 2 11:23). The audit disposition holds: the `.wl` is genuinely independent (Lagrange-Solve / matrix quadratic forms / `Eigenvalues` / quadratic-Solve, none mirroring the `.py`'s posit-and-verify), the assertions are non-tautological, and reconciliation is 14/14 MATCH with 0 misalignments. No regressions; `material_change: false`. No wrong-root path written.
