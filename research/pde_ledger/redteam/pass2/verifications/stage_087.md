---
unit_id: 087
batch: III.5
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 087

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
- SymPy (`...stage087..._sympy_audit.py`): the three `expect_close("rho_X vs stage-086", ...)` self-comparisons (former py:73-75) were replaced with structural ordering/gap asserts (py:75-87): `rho_suff < rho_fail` (py:75,80-81), `rho_fail < rho_max` (py:76,82-83), and `0 < rho_max - rho_fail < 1e-6` (py:77-78,85-87), each printing its boolean and raising `AssertionError` on violation. The literals are still printed (py:72-74). The overclaiming inline comment was reworded to drop the "cross-check against upstream stage-086" framing (py:69-71).
- Mathematica (`...stage087..._mathematica_audit.wl`): the three `expectApprox["rho_X vs stage-086", ...]` self-comparisons (former wl:61-63) were replaced with the same ordering/gap checks via `If[!TrueQ[...], fail[...], pass[...]]` (wl:64-70). The comment block (wl:53-55) was reworded the same way. The three `zeta_*` numeric checks (wl:80-82) were re-anchored from hardcoded 25-digit decimals to `rhoSuff - 1` / `rhoFail - 1` / `rhoMax - 1` with tol `10^-20`, genuinely testing the `epsBlk->0` substitution of `zetaReq`.

**Assessment:**
The fix correctly matches the directive's required change (option (a) applied symmetrically to both engines), and the new checks are genuinely can-fail rather than new tautologies:
- The ordering asserts compare three DISTINCT literals against each other, not a literal against a re-typed copy of itself. Mistyping any one literal can flip an ordering or the gap-magnitude outcome, so the assert can fail. This is the structural-drift guard the original comment falsely claimed.
- The gap assert `0 < rho_max - rho_fail < 1e-6` is doubly load-bearing: it would fail if rho_max ≤ rho_fail (sign/ordering error) OR if the gap exceeded 1e-6 (a coarse mistype). Output confirms the real gap is 9.67e-8, comfortably inside the band but not at either edge, so the check is meaningful, not vacuous.
- The Mathematica `zeta_*` re-anchoring is now non-tautological in the right way: `zetaSuff` is computed by substituting `rhoSuff, epsBlk->0` into the actual `zetaReq` rational function, while the target `rhoSuff - 1` is the independent algebraic prediction of the unblocked reduction. A bug in `zetaReq` that failed to vanish the epsBlk term at epsBlk=0 would make these diverge and the assert fail. (Both still print `diff = 0` because the reduction is correct — that is a passing genuine check, not a self-comparison.)

The load-bearing asserts are untouched and still pass: `unblocked zeta_req` (py:55, wl:44) and `d zeta_req exact formula` (wl:43) are byte-identical to pre-fix and PASS in both outputs (sympy txt:7; math txt:7-10). The three window literal VALUES are unchanged: py:58-60 and wl:57-59 remain `3.46622291347846 / 3.46752913273870 / 3.46752922945601`. No collateral edits beyond the two targeted blocks; the diff touches only the two script files.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `unblocked zeta_req = 0` (load-bearing reduction still passes)
- `rho_suff < rho_fail = True` / `rho_fail < rho_max = True`
- `0 < rho_max - rho_fail < 1e-6 = True ; gap = 0.0000000967173...` (real gap inside band, not at edge)

**Mathematica:** exit=0. Notable lines:
- `PASS: d zeta_req exact formula` and `PASS: unblocked zeta_req` (both load-bearing asserts intact)
- `PASS: rho_suff < rho_fail`, `PASS: rho_fail < rho_max`, `PASS: 0 < rho_max - rho_fail < 10^-6`
- `PASS: zeta_suff = rho_suff - 1` (re-anchored zeta check passes; diff = 0 because reduction is exact)

**Output freshness:** confirmed. Both committed `.txt` outputs mtime 2026-06-05 18:15:21; both scripts mtime 2026-06-05 18:07:30. Outputs are newer than scripts and their content matches the post-fix scripts (ordering/gap lines present, no `vs stage-086` self-comparison diff lines remain).

## Material-change assessment

`material_change`: false.

No deliverable value changed. The symbolic one-ratio finish line (`rho_alpha = alpha_req/alpha_mix`) and the three window literals (`3.46622291347846 / 3.46752913273870 / 3.46752922945601`) are all unchanged. The edit only swapped a hollow verification mechanism for a genuine one. No downstream unit can be affected by a stricter-but-passing self-check.

## Side observations (non-blocking)

- The SymPy top-of-file docstring (py:17-20) still says the script "cross-checks the Family-1 window literals against the upstream stage-086 quoted values to catch renumber or transcription drift." The F1 fix correctly reworded the directive-cited inline comment block (py:69-72) and the Mathematica docstring overclaim was not present in the same form, but this SymPy docstring sentence retains the original overclaiming framing. Minor residual prose that escaped the rewording; the actual mechanism is now genuine so this is cosmetic, non-blocking.
- Pre-existing stale cross-references survive (and are correctly out of scope per numbering policy): py:12 "former stage 65", py:13 "former stage 69 closure", py:92 / output line 21 "Stage-69 ratio window" (now stage 086). Not introduced by this fix; deferred to the dedicated numbering pass.

## Verdict justification

The single finding F1 is genuinely resolved in both engines. The hollow each-literal-against-its-own-copy comparisons were replaced with structural ordering/gap relations (three distinct literals) and a re-anchored zeta substitution check — all of which can fail under mistranscription and all of which pass with the correct literals. The load-bearing `unblocked zeta_req` and `d zeta_req exact formula` asserts are byte-untouched and still pass; the three window literal values are unchanged; the diff is confined to the two scripts (no paper/notes, no scope creep); and both scripts exit 0 with freshly regenerated outputs. No deliverable value changed, so material_change is false.
