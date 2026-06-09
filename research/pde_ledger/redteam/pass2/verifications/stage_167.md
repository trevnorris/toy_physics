---
unit_id: 167
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-08T20:05:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 167

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Codex added a self-contained "Independent numeric closure" block to the `.wl` only, inserted after the existing `expectZero["delta_perp", ...]` at line 87 and before the print-summary block (diff hunk `@@ -86,6 +86,59 @@`). The block defines `numericClosureResiduals[thetaN_, ksN_, kqN_, pN_, gN_, radiusN_]` (`wl:91-121`), a pure function that recomputes every primitive drift numerically (`drhoN, daN, dcsN, dZN, dcswN, dellN, dLWN`, then `dvN, dTN, dgqN, dgsN, dIsqN, dlamN`) from the ground-truth physical formulas and forms `rcN, rFrakN, gFrakN, chan1N, chan2N, deltaPerpN` from those locals. It then `Scan`s five tuples — `{1,0,0,0}`, `{0,1,0,0}`, `{0,0,1,0}`, `{0,0,0,1}`, and mixed `{2,-3,5,-7}` with `gstar->1/3`, `rstar->2` (`wl:123-140`) — and `expectZero`s all six residuals per tuple. No edit to the existing symbolic block; `.py` not in the diff.

**Assessment:**
Correct and matches the directive's primary (not minimal) option — the full multi-tuple numeric route. The block is genuinely ADDITIVE: every original symbolic check (A1–A15 carry-forwards, invariants, channels, delta_perp) remains intact and still prints PASS in the log. The numeric route is genuinely INDEPENDENT of the symbolic-assembly echo at the variable-binding level: `numericClosureResiduals` rebuilds its own local `drhoN…dlamN` and constructs `rcN/chan1N/chan2N/deltaPerpN` from those locals — it does NOT re-evaluate the module-level symbolic `chan1/chan2/deltaPerp/drc/dr/dg`, so the symbolic cancellation is exercised a second time via concrete-number substitution. The mixed tuple is non-trivial: a wrong assembly coefficient would leave a nonzero rational residual (auditor self-test confirms `dv=8`, `dlam=1`, yet `r_c=0`, `frak r=0`, `chan2=0`). Tuples 1–4 each light up exactly one drift, so the cancellation is not an origin artifact.

Honest caveat: this stage is pure first-order substitution, so the numeric route shares the SAME physical formulas as the symbolic route (it must — they are the Stage-165/166 carry-forward identities). The independence achieved is therefore "token-but-additive" / route-independence (symbolic-cancellation vs. concrete-substitution + recomputation), not algebraically-divergent re-derivation. A coefficient typo replicated identically in the numeric formulas would still escape. This limited achievable independence was explicitly user-accepted for this transliteration class, and the numeric block does add real standalone-evaluation value (catches errors that survive symbolic FullSimplify but fail on concrete numbers). Acceptable under that standard.

## Exec log assessment

**SymPy:** exit=n/a. No sympy script was changed (directive scoped the fix to the `.wl` only; `.py` is absent from the diff and last modified 2026-05-28), so no sympy re-run was captured. The original sympy assertions are unaffected.

**Mathematica:** exit=0. Notable lines:
- `PASS: delta_perp` (line 64) — original symbolic block intact.
- `STAGE 167 ... Independent numeric closure` banner (line 67) — new block ran.
- `delta ln r_c numeric (tuple 5) = 0` / `PASS: delta ln r_c numeric (tuple 5)` (lines 117–118) — mixed nonzero-drift tuple vanishes non-trivially.
- `delta_perp numeric (tuple 5) = 0` / PASS (lines 127–128) — full channel/normal-coordinate closure under concrete substitution.
- All 30 numeric residual lines print `= 0` with matching PASS; `Stage 167 Mathematica audit passed.` (line 139), `exit_code: 0`.

**Output freshness:** confirmed. `.wl` mtime 2026-06-08T16:07:21; committed `mathematica/output/...167...txt` mtime 2026-06-08T17:12:20 (newer than script), and contains 30 `numeric (tuple ...)` PASS lines — regenerated post-fix. SymPy `.py` (2026-05-28) and its `.txt` (2026-05-28) untouched, consistent with the no-change scope.

## Material-change assessment

`material_change`: false. The edit adds verification-only assertions (numeric `expectZero`s) and changes no derived constant, transport law, or invariant. Every primitive/branch/coupling value emitted is identical to pre-fix. No downstream unit (>167) consumes a value that this change altered.

## Side observations (non-blocking)

- The numeric block reuses the SAME closed-form primitive definitions as the symbolic path (unavoidable given the stage is pure substitution), so its independence is route-level not formula-level — see F1 caveat. Not a defect; called out for transparency.
- `FullSimplify` on already-concrete rationals (`wl:101-112`) is harmless overhead; the residuals are exact rationals that reduce to 0 regardless. No correctness impact.

## Verdict justification

The single low-severity `mathematica_transliteration` finding is resolved: Codex added an additive, route-independent numeric-substitution closure that recomputes all primitives from scratch and confirms every invariant/channel/δ⊥ vanishes for five sampled drift tuples (including a non-trivial mixed tuple), while leaving the entire original symbolic block — and the SymPy script — untouched. Mathematica exits 0 with all original and 30 new PASS lines; the committed output is fresh. Independence is genuine at the evaluation-route level but, given this is a pure-substitution stage, is the token-but-additive kind that was user-accepted for this class. No regressions, no material change. Verdict: verified.
