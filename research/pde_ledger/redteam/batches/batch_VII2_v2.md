---
batch: VII.2
range: 231-242
total_stages: 12
verified: 12
findings_count: 18
findings_resolved: 18
findings_blocked_legitimate: 0
material_change_count: 0
clean_stages: []
status_only: []
dirty_stages: [231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242]
checkpoints: [239, 242]
consult: none (238 F1/F4 reframed via Claude math-coverage resolution, non-conceptual, no paper edit)
audit_date: 2026-06-03
verify_date: 2026-06-03
status: closed
---

# Red-team batch VII.2 — Rigid-mouth orbit-lock / branch-dressing / twin-support

## Summary

12-stage audit unit for VII.2 (`Part VII — Mixed-bundle / resonance /
branch-packet`, continuing into the rigid-mouth orbit-lock / branch-dressing /
twin-support physics), the **third batch of Part VII**, forward first-pass under
the v2 paper-grounded auditor **WITH the dual-engine rule** in force. **Two
checkpoints this batch — 239 and 242** (both at the higher verify bar); no
status-only units. All 12 reached `verified`; **`material_change: false` on all
12**; **0 stop-cold, 0 ultimately-blocked, 0 needs_rework left open.** Every
change is a strengthening / route change / notes typo correction; no derived
value, constant, identity target, or PUBLISHED paper number on the SCRIPT side
moved, so no `upstream_stale` propagation.

**11 of 12 codex-invoke runs exited 0 on iteration 1; one iteration-2 (stage
238 — see the reframe narrative below).**

**Cumulative: 230/253 → 242/253 stages red-team verified (95.7%)**; the entire
range 001–242 is now paper-aligned at v2 depth. Stages 243–253 remain `pending`.

### Headline — full dual-engine coverage with zero sanctioned mirrors

VII.2 ran with the dual-engine rule (a Mathematica `.wl` is **REQUIRED wherever
Mathematica CAN independently verify** — the test is "is it possible," not "is
it necessary") in force. **Every stage now has an INDEPENDENT second engine; 0
sanctioned mirrors were accepted in VII.2.**

- **10 stages had NO `.wl`** and got a **NEW independent-route `.wl`**: 231, 232,
  233, 234, 235, 236, 237, 238, 240, 241.
- **2 checkpoint `.wl` were caught as line-by-line transliterations and
  RE-AUTHORED to genuinely independent routes**: 239 and 242.

Every new `.wl` was confirmed independent by a clean verify agent — native
primitives via a DIFFERENT decomposition than the SymPy `.py`. Representative
routes: **239** forward Jacobian of the boxed dependent vector + native
`PseudoInverse` left-inverse + `Reduce`/`Equivalent` orbit-lock (vs the `.py`'s
backward hardcoded `SrmDep`); **242** `Resolve[ForAll,Reals]` strict-inequality
certificate + `D[]` on real closed forms (vs the abstract-ζ device) +
`logDrift` total-log-differential (vs `Exp[t·d]`).

The labor split was strictly enforced: **Claude reviews** (audit + verify);
**Codex writes ALL script code** (designs and writes the new/re-authored `.wl`);
the directives stated only the requirement + acceptance criteria, never script
code.

## Checkpoint findings (the higher bar, no rubber-stamp)

Both VII.2 checkpoints had EXISTING `.wl` that were line-by-line
transliterations of their `.py`; both were RE-AUTHORED to independent routes and
both cleared the higher bar.

### 239 — rigid-mouth physical normal form / Cartesian orbit-lock

Cleared the higher checkpoint bar:

- **`.wl` RE-AUTHORED** from a transliteration to a genuinely independent route:
  a forward Jacobian of the boxed dependent vector + a native `PseudoInverse`
  left-inverse + a `Reduce`/`Equivalent` orbit-lock, versus the `.py`'s backward
  hardcoded `SrmDep`. The orbit-lock is now established by two structurally
  different mechanisms across engines.
- In-file stale-label cosmetics rode the fix loop: `.wl` banner `STAGE 222→239`
  and `Stage221→Stage238` suffixes (single-file, NOT a batch renumber).

**Checkpoint bar MET.**

### 242 — twin-support strict inclusion / branch-dressing drift

Cleared the higher checkpoint bar:

- **`.wl` RE-AUTHORED** from a transliteration to a genuinely independent route:
  a `Resolve[ForAll,Reals]` strict-inequality certificate + `D[]` on the real
  closed forms (versus the abstract-ζ device) + a `logDrift`
  total-log-differential (versus `Exp[t·d]`).
- The load-bearing **twin-window strict inclusion `C_mix < Pi_tr < 2 C_mix`** is
  now tested STRICTLY on BOTH engines (the `ρ_α`-style ratio 4/3 region).
- In-file stale-label cosmetic rode the fix loop: `.wl` banner `STAGE 225→242`
  (single-file, NOT a batch renumber).

**Checkpoint bar MET.**

## The dominant defect theme — the "variable-independence self-test trap"

The recurring VII.2 defect was support-blindness / independence "verified" by
differentiating an expression w.r.t. variables it never contains (vacuous). Hit
at **237-F2, 238-F1, 240-F1**. Fixes:

- **237 / 238** — negative control + leak-detector + structural exclusion. 237-F2
  was fixed with a live-channel negative control + a `Not[FreeQ[#,Derivative]]`
  guard; 238 followed the same shape (see the iter-2 narrative).
- **240-F1** — extract-the-weight-from-the-variable-bearing-object: the weights
  are now pulled from `Y_support`, which carries `Omega_Q` via its pole, so the
  independence claim is exercised against a variable-bearing object rather than a
  vacuous one.

### Stage 238 iter-2 (notable — the one iteration-2 this batch)

On iter-1, F1 (support-blindness) and F4 (its `.wl`) were correctly **BLOCKED**:
the notes define `M_tr = M_mix[1+ζ(1−ε)/(1−ζε)]` (§4) but give NO pre-reduction
observable form where `M_tr` cancels, and inventing one was forbidden — Codex
refused to fabricate. The **ORCHESTRATOR REFRAMED F1/F4** (Claude math-coverage
resolution, NON-conceptual, no paper edit) to a faithful non-vacuous form:

1. a negative control `∂_ζ M_tr ≠ 0`;
2. a leak detector `∂_ζ(Rtr·M_tr/M_mix) ≠ 0`;
3. exclusion of `{ζ, M_mix, M_tr}` from the reduced observables.

Codex applied on iter-2; the independent verify confirmed non-vacuity. (Cross-
check: 237-F2 had been fixed the same way — live-channel negative control +
`Not[FreeQ[#,Derivative]]` guard.)

## Notes-only paper_misalignment resolutions (3, user-resolved)

**3 notes-only `paper_misalignment` items, user-resolved (direction: correct the
notes to the script/canonical value; published paper cards + appendices
UNAFFECTED; each cross-engine-corroborated by the new `.wl`):**

- **231** — notes line 98 `dF/dξ` numerator coeffs `240·δ²ξ → 189·δ²ξ` and
  `189·ξ³ → 121·ξ³` (SymPy `Factor` + new `.wl` M1 both give 189/121).
- **232** — notes lines 153/157 figure-of-merit prefactor `168 → 100` (only 100
  reproduces the notes' own quoted decimals `Ξ_χ≈5.5548e5` / `Ξ_J≈1.2664e5`; new
  `.wl` M4 uses c=100). This mirrors the recurring stale "168" typo previously
  corrected at stage 148 (`168π² → 100π²`, IV.5).
- **241** — notes line 577 `ϱ_WΛ` upper bound `193/369 → 125/369` (matches the
  notes' own printed decimal 0.338753; new `.wl` M7 computes
  `ϱ_WΛ|_{β=2/11} = 125/369`).

All published cards/appendices were UNAFFECTED (they carry abstract forms); the
corrections were Codex-applied + Claude-reviewed.

## Mathematica mirror policy — full dual-engine, zero sanctioned mirrors

VII.2 ran with the dual-engine rule in force. **0 sanctioned mirrors were
accepted.** All 10 new `.wl` (231–238, 240, 241) are GENUINE INDEPENDENT routes;
the two checkpoint `.wl` (239, 242), both caught as transliterations, were
**RE-AUTHORED** to independent routes. All 12 are recorded in the Independent-
Mirror Set in `MATHEMATICA_MIRROR_POLICY.md`.

## Per-stage findings tally

| Stage | Status | Findings | Notes |
|-------|--------|----------|-------|
| 231 | dirty | 1 + `.wl` | New independent-route `.wl`. **paper_misalignment (notes-only):** line-98 `dF/dξ` numerator coeffs `240·δ²ξ→189·δ²ξ`, `189·ξ³→121·ξ³`; SymPy `Factor` + new `.wl` M1 both give 189/121 (cross-engine corroborated). Mathematica 39 PASS |
| 232 | dirty | 1 + `.wl` | New independent-route `.wl`. **paper_misalignment (notes-only):** lines 153/157 figure-of-merit prefactor `168→100` (only 100 reproduces the notes' own decimals Ξ_χ≈5.5548e5 / Ξ_J≈1.2664e5; M4 uses c=100); mirrors the recurring stale "168" typo (cf. 148 IV.5). Residual notes-TITLE drift "Stage 249" DEFERRED. Mathematica 38 PASS |
| 233 | dirty | + `.wl` | New independent-route `.wl`. In-file stale-label cosmetics fixed (`.wl` comments Stage 188/223/224→239/240/241, single-file). Mathematica 14 PASS |
| 234 | dirty | + `.wl` | New independent-route `.wl`. Residual notes-TITLE drift "Stage 251" DEFERRED. Mathematica 25 PASS |
| 235 | dirty | + `.wl` | New independent-route `.wl`. Residual notes-TITLE drift "Stage 251/252/253" DEFERRED. Mathematica 31 PASS |
| 236 | dirty | + `.wl` | New independent-route `.wl`. Residual notes-TITLE drift "Stage 253" DEFERRED. Mathematica 37 PASS |
| 237 | dirty | 1 + `.wl` | New independent-route `.wl`. **F2 (variable-independence self-test trap):** support-blindness "verified" by differentiating w.r.t. an absent variable (vacuous) → fixed with a live-channel negative control + `Not[FreeQ[#,Derivative]]` guard. Mathematica 24 PASS |
| 238 | dirty | 1 + `.wl` | New independent-route `.wl`. **ITER-2.** F1 (support-blindness) + F4 (its `.wl`) correctly BLOCKED on iter-1 (no pre-reduction observable form where `M_tr` cancels; inventing one forbidden, Codex refused to fabricate). ORCHESTRATOR REFRAMED F1/F4 (Claude math-coverage resolution, non-conceptual, no paper edit) to: negative control `∂_ζ M_tr≠0`; leak detector `∂_ζ(Rtr·M_tr/M_mix)≠0`; exclusion of {ζ,M_mix,M_tr} from the reduced observables. Codex applied iter-2; independent verify confirmed non-vacuity. Mathematica 36 PASS |
| 239 | dirty (ckpt) | `.wl` re-authored | **Checkpoint.** Transliteration `.wl` RE-AUTHORED to an independent route: forward Jacobian of the boxed dependent vector + native `PseudoInverse` left-inverse + `Reduce`/`Equivalent` orbit-lock (vs the `.py`'s backward hardcoded `SrmDep`). In-file stale-label cosmetics: `.wl` banner `STAGE 222→239` + `Stage221→Stage238` suffixes (single-file). Checkpoint bar MET. Mathematica 42 PASS (SymPy 52 [ok]) |
| 240 | dirty | 1 + `.wl` | New independent-route `.wl`. **F1 (variable-independence self-test trap):** weights now extracted from `Y_support` (which carries `Omega_Q` via its pole), exercising the independence claim against a variable-bearing object. Non-checkpoint constants now independently corroborated: `ρ_α=4/3`, `ζ_req=1/3`, `Pi_tr=(4/3)C_mix`. Mathematica 18 PASS |
| 241 | dirty | 1 + `.wl` | New independent-route `.wl`. **paper_misalignment (notes-only):** line-577 `ϱ_WΛ` upper bound `193/369→125/369` (matches the notes' own printed decimal 0.338753; M7 computes `ϱ_WΛ|_{β=2/11}=125/369`). `ϱ` windows now corroborated: `1/3, 125/369, 2/3, 250/441`. Mathematica 38 PASS (SymPy 32) |
| 242 | dirty (ckpt) | `.wl` re-authored | **Checkpoint.** Transliteration `.wl` RE-AUTHORED to an independent route: `Resolve[ForAll,Reals]` strict-inequality certificate + `D[]` on real closed forms (vs the abstract-ζ device) + `logDrift` total-log-differential (vs `Exp[t·d]`). The load-bearing twin-window strict inclusion `C_mix < Pi_tr < 2 C_mix` is now tested STRICTLY on both engines (ρ_α-style ratio 4/3 region). In-file stale-label cosmetic: `.wl` banner `STAGE 225→242` (single-file). Checkpoint bar MET. Mathematica 24 PASS (SymPy 36) |

**Totals:** 18 findings closed (the missing-`.wl` dual-engine gap on the 10
non-checkpoint stages + the two checkpoint transliterations re-authored; the
variable-independence self-test-trap fixes on 237/238/240; and the 3 notes-only
paper_misalignment typos on 231/232/241), 0 blocked, 0 status-only. **10 new
independent-route `.wl` (231–238, 240, 241)** + **2 re-authored checkpoint `.wl`
(239, 242)**.

## Paper / notes edits (Codex-applied, Claude-reviewed)

This batch APPLIED notes edits (per the file-ownership contract: Codex owns
`paper/*.tex` + `notes/stages/*.md` edits, Claude reviews). **ALL
paper_misalignment items were NOTES-ONLY; the PUBLISHED paper cards/appendices
were UNAFFECTED — they carry abstract forms.** Each correction was cross-engine
corroborated by the new `.wl` independently computing the corrected value. All
were orchestrator-reviewed: correct, isolated, no collateral.

**Numerical/coefficient typos:**

- **231:** notes line 98 `dF/dξ` numerator `240·δ²ξ → 189·δ²ξ` and `189·ξ³ →
  121·ξ³`.
- **232:** notes lines 153/157 figure-of-merit prefactor `168 → 100`.
- **241:** notes line 577 `ϱ_WΛ` upper bound `193/369 → 125/369`.

**In-file stale-label fixes (single-file, NOT a batch renumber, rode the fix
loop):**

- **233** `.wl` comments (Stage 188/223/224 → 239/240/241).
- **239** `.wl` banner `STAGE 222 → 239` + `Stage221 → Stage238` suffixes.
- **242** `.wl` banner `STAGE 225 → 242`.

## Out-of-scope residuals — LOGGED, not fixed this batch

Logged to `PAPER_CLEANUP_TRACKER.md` P4-54 (a continuation of the existing
P4-53 numbering-drift root-cause entry) so they are not lost. These are
notes-title cosmetics with **ZERO effect on verification**:

- **VII.2 notes-title renumber drift** — notes self-titles seen during audit:
  232 "Stage 249", 234 "Stage 251", 235 "Stage 251/252/253", 236 "Stage 253".
  These are the known EM-extension incomplete-renumber; canonical (paper-card /
  script / MANIFEST) = 232/234/235/236 is ground truth. **DEFER to the post-253
  stem-keyed reconciliation pass** per the project numbering-drift policy; NEVER
  offset-sweep (offsets are inconsistent across the realignment).

## Iteration-2 reworks

One: **stage 238** (the F1/F4 support-blindness reframe — see the iter-2
narrative above). The other 11 codex-invoke runs exited 0 on iteration 1.

## Consult

None escalated to the user. The 238 F1/F4 reframe was a Claude math-coverage
resolution (NON-conceptual, no paper edit — a faithful non-vacuous restatement
of support-blindness via negative-control + leak-detector + structural
exclusion). The 3 notes-only paper_misalignment numerical typos (231/232/241)
were USER-RESOLVED 2026-06-03 — corrected to the already-correct SymPy scripts,
each cross-engine corroborated by the new `.wl`.

## Verification

All 12 verification files under `redteam/verifications/stage_231.md` …
`stage_242.md`. Final verdicts:
- `verified` (12): 231–242.
- `needs_rework` → reworked → re-`verified`: 238 (iter-2).
- `blocked_unfixable` (0).

Per-stage Mathematica PASS/[ok] (all exit 0, 0 FAIL): 231=39, 232=38, 233=14,
234=25, 235=31, 236=37, 237=24, 238=36, 239=42, 240=18, 241=38, 242=24. SymPy
all exit 0 (e.g. 239=52 [ok], 241=32, 242=36).

Material change: **0** (`material_change: false` on all 12 — new/re-authored
second engine, de-tautologized / de-vacuized checks, and notes-only
numerical-typo corrections that align prose to the already-correct script; no
SCRIPT-side derived value, constant, identity target, or PUBLISHED paper number
moved).

## Cumulative

Range 001-242 paper-aligned at v2 depth. **242/253 stages red-team verified
(95.7%)** (was 230 after VII.1; VII.2 adds 12 across 231–242). **Third batch of
Part VII**; zero stop-cold, zero material change, full dual-engine coverage with
zero sanctioned mirrors, and BOTH checkpoints (239, 242) cleared the higher bar
(each via a transliteration `.wl` re-authored to an independent route).

Next batch (sequential-audit-chunks rule, awaits explicit user authorization):
**VII.3 onward = stages 243–253** (all currently `pending`). The planned full
end-to-end **second pass** remains a later cross-check, only after the first pass
reaches stage 253.
