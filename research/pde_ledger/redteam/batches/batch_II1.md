---
batch_id: II.1
label: Part II.1 — Overlap isotropy through continuum kernel
started: 2026-05-21
completed: 2026-05-22
stages_total: 13
stages_verified: 13
stages_blocked: 0
material_change_any: false
---

# Batch II.1 — Final summary

All 13 stages reached `verified` status. Both engines (SymPy + Mathematica) exit 0 on every unit; all 43 findings closed substantively. No `material_change` flags, so no downstream cascade — batch III.1 may proceed without re-auditing II.1.

## Per-stage outcomes

| Stage | Audit findings | Codex iterations | Verifier verdict | Notes |
|---|---|---|---|---|
| 024 | 1 (math_transliteration) | 1 | verified | checkpoint — Wick-pair `pairings` recursion replaced with direct sphere integration |
| 025 | 8 (taut×3, hardcoded, insuff_verif×3, math_transliteration) | 1 | verified | `54/5` target reanchored to Stage 023 derivation; non-tautological numerical sample-point checks added; `.wl` switched to `Factor`/`Apart`/`Limit`/`Reduce` |
| 026 | 4 (taut×3, math_transliteration) | 1 | verified | two algebraically-distinct overlap-law routes added; `_expected` self-substitution rebuilds removed |
| 027 | 1 (taut) | 1 | verified | hard-coded `kGeo` answer replaced with explicit `Integrate[chi*gEta, {s,0,l}]`; output canonical form `Cos[2*theta]` proves the integral was actually evaluated |
| 028 | 4 (taut×3, math_transliteration) | 1 | verified | `Eigenvalues[kEff]` sum/product cross-checks added; `Solve[detEff == 0, alpha]` |
| 029 | 3 (taut×2, math_transliteration) | 1 | verified | sequential elimination + `Eigensystem[keffAl]`; codex deviated from directive's `sigmaW` formula in a correctness-preserving way (verifier confirmed) |
| 030 | 3 (taut×2, math_transliteration) | 1 | verified | `Eigenvalues[mMat]` on explicit 2x2 wall block |
| 031 | 4 (insuff_verif, taut, hardcoded, math_transliteration) | 1 | verified | abstract `sp.Function` replaced with physical `s_-`/`lam_-`; `radcrit` derived from `T0^2*R^2`; codex used `Pow.replace` denester instead of `sp.radsimp` (acknowledged deviation) |
| 032 | 3 (taut, math_transliteration, insuff_verif) | 1 | verified | `LinearSolve` route + `Inverse` cross-check + `delta_kappa^2 + 4*Kprod = sigma^2` identity covering interior |
| 033 | 3 (taut×2, math_transliteration) | 1 | verified | gate denominator derived via `cancel(together(...))` + ratio guard; numerical cross-check at two rational rule sets |
| 034 | 1 (taut) | 1 | verified | linear-solve self-check replaced with lambda-form vs x-form cross-check |
| 035 | 2 (taut, math_transliteration) | 1 | verified | `expectZero` LHSs switched from `fTarget`/`alphaReqTarget` to engine-derived `f`/`alphaReq` so wrong coefficients would surface |
| 036 | 6 (taut×3, math_transliteration, insuff_verif) | 1 | verified | checkpoint — symbolic kappa-based `F`-`R_target` identity; `dGTarget`/`gMaxTarget`/`gSeriesTarget` derived natively; discriminant check `disc + 72*delta^2 == 0` |

## Findings breakdown

- `tautological_check`: ~22
- `mathematica_transliteration`: 13 (every single stage — biggest theme)
- `insufficient_verification`: ~4
- `hardcoded_result`: 2 (stage 025 `54/5`, stage 031 `radcrit`)

## Process observations

1. **Transliteration is the default state, not an exception.** Every stage in this batch carried a SymPy-translated Mathematica mirror. The MATHEMATICA_MIRROR_POLICY prose was updated mid-batch to make transliteration screening a named first-pass audit step rather than an opportunistic finding. Expectation for future batches: most stages will need `.wl` rewrites unless the file was authored as a native Mathematica derivation.

2. **Zero codex iter-2 needed.** Compare I.1 (1 iter-2) and I.2 (1 iter-2). Two contributing factors:
   - The `fix_loop.sh` sanity-exec patch (commit `3534b80`) landed before this batch; it now actually runs and writes refreshed transcripts to canonical paths, then greps for `FAIL`/`Traceback`/`AssertionError`/`$Failed`. Previously the sanity exec was silently no-op'ing against the multi-line YAML from `$RT paths`.
   - The `codex.md` "Common Mathematica pitfalls" section (also commit `3534b80`) documents the two engine-level quirks codex hit in I.2: the multi-line `lRed = ...` continuation defect and the `D[expr, f[t]] = 0` quirk. Codex avoided both in this batch.

3. **Auditor self-test step continues to pay off.** Zero directive-level math errors observed in this batch (compare I.1: 3 of 12, I.2: 0 of 11). One auditor caught its own discriminant assertion bug mid-write (stage 036) and corrected it before finalizing the directive.

4. **No engine-level surprises this batch.** No `lRed` continuation recurrence (although stage 029 had a comparable Schur-derivation complexity), no `D[..., qFun[t]] = 0` patterns, no race conditions.

5. **One controlled codex deviation per batch.** This batch's deviation: stage 031 used a `Pow.replace` denester instead of `sp.radsimp` (codex explicitly acknowledged in `## Applied: F2` that `sp.radsimp` "did not reduce in this environment"). Verifier confirmed the matcher only fires when the radical's base equals `root_crit**2`, so a sign error would still surface. Acceptable deviation.

## Prompt hardening landed this batch

None. The toolchain improvements were landed between I.2 and II.1 (commit `3534b80`), not mid-batch.

## Tracker updates landed in this commit

- `notes/MATHEMATICA_MIRROR_POLICY.md`: 13 new Independent-Mirror Set entries (024-036); policy prose updated to require transliteration screening as a named first-pass step (13/13 hit rate in II.1 is signal worth surfacing).
- `notes/CHECKPOINT_TRUST_AUDIT.md`: 024 and 036 rows extended with batch II.1 reverification notes; snapshot bumped to 2026-05-22.
- `notes/CHECKPOINT_CONSTANT_PROVENANCE.md`: no edits (no new constants surfaced at checkpoint stages; 025/031 hardcoded-result findings are at non-checkpoint stages and out of scope for that tracker).
- `notes/PAPER_CLEANUP_TRACKER.md`: new P4-33 (II.1 batch), P3-05 (paper-side propagation), two change-log entries dated 2026-05-22 (toolchain hardening + II.1 completion).
- `notes/EM_PROJECTED_INTEGRATION_TRACKER.md`: `Updated` bumped to 2026-05-22; new completed-checks bullet for batch II.1.
