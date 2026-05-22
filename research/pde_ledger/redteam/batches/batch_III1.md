---
batch_id: III.1
label: Part III.1 — Continuum kernel, generalized branch, rank-2
started: 2026-05-22
completed: 2026-05-22
stages_total: 12
stages_verified: 12
stages_blocked: 0
material_change_any: false
---

# Batch III.1 — Final summary

All 12 stages reached `verified` status. Both engines (SymPy + Mathematica) exit 0 on every unit; all 27 findings closed substantively across the 10 dirty stages; the remaining 2 (042 and 048) passed the auditor's transliteration / tautology / hardcoded probes on first read. No `material_change` flags, so no downstream cascade — batch III.2 may proceed without re-auditing III.1.

## Per-stage outcomes

| Stage | Audit findings | Codex iterations | Verifier verdict | Notes |
|---|---|---|---|---|
| 037 | 1 (math_transliteration) | 1 | verified | hand-supplied `xiTerm`/`alphaTerm`/`sigmaExpected` and `aExpected`/`deltaExpected` removed; `xi`/`alpha` now recovered from two `sigmaWall` entries with the third cross-checked; `A`/`delta` via `Together` numerator-denominator extraction |
| 038 | 2 (math_transliteration, insuff_verif) | 1 | verified | dropped the pre-baked `(cEtaU*cUW + cEtaW*kU)^2` substitution rule in `applyDimless`; added nine non-tautological sign assertions in both engines (multiply derivative by manifestly-positive template under transfer-branch assumption) |
| 039 | 3 (taut×2, math_transliteration) | 1 | verified | `deltaSplit`/`epsWSplit`/`dDir` now derived in `.wl` (SymPy postulate moved to RHS of `derived matches postulated` check); `z1/z0 - (kappa1/kappa0)*R_U` tautology replaced with explicit kappa-rho residual; flat-U baseline substitution checks added |
| 040 | 4 (insuff_verif, taut, math_transliteration, hardcoded) | 1 | verified | added genuine perturbed-matrix eigenvector residual against `M - alpha z z^T` (both rows = 0); two-path cross-check for `H_F` via `F_U` vs `F_general`; `alphaReq` via `Solve[Det[...] == 0, alpha]` and eigenvector via `NullSpace`; Stage 18/19 closed-form provenance comments added |
| 041 | 1 (tautological_check) | 1 | verified | source-tied `n_src` made non-tautological by deriving from general `n_expected` via `q -> t R_U, r -> t, t^2 -> lambda0` substitution in both engines |
| 042 | 0 (clean) | 0 | verified | rank-2 selected-mode mirror is structurally parallel to SymPy but cross-checks via independent canonical paths through `FullSimplify`/`Together`; not flagged because the claims are pure closed-form identities |
| 043 | 3 (math_transliteration, taut×2) | 1 | verified | five independent algebraic paths in `.wl` (`Det` for `dPhi`/`dPhiZ`, residue-ratio for `rPhi`, endpoint limits for `v.D_U.v`, `Series` for mismatch); `A_phi^eff` and `M_supp` self-comparisons replaced with minimal-overlap and split-vs-minimal ratio anchors plus mu-independence derivatives |
| 044 | 5 (math_transliteration, taut×2, insuff_verif, sym_assumption) | 1 | verified | independent `Solve` route for `xiPhys`; tautological renames replaced with `branch_eq` coefficient extraction; algebraically redundant tracking-total-loading assertion deleted; `Rphi=2` literal slice added; unused `sigma0` declarations and incorrect docstring corrected |
| 045 | 3 (taut×2, math_transliteration) | 1 | verified | polynomial-extraction route (`coupling_density` -> `coeff` -> `g_X_ext`) with four firewall assertions; enumerated `channels` list gives `M_tr_channel_sum`; `mTrReq` self-comparison replaced with `Solve[collapsedNum == 0, mTrSym]`; branch numerator via `Series[..., {rPhi, rU, 0}] // Normal` |
| 046 | 2 (math_transliteration, insuff_verif) | 1 | verified | hand-typed `pR`/`p1`/`p2`/`*Expected` literals removed; `.wl` uses `Together[D[...]]`, `Reduce[ForAll[...]]` sign claims, `PolynomialQuotientRemainder`; both engines gained boundary and three-point sign-sample assertions on `G_tr - G_flat` / `F_flat - F_tr` |
| 047 | 2 (taut, math_transliteration) | 1 | verified | `rho_0 - chi_0` and `sigma_0 - chi_0` tautologies (cancelling `lamW/lamW`, `lamphi/lamphi`) closed; `mSupp`/`sEnhance` rewritten via independent algebraic routes; cross-engine `PASS: S from ratio agrees with closed-form S` identity added |
| 048 | 0 (clean) | 0 | verified | support-compensation theorem mirror independently `Solve`s for `zeta_req` and adds two limit-coefficient checks (softening, pole) absent from SymPy; not a transliteration |

## Findings breakdown

- `mathematica_transliteration`: 10 (every dirty stage)
- `tautological_check`: ~11
- `insufficient_verification`: 4
- `hardcoded_result`: 1 (stage 040, Stage 18/19 closed-form readbacks; resolved by provenance comments, not numerical change)
- `symbol_assumption_error`: 1 (stage 044, unused `sigma0` + docstring mismatch)

## Process observations

1. **Transliteration default persists.** 10 of 12 stages had `.wl` files that were SymPy line-by-line ports, matching II.1's 13/13 rate. Stages 042 and 048 cleared the transliteration screen because their pre-existing structure already broke line-by-line correspondence: 042's claims are pure closed-form identities that `FullSimplify` canonicalizes through Mathematica's own engine path; 048's mirror does independent `Solve[zeta_req]` and adds two limit checks SymPy doesn't have. So the dominant pattern from II.1 holds — Mathematica mirrors authored via batch translation rather than independent re-derivation default to transliteration, and the policy doc's "first-pass screening" framing is correct.

2. **Zero codex iter-2 fixes.** Matches II.1 (vs 1 iter-2 each in I.1 and I.2). The `fix_loop.sh` sanity-exec patch (commit `3534b80`) and the `codex.md` "Common Mathematica pitfalls" section landed before this batch and continue to pay off: no `lRed` continuation defects surfaced, no `D[expr, f[t]] = 0` patterns, no race conditions.

3. **Auditor self-test step continues to pay off.** Zero directive-level math errors observed (compare I.1: 3 of 12, I.2: 0 of 11, II.1: 0 of 13). One stage (046) auditor caught its own would-be tautological-quotient check mid-self-test and revised F2 to use three-point sign sampling instead.

4. **One controlled codex deviation.** Stage 038's `applyDimless` rewrite uses `Factor[PowerExpand[...]]` as a generic canonicalization rather than the specific `(cEtaU*cUW + cEtaW*kU)^2 -> zW kEtaEff kWEff kU^2 (1+rho)^2` rule the directive asked to remove. Verifier confirmed the four closed-form reductions still PASS — Codex didn't bake the answer in via a different rule, the substitution chain genuinely reaches the result on its own. Acceptable deviation. (Stage 043 also rotated the `dPhiExpected` sign relative to the directive's wording, but verifier confirmed the chosen `Det` row order makes the sign self-consistent.)

5. **Two clean-on-first-read stages.** This is the first batch with stages that didn't need codex at all. They were marked `verified` directly from the audited state (audit report = clean, no directive written, no fix loop needed). Worth noting because it confirms the audit's bar isn't "find a finding no matter what" — when the script structure genuinely breaks the transliteration screen and the assertions exercise the claim, the auditor says so.

6. **Initial orchestration miss caught.** Wave 1 was sized to 10 (037-046) and wave 2 was intended to be (047-048), but 047 was dropped from the launch when wave 2 only fired 048. Caught on resumption; 047 was relaunched solo and completed clean. Process note: when batching parallel audits in waves, double-check the wave-2 dispatch enumerates every stage not in wave 1.

## Prompt hardening landed this batch

None. The toolchain improvements (commit `3534b80`) had landed before this batch; no further changes were necessary mid-batch.

## Tracker updates landed in this commit

- `notes/MATHEMATICA_MIRROR_POLICY.md`: 12 new Independent-Mirror Set entries (037-048); policy prose updated to record III.1's 10-of-12 transliteration rate alongside II.1's 13-of-13.
- `notes/CHECKPOINT_TRUST_AUDIT.md`: snapshot date bumped to "2026-05-22 (batch III.1 close)"; tier table unchanged because no checkpoints are in the 037-048 range.
- `notes/CHECKPOINT_CONSTANT_PROVENANCE.md`: no edits (no new constants surfaced; stage 040 hardcoded-result finding was resolved by adding provenance comments to existing Stage 18/19 closed-form references, not by changing values).
- `notes/PAPER_CLEANUP_TRACKER.md`: new P4-34 (III.1 batch), P3-06 (paper-side propagation), one change-log entry dated 2026-05-22.
- `notes/EM_PROJECTED_INTEGRATION_TRACKER.md`: `Updated` already at 2026-05-22; new completed-checks bullet for batch III.1 added.
