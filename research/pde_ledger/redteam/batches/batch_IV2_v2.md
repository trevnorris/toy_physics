---
batch_id: IV.2
range: 103-114
stage_count: 12
label: Part IV.2 — Outgoing DtN, deformation, robustness, robin
audit_pass: v2 (first-pass under paper-grounded auditor)
audit_date: 2026-05-27
verify_date: 2026-05-27
verified: 12
clean_on_first_read: 5   # 103, 104, 110, 113, 114
dirty: 7                  # 105, 106, 107, 108, 109, 111, 112
findings_count: 16
findings_resolved: 16
findings_blocked_legitimate: 0
user_redirection_rate: 0
material_change_count: 1  # 108 (Cluster B coverage extension; no derived value changed)
orchestrator_hotfixes: 0
codex_invocations: 0
checkpoints: 1            # 105 (chi_Q fix from outgoing DtN)
---

# Batch IV.2 v2 — Red-Team Audit Close

Snapshot date: `2026-05-27`. Batch closed.

## Why this matters

Batch IV.2 covers the **outgoing-DtN spherical-Hankel fingerprint and the deformation-class robustness ledger**: from the exact `l=2` outgoing Λ-form (104) through the checkpoint chi_Q=1 derivation (105 CHECKPOINT) to the reduced 2.5PN closure (106), the general isotropic-DtN deformation algebra (107), the β-parameterized robustness submanifold (108), the linearized branch-selection law (109), the explicit Robin outlet (110), the standalone mixed side-channel pole (111), the hybrid Robin-mixed compensation law (112), and the concrete two-channel core Schur outlet (114). Stages 103 and 113 are status-only consolidation cards. Checkpoint 105 (chi_Q=1) is the load-bearing pin that Parts IV.3+ consume.

This is the **first** IV.2 audit pass (no prior v1 for this batch — auditor used the v2 paper-grounded prompt directly).

## Coverage

| Stage | SymPy | Mathematica | Checkpoint | Pre-Audit | Post-Audit |
|---|---|---|---|---|---|
| 103 | absent (status-only) | absent (status-only) | — | pending | verified, 0 findings |
| 104 | present | present | — | pending | verified, 0 findings (banner sweep only) |
| 105 | present | present | ✓ | pending | verified, 2 findings (math_translit + paper_misalign) |
| 106 | present | present | — | pending | verified, 4 findings (paper_misalign + tautological + math_translit + insufficient) |
| 107 | present | present | — | pending | verified, 1 finding (insufficient_verification SymPy) |
| 108 | present | present | — | pending | verified, 3 findings (paper_misalign + insufficient parse-bug + stale_output) |
| 109 | present | present | — | pending | verified, 3 findings (tautological + math_translit + paper_misalign) |
| 110 | present | present | — | pending | verified, 0 findings (banner sweep only) |
| 111 | present | present | — | pending | verified, 1 finding (math_translit) |
| 112 | present | present | — | pending | verified, 2 findings (paper_misalign + math_translit) |
| 113 | absent (status-only) | absent (status-only) | — | pending | verified, 0 findings |
| 114 | present | present | — | pending | verified, 0 findings (banner sweep only) |

## Finding breakdown (16 substantive)

- `paper_misalignment`: 5 user-gated, all resolved (105 F2 banner-label; 106 F1 Checks-coverage; 108 F1 β-locus; 109 F3 Checks-coverage; 112 F1 banner-label). Resolved per `redteam/resolutions/batch_IV2_paper_alignment.md` Clusters A (106/109 — upstream/downstream Checks carry-forward docstrings), B (108 — substantive β-parameterized preservation submanifold), C (105/112 + 24-site banner sweep across all 10 scripted stages).
- `mathematica_transliteration`: 5 (105 F1 full `.wl` rewrite, 106 F3 full `.wl` rewrite, 109 F2 linearization rewrite, 111 F1 independent chi_Q^mix re-derivation, 112 F2 independent Stage-92 linearized cross-check).
- `tautological_check`: 2 (106 F2 `K4 - 4 K2^2/K0` → target-literal mutual consistency; 109 F1 SymPy a5 closed-form anchor).
- `insufficient_verification`: 3 (106 F4 Δ_Q first-order sensitivity, 107 F1 SymPy Σ2/Σ4 exact-formula parity, 108 F2 Math parse-bug fix).
- `stale_output` banner: 1 (108 F3).
- Banner sweep (Cluster C, not in finding count): 24 banner-relabel sites across all 10 IV.2 scripted stages (104/105/106/107/108/109/110/111/112/114, both engines).

## Resolution pattern

Per the III.4 stall lesson + III.5/IV.1 confirmation: **Codex bypassed entirely**; orchestrator-direct math-authority workflow held cleanly. The three user-gate clusters (A/B/C) were consolidated into `redteam/resolutions/batch_IV2_paper_alignment.md` with orchestrator-direct recommendations + destination-verification greps. The user approved all three as recommended.

**Eight consecutive zero-redirection batches** (II.1, III.1, III.2, III.3, III.4, III.5, IV.1, IV.2).

## User-gate resolutions

All three approved as orchestrator-recommended:

- **Q1 (Cluster A — 106 F1, 109 F3)** — (a) **script-side docstring carry-forward**. 106 docstring/comment names stage 102 (higher-odd irrelevance) and stage 104 (l=2 DtN fingerprint, chi_Q=1 source) as upstream verifiers of paper Checks (ii)/(iii); 109 docstring/comment names stages 110 (Robin), 111 (mixed-pole), 112 (hybrid compensation, even-coefficient preservation) as downstream verifiers of paper Checks. No new asserts, no paper edits.

- **Q2 (Cluster B — 108 F1)** — (a) **substantive script extension**. Added Class D (general β-parameterized preservation submanifold) to both engines. Builds `Lambda_gen(z) = S*Lambda_out(beta*z) + Sigma0 + Sigma2 z^2 + Sigma4 z^4 + I Sigma5 z^5`, re-solves (Σ₂,Σ₄) under canonical-even matching with solutions now β-dependent, then asserts the preservation locus `Sigma_5 = S(1-β^5)/9 - Sigma_0/27`. Existing β=1 reduction (Class C) is subsumed as sanity. Mirrors IV.1 Cluster B (stage 100 closure derivation) pattern.

- **Q3 (Cluster C — 105/112 banner labels + sweep)** — (a) **script-side 24-site banner sweep across all 10 IV.2 scripted stages; paper card display titles deferred to PAPER_CLEANUP_TRACKER P3-13** for a future paper-side pass. Mirrors III.5 / IV.1 Cluster C pattern.

## Stage 108 material_change note

Stage 108 was flagged `material_change: true` by the verifier because Cluster B added the general β-locus check, widening the verification surface. The β=1 reduction (Class C `Sigma_5 = -Sigma_0/27`) is **unchanged** — it's now framed as a sanity-reduction of the general locus. No downstream-visible numeric or symbolic export changes. The verifier did not propagate `upstream_stale: true` to stages > 108 (per IV.1 stage 100 precedent: verification-surface strengthening without value change is not a downstream-staling event).

## Notable items per stage

- **104** (clean) — Both engines already verified the boxed appendix fingerprint `Y_2^out(z) = 1 + z^2/9 + 4z^4/81 + i z^5/27 + ...` non-tautologically through ω^7. Only Cluster C banner edits applied.

- **105 CHECKPOINT** — Mathematica was a line-by-line port of SymPy with the same intermediate names. Full `.wl` re-author replaced this with a structurally distinct path: unfactored ratio `yQretRatio = (4 - 3 omega^2/polescl^2 - 3 I chiQ sigmaQcan omega^5)/(4 denomRet)`; `Apart` round-trip check confirming it matches the canonical `3/4 + 1/(4 denomRet)` decomposition; `SeriesCoefficient` operator-form coefficient extraction (vs `Coefficient[Normal[Series[...]]]`); `Reduce[..., chiQ, Reals]` for the chi_Q identification; polynomial inversion `Solve[Lambda_def · y_ansatz == -3]` for the deformed branch. Variable names changed throughout.

- **106** — F2 tautological fix: `K4 - 4 K2^2/K0 == 0` (forced by `K_n = K0/(4 OmegaQ^n)` construction) replaced with `K0_target K4_target - 4 K2_target^2 == 0` testing the four hardcoded literals' mutual consistency. F3 mathematica re-author: full retarded one-pole + ω^7 series-expansion path with no `nqGeneral/k0/k2/k4/gamma5` intermediates. F4 sensitivity: added Δ_Q = chi_Q - 1 first-order check, asserting linear slope = -2G/(5c^5) and zeroth coefficient = Gamma5_target. F1 paper alignment via Cluster A docstring.

- **107** — Single SymPy F1: added Σ2 and Σ4 exact-formula expect_zero asserts so SymPy matches Mathematica twin's coverage of paper notes' boxed closed forms `Σ_2 = -(3S β^2 - 3S + Σ_0)/9` and `Σ_4 = -(3S β^4 - 3S + Σ_0)/27`.

- **108** — F1 substantive β-locus (Cluster B). F2 parse-bug fix: `chiArg /. beta -> 1 - 1` had been parsing as `chiArg /. (beta -> 0)` (because Mathematica `Plus` binds tighter than `Rule`), silently passing `0^5 = 0` instead of testing χ_arg at β=1; fix `(chiArg /. beta -> 1) - 1`. F3 banner.

- **109** — F1 tautological: added `expect_zero("a5 preservation closed-form", a5_sol - (-5b/9 - a0/27))` between the existing condition print and the substitution check, anchoring the SymPy side to the paper's closed form (the Mathematica side already had this anchor). F2 mathematica rewrite: numerator/denominator-separate series + `1/denominator` inversion + multiply, solving `chiSeries - 1 == 0` directly for `a5` (no `coeff = (chiSeries - 1)/eps` intermediate). F3 paper alignment via Cluster A docstring.

- **110** (clean) — Both engines verify the 5 boxed Robin identities (c2, c4, c5, χ_Q^R, χ_Q^R linearized). Only Cluster C banner edits.

- **111** — F1 mathematica: added independent chi_Q^mix re-derivation block bypassing the L0/L5 extraction (computing directly from the geometric-series form of the pole alone) + new `expectZero["chi_Q^mix routes agree", chiMix - chiMixAlt]` cross-checking the two routes.

- **112** — F1 banner labels via Cluster C (script docstring/print/banner all updated to "Stage 112"; provenance comments updated to reference current stage 104 instead of stale 087/088). F2 mathematica: added independent Stage-92 linearized cross-check block extracting `(a_0, a_5) = (3 sigma, -sigma gamma)` from solB's deformation by reading off constant-piece deviation and imaginary z^5 deviation, then solving the preservation condition `a_0/3 + 9 a_5 = 0` for `gamma_W = 1/9`. Algebraically distinct route from the chi_Q-based solve.

- **114** (clean) — Both engines verify the boxed two-channel Schur outlet identity. Only Cluster C banner edits.

## Pitfall notes

- **Pitfall #11 re-confirmed (PASS-line discipline + Mathematica precedence)**: 108 F2 caught the buggy Mathematica `chiArg /. beta -> 1 - 1` line that had been silently passing on every prior run. The Mathematica parser binds `+`/`-` tighter than `Rule` (`->`), so `beta -> 1 - 1` is `beta -> (1 - 1)` = `beta -> 0`. Combined with `chiArg = beta^5`, the substituted residual is `0^5 = 0`, passing the buggy assertion regardless of whether `chi_arg(beta=1)` actually equals 1. Augmentation: directives instructing Codex to write `expr /. var -> value <op> rest` in Mathematica must use explicit parentheses around the `Rule`. The auditor's structural read caught it where rc=0 alone never would.

- **No new pitfall candidates** from IV.2 — all orchestrator-direct edits applied first-attempt clean.

## Status update on existing pitfalls

- Pitfall #6 (`Dt[..., Constants]` residuals) — not exercised in IV.2.
- Pitfall #7 (primitive-vs-derived substitution) — not exercised in IV.2.
- Pitfall #8 (heavy BVP `dsolve`) — not exercised in IV.2.
- Pitfall #9 (Mathematica `Integrate[]` constant factoring) — not exercised in IV.2.
- Pitfall #10 (SymPy `nsolve` near singularities) — not exercised in IV.2.
- Pitfall #11 (PASS-line discipline + Mathematica precedence inside `/.`) — **re-confirmed prominently** at 108 F2 with the new precedence-around-Rule augmentation.

## See also

- `redteam/resolutions/batch_IV2_paper_alignment.md` — 3-cluster user-gate consolidation with destination-verification greps.
- `redteam/reports/stage_{103..114}.md` — per-stage auditor reports.
- `redteam/directives/stage_{105..109,111,112}.md` — per-stage Codex directives (10 stages had no findings → no directive).
- `redteam/verifications/stage_{103..114}.md` — per-stage verifier outcomes.
- `notes/STAGE_VERIFICATION_COVERAGE.md` — top-level batch table updated (126/253 verified, range 001-114 paper-aligned at v2 depth).
- `notes/MATHEMATICA_MIRROR_POLICY.md` — IV.2 transliteration summary (5/7 dirty rewrites).
- `notes/CHECKPOINT_TRUST_AUDIT.md` — 105 checkpoint entry (third paper-grounded v2 checkpoint at higher-bar standard).
- `notes/CHECKPOINT_CONSTANT_PROVENANCE.md` — 105 provenance entry (zero unexplained literals).
- `notes/PAPER_CLEANUP_TRACKER.md` — P3-13 deferred paper card display title rewrite for stages 103-114.
- `notes/EM_PROJECTED_INTEGRATION_TRACKER.md` — range updated to 001-114 v2-aligned.
