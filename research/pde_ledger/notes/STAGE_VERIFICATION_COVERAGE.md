# Stage Verification Coverage Baseline

This document is the stage-coverage control sheet for the PDE ledger archive.

Use it together with `STAGE_PROVENANCE_INDEX.md`:

- `STAGE_PROVENANCE_INDEX.md` holds the exact raw note and audit artifact paths.
- `STAGE_VERIFICATION_COVERAGE.md` summarizes the current verification surface,
  exposes the main gaps, and gives us a stable baseline for audit planning.

Snapshot date: `2026-05-27` (batch IV.3 close — first-pass paper-grounded audit under v2 prompt for stages 115-126, no checkpoints; one material_change at stage 118 — script-side λ sign flip after auditor caught internal section IV vs section V disagreement)

## Scope

- Canonical stage range: `001--253`
- Canonical stage source model:
  - Parts `I--VII`: compact stage cards
  - Part `VIII`: full derivation-stage files

## Coverage Totals

| Metric | Count |
|---|---:|
| Total stages in archive | 253 |
| SymPy audits present | 240 |
| Mathematica audits present | 165 |
| Numerical stress artifacts present | 15 |
| Stage-specific review notes present | 178 |
| Stages with no SymPy audit | 13 |
| Stages with no executable audit | 11 |
| Stages with Mathematica but no SymPy | 2 |

Mathematica counts above are presence counts only. They are **not**
independence counts. See `MATHEMATICA_MIRROR_POLICY.md` for the current rule
that distinguishes secondary replay coverage from genuinely independent
Mathematica mirrors.

The SymPy runner summary reports `241` passing files because it includes the
repo-level `moving_throat_pde_master_sympy_audit.py` in addition to the `240`
stage-level SymPy audits counted above.

### Red-team verification status

Distinct from "audit file present" (counted above) is "red-team verified" —
the stage has passed the `audit → directive → codex → verifier` pipeline run
out of `redteam/`, with both engines independently checking load-bearing
claims and a clean-context verifier agent confirming the directive's intent
was honored. See `redteam/BATCHES.md` for the live batch table.

As of `2026-05-27`: **138** of 253 stages red-team verified. With IV.3 closed (first audit pass for these 12 stages under v2), the entire range 001–126 is now paper-aligned at v2 depth. No checkpoints in IV.3; stages 120 and 124 were status-only-clean (consolidation cards). One material_change at stage 118 (λ sign flip — internal script inconsistency between sections IV and V; resolved by aligning section V with both section IV's bilinear derivation and the notes' boxed form). Three stages (121, 122, 123) have no Mathematica mirror by design.

| Batch | Range | Stages | Verified | Date |
|---|---|---:|---:|---|
| I.1 | `001--012` | 12 | 12 | 2026-05-21 (v1) / 2026-05-25 (v2 paper-grounded) |
| I.2 | `013--023` | 11 | 11 | 2026-05-21 (v1) / 2026-05-25 (v2 paper-grounded) |
| II.1 | `024--036` | 13 | 13 | 2026-05-22 (v1) / 2026-05-26 (v2 paper-grounded) |
| III.1 | `037--048` | 12 | 12 | 2026-05-22 (v1) / 2026-05-26 (v2 paper-grounded) |
| III.2 | `049--060` | 12 | 12 | 2026-05-22 (v1) / 2026-05-26 (v2 paper-grounded) |
| III.3 | `061--072` | 12 | 12 | 2026-05-22 (v1) / 2026-05-26 (v2 paper-grounded) |
| III.4 | `073--084` | 12 | 12 | 2026-05-25 (v1) / 2026-05-27 (v2 paper-grounded) |
| III.5 | `085--090` | 6 | 6 | 2026-05-27 (v2 paper-grounded, first-pass) |
| IV.1 | `091--102` | 12 | 12 | 2026-05-27 (v2 paper-grounded — first pass) |
| IV.2 | `103--114` | 12 | 12 | 2026-05-27 (v2 paper-grounded — first pass) |
| IV.3 | `115--126` | 12 | 12 | 2026-05-27 (v2 paper-grounded — first pass) |
| IV.4 onward | `127--253` | 127 | 0 | pending |

Cumulative findings closed: ~398 (~219 v1 + 10 v2 from I.1 + 10 v2 from I.2 + 18 v2 from II.1 + 13 v2 from III.1 + 16 v2 from III.2 + 13 v2 from III.3 + 14 v2 from III.4 + 15 from III.5 + 27 from IV.1 + 16 from IV.2 + **27 from IV.3**, plus 1 blocked-legitimate from IV.1). Of the 27 IV.3 closes, 7 were `paper_misalignment` (Cluster A 4 notes-side numerical typos at 121/122/123/126 — `168π²→100π²` for three stages and `228→160` for stage 123; Cluster B λ sign flip at 118 — internal script inconsistency; Cluster C integral inequality coverage gap at 125; **Cluster D 117 transliteration + 2 tautological resolved via cite-upstream-and-downgrade**), 3 `mathematica_transliteration` (115, 117 F1, 125 F2 — independent-derivation insertions), 6 `tautological_check` (116 F1, 116 F2, 117 F2, 117 F3, 117 F4, 119 F1), 5 `insufficient_verification` (116 F4, 119 F2, 121 F2, 122 F2, 126 F3), 2 `hardcoded_result` (116 F3, included in F3 count), 1 `script_missing_paper_claim` (121 F3 Ω_W), 1 `paper_missing_script_claim` (126 F2 positivity), 1 `stale_output` banner (126 F4). Plus a 10-site banner-relabel sweep across all IV.3 scripted stages (115, 116, 117, 121, 122, 123, 125, 126 — orchestrator-direct, matching IV.2's pattern).
`tautological_check` dominant overall, `mathematica_transliteration` second.
`hardcoded_result` rose sharply in III.4 to 12 because the Family-1 numerology
cluster 075-084 packs many literal constants; III.5 quieted again (1 hardcoded
in 089 F4). v2 added `paper_misalignment` as the 10th category — **31** items
total across the eight v2 batches (7 in I.1, 3 in I.2, 3 in II.1, 3 in III.1,
2 in III.2, 4 in III.3 — 2 substantive + 2 banner relabels, 7 in III.4 — 4
substantive + 3 audit-flagged banner relabels, plus an 8-stage
orchestrator-direct banner-relabel sweep when the global-renumber leftover
turned out to be pervasive across III.4; **2 in III.5 — both substantive (087 F1
status/checkpoint consolidation, 089 F1 Pe_req=0 chain closure), plus a 12-script
orchestrator-direct banner-relabel sweep**);
zero user redirections in II.1, III.1, III.2, III.3, III.4, III.5, IV.1, IV.2, **IV.3** (9 consecutive
batches — Codex was bypassed in III.5, IV.1, IV.2, and IV.3 per the III.4 availability lesson; orchestrator-direct
math-authority worked cleanly because the audit + grep evidence was conclusive).
v2 surfaces `insufficient_verification` prominently — 8 in II.1, 5 in III.1, 8 in III.2,
4 in III.3, 1 in III.4, 1 in III.5, 7 in IV.1, 2 in IV.2, 5 in IV.3 = **41** cumulative.
Stage 060 (v1 `material_change: true`) returned **clean (0 findings)** under III.2 v2.
**Stage 068 (v1 `material_change: true`) returned clean at v2**.
Stages now carrying `material_change: true`: 001, 004
(I.1 v2); 013, 014, 015, 018 (I.2 v2); 045 (III.1 v2 — structural-only, F_tr
export value unchanged); 060 (v1, clean at v2); 068 (v1, clean at v2). II.1,
III.1, III.2, III.3, III.4, III.5 v2 each added **zero** value-changing material_change. **IV.1 added one structural material_change at stage 100** (closure derivation strengthened from tautological cross-check to substantive `mhat_0^2 Gamma_5 = Gamma_5_target` imposition; no derived value changed; downstream stages > 100 not marked `upstream_stale`). **IV.2 added one structural material_change at stage 108** (Cluster B β-parameterized preservation submanifold added; the β=1 reduction already verified previously is unchanged, only the verification surface widened; downstream stages > 108 not marked `upstream_stale`). **IV.3 added one substantive material_change at stage 118** (λ sign flip from + to − in section V's `lam_uniform = qstar*v0*I_sq` closure; internal script inconsistency where section IV's bilinear derivation and the notes' boxed form both had `−` but section V dropped the minus during integration. Resolved by aligning section V with section IV. Downstream Schur reductions use `K_s K_q + λ²` and `(K_s g_q − λ g_s)²` — both sign-invariant under squaring — so no numerical downstream propagation; upstream_stale NOT flagged. Future-batch caveat: revisit if a downstream stage uses λ in an unsquared cross-term).
III.3 v2 introduced one orchestrator hot-fix on stage 064 Mathematica
(`Integrate[]` with symbolic functions does not factor constants — verify
integrands first; pitfall #9 candidate). III.4 v2 introduced one orchestrator
fix on stage 082 SymPy (`sp.nsolve` is unstable for `y tan y = 37` near
`pi/2` and jumps to far-away roots — use `mpmath.findroot(..., solver="bisect")`
instead; pitfall #10 candidate). **III.5 introduced two orchestrator hot-fixes:
(a) on stage 088 SymPy, `Y_rho.subs(omega**2/Omega_Q**2, u)` failed silently
because `sp.simplify` reshapes the denominator into `(Omega_Q**2 - omega**2)`
form and the combined ratio is no longer a syntactic subexpression — fix:
substitute `omega**2 -> u * Omega_Q**2` then `sp.simplify`. (b) On stage 088
Mathematica, a comment containing the substring `stage085_*)` was prematurely
closed by the embedded `*)`, causing `Syntax::sntx` and silently skipping the
F1 assertion and regime trichotomy while still reaching `Exit[0]` (rc=0
masking a partial run) — fix: reword to avoid `*)` substrings in comment
text. New pitfall #11 candidate.**
Pitfall #8 was promoted from candidate to
documented in `codex.md` "Common cross-engine pitfalls" item #1 before
III.3 launched. Pitfalls #6, #7 remain candidates; #9 (Mathematica
`Integrate[]` constant factoring), #10 (SymPy `nsolve` near
singularities), and #11 (Mathematica `*)` substring inside comment body
closes prematurely; verifier must check that all expected PASS lines appear,
not just `rc=0`) added in III.5. **IV.1 added no new pitfall candidates** —
all orchestrator-direct edits applied first-attempt clean (Cluster A docstring
carry-forwards + Cluster B 100 closure derivation + Cluster C 23-site banner
sweep), with verifier PASS-line counting confirming all expected substantive
checks across the 12 stages. **IV.2 added no new pitfall candidates** — all
orchestrator-direct edits applied first-attempt clean (Cluster A 106/109
docstring carry-forwards + Cluster B 108 β-locus extension + Cluster C 24-site
banner sweep), plus a Mathematica parse-bug correction (`chiArg /. beta -> 1 - 1`
parsed as `beta -> 0` at 108 F2; fix: `(chiArg /. beta -> 1) - 1`) which
re-confirms pitfall #11 PASS-line discipline (the buggy line passed by accident).
**IV.3 added pitfall #12 candidate** (Mathematica `Solve[expr == 0, frakG]`
fails with "frakG is not a valid variable" when `frakG` is bound to its definition
`= gQ*Sqrt[kS]/(gS*Sqrt[kQ])` — always introduce a fresh symbol for Solve's
target variable, then substitute back). Also re-confirmed: Mathematica `Minimize[]`
on `cos` over an interval frequently returns the input unevaluated; use boundary-value
checks instead under monotone-decreasing assumptions. Directive-correction during
edit: (a) stage 116 directive's `gamma0_from_D = -I·coeff(z,5)` produced wrong sign;
orchestrator corrected to `+I·coeff(z,5)` (analytic check confirms `+I` is right);
(b) stage 115 directive's `parentFamilyResidual - balanceEq*(kS*kQ)/(gS^2)`
multiplicative factor was wrong because `balanceEq` carries denominator
`kS*(kQ*kS + lam^2)`; correct factor is `(kS*kQ + lam^2)/(gS²*kQ)`. Both
corrections caught by re-running the scripts; sub-agents would otherwise have
silently produced FAIL residuals. See per-batch summaries in
`redteam/batches/batch_<ID>.md`.

### Linear projected-EM update

Snapshot addendum: `2026-05-11`

Stages `004--021` are now canonical linear stages for the projection-first
Maxwell integration, parent-action packet, and retained reduced one-port normal
form.  They are no longer counted as Stage `004` substages.  The old compact
Stage `004` reduced Maxwell/mixed calculation is retained as Stage `021`.

Stages `004--020` have file-for-file SymPy migrations from the derivation-only
`notes/em_projected` scripts through
`step_18_parent_throat_action_weak_axisym_packet`.  The `step_19_*`
branch-export packet and `step_20+` computational/runtime diagnostics remain
excluded from paper coverage.  File-for-file Mathematica mirrors for
Stages `004--020` have not been independently derived yet; Stage `021` retains
the legacy reduced Mathematica audit.

## Coverage By Part

| Part | Stage Range | Total | SymPy | Mathematica | Numerical | Review |
|---|---|---:|---:|---:|---:|---:|
| I | `001--023` | 23 | 23 | 6 | 2 | 6 |
| II | `024--036` | 13 | 13 | 13 | 0 | 13 |
| III | `037--090` | 54 | 53 | 54 | 2 | 54 |
| IV | `091--163` | 73 | 61 | 59 | 8 | 69 |
| V | `164--200` | 37 | 37 | 25 | 1 | 25 |
| VI | `201--218` | 18 | 18 | 2 | 0 | 2 |
| VII | `219--242` | 24 | 24 | 3 | 0 | 3 |
| VIII | `243--253` | 11 | 11 | 3 | 2 | 3 |

## Coverage Classes

| Coverage class | Count | Stage ranges |
|---|---:|---|
| SymPy + Mathematica | 163 | See `STAGE_PROVENANCE_INDEX.md` for the file-by-file list. |
| SymPy only | 77 | Includes projected-EM Stages `004--020` and later stages whose Mathematica mirrors are not yet present. |
| Mathematica only | 2 | `084`, `093` |
| No executable audit | 11 | `103`, `113`, `120`, `124`, `128`, `132`, `136`, `141`, `145`, `149`, `153` |

## Constant Provenance Rule

Coverage counts are not trust grades.

Likewise, `Mathematica present` is not the same thing as `independent second
CAS derivation`. Repo-wide counts track coverage breadth; independence is now a
separate policy classification in `MATHEMATICA_MIRROR_POLICY.md`.

For this archive, an audit should be treated as insufficient until every
nontrivial constant used in it is classified as one of:

- `derived in audit`
- `carried forward with source anchor`
- `probe-only numeric value labeled`

Any unexplained literal should be treated as a verification defect, not a style
issue.

## Immediate Gaps

### 1. No executable audit

These stages currently have neither SymPy nor Mathematica nor numerical stress:

`103`, `113`, `120`, `124`, `128`, `132`, `136`, `141`, `145`, `149`, `153`

These are the first places where the archive has mathematical content without an
executable backstop.

### 2. Mathematica without SymPy

These stages have a Mathematica artifact but no SymPy mirror:

`084`, `093`

These need reconciliation before we can make strong claims about dual-CAS
coverage.

### 3. SymPy-only frontier

The current SymPy-only region is:

`121--123`, `188--199`, `201--202`, `204--217`, `219--220`, `222--238`, `240--241`, `244--247`, `249--252`

Operationally, the main contiguous Mathematica mirror still ends at Stage `187`,
with isolated later checkpoints now hardened at Stages `200`, `203`, `218`,
`221`, `239`, `242`, `243`, `248`, and `253`.

### 4. Review-coverage gap

Stage-specific review notes are missing for:

`121--124`, `188--199`, `201--202`, `204--217`, `219--220`, `222--238`, `240--241`, `244--247`, `249--252`

That means the late-stage archive is not only missing Mathematica coverage, but
also lacks the earlier review-note pattern that exists for most of Stages
`001--187`.

### 5. Numerical stress remains sparse

Numerical stress coverage exists only for the following stage families:

`003/021`, `045--048`, `058/170`, `106`, `125`, `142--144`,
`146--148`, `150--152`, `154`, `155--156`, `157`, `185--187`,
`248`, `253`

This is a narrow slice relative to the symbolic verification surface and should
be treated as targeted spot-checking, not broad numerical validation.

## How To Use This Baseline

1. Use `STAGE_PROVENANCE_INDEX.md` when you need the exact note or script path
   for a stage.
2. Use this baseline when deciding which verification gaps are structural and
   which are only missing a second engine or review note.
3. Use the part-level counts to prioritize the next audit wave.
4. Use the coverage classes to define the future-paper citation support set.
5. Use `CITATION_SUPPORT_SET.md` when deciding which gaps are most important for
   downstream citation hardening.
6. Use `CHECKPOINT_CONSTANT_PROVENANCE.md` for the growing no-magic-numbers log
   on the checkpoint subset.
7. Use `CHECKPOINT_TRUST_AUDIT.md` for the current stage-level trust baseline
   on that checkpoint subset.

## Recommended Next Verification Sequence

1. Use `CHECKPOINT_TRUST_AUDIT.md` as the current checkpoint trust baseline.
2. Reconcile the remaining repo-wide Mathematica-only outliers:
   `084`, `093`, `103`, `113` (last two are IV.2 status-only).
3. Use `MATHEMATICA_MIRROR_POLICY.md` when deciding whether a mirror gap is an
   execution-coverage gap or an independence gap.
4. Then widen the audit wave to the remaining repo-wide gaps.
5. Backfill executable audits for the remaining no-executable stages:
   `128`, `132`, `136`, `141`, `145`, `149`, `153`. (103, 113, 120, 124
   cleared as status-only consolidation cards under IV.2/IV.3.)
6. Use the new `248` / `253` harnesses as the template if numerical-stress
   coverage is widened beyond the current spot-check set.
