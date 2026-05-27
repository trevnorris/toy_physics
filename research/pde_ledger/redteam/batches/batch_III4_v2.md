---
batch_id: III.4_v2
label: Part III.4 v2 — Family-1 geometry, thresholds, quadrupole (paper-grounded re-audit)
started: 2026-05-27
completed: 2026-05-27
stages_total: 12
stages_verified: 12
stages_blocked: 0
material_change_any: false
audit_version: v2 (paper-grounded)
predecessor: III.4 (v1, completed 2026-05-25)
---

# Batch III.4 v2 — Final summary

All 12 stages 073-084 re-verified at v2 paper-grounded depth. Five stages came back clean (073, 077, 079, 080, 083); seven had findings (074, 075, 076, 078, 081, 082, 084 banner-only). Cumulative v2 findings: 14 audit-flagged + 11-stage orchestrator-direct banner-relabel sweep.

## Per-stage outcomes

| Stage | v2 findings | Codex iterations | Verifier verdict | Notes |
|---|---|---|---|---|
| 073 | 0 | n/a | verified | clean (with orchestrator-direct banner relabel post-verifier) |
| 074 | 1 paper_misalignment (value_mismatch) | n/a (orchestrator-direct) | verified | Q1 (a) — paper:31 `128/sqrt(5)` → `111/sqrt(5)`, notes:117 + notes:075:63 `179/sqrt(5)` → `111/sqrt(5)`; added `alpha_ref - 111/sqrt(5) == 0` assertion both engines |
| 075 | 1 high tautological_check + 1 high paper_misalignment + 1 mathematica_transliteration + 1 new F4 lock | n/a (orchestrator-direct) | verified | F1: replaced tautological round-trip and definitional identity with asymptotic limits `alpha * Delta_inf → 1` (large alpha) and `Delta_0 → 1/2` (small alpha) in both engines. F3: subsumed by F1 (comment cleanup). F2 (Q2 a): paper:7,24 `117 Theta_w` → `100 Theta_w`; notes:108,116 `168 Theta_w` → `100 Theta_w`; notes:124-128 `/168` → `/100`. F4 (new): `alpha_r^2 == 100` script-side lock both engines |
| 076 | 1 tautological_check + 1 hardcoded_result | n/a (orchestrator-direct) | verified | F1: trivial `(2x)^2 = 4 x^2` replaced with closed-form target comparison exercising enthalpy-lock `1/4` factor. F2: TODO(provenance) replaced with direct citation to `notes/.../stage076_n5_wall_depth_lock.md` section 4 |
| 077 | 0 | n/a | verified | clean (with orchestrator-direct banner relabel post-verifier) |
| 078 | 1 mathematica_transliteration + 1 paper_misalignment (notes_contradicts_script, banner) | n/a (orchestrator-direct) | verified | F1: `thetaSuffSym = thetaFailSym * decimal_ratio` replaced with explicit Stage-75 closed form `-(45 cosh(α) + 27√5 sinh(α))/(2500 − 2500 cosh(α))` at `α = 111√5/5`; misleading "verify ratio" comment removed. F2: banner `Stage 61` / `STAGE 061` → `Stage 078` / `STAGE 078` both .py and .wl. One Mathematica syntax fix during apply (line-continuation broke at `)` — wrapped outer parens) |
| 079 | 0 | n/a | verified | clean (with orchestrator-direct banner relabel post-verifier) |
| 080 | 0 | n/a | verified | clean (with orchestrator-direct banner relabel post-verifier) |
| 081 | 1 paper_misalignment (target_mismatch, banner) | n/a (orchestrator-direct) | verified | F1: banner `Stage 64` / `STAGE 64` → `Stage 081` / `STAGE 081` on .py docstring + banner (.wl was already correct) |
| 082 | 2 paper_misalignment (script_missing_paper_claim, ×2) + 1 insufficient_verification | n/a (orchestrator-direct) | verified | F1+F2 (Q3 a): closed-form `zeta_phys = Omega_Pe^2 (kappa + pi^2/4)/(kappa + y^2)` pin block added in both engines; verified `Omega_Pe → pi/2` as `Pe → oo`; computed `y_F1 ≈ 1.52948` via `mpmath.findroot(..., solver='bisect')` (sp.nsolve unstable near pi/2 — pitfall #10 candidate); matched `zeta_phys(Pe→oo, kappa_F1, y_F1) ≈ 2.4675292` against upstream `zeta_max^(F1)` to `1.77e-13`. F3: replaced tautological derivative checks `dR_quad/dzeta_phys + 1 == 0` and `dR_quad/dPi_tr − dzeta_req/dPi_tr == 0` with numerator/denominator factorization `numerator − C_mix(1−eps_blk) == 0` and `denominator − (C_mix − eps_blk(2 C_mix − Pi_tr))^2 == 0` exercising notes §4 strict-positivity content |
| 083 | 0 | n/a | verified | clean (with orchestrator-direct banner relabel post-verifier) |
| 084 | 0 audit findings + 1 informal banner-relabel | n/a (orchestrator-direct) | verified | clean at v2; banner `STAGE 067` → `STAGE 084` applied orchestrator-direct after auditor's side observation |

## Findings breakdown (v2 audit-flagged only)

- `paper_misalignment` (substantive): 4 (074 F1 alpha, 075 F2 Upsilon_w, 082 F1 zeta_phys, 082 F2 Family-1 pair)
- `paper_misalignment` (banner relabel): 3 audit-flagged (078 F2, 081 F1, 084 noted) + 8-stage orchestrator-direct sweep (073, 075, 076, 077, 079, 080, 082, 083) = 11 total banner items
- `tautological_check`: 2 (075 F1, 076 F1)
- `hardcoded_result`: 1 (076 F2)
- `mathematica_transliteration`: 2 (075 F3 subsumed by F1, 078 F1)
- `insufficient_verification`: 1 (082 F3)
- new F4 added by orchestrator post-Q2 resolution: 1 (075 F4 — script-side `alpha_r^2 == 100` lock)

Total v2 audit-flagged: 14 findings closed (4 substantive paper_misalignment + 3 banner paper_misalignment + 7 other script-side). Plus 8-stage orchestrator-direct banner sweep (informal cleanup; no formal audit finding).

Pattern: substantive paper_misalignment dominated (4 of 14 audit-flagged), reflecting the v2 auditor's stronger paper-grounded comparisons catching value mismatches that v1 missed (alpha = 128 vs 111, Upsilon_w = 117 vs 100, both with engine evidence). Banner relabels surfaced ubiquitously from the global stage renumber.

## Substantive paper-misalignment cluster (4 items, all (a) approved)

Three substantive user-gate questions, all resolved as direction (a):

- **Q1 — Stage 074 `alpha = sqrt(kappa)` value:** paper `128/sqrt(5)`, notes `179/sqrt(5)` (twice — in 074 and 075 notes), engines both compute `111/sqrt(5)` (since `111² = 12321` exactly). The numerical value `~49.6407091` cited in both notes files matches `111/sqrt(5)`, not the other two. Direction (a) approved: paper-side and notes-side typo fix to `111`. Script-side assertion `alpha_ref - 111/sqrt(5) == 0` locks the value going forward.

- **Q2 — Stage 075 `Upsilon_w` conversion factor:** paper Inputs `117 Theta_w`, notes §3 `168 Theta_w`, script `100 Theta_w`. The paper's own boxed `Theta_fail = 3.626e-4` and the `Xi_F1 = 1369·Upsilon_w = 136900·Theta_w` carry-forward in notes/stages 082/083/084 are mathematically consistent ONLY with 100. Direction (a) approved: paper Inputs and notes §3 prose updated to `100 Theta_w` (and notes' `/168` arithmetic updated to `/100`). Script-side `alpha_r^2 == 100` assertion lock added.

- **Q3 — Stage 082 F1+F2 extend scripts:** paper-claim closed form `zeta_phys = Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2)` was not exercised in stage 082 scripts (stage 084 already pinned it). Direction (a) approved: extend both stage 082 scripts to instantiate the closed form, verify `Omega_Pe → pi/2` as `Pe → oo`, compute `y_F1` (smallest positive root of `y tan y = 37` in `(0, pi/2)`), and match `zeta_phys(Pe→oo, kappa_F1, y_F1)` against upstream `zeta_max^(F1)`. F2's Upsilon_w cross-reference resolved automatically by Q2.

User-redirection rate this batch: **0** (fifth consecutive batch — Codex stalled mid-consultation but the audit + grep evidence was already conclusive; orchestrator-direct apply substituted and produced the same recommendations the user approved).

## Codex stall and orchestrator-direct apply

The first codex-chat consultation (Q1 on stage 074 alpha) stalled with no session log written to `~/.codex/sessions/2026/05/27/`. Codex processes were running (PIDs 2260200/2260220/2260223/2260241) but produced no output. Orchestrator killed the processes after ~5min.

Pivoted to orchestrator-direct apply for the substantive items because:
1. The audit + grep evidence on all 3 questions was already conclusive (no engineering judgment ambiguity).
2. The math in each directive was fully specified by the auditors.
3. The verifier wave catches any error afterward.

Verifier wave confirmed all 7 dirty stages `verified`. Subsequent stage 076 verifier note prompted the comprehensive III.4 banner-relabel sweep, which was applied orchestrator-direct for consistency with the audit-flagged banner relabels.

**Recommend cross-checking codex availability before III.5 launches** — the math-authority delegation should resume if Codex is back online.

## Orchestrator-applied math fixes

**Stage 082 SymPy `nsolve` instability (new pitfall #10 candidate):**

The closed-form pin block for stage 082 needed `y_F1 = smallest positive root of y tan y = 37 in (0, pi/2)`. Initial attempt:

```python
y_F1 = sp.nsolve(y_sym * sp.tan(y_sym) - 37, y_sym, sp.Float("1.527"), prec=30)
```

returned `29.1773926972360908276236847979` — Newton iteration overshot from initial guess `1.527` because `tan'(y) = sec^2(y)` blows up near `pi/2`, sending the next iterate to a far root.

Replaced with bracketing bisection:

```python
import mpmath
mpmath.mp.dps = 30
y_F1 = sp.Float(
    mpmath.findroot(lambda yv: yv * mpmath.tan(yv) - 37, (1.5, 1.55), solver="bisect"),
    30,
)
```

which returns `1.52948248371469964992710762240` stably. The Mathematica counterpart `FindRoot[y*Tan[y] - 37, {y, 1.527}, WorkingPrecision -> 30]` is stable because Mathematica's hybrid solver detects the near-singularity and falls back.

**Defense:** when finding roots near a derivative singularity, prefer bracketing solvers (`mpmath.findroot(..., solver="bisect"|"anderson")`) over `sp.nsolve` (Newton). Promote to `codex.md` if recurs.

**Stage 078 Mathematica line-continuation fix:**

The directive's F1 prescribed multi-line `thetaSuffSym = -(...) / (...)`. Mathematica's auto-continuation broke at the closing `)` of the first line (a complete expression), and the next line's leading `/` was a parse error. Fixed by wrapping the whole expression in outer parens: `thetaSuffSym = (-(...) / (...));`.

## Banner relabel sweep (orchestrator-direct, 23 edits across 11 stages)

After verifier on 076 flagged stale `STAGE 59/059` as non-blocking observation, investigation revealed pervasive stale banners across all III.4 stages from commit `0d09ef6` (global renumber). Applied uniform relabel to align all III.4 self-banners with post-renumber numbering:

| Stage | Old (.py docstring + banner) | Old (.wl banner) | New |
|---|---|---|---|
| 073 | "Stage 56" / "STAGE 56" | "STAGE 056" | "Stage 073" / "STAGE 073" |
| 074 | "Stage 57" / "STAGE 57" | "STAGE 057" | "Stage 074" / "STAGE 074" |
| 075 | "Stage 58" / "STAGE 58" | "STAGE 058" | "Stage 075" / "STAGE 075" |
| 076 | "Stage 59" / "STAGE 59" | "STAGE 059" | "Stage 076" / "STAGE 076" |
| 077 | "Stage 60" / "STAGE 60" | "STAGE 060" | "Stage 077" / "STAGE 077" |
| 078 | "Stage 61" / "STAGE 61" | "STAGE 061" | "Stage 078" / "STAGE 078" |
| 079 | "Stage 62" / "STAGE 62" | "STAGE 062" | "Stage 079" / "STAGE 079" |
| 080 | "Stage 63" / "STAGE 63" | "STAGE 063" | "Stage 080" / "STAGE 080" |
| 081 | "Stage 64" / "STAGE 64" | (already correct) | "Stage 081" / "STAGE 081" |
| 082 | "Stage 65" / "STAGE 65" | "STAGE 065" | "Stage 082" / "STAGE 082" |
| 083 | "Stage 66" / "STAGE 66" | "STAGE 066" | "Stage 083" / "STAGE 083" |
| 084 | (no .py) | "STAGE 067" | "STAGE 084" |

All 23 edits are pure label changes (no math content affected). All 11 stages re-run with exit 0 and correct new banners visible in transcripts.

## Process observations

1. **Zero `material_change: true` this batch** (matching III.2/III.3 v2; II.1 v2 also zero; only III.1 v2's stage 045 had structural-only material_change with unchanged export value). The paper-side prose updates align text to math the scripts already computed correctly; the stage 082 closed-form pin block produces only verification (matches upstream `zeta_max^(F1)` to 1.77e-13), introducing no new numeric constants.

2. **First batch where Codex stalled.** First five v2 batches saw Codex deliver math-authority recommendations cleanly; III.4's first codex-chat consultation produced no session log despite spawning processes. Substituted with orchestrator-direct apply leveraging the auditors' fully-specified directives plus local grep evidence. Outcome: same recommendations the user would have approved, but bypassing Codex iteration. Net effect on quality: none (verifier wave confirmed all 7 dirty stages). Net effect on robustness: a single point of failure surfaced in the math-authority delegation flow — recommend monitoring Codex availability before each batch's user-gate kicks off.

3. **Banner-relabel sweep was the largest cosmetic delta.** Although the substantive math fixes are 4 paper-side typo/value updates plus the stage 082 closed-form pin, the diff size is dominated by 23 banner relabels across 11 stages — none of which change any math. This is the first batch where the orchestrator broadened scope beyond audit-flagged items to clean up the pervasive renumber-leftover, justified because (a) verifier on 076 prompted investigation, (b) the pattern is identical to audit-flagged 078/081/084 relabels, and (c) leaving III.4 in a half-relabeled state would create confusion. Future batches should consider running a `grep -rE 'STAGE [0-9]+|Stage [0-9]+' scripts/ mathematica/` pre-audit pass to surface the same class of stale text upfront.

4. **Pitfall #10 candidate** (SymPy `nsolve` instability near singularities) joins pitfall #9 (Mathematica `Integrate[]` constant factoring) as the second candidate in the v2 sweep that emerged from substantive math content rather than from style/format. Both promote to `codex.md` if recurs.

5. **`mathematica_transliteration` share dropped sharply.** Only 1 of 7 dirty stages (078) had a transliteration finding at v2, vs 9/12 in III.3 v1 and 7/12 in III.4 v1. The v2 paper-grounded reading naturally surfaces more value-mismatch issues than line-by-line port issues.

## Tracker updates landed in this commit

- `notes/STAGE_VERIFICATION_COVERAGE.md`: snapshot date `2026-05-27`, cumulative findings ~313, III.4 row v2 date added, paper_misalignment cumulative count to 29, user-redirection streak to 5 consecutive, pitfall #10 candidate noted.
- `notes/MATHEMATICA_MIRROR_POLICY.md`: snapshot bumped; III.4 v2 entry added with the new SymPy nsolve instability defense note (pitfall #10 candidate).
- `notes/CHECKPOINT_TRUST_AUDIT.md`: snapshot date bumped (no checkpoint stages in III.4).
- `notes/CHECKPOINT_CONSTANT_PROVENANCE.md`: snapshot date bumped (no new hardcoded constants surfaced; existing literals in 075-084 numerology cluster all have either derivation or paper-anchor provenance).
- `notes/PAPER_CLEANUP_TRACKER.md`: new P4-40 row (III.4 v2 verification status) and P3-12 row (optional paper-stage-card v2 verification stamps); cumulative change-log entry dated 2026-05-27.
- `notes/EM_PROJECTED_INTEGRATION_TRACKER.md`: snapshot date bumped to 2026-05-27 with III.4 v2 close narrative added (III.4 is out of the linear projected-EM core range 004-021).
