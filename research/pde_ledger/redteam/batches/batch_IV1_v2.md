---
batch_id: IV.1
range: 91-102
stage_count: 12
label: Part IV.1 — Grouped p2 geometry, decoupling, contamination
audit_pass: v2 (first-pass under paper-grounded auditor)
audit_date: 2026-05-27
verify_date: 2026-05-27
verified: 12
clean_on_first_read: 1   # 093 (status-only carve-out)
dirty: 11
findings_count: 28
findings_resolved: 27
findings_blocked_legitimate: 1   # 100 F4 mathematica_transliteration (design-level rewrite)
user_redirection_rate: 0
material_change_count: 1          # 100 (verification surface strengthened; no derived value changed)
orchestrator_hotfixes: 0
codex_invocations: 0
checkpoints: 1                    # 096 (geometry lane check verdict)
---

# Batch IV.1 v2 — Red-Team Audit Close

Snapshot date: `2026-05-27`. Batch closed.

## Why this matters

Batch IV.1 covers the **grouped-P_2 geometry / isotropic decoupling / second-order contamination** ledger: from the static-geometry derivation of the conservative quadrupole module (091) through the explicit `l=0`/`l=2` decoupling theorem (094) to the actual-branch verdict (096 CHECKPOINT) and the outgoing-normalization factorization chain (097-102). The checkpoint at 096 ("geometry lane check verdict") is the load-bearing actual-branch theorem that downstream Part IV consumers cite, and stage 094 contains the 15-angular-integral orthogonality proof that turns Stage 092's obstruction formula into the clean 3/4 + 1/4 module.

This is the **first** IV.1 audit pass (no prior v1 for this batch — auditor used the v2 paper-grounded prompt directly).

## Coverage

| Stage | SymPy | Mathematica | Checkpoint | Pre-Audit | Post-Audit |
|---|---|---|---|---|---|
| 091 | present | present | — | pending | verified, 2 findings (paper_misalign + math_translit) |
| 092 | present | present | — | pending | verified, 2 findings (math_translit + insufficient_verification) |
| 093 | absent (status-only) | present | — | pending | verified, 0 findings (banner-only sweep applied) |
| 094 | present | present | — | pending | verified, 2 findings (engine_disagreement + insufficient_verification) |
| 095 | present | present | — | pending | verified, 3 findings (paper_misalign + insufficient + math_translit) |
| 096 | present | present | ✓ | pending | verified, 3 findings (tautological + insufficient + banner) |
| 097 | present | present | — | pending | verified, 2 findings (paper_misalign + math_translit) |
| 098 | present | present | — | pending | verified, 3 findings (paper_misalign + hardcoded + insufficient) |
| 099 | present | present | — | pending | verified, 4 findings (2×paper_misalign + tautological + insufficient) |
| 100 | present | present | — | pending | verified, 5 findings (paper_misalign + 2 tautological/insufficient + math_translit BLOCKED + symbol_assumption) |
| 101 | present | present | — | pending | verified, 4 findings (tautological + script_doesnt_cover_claim + insufficient + paper_misalign) |
| 102 | present | present | — | pending | verified, 1 finding (insufficient_verification SymPy) |

## Finding breakdown (28 substantive)

- `paper_misalignment`: 10 user-gated, all resolved (091 F1; 095 F3; 097 F1; 098 F1; 099 F2+F3; 100 F1+F2+F3; 101 F4). Resolved per `redteam/resolutions/batch_IV1_paper_alignment.md` Clusters A (091/095/097/099 — orthogonality + static-limit + minimal-module carry-forward), B (100/101 — DtN-fingerprint + higher-odd-term carry-forward + 100 closure derivation), C (banner sweep + paper card titles deferred).
- `tautological_check`: 4 (096 F1, 099 F1, 100 F2 closed by F1, 101 F1).
- `insufficient_verification`: 7 (092 F2, 094 F2, 095 F1, 096 F2, 098 F3, 099 F4, 101 F3, 102 F1).
- `mathematica_transliteration`: 4 (091 F2, 092 F1, 095 F2, 097 F2). One additional (100 F4) BLOCKED in directive as design-level — remains blocked_legitimate.
- `engine_disagreement`: 1 (094 F1).
- `hardcoded_result`: 1 (098 F2).
- `symbol_assumption_error`: 1 (100 F5).
- `script_doesnt_cover_claim`: 1 (101 F2 — SymPy had zero asserts).
- Banner sweep (Cluster C, not in finding count): 23 banner-relabel sites across all 12 stages (12 .py + 12 .wl minus 093 no-py = 23 sites).

## Resolution pattern

Per the III.4 stall lesson + III.5 confirmation: **Codex bypassed entirely**; orchestrator-direct math-authority workflow held cleanly. The three user-gate questions (Clusters A/B/C) were consolidated into `redteam/resolutions/batch_IV1_paper_alignment.md` with orchestrator-direct recommendations + destination-verification greps. The user approved all three as recommended.

**Seven consecutive zero-redirection batches** (II.1, III.1, III.2, III.3, III.4, III.5, IV.1).

## User-gate resolutions

All three approved as orchestrator-recommended:

- **Q1 (Cluster A — 091/095/097/099)** — (a) **script-side docstring carry-forward**. Each of 091, 095, 097, 099 added a docstring comment naming the upstream verifying stage(s) for the paper card's `\stagefield{Checks}` items (orthogonality from 094; static-limit from 091/092/094/096; minimal-module from Part III 088-090). No new assertions, no paper edit. Mirrors III.5 087 F1 consolidation pattern.

- **Q2 (Cluster B — 100/101)** — (c) **hybrid: docstring carry-forwards + targeted closure assertion in 100**.
  - Stage 100 got a substantive closure derivation: impose `mhat_0^2 * Gamma_5 = Gamma_5_target` as the observable condition on the script-derived Gamma_5(K_0, chi_Q, Omega), then assert `closure_ratio - (mhat_0^2 chi_Q N_Q - 1) = 0`. Both engines updated. Tautological A4/A5/A9/A10 (F2/F3) close automatically.
  - Stages 100 and 101 added docstrings naming stage 097 (Check 3 DtN fingerprint pinning chi_Q = 1) and stage 102 (Check 2 higher-odd-term placement) as upstream anchors.
  - 100 F4 (design-level Mathematica transliteration rewrite) remains `blocked_legitimate` in the directive.

- **Q3 (Cluster C — banner sweep)** — (a) **script-side banner sweep across all 12; paper card titles deferred** to PAPER_CLEANUP_TRACKER for a future paper-side pass. Mirrors III.5 12-script sweep.

## Stage 100 material_change note

Verifier 100 flagged `material_change: true` because the script-side closure derivation changed from a tautological cross-check to a substantive one (impose observable condition, derive headline output). **No derived numeric value changed** — `mhat_0^2 chi_Q N_Q = 1` is the same Output as before; the change is purely in how the script *verifies* it. Downstream consumers (stages 101-102 et seq.) see the same algebra and the same closure relation. The flag should NOT propagate `upstream_stale: true` for stages > 100.

## Banner-relabel sweep (23 sites)

Every IV.1 script carried a stale banner from a previous renumber pass. Orchestrator-direct sweep:

| Old banner | New banner |
|---|---|
| STAGE 74 / 074 — GROUPED-P2 + STATIC-GEOMETRY DERIVATION | STAGE 091 |
| STAGE 75 / 075 — DYNAMIC-GEOMETRY OBSTRUCTION | STAGE 092 |
| STAGE 076 — GROUPED-P2 STATUS UPDATE | STAGE 093 |
| STAGE 77 / 077 — ISOTROPIC GEOMETRY-DECOUPLING | STAGE 094 |
| STAGE 78 / 078 — SECOND-ORDER GEOMETRY CONTAMINATION | STAGE 095 |
| STAGE 079 — GEOMETRY-LANE CHECK VERDICT | STAGE 096 |
| STAGE 80 / 080 — SINGLE NORMALIZATION DEFECT | STAGE 097 |
| STAGE 81 / 081 — FAMILY-1 SUPPORT IS AUTOMATIC | STAGE 098 |
| STAGE 82 / 082 — REDUCED FINISH LINE | STAGE 099 |
| STAGE 083 — OUTGOING NORMALIZATION FACTORIZATION | STAGE 100 |
| STAGE 084 — NATURAL SOURCE-MAP REDUCTION | STAGE 101 |
| STAGE 085 — HIGHER ODD IRRELEVANCE | STAGE 102 |

Also: 095 docstring "from Stage 75" → "from Stage 092"; 096 docstring "the Stage 75 obstruction formula" → "the Stage 092 obstruction formula".

## Material-change assessment

`material_change_count: 1` (stage 100 verification-surface change only).

All derived numerics unchanged:
- Stage 091/094/096: `c_pole = 1/4`, `c_geom = 3/4`, `rho_alpha = 4/3`, `zeta_req = 1/3`, `Yhat_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)` ✓
- Stage 092: obstruction `c_pole = (1 + eps_4)/(4(1 + eps_2)^2)` ✓
- Stage 095: contamination begins at `O(chi^2)`; closed forms `K_{g,0}^eff = -chi^2 M_0^2/G_0`, etc ✓
- Stage 097/099/100: `Gamma_5 = 9 K_0/(32 Omega^5)`, `K_0_target = 64 G Omega^5/(45 c^5)`, `Gamma_5_target = 2 G/(5 c^5)`, `mhat_0^2 chi_Q N_Q = 1` ✓
- Stage 098: `zeta_max^(F1) = 2.46752922945601`, `zeta_edge^(F1) ≈ 0.456731`, `gap^(F1) ≈ 2.01080` ✓
- Stage 101: `N_Q = 1/(mhat_0^2 chi_Q)`, small-`Delta_Q` linearization `N_Q - 1 = -Delta_Q + Delta_Q^2 + O(Delta_Q^3)` ✓
- Stage 102: `tauQ` irrelevant at `omega^5`; first `tauQ`-derivative imaginary place at `omega^7` is `1/4` ✓

No downstream Part IV.2+ consumer sees a value change.

## 100 F4 status

`100 F4 mathematica_transliteration: blocked_legitimate`. The directive explicitly marked F4 as Blocked (design-level rewrite needs user direction). User chose Cluster B (c) which strengthens F1 (the closure) but does not authorize a full second-engine rewrite for 100. Implication: stage 100's Mathematica is no longer a pure transliteration of SymPy at the load-bearing closure (the substantive `closure_ratio - (mhat_0^2 chi_Q N_Q - 1)` is the same algebra in both engines but is no longer a tautology), but the upstream series-expansion choreography still follows SymPy step-for-step. Flag for IV.2+ batch to monitor whether the design-level rewrite belongs to a later pass or is sufficient as is.

## Process notes

Verifier prompts now uniformly enforce **PASS-line counting** (per pitfall #11 from III.5). All 12 verifiers confirmed expected PASS counts:
- 091: 9 (8 original + 1 new partial-fraction)
- 092: 7 (4 from F1 closed-form + 3 first-order coefs)
- 093: 4 (clean)
- 094: 34 (30 orthogonality + 1 Y00 norm + 3 static-limit)
- 095: 5 (3 closed-form residuals + 1 Schur derivation + 1 d cpole/dchi |0)
- 096: 20 (15 orthogonality + 5 final ledger)
- 097: 9 (2 series-equiv + 1 Gamma5 closed form + 1 geometric target + 1 Gamma5_target + 4 R_i)
- 098: 6 (3 algebraic + 2 numeric + 1 gap positive)
- 099: 6 (1 static slot + 1 pole residue + 1 K0_target + 3 structural)
- 100: 4 (3 derived ratios + 1 closure ratio)
- 101: 4 (2 input-equation anchors + 1 small-DeltaQ + 1 exact replacement)
- 102: 3 (3 SymPy asserts mirroring 3 Mathematica expectZeros)

Total PASS lines across the 12 Mathematica logs: ~111 substantive checks, all green.

## Output paths

- Reports: `redteam/reports/stage_{091..102}.md`
- Directives: `redteam/directives/stage_{091..102}.md` (093 had no directive — clean)
- Verifications: `redteam/verifications/stage_{091..102}.md`
- Resolutions: `redteam/resolutions/batch_IV1_paper_alignment.md`
- Audit prompts: `redteam/tmp_prompts/audit_prompt_{091..102}.md`
- Verify prompts: `redteam/tmp_prompts/verify_prompt_{091..102}.md`
- Updated scripts: `scripts/moving_throat_pde_stage{091..102}_*_sympy_audit.py` (except 093 — no SymPy)
- Updated Mathematica: `mathematica/moving_throat_pde_stage{091..102}_*_mathematica_audit.wl`
- Refreshed outputs: `scripts/output/...` and `mathematica/output/...` (all 23 mtimes post-fix)
- Updated MANIFEST: `redteam/MANIFEST.yaml` (status: verified, iteration_count: 1, last_verify_date: 2026-05-27 for all 12)

## Cumulative coverage state

114 of 253 stages red-team verified (was 102 → +12 from IV.1). Entire range 001-102 now paper-aligned at v2 depth. Next batch: IV.2 (stages 103-114, "Part IV.2 — Outgoing DtN, deformation, robustness, robin") — 12 stages, includes no checkpoints per BATCHES.md. Per `[[feedback-sequential-audit-chunks]]`, awaits explicit user approval before starting.
