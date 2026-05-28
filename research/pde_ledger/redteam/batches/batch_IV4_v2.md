---
batch: IV.4
label: Part IV.4 — Penetration, mouth boundary, fixedpoint
range: 127-138
stage_count: 12
verified: 12
findings_count: 22
material_change_count: 0
material_change_stages: []
clean_stages: [128, 129]
checkpoints: []
status_only: [128, 132, 136]
no_mathematica_only: []
created_at: 2026-05-27
applied: true
verification_status: complete
---

# Red-team batch summary — IV.4

12 stages audited, 12 verified. Two clean (128, 129 — 128 is status-only, 129 had only a cosmetic banner note not raised as a formal finding). Zero checkpoints in this range. 10 stages produced findings; 22 findings total. **Ninth ↓ tenth consecutive zero-redirection batch** (II.1 → IV.4) — Codex bypassed entirely.

Three status-only stages in this batch (128, 132, 136). Zero material changes — every fix was either a renumbering, a notes re-attribution, a paper-card downgrade, or a script-substance addition; no derived numerical constants moved.

## Cluster resolutions (3 user-gated, all "Recommended")

- **Cluster A — Mass renumbering (mechanical, 19 edits).** Two uniform stage-number shifts surfaced across the batch:
  - 9 `.wl` banners had `STAGE N` where N = actual − 17 (127→110, 129→112, 130→113, 131→114, 133→116, 134→117, 135→118, 137→120, 138→121).
  - 10 notes/ H1 lines (stages 127-136) had `# Moving-Throat PDE — Stage M` where M = actual + 102 (127→229, 128→230, …, 136→238). 137 and 138 H1 were already correct.
  - Resolution applied: 9 `.wl` banner replacements + 4 `.py` banner replacements (5 banner occurrences; stage 133 has two) + 10 notes H1 fixes, via `/tmp/iv4_mass_renumber.py`. All shifts uniform.

- **Cluster B — Status-only carry-forward re-attribution (2 stages).** Stages 132 and 136 notes attributed load-bearing constants to *downstream* stages (132 cited 180-182; 136 cited 184-186), which can't form a valid carry-forward chain. The constants (`Π_*`, `g_Π`, `(M_s, M_q)`, `Σ_m^*`) actually originate within IV.4 itself.
  - 132 notes line 6: "Stages 180–182" → "Stages 129–131" (covers σ_Π at 129, g_Π at 130, Π_* at 131; verifier noted the broader citation is mathematically more accurate than the resolution doc's draft 130-131).
  - 136 notes line 6: "Stages 184–186" → "Stages 133–135" (covers coupled fixed-point law at 133, F1 first explicit branch at 134, Σ_m^* and outlet-consistent closure at 135).

- **Cluster C — Stage 134 paper-card Checks downgrade.** Card's `\stagefield{Checks}` listed two items neither script exercised: outlet consistency (subject of stage 135) and susceptibility closure (runs through 137-138). Edit applied to `paper/stages/stage_134.tex:21-25`: items 1 and 2 downgraded to carry-forward citations of `\ref{stage:135}` and `\ref{stage:137}` respectively; item 3 (numerical fixed point recording) unchanged. No script change for these Checks. This resolves stage 134's F3 paper_misalignment.

## Mechanical edits applied (no user gate)

- **127** F2: Mathematica independent `Integrate[...]` derivations for `g_slab` and `g_exp`; two new `pass[]` lines (slab/exp closed-form matches source integral). SymPy untouched (transliteration fix is one-sided per directive).
- **130** F1+F2: SymPy + Mathematica — 6-point `dg/dPi > 0` monotonicity sweep (covers Π_*) + `gPi == gPi_boxed` closed-form-vs-paper-form assertion. Three new SymPy raises + three new `expect*` lines per engine.
- **131** F1+F2+F4: SymPy — 4 anchored asserts (Π_*, slope at Π_*, threshold identity, lower-branch discrimination) + closed-form g_minus derivation. Mathematica — banner fix, tautological `expectApprox` replaced with same 4 anchored checks (rewritten with ASCII labels because Mathematica's parser chokes on the `_*)` substring in `Pi_*)` near a comment terminator), closed-form `g_-^F1` derivation with literal-vs-closed-form PASS check. Anchor 3's `Simplify === 0` test required `Chop[..., 10^-30]` because FindRoot's numeric `piStar` produces precision-79 near-zero residue rather than exact symbolic zero (new pitfall #13 candidate).
- **133** F1: Mathematica hand-ansatz block replaced with `DSolveValue[{ODE, BCs}, uFun[x], x]`; `cCoeff` and `aCoeff` removed. Four downstream `expectZero` assertions (ODE residual, u(0), u'(1), mouth derivative kernel) now serve as independent cross-checks against the `DSolve` answer.
- **134** F1+F2+F4: SymPy + Mathematica — replaced with substantive multi-anchor checks against external high-precision targets. The directive's target literals were WRONG (auditor fabricated them); orchestrator recomputed independently via mpmath at 50 digits and used the verified values `S_q(1/2)=0.608336415687717…`, `S_q(1)=0.633127670034487…`, `S_q(2)=0.681366857005321…`. The `static shell channel` and `Pi_*` checks retained. Cluster C handled F3.
- **135** F1+F2: SymPy combined edit — five new substantive asserts (outlet substitution identity, `0 < S_q(Π_*) < 1`, three numerical anchors `Σ_m^*`/`M_s^*`/`M_q^*`, mixed-lane correction `M_q^* * S_q(Π_*) ≈ -0.297111597463199`). The tautological conditional raise removed; residual print retained as transcript probe. Mathematica untouched (already had all six checks).
- **137** F1+F2: SymPy + Mathematica — three new anchored asserts in both engines: paper-card closed-form `M_s`/`M_q` (independent route via `gs^2/(Ks*Theta)` and `(Ks*gq-lam*gs)^2/(Ks*(Ks*Kq+lam^2)*Theta)`), Schur static limit (SymPy uses `sp.limit`, Mathematica uses `Normal[Series[…, {z, 0, 0}]]` — distinct algorithmic routes), outlet consistency at `S_q = 0`. F3 (matrix-Schur derivation) was `blocked_legitimate` per the directive's own caveat — the prescribed block was acknowledged tautological without sign-convention bookkeeping; F1's anchors cover the substantive claim.
- **138** F1: banner mass-fix only (no math content change).

## Notes-side cleanup applied this batch

| File | Lines | Change |
|---|---|---|
| `notes/stages/moving_throat_pde_stage127_penetration_families.md` | 2 | H1 `Stage 229` → `Stage 127` |
| `notes/stages/moving_throat_pde_stage128_mouth_source_bias_status.md` | 2 | H1 `Stage 230` → `Stage 128` |
| `notes/stages/moving_throat_pde_stage129_mouth_boundary_layer.md` | 2 | H1 `Stage 231` → `Stage 129` |
| `notes/stages/moving_throat_pde_stage130_mouth_bias_map.md` | 2 | H1 `Stage 232` → `Stage 130` |
| `notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md` | 2 | H1 `Stage 233` → `Stage 131` |
| `notes/stages/moving_throat_pde_stage132_mouth_boundary_layer_status.md` | 2, 6 | H1 `Stage 234` → `Stage 132`; carry-forward `180–182` → `129–131` |
| `notes/stages/moving_throat_pde_stage133_coupled_mouth_fixedpoint.md` | 2 | H1 `Stage 235` → `Stage 133` |
| `notes/stages/moving_throat_pde_stage134_family1_mouth_fixedpoint.md` | 2, 140 | H1 `Stage 236` → `Stage 134`; gain-line surd `605429` → `605428` (consistency with boxed forms at 86/92 and with scripts) |
| `notes/stages/moving_throat_pde_stage135_outlet_consistent_mouth_closure.md` | 2 | H1 `Stage 237` → `Stage 135` |
| `notes/stages/moving_throat_pde_stage136_coupled_mouth_status.md` | 2, 6 | H1 `Stage 238` → `Stage 136`; carry-forward `184–186` → `133–135` |
| `paper/stages/stage_134.tex` | 21-25 | Checks items 1 and 2 downgraded to `\ref{stage:135}` / `\ref{stage:137}` carry-forward citations |

## Directive-correction catches (orchestrator math-authority)

- **Stage 134 F1/F2 target literals.** The auditor agent fabricated literal targets for `S_q(1/2)`, `S_q(1)`, `S_q(2)` — values disagreed with my independent mpmath calculation by ~0.3-0.4 (not float noise; these were just wrong). Orchestrator detected via independent re-derivation and substituted the mpmath-verified values before applying the directive. The script now passes against the correct externally-sourced targets.
- **Stage 131 Mathematica comment quirk.** Mathematica's parser fails on the comment substring `g'(Pi_*)` adjacent to `*)` (a `(2): g'(Pi_*) ...` pattern parses as malformed). Orchestrator rewrote the new block with ASCII labels (`piStar`, `slope at piStar`) and verified it parses cleanly. This becomes a pitfall #13 candidate.
- **Stage 131 Anchor 3 numerical residue.** `Simplify[thresholdAtStar - expectedForm] === 0` failed because the FindRoot-derived `piStar` produces precision-79 near-zero residues, not exact zero. Added `Chop[..., 10^-30]` to handle.

## Stage-level table

| Stage | Findings | Verdict | Material | Notes |
|---:|---:|---|:-:|---|
| 127 | 2 | verified | no | mathematica_transliteration + banner |
| 128 | 0 | clean (status-only) | no | only Cluster A H1 |
| 129 | 0 | clean | no | only Cluster A banner + H1 |
| 130 | 2 | verified | no | insufficient_verification × 2 |
| 131 | 4 | verified | no | missing-script + tautological + insufficient + hardcoded |
| 132 | 1 | verified | no | paper_misalignment, Cluster B carry-forward |
| 133 | 1 | verified | no | mathematica_transliteration → DSolveValue |
| 134 | 4 | verified | no | missing + tautological + paper_misalignment (Cluster C) + transliteration |
| 135 | 2 | verified | no | tautological + insufficient (combined SymPy edit) |
| 136 | 1 | verified | no | paper_misalignment, Cluster B carry-forward |
| 137 | 4 | 3 resolved + 1 blocked_legitimate | no | insufficient + transliteration + hardcoded (blocked) + banner |
| 138 | 1 | verified | no | paper_misalignment (banner only) |

## Cumulative state after IV.4

150/253 stages red-team verified (59.3%). Range **001-138 paper-aligned at v2 depth**.
