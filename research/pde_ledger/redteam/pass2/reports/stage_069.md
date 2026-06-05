---
unit_id: 069
batch: III.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage069_final_reduced_verdict.md]
  paper_appendix: present
---

# Audit unit 069 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_069.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage069_final_reduced_verdict.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows 116, 256, 318 reference stage 069)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.txt`

## What the paper claims

Stage 069 is a consolidation/checkpoint that states the final pre-Family-1 support/source verdict. `\stagefield{Output}` reads: "The final pre-Family--1 support/source verdict \eqref{eq:app-stage069-three-zone}." The boxed three-zone result is: (Zone A) `W_wall < Pe_req/Delta_inf` ⇒ universal fail; (Zone B) `W_wall > Pe_req/Delta_0` ⇒ universal matched success; (Zone C) `Pe_req/Delta_inf ≤ W_wall ≤ Pe_req/Delta_0` ⇒ profile-sensitive band. The notes add the refinement that the explicit sech–Gaussian benchmark only narrows the profile-sensitive region to two thin side-bands of relative width `P_res - 1 = (1 - C_res^2)/C_res^2 ≈ 0.56%`, with `C_res^2 = 0.994418836451529...` and `P_res = 1.005612487760576...` carried in from Stages 067/068. These numeric values are upstream deliverables (Stages 067/068), reported in the notes' preamble and the appendix row for stage 068 (`P_res ≃ 1.005612`); the stage-069 card itself states only the symbolic three-zone structure. No numeric constant is a stage-069 deliverable; the stage's job is to assemble the symbolic verdict and the exact side-band geometry from the upstream forms.

## What the script claims to verify

Both scripts construct the matched thresholds `Wfail_match = Pe_req/Delta_inf`, `Wsuff_match = Pe_req/Delta_0` (with `Delta_inf = Delta_0 + Delta_gap`) and the resonance-shifted thresholds, then assert: (i) the matched window edges come from the generating function `W_match(Delta_eff) = Pe_req/Delta_eff` and `W_match` is monotone decreasing in `Delta_eff`; (ii) `P_res = 1/C_res^2`; (iii) the matched-window width equals `Pe_req(Delta_inf-Delta_0)/(Delta_0 Delta_inf)`; (iv) the resonance thresholds equal `Pe_req/(C_res^2 Delta_inf)` and `Pe_req/(C_res^2 Delta_0)`; (v) the two profile-sensitive side-bands have exact widths `Pe_req(1-C_res^2)/(C_res^2 Delta_inf)` and `.../(C_res^2 Delta_0)`, with relative width `(1-C_res^2)/C_res^2 = P_res - 1`; (vi) ordering/positivity that places an interior point (parameterized by `v_fail`,`v_succ` via `u = v/(1+v)`) strictly inside each side-band. Everything is symbolic over free positive parameters (`Pe_req, Delta_0, Delta_gap, Pres_gap` in SymPy; `Cres2Prim` as the primitive in Mathematica). No numeric constants are pinned.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Zone A/B thresholds `Pe_req/Delta_inf`, `Pe_req/Delta_0` | `matched fail/success edge from W_match(...)` (sympy 99–106; wl 103–110) | match |
| `W_match` monotone (so thresholds order correctly) | `W_match decreasing in Delta_eff` (sympy 107–110; wl 111–114) | match |
| Profile-sensitive band = matched window `[Pe_req/Delta_inf, Pe_req/Delta_0]` | `matched window width` + ordering (sympy 114–119; wl 128–132) | match |
| Resonance refinement `P_res = 1/C_res^2`, thresholds `/(C_res^2 Δ)` | sympy 113, 133–140; wl 86, 88–95 | match |
| Side-band relative width `P_res-1 = (1-C_res^2)/C_res^2 ≈ 0.56%` | side-band width + `P_res-1-(1-C_res^2)/C_res^2` (sympy 151–162; wl 140–151) | match |
| Side-bands strictly inside `[matched, resonance]` edges | `u`-parameterized interior-point positivity (sympy 167–173; wl 155–161) | match |

`paper_alignment: aligned`. Every paper-side symbolic deliverable has a faithful script-side check. The numeric `C_res^2`/`P_res` values are upstream deliverables (Stages 067/068), correctly carried as symbols here; not a stage-069 deliverable, so no MISSING-DELIVERABLE.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 99–102 | `W_match(Delta_inf) - Wfail_match == 0` | Zone A threshold | yes |
| A2 | sympy | 103–106 | `W_match(Delta_0) - Wsuff_match == 0` | Zone B threshold | yes |
| A3 | sympy | 107–110 | `-dW/dDelta_eff·Δ²/Pe_req > 0` | monotone ordering | yes |
| A4 | sympy | 113 | `Pres - 1/Cres2 == 0` | `P_res=1/C_res^2` | partial (Pres,Cres2 both from Pres_gap; identity by construction but exercises algebra) |
| A5 | sympy | 114–117 | `matched window width == 0` | band = matched window | yes |
| A6 | sympy | 124–127 | `Pres_from_ratio - (1+Pres_gap) == 0` | penalty=band ratio | no (tautological in sympy; see Notes) |
| A7 | sympy | 128–132 | fail-ratio == success-ratio | penalty consistency | partial |
| A8 | sympy | 133–140 | resonance thresholds `=Pe_req/(Cres2·Δ)` | resonance thresholds | yes |
| A9 | sympy | 151–162 | side-band widths + `P_res-1=(1-C²)/C²` | side-band geometry | yes |
| A10 | sympy | 167–173 | interior point strictly inside side-bands | band ordering | yes |
| B1 | mathematica | 86–95 | `P_res`, resonance thresholds via `Cres2` primitive | resonance thresholds | yes |
| B2 | mathematica | 48–57, 97 | `Pres-PresGap consistency via Solve` | Solve-route consistency | partial |
| B3 | mathematica | 103–114 | matched edges + monotone | Zone A/B + ordering | yes |
| B4 | mathematica | 117–151 | window width + side-band widths | band/side-band geometry | yes |
| B5 | mathematica | 155–161 | interior-point positivity | band ordering | yes |

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.txt:3`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.txt:3`

**What's wrong:**
Both committed transcripts predate the scripts (output mtime `May 22 20:05`; both scripts mtime `Jun 3 15:59`). The content disagreement is in the banner self-label: the current scripts print `STAGE 069 — FINAL REDUCED SUPPORT/SOURCE VERDICT` (sympy `.py:65`, wl `.wl:33`), but both saved `.txt` files print `STAGE 052 — FINAL REDUCED SUPPORT/SOURCE VERDICT` (line 3 of each). The body of each transcript otherwise matches the current scripts' printed expressions line-for-line (every `expect_zero`/`expectZero` result and printed threshold form agrees), so the staleness is confined to the stage-number self-label — the known pass-2 "committed `.txt` predates the banner fix" pattern.

**Why this matters:**
The transcript is the committed evidence the stage card cites (`\stagefield{Verification}` points at both `.txt` files). A transcript self-labelled "STAGE 052" for a stage-069 checkpoint is a numbering self-label error in the script/output band and contradicts the current script. On a checkpoint (higher bar), the committed evidence must match the script it certifies.

**Required change:**
Re-run both scripts and overwrite the committed transcripts so the banner reads `STAGE 069`. No source edit is required — the scripts already emit the correct banner; only the stale `.txt` files need refreshing.

**Verification:**
After refresh, line 3 of both `.txt` files reads `STAGE 069 — FINAL REDUCED SUPPORT/SOURCE VERDICT`; all `PASS:`/`= 0` lines unchanged; both scripts exit 0.

## Independent-derivation check (Mathematica)

Not a transliteration. The two engines use genuinely different primitives and routes:
- SymPy parameterizes the penalty through `Pres_gap` directly: `Pres = 1 + Pres_gap`, `Cres2 = 1/Pres` (`.py:71–72`), and builds resonance thresholds as `Pres * Wfail_match` (`.py:78`).
- Mathematica takes `Cres2Prim` as the primitive: `Cres2 = Cres2Prim`, `Pres = 1/Cres2` (`.wl:46–47`), recovers `PresGap` by `Solve[Pres == 1 + presGapFree, presGapFree]` (`.wl:48–53`), and builds resonance thresholds as `WfailMatch/Cres2` (`.wl:66`).
The Mathematica `presGapConsistency` check (`.wl:57,97`) has no SymPy counterpart and verifies the Solve inversion. The choreography differs (Pres-first vs Cres2-first, multiplication vs division, Solve-derived gap vs free gap), so this is an independent re-derivation of the same identities, not a line-by-line port.

## Engine cross-check

Outputs agree on every shared deliverable (modulo the symbol naming `Delta_gap`/`DeltaGap` and the primitive choice). Examples: matched fail threshold sympy `Pe_req/(Delta_0 + Delta_gap)` vs wl `PeReq/(Delta0 + DeltaGap)`; resonance fail threshold sympy `Pe_req*(Pres_gap+1)/(Delta_0+Delta_gap)` vs wl `PeReq/(Cres2Prim*(Delta0+DeltaGap))` — equal since `Pres_gap+1 = 1/Cres2Prim`; `matched window width = 0`, all side-band-width residuals `= 0`, all positivity `PASS` in both. No `engine_disagreement`.

## Verdict justification

verdict = findings, driven solely by F1 (`stale_output`, low). On substance the checkpoint holds: the scripts assemble the three-zone verdict and the exact side-band geometry as genuine symbolic identities (window width, side-band widths, `P_res-1=(1-C²)/C²`, and strict interior-point ordering) that would fail if the threshold forms were wrong; both engines are present and independent; paper alignment is exact. Attacks tried that failed: (a) tautology hunt — the load-bearing checks (A5, A8, A9, A10 / B4, B5) are real algebraic identities, not re-assertions; only the secondary `Pres_from_ratio` check (A6) is tautological-by-construction in SymPy, but it is not load-bearing and the Mathematica route gives the same identity genuine cross-route content (B1/B2). (b) hardcoded checkpoint constant — none; `C_res^2`/`P_res` are kept symbolic (`Pres_gap`/`Cres2Prim`), so no numeric is re-asserted against itself. (c) symbol-domain attack — all symbols positive/real as the physics requires; `expect_positive` results (`Delta_gap`, `Pres_gap`, etc.) are genuinely positive under the declared assumptions. (d) paper-misalignment — every card/notes deliverable maps to a matching check. The only defect is the committed `.txt` banner self-label "STAGE 052", which the refresh corrects to "STAGE 069".

## Self-test notes

Checked traps: (1) variable independence — the single derivative `dW_match/dDelta_eff` (sympy 109 / wl 113) genuinely depends on `Delta_eff`, so it is nonzero and the monotonicity check is real, not vacuous. (2) No integrals/parity to check. (3) Trivial-case: substituting `Pres_gap→0` (no penalty) collapses both side-bands to zero width, consistent with `delta_fail,delta_succ ∝ Pres_gap`; nonzero `Pres_gap` gives strictly positive widths — matches outputs. (4) F1 is output-refresh only (no source edit), no path-spec risk. (5) Paper round-trip: the fix changes only the stale banner self-label, introducing no new constant or claim — alignment stays exact.

## Value Reconciliation (pass-2 augmentation)

The scripts emit only **symbolic** deliverables (no pinned numeric constants — `C_res^2`/`P_res` are carried as the free symbols `Cres2Prim`/`Pres_gap`). The deliverable-level values are the threshold and side-band forms.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Matched fail threshold `Pe_req/Delta_inf` | py:76, wl:61; txt sympy:5, mma:5 | `.tex:18` (eq box), `.md:34,68` | MATCH |
| Matched success threshold `Pe_req/Delta_0` | py:77, wl:62; txt sympy:6, mma:6 | `.tex:19` (eq box), `.md:35,79` | MATCH |
| Resonance fail threshold `Pe_req/(C_res^2 Delta_inf)` | py:134–135, wl:88–91; txt sympy:18, mma:7 | `.md:96` (`P_res·Pe_req/Delta_inf`) | MATCH |
| Resonance success threshold `Pe_req/(C_res^2 Delta_0)` | py:138–139, wl:92–95; txt sympy:19, mma:8 | `.md:100` | MATCH |
| Matched window width `Pe_req(Δ_inf-Δ_0)/(Δ_0 Δ_inf)` | py:114–117, wl:128–131; txt sympy:13, mma:33 | `.tex:20` (band), `.md:88` | MATCH |
| Side-band relative width `P_res-1=(1-C_res^2)/C_res^2` | py:159–162, wl:148–151; txt sympy:26, mma:43 | `.md:92,104,151` (`≈0.56%`, `P_res=1.005612`) | MATCH |
| Failure-side abs width `Pe_req(1-C²)/(C² Δ_inf)` | py:151–153; txt sympy:24 | `.md:96` (side-band def) | MATCH |
| Success-side abs width `Pe_req(1-C²)/(C² Δ_0)` | py:155–157; txt sympy:25 | `.md:100` (side-band def) | MATCH |
| Three-zone verdict (fail/success/band) | FINAL LEDGER (py:175–186, wl:163–174) | `.tex:14–22` boxed, `.md:64–104` | MATCH |

INTERNAL scaffolding (accounted for, no finding): `u_fail = v_fail/(1+v_fail)`, `u_succ` (interior-point parameterization); `W_failure_band`, `W_success_band` (interior probe points); `Delta_inf - Delta_0 = Delta_gap`, `Delta_gap`-positivity residuals; `presGapConsistency` (Solve self-consistency); `PresFromSolve`, `PresGapFromSolve` (Mathematica intermediate); all `PASS`/`= 0` flags.

reconciliation: complete; 9 deliverable values checked, 0 misaligned. The numeric `C_res^2 = 0.994418836451529...` / `P_res = 1.005612487760576...` are Stage 067/068 deliverables (notes `.md:19–21,49–52,92`); they are not pinned numerically in the stage-069 scripts (kept symbolic), so they are not stage-069 reconciliation items — no MISMATCH/MISSING-DELIVERABLE.
