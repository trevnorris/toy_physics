---
unit_id: 229
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 229 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_229.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at lines 70, 809-837 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_mathematica_audit.txt`

## What the paper claims

Stage 229 pulls the selected-branch normalization product `N_-(x)` onto dimensionless softening-depth variables `ξ=x/A`, `δ=ΔK_ax/A` and classifies which of the two rigid Stage-228 signatures (numerator-rigid vs. denominator-rigid) the real selected branch resembles. The `\stagefield{Output}` reads verbatim: "Selected-branch classifier theorem: the pure-transfer defect is denominator-like near softening, with crossover controlled by the exact threshold $\delta=8/9$." The derivation ledger states the classifier threshold `δ=8/9` and the near-softening denominator-like signature. The notes (the authoritative detail source, since the `.tex` is terse) enumerate the distinct deliverables: (1) the exact factorization `F = F_num · F_den`; (2) the exact log-slope classifier `R_ND` and its component log-slopes `L_num`, `L_den`; (3) the onset law `R_ND(0,δ)=8/(9δ)`; (4) the near-softening limit `lim_{ξ→1⁻} R_ND = 0` (with `L_num` staying finite); (5) the exact crossover cubic `P(ξ,δ)=121ξ³+297δξ²+333δ²ξ+81δ³−72δ²=0`, its derivative `∂_ξP=363ξ²+594δξ+333δ²>0` (strict monotonicity), and `P(0,δ)=9δ²(9δ−8)`; (6) the `δ=8/9` threshold separating always-denominator from mixed regimes; (7) the three sample crossover depths for `δ=1/4,1/2,3/4`. The appendix narrative (lines 809-846) restates the classifier and the `δ=8/9` threshold. The card explicitly states "Mathematica audit: none yet" — a stale `\stagefield{Verification}` line, since a `.wl` now exists (noted below, not a finding requiring user resolution since the appendix/notes do not claim single-engine, and a richer-than-`.py` `.wl` is present).

## What the script claims to verify

The SymPy script derives `s_-`, `N_-`, reduces `N_-(Aξ, Aδ)` to `(8β₀/π²A)·F` and asserts the residual is exactly 0 (line 49); asserts the factorization `F = F_num·F_den` (line 53); computes `L_num`, `L_den`, `R_ND` by symbolic differentiation of `log F_num`, `log F_den` and asserts each matches the boxed closed form (lines 71-73); asserts onset `R_ND(0,δ)=8/(9δ)` (line 76), softening limit `=0` (line 79), and finite `L_num` softening limit (line 83); derives the crossover cubic numerator of `R_ND−1` and asserts it equals `P` with leading coeff `121ξ³` (lines 95-97), asserts `∂_ξP` form (lines 99-101) and `P(0,δ)=9δ²(9δ−8)` (line 104); then numerically solves `P=0` for the three sample `δ` and asserts the roots plus left/right side classification (lines 122-135), and a `δ=1` always-denominator slice (lines 141-144). The Mathematica script (M1-M10) verifies the same deliverable list by independent operations and additionally proves `∂_ξP>0` universally via `Resolve[ForAll]` (lines 170-180).

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| Factorization `F=F_num·F_den`, `F_den=1/(1−ξ)` | py:49,53 / wl:M1,M2 | match |
| Log-slope classifier `R_ND` + `L_num`,`L_den` | py:71-73 / wl:M3,M4 | match |
| Onset `R_ND(0,δ)=8/(9δ)` | py:76 / wl:M5 | match |
| Near-softening `lim R_ND=0`, finite `L_num` | py:79,83 / wl:M6 | match |
| Crossover cubic `P` (`121ξ³…`) + `P(0,δ)=9δ²(9δ−8)` | py:95-97,104 / wl:M7,M9 | match |
| Strict monotonicity `∂_ξP>0` | py:99-101 (form only) / wl:M8 (form + `Resolve[ForAll]` proof) | match (wl stronger) |
| Threshold `δ=8/9` | py:110 narrative + P(0,δ) sign / wl:M9 `threshold−8/9` | match |
| Sample depths `δ=1/4,1/2,3/4` | py:122-135 / wl:M10 | match |
| Always-denominator `δ≥8/9` regime | py:141-144 (`δ=1` slice) / wl:M10 (`δ=1` slice) | match |

`paper_alignment: aligned`. Every notes/appendix deliverable has a non-tautological script-side check on both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 49 | `reduced == 0` | factorization (deliverable 1) | yes |
| A2 | sympy | 53 | `simplify(F − F_num·F_den)==0` | factorization | yes |
| A3 | sympy | 71-73 | classifier residuals `==0` | classifier (2) | yes |
| A4 | sympy | 76 | `onset − 8/(9δ)==0` | onset (3) | yes |
| A5 | sympy | 79,83 | softening limits | near-softening (4) | yes |
| A6 | sympy | 97 | `numerator − P == 0` (derives P from R_ND−1) | crossover cubic (5) | yes |
| A7 | sympy | 101,104 | `∂_ξP` form, `P(0,δ)` | monotonicity/threshold (5,6) | yes (form) |
| A8 | sympy | 125,131-134 | root + side classification | sample depths (7) | yes |
| M1 | mathematica | 82-83 | `selectedN − scaleN·FClaim`, `F − FClaim` zero | factorization (1) | yes |
| M2 | mathematica | 90 | `F − FNum·FDen` zero | factorization | yes |
| M3 | mathematica | 100-101 | log-derivative residuals zero | classifier (2) | yes |
| M4 | mathematica | 109 | `RND − RNDClaim` zero | classifier | yes |
| M5 | mathematica | 112 | onset zero | onset (3) | yes |
| M6 | mathematica | 130-135 | softening limits + pole reciprocal | near-softening (4) | yes |
| M7 | mathematica | 152-157 | CoefficientList of derived P vs target (incl. `121`) | crossover cubic (5) | yes |
| M8 | mathematica | 163-180 | `∂_ξP` coeffs + `Resolve[ForAll] dP>0` | monotonicity (5) | yes (proves positivity) |
| M9 | mathematica | 185-186 | `P(0,δ)−9δ²(9δ−8)`, `threshold−8/9` | threshold (6) | yes |
| M10 | mathematica | 197-221 | roots + side classification + δ=1 slice | sample depths / always-denom (7) | yes |

No tautological rows: in every case the asserted closed form is compared against an *independently derived* expression (substitution+simplify, or symbolic differentiation, or numerator-of-`R_ND−1`), not against itself.

## Findings

None. (Zero findings; no directive written.)

## Independent-derivation check (Mathematica)

The `.wl` is GENUINELY INDEPENDENT, not a transliteration. Evidence by corresponding section:

1. **F reduction.** SymPy: `sp.simplify(N_minus.subs(...) - (8β₀/π²A)·F)` then `assert == 0` (py:45-49). Mathematica: `selectedN = Cancel[Together[nMinus /. {x->Aξ,...}]]`, `F = Cancel[Together[selectedN/scaleN]]`, then `expectZero[selectedN - scaleN·FClaim]` and `expectZero[F - FClaim]` (wl:72-83). Different simplification engines (`simplify` vs `Cancel∘Together`) and the `.wl` additionally re-derives `F` from the constants and checks it against the claim, a check the `.py` folds into a single residual.
2. **Log-slope.** SymPy `sp.diff(sp.log(F_num), xi)` + `sp.simplify` (py:63-65). Mathematica `FullSimplify[D[Log[FNum], xi]]` (wl:92-93). Independent differentiation; the printed forms even differ in presentation (py output line 8 expands the denominator to `81δ³+261δ²ξ+297δξ²+121ξ³`, the `.wl` keeps the factored `(9δ+11ξ)(9δ²+18δξ+11ξ²)` form, M3 output line 10).
3. **Crossover cubic + positivity.** This is the decisive divergence. The `.py` only checks the *coefficient form* of `∂_ξP` (py:99-101) and never proves `∂_ξP>0`. The `.wl` adds a genuine universal-quantifier decision via `Resolve[ForAll[{xiPos,deltaPos}, Implies[xiPos>=0 && deltaPos>0, (dP/.…)>0]], Reals]` (wl:170-180), output line 46 `M8 derivative positivity = True`. This is a CAD-based proof operation with no SymPy counterpart — the opposite of a transliteration; the second engine does *more* than the first. Root-finding also differs (`sp.nroots` vs `NSolve[…, Reals]`).

The shared M1-M10 ordering tracks the stage's own deliverable list (notes section 9), not the `.py`'s internal choreography; that shared ordering is intrinsic to the physics, not evidence of porting. INDEPENDENCE CALL: **independent**.

## Engine cross-check

Both engines emit the same load-bearing results and agree:
- `F = (9δ+11ξ)⁴/(81(1−ξ)(9δ²+18δξ+11ξ²)²)` — py output line 5, wl output line 3 (same up to sign-normalized `(−1+ξ)` form). Match.
- `R_ND(0,δ)=8/(9δ)`, `lim_{ξ→1⁻}R_ND=0` — py lines 11-12, wl lines 19-22. Match.
- Crossover cubic leading coeff `121` — py line 15 (`121*xi**3`), wl line 30 (CoefficientList `…, 121`). Match.
- `P(0,δ)=δ²(81δ−72)=9δ²(9δ−8)` — py line 17, wl line 48. Match.
- Sample roots `0.107223051105697`, `0.081847937860074`, `0.032505121082825` — py lines 19-21, wl lines 53-66, errors `~1e-16`. Match.
- `δ=1` always-denominator slice all `<1` — py line 23, wl lines 68-70. Match.

`engines_agree: true`. Both outputs are fresh (`.txt` mtime Jun 2 17:40, both newer than their respective sources). `outputs_fresh: true`.

## Crossover-cubic reconfirmation (re-audit mandate)

CONFIRMED. The corrected value `121ξ³` holds on all three carriers:
- SymPy derives `P` from the numerator of `R_ND−1` (py:95) and asserts it equals `121ξ³+297δξ²+333δ²ξ+81δ³−72δ²` (py:96-97); output line 15 emits `…+ 121*xi**3`.
- Mathematica derives `P` from `Numerator[Together[RND−1]]` (wl:137-139), forms the `CoefficientList`, and matches it against `targetCoeffs = {…, 121}` for the `ξ³` slot (wl:145, M7 coefficient 3 residual 0, output line 30/37).
- Notes line 292 reads `121\xi^3+297\delta\xi^2+333\delta^2\xi+81\delta^3-72\delta^2`.
- NO surviving `189` anywhere in the stage's notes, `.tex` card, or the Part VII appendix (grep over all three returned nothing). The leading coefficient is `121ξ³` on every side; the −68 correction (189→121) is fully landed and cross-engine corroborated.

## Verdict justification

`clean`. Attacks attempted that failed to find a defect: (1) checked whether the crossover cubic is hardcoded — it is *derived* from `R_ND−1` on both engines, not pinned; the asserted `121` form is compared to the derived numerator; (2) checked the `∂_ξP>0` monotonicity claim for a missing-branch/insufficient-verification gap — the `.py` only checks its coefficient form but the `.wl` supplies a genuine `Resolve[ForAll]` positivity proof, so the deliverable is covered; (3) checked the symbol domains (`positive=True` / `xi>0,xi<1,delta>0`) against the physical setup (`0≤ξ<1`, `δ>0`) — consistent, and the `xi→1⁻` directional limits respect the open upper bound; (4) checked the sample-root branch selection (`select_stable_real_root` / `stableRoot` both require a *unique* root in `(0,1)`, matching the strict-monotonicity theorem) — non-degenerate; (5) checked engine agreement on every emitted value — all match; (6) re-confirmed the 189→121 correction holds on script, notes, `.tex`, and appendix. I read the `.tex` card, the single notes file, and the Part VII appendix rows; the script's verified claim matches the paper's stated claim. One non-blocking observation: the card's `\stagefield{Verification}` line still says "Mathematica audit: none yet" while a substantive `.wl` now exists — this is a stale prose line in the card (deferred prose-drift class, like the residual stage-number labels), not a math defect and not a value mismatch, so it does not change the verdict; flagged here for the doc-pointer sync rather than as a user-gated `paper_misalignment`.

## Self-test notes

Trap 1 (variable independence): the only derivatives are `D[Log[F_num/F_den], ξ]` and `D[P, ξ]`; `F_num`,`F_den`,`P` all genuinely depend on `ξ`, so no identically-zero derivative — the `assert == claim` checks are non-trivial. Trap 2 (parity): no unbounded-domain integrals in this stage; n/a. Trap 3 (trivial-case): substituting `δ=1` into `P(0,δ)=9δ²(9δ−8)=9>0` (always-denom, consistent with `δ≥8/9`) and `δ=1/4` into `P(0,δ)=9/16·(9/4−8)<0` (mixed, root exists) both reduce correctly, matching the regime split. Trap 5 (paper round-trip): no script-side fix prescribed, and the verified `121ξ³` exactly matches the corrected notes — no new misalignment introduced.

## Value Reconciliation (pass-2 augmentation)

Every emitted deliverable value was located in the notes `.md` (authoritative detail carrier) and, for load-bearing items, the `.tex` card / Part VII appendix.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `F=(9δ+11ξ)⁴/(81(1−ξ)(9δ²+18δξ+11ξ²)²)` | py:44 / wl:75-77, py out:5, wl out:3 | notes:135-141 (boxed); appendix eq via N_-/s_- :811-824 | MATCH |
| `F_num=(9δ+11ξ)⁴/(81(9δ²+18δξ+11ξ²)²)`, `F_den=1/(1−ξ)` | py:51-52 / wl:85-89 | notes:152-162 (boxed) | MATCH |
| `L_num=72δ²/((9δ+11ξ)(9δ²+18δξ+11ξ²))` | py:67, py out:8 / wl:94-95, wl out:10 | notes:191-197 (boxed) | MATCH |
| `L_den=1/(1−ξ)` | py:68 / wl:97 | notes:199-203 (boxed) | MATCH |
| `R_ND=72δ²(1−ξ)/((9δ+11ξ)(9δ²+18δξ+11ξ²))` | py:69, py out:10 / wl:104-107, wl out:16 | notes:206-214, 447-453 (boxed) | MATCH |
| onset `R_ND(0,δ)=8/(9δ)` | py:76, py out:11 / wl:111-112, wl out:19 | notes:239, 473 | MATCH |
| `lim_{ξ→1⁻}R_ND=0` | py:79, py out:12 / wl:130, wl out:21 | notes:263-267 | MATCH |
| `L_num` softening limit `=72δ²/((9δ+11)(9δ²+18δ+11))` | py:82 / wl:129 | notes:255-260 | MATCH |
| crossover cubic `P=121ξ³+297δξ²+333δ²ξ+81δ³−72δ²` | py:96, py out:15 / wl:141-146, wl out:29-30 | notes:288-295 (boxed) | MATCH (121, not 189) |
| `∂_ξP=363ξ²+594δξ+333δ²` | py:100, py out:16 / wl:159-161, wl out:39 | notes:298-305 (boxed) | MATCH |
| `P(0,δ)=9δ²(9δ−8)` | py:104, py out:17 / wl:185, wl out:48 | notes:316-319 | MATCH |
| threshold `δ=8/9` | py:110 / wl:183,186, wl out:51 | tex:13,15 (`δ=8/9`); appendix:835; notes:60-63,455 | MATCH |
| sample `δ=1/4 → ξ_*≈0.107223051105697` | py:117, py out:19 / wl:189, wl out:53-54 | notes:385-387 | MATCH |
| sample `δ=1/2 → ξ_*≈0.081847937860074` | py:118, py out:20 / wl:190, wl out:58-59 | notes:389-391 | MATCH |
| sample `δ=3/4 → ξ_*≈0.032505121082825` | py:119, py out:21 / wl:191, wl out:63-64 | notes:393-395 | MATCH |
| D/N constants `κ₀²=8/π²`, `κ₁²=16/(9π²)` | py:35-36 / wl:61-62 | notes:100-104 | MATCH |
| `Output` classifier theorem + `δ=8/9` | (overall) | tex:15 verbatim | MATCH |

INTERNAL scaffolding (no finding): `assert_close`/`expectNumeric` tolerances (`1e-12`, `5e-13`); the `δ=1` always-denominator probe set `{1/100,1/5,3/5,9/10}` and their numeric `R_ND` values (a slice-check, not a stated deliverable); left/right side-probe values in M10 (classification scaffolding); residual `= 0` lines; `select_stable_real_root`/`stableRoot` helper logic.

reconciliation: complete; 17 deliverable values checked, 0 misaligned.
