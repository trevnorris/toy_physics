---
unit_id: 040
batch: III.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T06:04:06Z
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
  notes_stage_files: [moving_throat_pde_stage040_generalized_selected_branch.md]
  paper_appendix: present
---

# Audit unit 040 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_040.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage040_generalized_selected_branch.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 58)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.txt`

## What the paper claims

Stage 040 generalizes the flat one-vector selected branch (Stages 035/036) to a **source/loading mismatch** law: the loading vector `z = z0(1,q)^T` and source vector `s = s0(1,t)^T` are no longer collinear. `\stagefield{Output}{The generalized selected functions \eqref{eq:app-stage040-Fqeta}--\eqref{eq:app-stage040-Gq} and their split--U specialization \eqref{eq:app-stage040-FU}--\eqref{eq:app-stage040-GU}.}` The distinct deliverables are: (1) required loading `α_req = A0 ξ(ξ+δ)/(z0^2[δ+(1+q^2)ξ])`; (2) eigenvector ratio `e1/e0 = qξ/(δ+ξ)`; (3) the boxed two-vector normalization `F_{q,η_s} = [δ+(1+q^2)ξ]^2[δ+(1+η_s)ξ]^2 / {(1−ξ)[(δ+ξ)^2+q^2ξ^2]^2}`; (4) the boxed loading `G_q = ξ(ξ+δ)/[δ+(1+q^2)ξ]`; (5) the split-U specialization with `q = −(√2/3)R_U`, `η_s = (2/9)R_U`, `λ0 = 2/9`, giving boxed `F_U` and `G_U`; (6) recovery of the Stage-035/036 flat-U branch at `R_U=1`. The notes additionally enumerate (Section 5) the first-order deformation coefficients `H_F = 4ξ(27δ^2+36δξ+11ξ^2)/[(9δ+11ξ)(9δ^2+18δξ+11ξ^2)]` and `H_G = −4ξ/(9δ+11ξ)`. The appendix row (58) summarizes: "Source/loading mismatch functions F_{q,η_s}, G_q, and split--U functions F_U, G_U."

## What the script claims to verify

The SymPy script (docstring lines 3-16) verifies five things matching the paper deliverables: (1) closed-form `α_req` and eigenvector for a diagonal 2×2 baseline plus rank-1 loading `z`; (2) the two-vector `F_{q,η}` deformation when `s` is non-collinear with `z`; (3) the split-U specialization `F_U`, `G_U`; (4) `R_U=1` recovers Stage-035/036; (5) the first-order deformation `H_F`, `H_G` is exact, cross-checked two independent ways (differentiate `F_U` vs. differentiate `F_general` under the eps-parametrized `(q,η)`). The Mathematica script verifies the same chain but derives `α_req` by `Solve[Det[M−αzz^T−λI]==0]`, the eigenvector by `NullSpace`, and the overlaps by dotting the NullSpace eigenvector with explicit `z` and `s` vectors — a genuinely independent route. Every `expect_zero`/`expectZero` in both engines reports `= 0`/PASS in the saved outputs.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) `α_req` form | sympy L59-61 / wl L52-55 (`Solve[Det=0]`); output `A0 xi(δ+xi)/(z0^2(δ+xi+q^2 xi))` = paper `δ+(1+q^2)ξ` | match |
| (2) `e1/e0 = qξ/(δ+ξ)` | sympy L64-66, residual L73-79 / wl L58-65, residual L70-77 | match |
| (3) `F_{q,η_s}` boxed | sympy L89,95,98 (`F_general − F_expected == 0`) / wl L90,99-105 | match |
| (4) `G_q` boxed | sympy L90,96,99 / wl L91,104,106 | match |
| (5a) split-U `q=−(√2/3)R_U`, `η_s=(2/9)R_U`, `λ0=2/9` | sympy L51,105-106 / wl L40,110-111 | match |
| (5b) `F_U`, `G_U` boxed | sympy L107-111 / wl L112-113,126-127 | match |
| (6) `R_U=1` recovers Stage035/036 F,G | sympy L118-122 / wl L120-129 | match (anchored to stage035 L54, stage036 L53-54) |
| notes §5 `H_F`, `H_G` | sympy L132-149 / wl L136-150 (dual-route cross-check) | match |

`paper_alignment: aligned`. Every paper-side boxed result and the notes' first-order coefficients have a non-tautological script-side check, and the constants (`λ0=2/9`, the `q,η` map, the `11`-coefficients) match the paper verbatim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 66 | `expect_zero(r − xiq/(δ+ξ))` | claim 2 (e1/e0) | yes |
| A2 | sympy | 78-79 | `expect_zero(eig_residual[0]/[1])` | claim 1+2 (α_req+evec satisfy M e=λe) | yes |
| A3 | sympy | 98 | `F_general − F_expected == 0` | claim 3 (F_{q,η}) | yes |
| A4 | sympy | 99 | `G_general − G_expected == 0` | claim 4 (G_q) | yes |
| A5 | sympy | 121 | `F_U(R_U=1) − F_stage18 == 0` | claim 6 (recovery) | yes |
| A6 | sympy | 122 | `G_U(R_U=1) − G_stage19 == 0` | claim 6 (recovery) | yes |
| A7 | sympy | 148 | `HF − HF_direct == 0` | notes §5 (H_F dual-route) | yes |
| A8 | sympy | 149 | `HG − HG_direct == 0` | notes §5 (H_G dual-route) | yes |
| B1 | mathematica | 65 | `expectZero(r − xi q/(δ+ξ))` | claim 2 (e1/e0 from NullSpace) | yes |
| B2 | mathematica | 76-77 | `expectZero(eigResidual[[1]]/[[2]])` | claim 1+2 | yes |
| B3 | mathematica | 105 | `expectZero(fGeneral − fExpected)` | claim 3 | yes |
| B4 | mathematica | 106 | `expectZero(gGeneral − gExpected)` | claim 4 | yes |
| B5 | mathematica | 128 | `expectZero(fU(rU=1) − fStage18)` | claim 6 | yes |
| B6 | mathematica | 129 | `expectZero(gU(rU=1) − gStage19)` | claim 6 | yes |
| B7 | mathematica | 149 | `expectZero(hF − hFDirect)` | notes §5 | yes |
| B8 | mathematica | 150 | `expectZero(hG − hGDirect)` | notes §5 | yes |

No row is tautological or unanchored. A2/B2 are the load-bearing physics checks: they substitute the *derived* `α_req` and the *claimed* closed-form eigenvector into `(M − αzz^T)e − λe` and require it to vanish — this fails if either `α_req` or the eigenvector closed form is wrong. A3/A4 require the overlap-built `F_general`/`G_general` to equal the boxed `F_expected`/`G_expected`; the SymPy `F_expected` is the paper's boxed form transcribed (L95-96), so A3/A4 verify the overlap construction reproduces the boxed claim — non-tautological because `F_general` is built independently from `z_overlap_sq · s_overlap_sq · A0/λ_minus`.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.txt` (mtime 2026-05-22 12:32:20)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.txt` (mtime 2026-05-22 12:32:29)
- vs. sympy `.py` mtime 2026-06-03 15:59:11 and mathematica `.wl` mtime 2026-06-03 15:59:11

**What's wrong:**
Both saved transcripts predate both scripts by ~12 days, and their content shows the pre-renumber stage label. The current SymPy banner is `banner("STAGE 40 — ...")` (`.py:53`) and the ledger header is `"STAGE 40 THEOREM LEDGER"` (`.py:151`), but the saved output reads `STAGE 23 — ...` (`.txt:3`) and `STAGE 23 THEOREM LEDGER` (`.txt:44`). The current Mathematica banner is `banner["STAGE 040 — ..."]` (`.wl:33`) but the saved output reads `STAGE 023 — ...` (`.txt:3`) while its closing line `Stage 040 Mathematica audit passed.` (`.txt:52`) does match `.wl:153` — a mid-renumber transcript. The SymPy docstring (`.py:3` "Stage 23", and "Stage-18/19/22" references at `.py:7-15`) and the subbanner labels (`.py:55,81,101,124` "23.1"–"23.4") are likewise the pre-renumber numbering, reflected verbatim in the saved output. **All numeric/symbolic results in the transcripts are identical to what the current scripts compute** (α_req, F_{q,η}, G_q, F_U, G_U, H_F, H_G all match the script source); only the banner/docstring stage *labels* are stale.

**Why this matters:**
The transcript is informational and not blocking — the math is unchanged — but the committed `.txt` does not reflect the current banner. This is the known SCRIPT/OUTPUT-band numbering drift (committed outputs predate the banner fix; bare docstring/subbanner self-labels were not part of the prior notes-only renumber sweep). The verifier's independent re-run will refresh the transcript.

**Required change:**
None applied by Codex for the stale label itself — this is the deferred SCRIPT/OUTPUT-band cleanup (`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`, user-decided to defer). The verifier's fresh re-run of both scripts will regenerate the transcripts with the current `STAGE 40`/`STAGE 040` banners. No math edit is warranted.

**Verification:**
After the orchestrator's independent exec re-run, the refreshed `.txt` headers should read `STAGE 40 — ...` (sympy) and `STAGE 040 — ...` (mathematica), and all result lines (α_req, F, G, H_F, H_G) should be byte-identical to the current saved values. If any *result* line changes on re-run, escalate; if only the banner changes, this finding is fully resolved.

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration. Corresponding sections take materially different routes:

1. **α_req.** SymPy types the closed form and simplifies (`.py:59` `alpha_req = sp.simplify(A0*xi*(delta+xi)/(z0**2*(delta+xi)+(q*z0)**2*xi))`). Mathematica derives it: `charEq = Det[mPert[alpha] − lamMinus IdentityMatrix[2]] == 0; alphaSol = Solve[charEq, alpha]` (`.wl:52-54`). Different mechanism (posited-then-validated vs. solved-from-determinant).
2. **Eigenvector.** SymPy computes `r` from row-0 algebra `(A0 xi − α z0^2)/(α z0(q z0))` (`.py:64`). Mathematica uses `NullSpace[mPert[alphaReq] − lamMinus I]` then normalizes the first component to 1 (`.wl:58-61`).
3. **Overlaps.** SymPy pre-substitutes the η relation and writes `s_overlap = (1 + eta*xi/(delta+xi))^2/(1+r^2)` (`.py:84`). Mathematica builds `sVec = {1, eta/q}`, `sDotE = sVec.eMinus`, `normESq = eMinus.eMinus` using the *NullSpace* eigenvector (`.wl:84-89`) — a fuller, more independent construction.

The two engines arrive at the same boxed `F`, `G`, `F_U`, `G_U`, `H_F`, `H_G`. No `mathematica_transliteration` finding.

## Engine cross-check

Both transcripts agree at the level they claim (all `expectZero`/`expect_zero` report 0 / PASS):

- `α_req`: sympy `A0*xi*(delta+xi)/(z0**2*(delta + q**2*xi + xi))` (out L9) ≡ mathematica `(a0*xi*(delta+xi))/((delta+xi+q^2*xi)*z0^2)` (out L9). Same.
- `F_(q,η)`: sympy `(delta+xi+eta*xi)^2*(delta+xi+q^2*xi)^2/((1−xi)*(...)^2)` (out L23) ≡ mathematica `−(delta+eta*xi+xi)^2*(delta+q^2*xi+xi)^2/((xi−1)*(...)^2)` (out L20). The Mathematica `−.../(xi−1)` form is algebraically identical to the `1/(1−xi)` SymPy form; the `expectZero["F_general − expected"]` PASS confirms equality.
- `F_U`, `G_U`, `H_F`, `H_G`: identical between the two transcripts (sympy out L33-46, mathematica out L33-49). `H_F = 4*xi*(27*delta^2+36*delta*xi+11*xi^2)/((9*delta+11*xi)*(9*delta^2+18*delta*xi+11*xi^2))` and `H_G = −4*xi/(9*delta+11*xi)` in both, matching notes §5 verbatim.

`engines_agree: true`.

## Verdict justification

`findings` with a single low-severity `stale_output` (deferred SCRIPT/OUTPUT-band numbering label). The math is sound and **paper-aligned**: every boxed paper deliverable (`α_req`, `e1/e0`, `F_{q,η_s}`, `G_q`, `F_U`, `G_U`, the `R_U=1` recovery) and both notes-§5 first-order coefficients have a non-tautological, well-anchored check in both engines, and the two engines derive the chain by independent routes and agree. Attacks tried that failed: (a) checked the `α_req` denominator `z0^2(δ+xi)+q^2 z0^2 xi` actually equals the paper's `z0^2[δ+(1+q^2)ξ]` — it does; (b) checked the `R_U=1` recovery is non-trivial (needs `λ0=2/9` so `9(1+q^2)=9(1+2/9)=11` and `9(1+η)=11`) — it is, the `11`-coefficients only appear after correct specialization; (c) checked the eigenvector residual A2/B2 can fail if `α_req` or the closed-form eigenvector is wrong — it can; (d) checked the upstream `F_stage18`/`G_stage19` literals against stage035 (`:54`) and stage036 (`:53-54`) — they match, so the recovery target is genuinely anchored, not a free hardcoded value; (e) checked the Mathematica is not a SymPy port — it solves the determinant and uses NullSpace, materially independent. I read the stage card, the notes, and the appendix row; the script's claim matches the paper's claim. The only defect is the stale transcript label, which the verifier's re-run resolves.

## Self-test notes

Checked: (1) Variable-independence of the `D[·,eps]` deformation derivatives — `F_U`/`F_general` genuinely depend on `eps` after the `R_U→1+eps` / `(q,η)→eps`-parametrized substitution (both routes give the same nonzero `H_F`, `H_G`), so the derivative is not identically zero and the cross-check `HF − HF_direct == 0` is substantive, not trivially-passing. (2) No unbounded-domain integrals in this unit, so parity traps are N/A. (3) Trivial-case: at `R_U=1` the split-U numerator coefficients collapse to `11` (`9+2=11`), recovering the Stage-035/036 quartic — confirmed the recovery checks reduce correctly. No directive-applied math change is prescribed (the sole finding is a deferred-band stale label), so the paper round-trip self-test is moot.

## Value Reconciliation (pass-2 augmentation)

Enumerating every RESULT/deliverable value the scripts emit (from `.py`/`.wl` source + saved `.txt`), and locating each in the `.tex` card and `.md` notes:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `α_req = A0 ξ(δ+ξ)/(z0^2[δ+(1+q^2)ξ])` | py L59 / wl L52-55; out L9 | `.tex:16-18` eq:app-stage040-alpha-req; `.md:53` | MATCH |
| `e1/e0 = qξ/(δ+ξ)` | py L66 / wl L65; out L11/L10 | `.tex:22-24` eq:app-stage040-evec; `.md:63` | MATCH |
| `F_{q,η_s} = [δ+(1+q^2)ξ]^2[δ+(1+η_s)ξ]^2/{(1−ξ)[(δ+ξ)^2+q^2ξ^2]^2}` | py L92,95 / wl L95,99-103; out L20/L23 | `.tex:28-31` boxed eq:app-stage040-Fqeta; `.md:93-95` | MATCH |
| `G_q = ξ(ξ+δ)/[δ+(1+q^2)ξ]` | py L93,96 / wl L96,104; out L21/L24 | `.tex:35-36` boxed eq:app-stage040-Gq; `.md:99` | MATCH |
| split-U map `q = −(√2/3)R_U`, `η_s = (2/9)R_U`, `λ0 = 2/9` | py L51,105-106 / wl L40,110-111 | `.tex:38,41-43` eq:app-stage040-qeta-RU; `.md:117,121,125` | MATCH |
| `F_U = [9δ+(9+2R_U^2)ξ]^2[9δ+(9+2R_U)ξ]^2/{81(1−ξ)[9δ^2+18δξ+(9+2R_U^2)ξ^2]^2}` | py L110 / wl L126; out L28/L33 | `.tex:48-51` boxed eq:app-stage040-FU; `.md:131-133` | MATCH |
| `G_U = 9ξ(ξ+δ)/[9δ+(9+2R_U^2)ξ]` | py L111 / wl L127; out L29/L34 | `.tex:55-56` boxed eq:app-stage040-GU; `.md:135-136` | MATCH |
| recovery `F_U(R_U=1)=(9δ+11ξ)^4/[81(1−ξ)(9δ^2+18δξ+11ξ^2)^2]`, `G_U(R_U=1)=9ξ(ξ+δ)/(9δ+11ξ)` | py L118-122 / wl L120-129; out L30-31/L35-38 | `.tex:58` "flat branch recovered at R_U=1"; `.md:147-151` (explicit forms) | MATCH |
| `H_F = 4ξ(27δ^2+36δξ+11ξ^2)/[(9δ+11ξ)(9δ^2+18δξ+11ξ^2)]` | py L143-144 / wl L144-145; out L36-37/L43-44 | not in `.tex` (terse card); `.md:173-175` §5 | MATCH (notes carrier) |
| `H_G = −4ξ/(9δ+11ξ)` | py L145-146 / wl L146-147; out L38-39/L45-46 | not in `.tex`; `.md:177` §5 | MATCH (notes carrier) |

INTERNAL scaffolding (accounted for, no finding): `lam_minus = A0(1−ξ)`, `M_base/mBase`, `z_vec/zVec`, `M_perturbed`, `e_minus`, `eig_residual`, intermediate `z_overlap_sq`/`s_overlap_sq` (these feed the boxed `F` and are surfaced as the paper's eq:app-stage040-evec-derived overlaps but are not standalone deliverables), `q_U_eps`/`eta_U_eps`, the two `HF`/`HF_direct` route variables (cross-check pair), all `expect_zero`/`expectZero` residual flags and PASS lines.

The `H_F`/`H_G` first-order coefficients are deliverables enumerated in the notes (§5) but legitimately omitted from the terse `.tex` card; per the augmentation guards, a value living correctly in the `.md` notes is a MATCH, not MISSING. All other deliverables appear in both `.tex` (boxed) and `.md`.

reconciliation: complete; 10 deliverable values checked, 0 misaligned
