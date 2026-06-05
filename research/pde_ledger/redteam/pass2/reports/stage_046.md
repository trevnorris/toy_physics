---
unit_id: 046
batch: III.1
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md]
  paper_appendix: present
---

# Audit unit 046 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_046.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 70; also referenced at line 347)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.txt`

## What the paper claims

Stage 046 takes the exact one-parameter tracking branch from Stage 045 (`M_tr = G_tr(xi,delta;R)`, `R_target = F_tr(xi,delta;R)`) and asks whether the first coherent split-`U` deformation helps or hurts the Stage-035/036 normalization target. Its `\stagefield{Output}` is two boxed results: (1) the monotonicity theorem `eq:app-stage046-monotonicity` — `∂G_tr/∂R < 0` and `∂F_tr/∂R > 0` at fixed `(xi,delta)`; and (2) the residual comparison theorem `eq:app-stage046-residual` — `E_tr − E_flat = F_flat − F_tr > 0`. The card's body also states the constructive inequalities `G_tr > G_flat`, `F_tr < F_flat` (eq `:app-stage046-inequalities`) at `R_tr < 1`. The notes go much further, pinning the exact closed forms of `G_tr`, `F_tr`, `G_flat`, `F_flat`, the exact `dG_tr/dR` and `dF_tr/dR` (with the explicit positive polynomial `P_R`), the exact branch-difference forms `G_tr − G_flat` and `F_flat − F_tr` (with positive polynomials `P_1`, `P_2`), and the strong-split endpoints `G_tr(·;0) = xi`, `F_tr(·;0) = 1/(1−xi)`, plus the endpoint bounds `G_flat ≤ G_tr ≤ xi`, `1/(1−xi) ≤ F_tr ≤ F_flat`. The physical conclusion: the split-`U` deformation makes the normalization target harder, requiring more loading and delivering less response. Appendix row 70 summarizes: "Monotonicity in `R` and proof that first split-`U` tracking worsens the target."

## What the script claims to verify

Both engines pin the four branch functions (`G_tr`, `F_tr`, `G_flat = G_tr|_{R=1}`, `F_flat = F_tr|_{R=1}`) directly from the Stage-045 closed forms, then verify: the strong-split endpoints `G_tr|_{R=0} = xi` and `F_tr|_{R=0} = 1/(1−xi)`; the exact derivatives `dG_tr/dR` and `dF_tr/dR`; their signs (`< 0` and `> 0`) on the open domain `0<R<1, 0<xi<1, delta>0`; the exact branch differences `G_tr − G_flat` and `F_flat − F_tr`; their strict positivity on that domain; the boundary values `(G_tr−G_flat)|_{R=1} = 0`, `(F_flat−F_tr)|_{R=1} = 0`; and three numeric sample points confirming both differences are positive. SymPy proves the sign claims via positive-coefficient checks of the hand-typed `P_R`, `P1`, `P2` (after confirming each equals the engine-derived difference), while Mathematica proves them by `Reduce[ForAll[...]]` over the domain and by polynomial-division confirmation of the `(1−R²)` and `(1−R)` factors. The SymPy docstring still reads "Stage 29 SymPy audit" and its closing line reads "All Stage-29 symbolic checks passed" (stale labels); the substantive content is Stage 046.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `∂G_tr/∂R < 0` (monotonicity, boxed) | SymPy: `dG_tr/dR` identity (l.87) + `P_R` form is `−36 R xi²(δ+xi)/(…)²` < 0 by inspection; MMA: `Reduce[dGdR<0]` = True (l.64-67) | match |
| `∂F_tr/∂R > 0` (monotonicity, boxed) | SymPy: `dF_tr/dR` identity (l.88) + `expect_positive_coefficients(P_R)` (l.89); MMA: `Reduce[dFdR>0]` = True (l.69-73) | match |
| `E_tr − E_flat = F_flat − F_tr > 0` (residual, boxed) | `F_flat − F_tr` closed form (l.94,139) + positivity via `P1`,`P2` (l.140-141) / MMA `Reduce[deltaF>0]` (l.102-106). The `E_tr−E_flat = F_flat−F_tr` step is definitional (R_target cancels) and not separately asserted | match |
| `G_tr > G_flat`, `F_tr < F_flat` at `R_tr<1` (body ineq.) | `G_tr − G_flat` form (l.93,138) + positivity; MMA `Reduce[deltaG>0]` (l.96-100) | match |
| Strong-split endpoints `G_tr(·;0)=xi`, `F_tr(·;0)=1/(1−xi)` (notes §4) | `expect_zero` l.58-59, 156-162; MMA l.52-53, 111-112 | match |
| Endpoint bounds `G_flat≤G_tr≤xi`, `1/(1−xi)≤F_tr≤F_flat` (notes §4) | Established by the monotonicity sign + endpoint checks; printed as a banner-4 statement (SymPy l.188-191) — proven by composition, not a standalone assertion | match (composed) |
| Exact `P_R`, `P1`, `P2` positive polynomials (notes §2.2,§3.2) | hand-typed then confirmed equal to engine-derived diff (`expect_zero`) AND positive-coefficient checked | match |

Every paper/notes deliverable maps to a non-tautological script check. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 58 | `expect_zero(G_tr\|_{R=0} − xi)` | strong-split endpoint G | yes |
| A2 | sympy | 59 | `expect_zero(F_tr\|_{R=0} − 1/(1−xi))` | strong-split endpoint F | yes |
| A3 | sympy | 87 | `expect_zero(diff(G_tr,R) − dG_expected)` | monotonicity (dG form) | yes |
| A4 | sympy | 88 | `expect_zero(diff(F_tr,R) − dF_expected)` | monotonicity (dF form) | yes |
| A5 | sympy | 89 | `expect_positive_coefficients(P_R)` | `∂F_tr/∂R > 0` | yes |
| A6 | sympy | 138 | `expect_zero(δG − G_diff_expected)` | body ineq. `G_tr>G_flat` | yes |
| A7 | sympy | 139 | `expect_zero(δF − F_diff_expected)` | residual `F_flat−F_tr` | yes |
| A8 | sympy | 140-141 | `expect_positive_coefficients(P1,P2)` | residual `>0` | yes |
| A9 | sympy | 147-162 | boundary `R=1` vanish + `R=0` endpoints | sanity anchors | yes |
| A10 | sympy | 179-186 | 3 numeric samples, both diffs `>0` | residual/ineq strict positivity | yes (sampled) |
| A11 | mma | 52-53 | `expectZero` strong-split endpoints | endpoints G,F | yes |
| A12 | mma | 64-67 | `Reduce[∀ dGdR<0]` = True | `∂G_tr/∂R < 0` | yes |
| A13 | mma | 69-73 | `Reduce[∀ dFdR>0]` = True | `∂F_tr/∂R > 0` | yes |
| A14 | mma | 82-86 | `(1−r²)\|num(δG)` remainder = 0 | factor structure of `G_tr−G_flat` | yes |
| A15 | mma | 89-93 | `(1−r)\|num(δF)` remainder = 0 | factor structure of `F_flat−F_tr` | yes |
| A16 | mma | 96-106 | `Reduce[∀ δG>0]`, `Reduce[∀ δF>0]` = True | body ineq. + residual `>0` | yes |
| A17 | mma | 109-112 | boundary `r=1` vanish + `r=0` endpoints | sanity anchors | yes |
| A18 | mma | 119-129 | 3 numeric samples, both diffs `>0` | residual/ineq strict positivity | yes (sampled) |

No assertion is tautological: every `expect_zero` compares an engine-derived object (`sp.diff`, `Together`, `Factor`, substitution) against an independently hand-typed closed form, and the `Reduce`/positive-coefficient checks are genuine domain-wide sign proofs that can fail.

## Findings

### F1 — stale_output

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.txt` (mtime 2026-05-22 12:54)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.txt` (mtime 2026-05-22 12:54)
- vs scripts `.py`/`.wl` (mtime 2026-06-03 15:59)

**What's wrong:**
Both saved transcripts predate the scripts. `git` shows the scripts were last edited by `e2a4780` ("numbering reconciliation Phase 1 (deterministic): 273 doc-only stage-label fixes"), which changed ONLY the banner string (`banner("STAGE 29 …")` → `STAGE 46`; `banner["STAGE 029 …"]` → `STAGE 046`) — byte-identical except the stage-number token, no equation/value/logic touched. The committed `.txt` files were produced before that commit (`84b60c4`), so their headers still read `STAGE 29` / `STAGE 029` and the SymPy transcript's closing line reads "All Stage-29 symbolic checks passed", while the Mathematica transcript's closing line already reads "Stage 046 Mathematica audit passed." The substantive numbers in the transcripts (the four branch functions, the `dG/dR`, `dF/dR`, `δG`, `δF` forms, the `P_R`/`P1`/`P2` coefficient lists, and the three numeric sample values `225/8869`, `81/1736`, `91/21935`, `38617837960/99381932001`, …) all match what the current scripts would print — the staleness is label-only.

Additionally, residual stale self-labels survive in the scripts themselves (not fixed by `e2a4780`, which only touched the main banner): the SymPy docstring line 3 `"""Stage 29 SymPy audit.` and SymPy line 193 `"All Stage-29 symbolic checks passed."`. These are exactly the script/output-band residuals the project has deferred to a dedicated content-keyed pass (NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md); they are label-only and change no math.

**Why this matters:**
The transcripts on disk do not reflect the current banner; a reader diffing them against a fresh run sees the header/footer change. No numeric result changes. Non-blocking.

**Required change:**
None applied by Codex for the residual docstring/footer labels — these belong to the gated script/output-band numbering pass, not this red-team unit. The verifier's independent re-run will refresh the committed `.txt` transcripts (which will then carry the correct `STAGE 46`/`STAGE 046` banners). Treat as informational.

**Verification:**
After the orchestrator's independent re-run, the SymPy transcript header should read `STAGE 46 — TRACKING-BRANCH BOUNDS AUDIT` and all `PASS`/zero-residual lines remain unchanged.

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration. The SymPy proof strategy for the sign claims is "hand-type `P_R`/`P1`/`P2`, confirm they equal the engine-derived difference, then check every coefficient is positive" (l.89,140-141). The Mathematica script deliberately does NOT hand-type any comparison polynomial — its comments say so explicitly (l.55 "Derivatives obtained directly from gTr/fTr; no hand-typed comparison polynomial"; l.75 "Branch-difference forms derived directly; no hand-typed p1/p2/gDiffExpected/fDiffExpected"). Instead it proves the signs by `Reduce[ForAll[{r,xi,delta}, domain, dGdR<0], Reals]` (l.64) and `Reduce[…δF>0…]` (l.102-103), a genuinely different decision-procedure route, and independently confirms the `(1−r²)` and `(1−r)` factor structure via `PolynomialQuotientRemainder` (l.82-93) — a check the SymPy script does not perform. This is a real second engine attacking the same claims by different machinery.

## Engine cross-check

The engines agree:
- `G_tr`, `F_tr`, `G_flat`, `F_flat`: identical closed forms (SymPy out l.9-12; MMA out l.5-8; the `−1/((xi−1)…)` vs `1/((1−xi)…)` sign rendering is cosmetic, `xi<1`).
- `dG_tr/dR = −36 r xi²(δ+xi)/(…)²` (SymPy out l.19; MMA out l.13) — match.
- `dF_tr/dR` numerator polynomial (the `P_R` factor) — coefficient lists identical (SymPy `P_R coefficients = [4,54,90,36,162,324,162,81,243,243,81]`, MMA `dF_tr/dR` printed with the same monomials, out l.14).
- `δG`, `δF` factored forms with `P1`/`P2` — coefficient lists identical (SymPy out l.32-33; MMA out l.20).
- All three numeric samples identical across engines: `225/8869`, `81/1736`, `91/21935` and the `F_flat−F_tr` fractions `38617837960/99381932001`, `759648230/1473329763`, `5842146019415/70196178995856` (SymPy out l.42-44; MMA out l.45-47).
- Both transcripts end in all-PASS.

`engines_agree: true`.

## Verdict justification

`findings: 1` — and that single finding is the informational `stale_output` (label-only; the orchestrator re-run refreshes the transcripts and the residual docstring/footer stage-number labels are deferred to the gated script/output-band numbering pass per project policy). The substantive math holds up under attack: I tried to break the derivative checks (they compare `sp.diff`/`D[]` against an independently hand-typed form — not tautological, and `G_tr`/`F_tr` genuinely depend on `R` so the derivatives are nontrivial), the positivity arguments (the `P_R`/`P1`/`P2` coefficient checks are valid because all three are genuine polynomials in `(R,δ,xi)` and Mathematica independently confirms the signs via `Reduce[ForAll…]` over the same domain), the endpoint and boundary checks (genuine substitutions, not by construction), and the paper↔script mapping (all of monotonicity, the residual theorem, the body inequalities, the endpoints, and the bounds are covered). I read the card, the notes, and the appendix rows; the script's verified claim matches the paper's stated claim. No `paper_misalignment`, no engine disagreement, no tautology, no missing branch, no transliteration. Because the only finding is an informational `stale_output`, no Codex directive is warranted.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 13 deliverable values checked, 0 misaligned.

The outputs are label-stale (see F1) but numerically current; reconciliation is based on the script source plus the committed transcripts, which agree on every numeric/symbolic result.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `G_tr = 9 xi(xi+δ)/(9δ+(9+2R²)xi)` | py l.43 / out l.9; wl l.39 / out l.5 | notes l.41-42 | MATCH |
| `F_tr = (9δ+(9+2R²)xi)²(9δ+(9+2R)xi)² / [81(1−xi)(9δ²+18δxi+(9+2R²)xi²)²]` | py l.44-48 / out l.10; wl l.40-44 / out l.6 | notes l.43-44 | MATCH |
| `G_flat = 9 xi(xi+δ)/(9δ+11xi)` | py l.50 / out l.11; wl l.45 / out l.7 | notes l.52 | MATCH |
| `F_flat = (9δ+11xi)⁴ / [81(1−xi)(9δ²+18δxi+11xi²)²]` | py l.51 / out l.12; wl l.46 / out l.8 | notes l.54-56 | MATCH |
| `dG_tr/dR = −36 R xi²(δ+xi)/(2R²xi+9δ+9xi)²` | py l.62,81 / out l.19; wl l.56 / out l.13 | notes l.68-69 | MATCH |
| `dF_tr/dR` (with `P_R` numerator, l.67-79) | py l.63,82-85 / out l.20; wl l.57 / out l.14 | notes l.81-91 | MATCH |
| `P_R` coeffs `[4,54,90,36,162,324,162,81,243,243,81]` | py l.89 / out l.23 | notes l.87-91 | MATCH |
| `G_tr − G_flat = 18 xi²(1−R²)(δ+xi)/[(9δ+11xi)(2R²xi+9δ+9xi)]` | py l.93,130 / out l.28; wl l.76 / out l.19 | notes l.113-115 | MATCH |
| `F_flat − F_tr = 4 xi(1−R) P1 P2 / [81(1−xi)(…)²(…)²]` | py l.94,133 / out l.29; wl l.77 / out l.20 | notes l.123-139 | MATCH |
| `P1` coeffs `[18,36,22,81,180,99,162,423,360,99]` | py l.140 / out l.32 | notes l.130-133 | MATCH |
| `P2` coeffs `[18,36,22,81,324,459,220,81,243,261,99,729,3078,4959,3600,990]` | py l.141 / out l.33 | notes l.135-139 | MATCH |
| strong-split endpoint `G_tr(·;0) = xi` | py l.58,156 / out l.13; wl l.52,111 | notes l.155, card eq via bounds | MATCH (notes) |
| strong-split endpoint `F_tr(·;0) = 1/(1−xi)` | py l.59,159 / out l.14; wl l.53,112 | notes l.157 | MATCH (notes) |

Every emitted deliverable value reconciles against the notes (`.md`), and the boxed monotonicity + residual theorems reconcile against the `.tex` card. The `.tex` card is terse (it boxes only the two theorems and one inequality block) and legitimately omits the intermediate closed forms and `P_R`/`P1`/`P2` coefficient lists — those all live correctly in the notes, so they are MATCH, not MISSING. No MISMATCH or MISSING-DELIVERABLE.

INTERNAL (scaffolding, no finding expected in prose): the three numeric sample residuals `225/8869`, `81/1736`, `91/21935`, `38617837960/99381932001`, `759648230/1473329763`, `5842146019415/70196178995856` (sign-test probes); the `(1−r²)` / `(1−r)` polynomial-division quotients/remainders; all `PASS`/`= 0` residual flags.

## Self-test notes

- **Variable independence:** Confirmed `G_tr` and `F_tr` both genuinely depend on `R` (l.43-48), so `sp.diff(·,R)` / `D[·,r]` are nonzero — the derivative `expect_zero` checks compare a real derivative against a hand-typed form, not 0 vs 0.
- **Symmetry/parity:** No unbounded-domain integrals; positivity is decided by polynomial-coefficient sign (SymPy) and `Reduce[ForAll…]` over the bounded open domain (Mathematica), both valid.
- **Trivial-case / boundary:** Verified the `R=1` differences vanish (branches coincide) and `R=0` hits the strong-split endpoints, consistent with the bound chain; the three concrete samples give strictly positive literals as printed.
- **Symbol domains:** `xi,δ,R` declared positive/real with `0<xi<1` (MMA) — matches the physical `0<xi<1, δ>0, 0<R<1` setup; the `1/(1−xi)` and `(1−xi)` denominators are safe.
- No directive written: the sole finding is informational `stale_output` (label-only, refreshed by the orchestrator re-run; residual docstring/footer labels deferred to the gated numbering pass).
