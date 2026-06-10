---
unit_id: 246
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-10T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 246 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_246.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (row 90 + result-anchor row 9 + checklist row 117)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_mathematica_audit.txt`

## What the paper claims

The stage compiles the compensated source lane declared at Stage 243. It claims: (i) the mean-preserving two-mode family `σ_{a,b}(x)=1+a cos(πx)+b cos(2πx)` integrates to 1 on `[0,1]`; (ii) the exact piecewise minimum `σ_min(a,b)` (boundary case `1+b-|a|` when `b≤0` or `|a|≥4b`; vertex case `1-b-a²/(8b)` otherwise) and thus the sign-change test `σ_min<0`; (iii) the two carried source moments `g(a,b)=(2/π)(1+a/3-b/15)` and `S(a,b)=(2tanh(π/2)/π)(1+a/5+b/17)`; (iv) the normalized two-moment matrix with `det M_src=14/425` and its inverse `a=(85/42)S̃+(25/14)g̃`, `b=(425/42)S̃-(85/14)g̃`; (v) the loading ratio `R(a,b)=(g-r_F1)²/(1+r_F1²)` with `r_F1=√((12/π²)(37/20)²-1)`; (vi) the quarter-ratio theorem `R(g_±)=1/4`; (vii) compensation line `b=5a+15-(15π/2)g_c`; (viii) transported closure `s(r)=r_σ²/(r²+r_σ²)`, transported moments, and sign-change threshold `r<r_σ√(a₀-b₀-1)`; (ix) Session-I readback `g≈0.828236674`, `S≈0.675847711`, `R≈0.216770372`, `σ_min≈-0.089795454`. The card lists `\stagefield{Verification}{...Mathematica audit: none yet.}` — but a `.wl` now exists (see F1).

## What the script claims to verify

The SymPy script verifies all nine deliverables: mean preservation, the quadratic rewrite and vertex/boundary candidates, the piecewise `σ_min` checked against an *independently* computed `Min` of the boundary+in-range-vertex candidates at three representative `(a,b)` points, the two closed-form moments via direct `integrate`, `det M_src=14/425` and the inverse via `M_src·[a_inv,b_inv]-[g̃,S̃]=0`, the quarter-ratio at both `g_±`, the compensation line via `solve`, the transported moments/threshold (with an orientation self-test that the `b₀<0` orientation lands on the boundary branch of the piecewise), and the numeric Session-I readback. The `.wl` mirrors the same deliverables (M1–M9) via independent routes (`MinValue`, `Inverse`, `Reduce`).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| mean preservation =1 | py L53 / wl M1 | match |
| σ_min piecewise + sign-change | py L64-96 / wl M4 (MinValue) | match |
| g, S moments | py L105-115 / wl M2,M3 | match |
| det M_src=14/425, inverse | py L123-139 / wl M5 | match |
| R(a,b), r_F1 | py L145-146 / wl M6 | match |
| quarter-ratio R(g_±)=1/4 | py L167-168 / wl M6 | match |
| compensation line | py L169 / wl M7 | match |
| transported s(r), g(r), S(r), R(r) | py L186-197 | match (sympy only; wl recomputes numerically in M9) |
| sign-change threshold | py L204,209,216-220 / wl M8 | match |
| Session-I readback numerics | py L254-259 / wl M9 | match |
| card "Mathematica audit: none yet" | `.wl` now present | mismatch (F1, prose) |

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `sigma_mean == 1` | (i) mean preservation | yes |
| A2 | sympy | 56 | `simplify(sigma - quadratic in cos) == 0` | (i) quadratic rewrite (real subst) | yes |
| A3 | sympy | 57 | `vertex - (1-b-a²/8b)==0` | (ii) vertex | yes |
| A4 | sympy | 95 | `sigma_min_test - sigma_min_true==0` (Min of candidates) | (ii) piecewise min, independent | yes |
| A5 | sympy | 96 | `sigma_min_test - sigma_min_expected==0` | (ii) — round-trip vs piecewise logic | no (redundant, see F-note) |
| A6 | sympy | 114-115 | `g_sigma-g_exp==0`, `S_sigma-S_exp==0` (from integrate) | (iii) | yes |
| A7 | sympy | 138-139 | `det==14/425`, `M·inv-vec==0` | (iv) | yes |
| A8 | sympy | 167-169 | `R_±-1/4==0`, comp-line | (vi),(vii) | yes |
| A9 | sympy | 209,220 | transported min + orientation self-test | (viii) | yes |
| A10 | sympy | 254-259 | numeric readback tolerances | (ix) | yes |
| M1 | math | 72 | `expectZero[meanMoment-1]` | (i) | yes |
| M2-3 | math | 73-80 | `expectTrue[moment-expected==0]` | (iii) | yes |
| M4 | math | 90-97 | `MinValue` vs target (independent route) | (ii) | yes |
| M5 | math | 111-115 | `Det-14/425`, `Inverse.{gt,St}==expected` | (iv) | yes |
| M6-7 | math | 123-127 | quarter-ratio + comp-line | (vi),(vii) | yes |
| M8 | math | 148-156 | boundary-min + `Reduce` threshold | (viii) | yes |
| M9 | math | 171-176 | numeric readback | (ix) | yes |

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim (stale Verification field)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_246.tex:4` quote: "Mathematica audit: none yet."
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_mathematica_audit.wl` (full file present, M1–M9, output PASS)

**What's wrong:**
The card's `\stagefield{Verification}` line states the Mathematica audit does not yet exist. Pass-1 added a full, independent `.wl` second engine that runs and passes (output file present, dated 2026-06-03). The prose card is now stale relative to the script set.

**Why this matters:**
Card under-reports the verification coverage; a reader auditing dual-engine status would conclude this stage is single-engine when it is not. Pure prose lag, no math impact.

**Required change:** None for Codex. Paper-side prose — route to user. The Verification line should be updated to cite the `.wl` audit; direction is the user's call (Codex never edits paper/).

**Verification:** Card line 4 updated to reference the Mathematica audit file once the user authorizes.

## Independent-derivation check (Mathematica)

Genuinely independent, not a transliteration. Three corroborating contrasts:
1. **Minimum** — `.py` derives `σ_min` analytically (vertex `y_*=-a/4b`, boundary candidates) and validates via `Min[...]`; the `.wl` uses `MinValue[{source /. rules, 0<=x<=1}, x]` — a black-box optimizer over the original `cos`-form `source`, never touching the quadratic/vertex algebra (wl L93).
2. **Inverse map** — `.py` retypes the inverse coefficients and checks `M·inv-vec==0` (round-trip on the *forward* matrix); the `.wl` computes `Inverse[Msrc].{gt,St}` directly and compares to the expected coefficients (wl L102-115) — different operation.
3. **Threshold** — `.py` proves the orientation via `piecewise_fold`+`simplify`; the `.wl` uses `Reduce[orientation && a0-b0>1 && transportMin<0, r, Reals]` and checks equivalence to `0<r<rsig√(a0-b0-1)` (wl L139-156). Distinct machinery. No transliteration finding.

## Engine cross-check

Engines agree. Closed forms: both give `g=(2/π)(1+a/3-b/15)` and `S=(2tanh(π/2)/π)(1+a/5+b/17)`, `det=14/425`, quarter-ratio `=1/4`, threshold `r_σ√(a₀-b₀-1)`. Numerics agree to reported precision: g delta `4.08e-9`, S `1.47e-9`, R `1.56e-9`, σ_min `3.89e-9` — all within the `5e-9` tolerance and matching the `.py` readback (`g=0.8282366740792915`, etc.). No disagreement.

## Verdict justification

`findings` (one low-severity prose-lag `paper_misalignment`). The math holds against the paper on all nine deliverables in both engines. Attacks tried and failed: (a) the pass-1 tautology — the round-trip `sigma_min_expected` (A5) is still present but is now *backed* by the independent `sigma_min_true=Min[...]` check (A4) computed straight from `sigma_y` without the piecewise branch logic, so a wrong piecewise branch would be caught; the fix holds (A5 alone is redundant, not the load-bearing check, so no separate finding). (b) Variable-independence/self-test trap in §8 — the orientation self-test substitutes `a0p>0, b0n<0` into the *actual* piecewise and folds it; it genuinely selects the boundary branch (not a hand-asserted form), no zero-derivative trap. (c) Numeric readback uses real session targets `0.82823667` etc. with `5e-9` tolerances, not retyped script outputs. (d) `r_F1` is derived symbolically in both engines, not hardcoded. I read the paper card, notes, and appendix; the script's verified claims match the paper, with the sole exception of the stale "Mathematica audit: none yet" line.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 17 deliverable values checked, 0 misaligned (1 prose-only stale-field flagged separately as F1).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `∫σ=1` | py out L6 / wl out L9 | tex L161-163 / md L160-163 | MATCH |
| `g=(2/π)(1+a/3-b/15)` | py out L30 / wl L69 | tex eq L108 / md L233-237 | MATCH |
| `S=(2tanh(π/2)/π)(1+a/5+b/17)` | py out L31 / wl L70 | tex eq L110-112 / md L252-257 | MATCH |
| `det M_src=14/425` | py out L40 / wl L111 | tex L71 / md L294-296 | MATCH |
| inverse `a=85/42 S̃+25/14 g̃`, `b=425/42 S̃-85/14 g̃` | py out L41 / wl L106-109 | tex L76-78 / md L302-307 | MATCH |
| `r_F1≈1.777993535474978` | py out L46 | md L484 (`≈1.77799353547498`) | MATCH |
| quarter-ratio `R(g_±)=1/4` | py out L50-51 / wl L123-124 | tex L382 / md L381-383 | MATCH |
| compensation line `b=5a+15-15π g_c/2` | py out L52 / wl L126-127 | tex (line form) / md L348-353 | MATCH |
| `σ_min(r)=1-(a₀-b₀)s(r)` | py out L70 | tex L117-119 / md L451-453 | MATCH |
| threshold `r_σ√(a₀-b₀-1)` | py out L71 / wl L169 | tex L119 / md L458-463 | MATCH |
| `s(r_eval)≈0.38921266` | py out L76 | md L489-492 | MATCH |
| `a(r_eval)≈0.85626786` | py out L77 | md L497 | MATCH |
| `b(r_eval)≈-0.23352760` | py out L78 | md L498-499 | MATCH |
| `g(r_eval)≈0.82823667` | py out L79 / wl out L55-56 | tex (Output packet) / md L505 | MATCH |
| `S(r_eval)≈0.67584771` | py out L80 / wl out L57-58 | md L511 | MATCH |
| `R(r_eval)≈0.21677037` | py out L81 / wl out L59-60 | md L517 | MATCH |
| `σ_min(r_eval)≈-0.08979545` | py out L82 / wl out L61-62 | md L523 | MATCH |
| `g(r=0)≈1.12893906>1` | py out L84 | md L539-544 | MATCH |

Session-I inputs cross-checked: `a₀=2.2=11/5`, `b₀=-0.6=-3/5`, `r_σ=0.8=4/5`, `r_eval=1.00217028=25054257/25000000` — both engines and notes agree.

INTERNAL (no finding): pass/FAIL flags, residual deltas, tolerances (`5e-9`, `10^-20`), `y_*=-a/4b`, boundary/vertex candidate intermediates, `Pi_sigma` packet print (intermediate display), `g̃`,`S̃` symbol scaffolding.

## Self-test notes

Checked: (1) variable independence — §8 orientation self-test substitutes signed surrogates into the real piecewise and folds; no `diff` against an absent variable; no zero-derivative trap. (2) Parity/symmetry — moment integrals are over `[0,1]` (not symmetric/unbounded); not applicable, integrals computed directly. (3) Trivial-case — `a=b=0` gives `σ=1`, `σ_min=1`, `R=(g-r_F1)²/(1+r_F1²)>0`; consistent. (4) Pass-1 tautology re-check — round-trip `sigma_min_expected` (A5) is redundant but the independent `Min` check (A4) is load-bearing and non-tautological; no surviving tautology/round-trip that controls the verdict. (5) Paper round-trip — only finding is prose-lag F1; no fix prescribed that could introduce a new misalignment.
