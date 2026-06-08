---
unit_id: 155
batch: IV.6
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage155_frozen_traction_fixedpoint.md]
  paper_appendix: present
---

# Audit unit 155 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_155.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage155_frozen_traction_fixedpoint.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (frozen-traction values block at lines 944–952; row at line 1344 `\input{stages/stage_155}`)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage155_frozen_traction_fixedpoint_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage155_frozen_traction_fixedpoint_mathematica_audit.txt`

## What the paper claims

The stage card's bottom-line claim (quoted block, line 16) is:
> "At old traction, fixed point survives but moves off exact compensation: \(\mathcal R_{\rm fp}\approx0.282714\)."

The notes elaborate the full deliverable set: solve the exact nonlinear fixed-point map
\(\Sigma\mapsto e^{-\Phi_{\Sigma_0^*}[\Sigma]}/\int_0^1 e^{-\Phi}\) at frozen canonical traction
\(\Sigma_0^*\approx1.80594111095636\), and report the converged moments
\(\mathfrak g_{\rm fp}\approx0.693352419668063\),
\(\mathcal S_{\rm fp}\approx0.6216013167514007\),
\(\mathcal R_{\rm fp}\approx0.2827139049082381\),
\(\Pi_{\rm fp}\approx1.4885734438300713\); the broadening
\(\delta\mathfrak g_{\rm fp}\approx-0.0646826592766000<0\); the drift
\(\delta\mathcal R_{\rm fp}=\mathcal R_{\rm fp}-1/4\approx0.0327139049082381\); and the
verification that this drift matches the exact first-order transport-law prediction
\(\delta\mathcal R_{\rm pred}=-\delta\mathfrak g/\sqrt{1+\mathfrak r_{F1}^2}+(\delta\mathfrak g)^2/(1+\mathfrak r_{F1}^2)\).
The card bottom-lines on the R value; the notes elevate the moment quartet and the
transport-law identity to load-bearing deliverables. The appendix block
(`eq:app-part04-frozen-traction-values`, lines 947–951) independently carries
\(\mathfrak g_{\rm fp}\), \(\mathcal S_{\rm fp}\), \(\mathcal R_{\rm fp}\).

## What the script claims to verify

Both scripts solve the same exact nonlinear Picard fixed-point iteration (grid `N=2401`,
trapezoidal weights, kernels `c=cos(πx/2)`, `K_q=cosh(κ(1−x))/cosh(κ)` with `κ=π/2`, operators
`T_s`, `T_q`), seeded from `normalize(Π_* e^{-Π_* x})`, to tolerance `1e-13`. The SymPy script
(docstring lines 5–11) reports `g_fp, S_fp, R_fp, Pi_fp` and now (lines 94–97) asserts all four
against the paper-stated values with `1e-12` tolerance, plus the transport-law consistency check
(lines 105–106, `raise AssertionError` if `|pred−direct|>1e-10`). The Mathematica script asserts
the same five quantities via `expectApprox` (lines 101–105). Both then declare the qualitative
conclusion: at frozen traction the fixed point stays close in `Π` but drifts to `R>1/4`.

## Paper ↔ script cross-check

| Paper-side deliverable | SymPy coverage | Mathematica coverage |
|---|---|---|
| \(\mathcal R_{\rm fp}\approx0.2827139049082381\) (card line 16 + notes §2 + appendix line 951) | match — `assert` line 96 | match — `expectApprox` line 103 |
| \(\mathfrak g_{\rm fp}\approx0.693352419668063\) (notes §1 + appendix line 947) | match — `assert` line 94 | match — `expectApprox` line 101 |
| \(\mathcal S_{\rm fp}\approx0.6216013167514007\) (notes §1 + appendix line 949) | match — `assert` line 95 | match — `expectApprox` line 102 |
| \(\Pi_{\rm fp}\approx1.4885734438300713\) (notes §2) | match — `assert` line 97 | match — `expectApprox` line 104 |
| Transport-law identity \(\delta\mathcal R_{\rm pred}=\delta\mathcal R_{\rm direct}\) (notes §3) | match — `AssertionError` line 105–106 | match — `expectApprox` line 105 |

Dominant pattern: every paper-side deliverable now has an anchored script-side check on BOTH
engines. `paper_alignment: aligned`. (This is materially better than the first-pass snapshot,
whose two findings — SymPy lacking the four moment asserts, and a "STAGE 138" banner — have both
been remediated: the four `assert` lines are present at 94–97 and both banners read
"STAGE 155".)

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 94 | `assert abs(g_fp - 0.693352419668063) < 1e-12` | g_fp (notes §1) | yes |
| A2 | sympy | 95 | `assert abs(S_fp - 0.6216013167514007) < 1e-12` | S_fp (notes §1) | yes |
| A3 | sympy | 96 | `assert abs(R_fp - 0.2827139049082381) < 1e-12` | R_fp (card bottom line) | yes |
| A4 | sympy | 97 | `assert abs(Pi_fp - 1.4885734438300713) < 1e-12` | Pi_fp (notes §2) | yes |
| A5 | sympy | 105–106 | `if abs(pred-direct)>1e-10: raise` | transport-law identity (notes §3) | yes — ties g_star, rF1, g_fp, R_fp |
| A6 | mathematica | 101 | `expectApprox["g_fp...", gFp, 0.693352419668063, 10^-12]` | g_fp | yes |
| A7 | mathematica | 102 | `expectApprox["S_fp...", sFp, 0.6216013167514007, 10^-12]` | S_fp | yes |
| A8 | mathematica | 103 | `expectApprox["R_fp...", rFp, 0.2827139049082381, 10^-12]` | R_fp | yes |
| A9 | mathematica | 104 | `expectApprox["Pi_fp...", piFp, 1.4885734438300713, 10^-12]` | Pi_fp | yes |
| A10 | mathematica | 105 | `expectApprox["transport law consistency", deltaRPred, deltaRDirect, 10^-10]` | transport-law identity | yes |

The transport-law check (A5/A10) is non-tautological despite both `g_fp` and `R_fp` coming from
the same solve: it requires `g_* = r − √(1+r²)/2` (the compensation relation tying `g_star` to
`rF1`) for the linear coefficient `−1/√(1+r²)` to be correct. I verified the pinned constants
satisfy `r − √(1+r²)/2 = 0.7580350789446639` vs `g_star = 0.758035078944663` (|Δ|≈9e-16), so the
identity is genuinely anchored to the input constants, not vacuous.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` shares the SymPy fixed-point skeleton: identical grid (`Subdivide[0.0,1.0,n-1]`,
`n=2401`), identical trapezoidal weights, identical kernels (`cGrid=Cos[(Pi*xGrid)/2]`,
`kqGrid=Cosh[kappa*(1-xGrid)]/Cosh[kappa]`), identical `tsOperator`/`tqOperator` `Accumulate`
choreography mirroring `Ts`/`Tq`, identical seed `normalize[piStar*Exp[-piStar*xGrid]]`, identical
tolerance `1e-13` and `maxIt=400`, and Mathematica function names that are direct renderings of
the SymPy ones (`tsOperator↔Ts`, `tqOperator↔Tq`, `gMoment↔g`, `rMoment↔R`, `nextSigma↔next_sigma`,
`solveFixedPoint↔solve_fixed_point`).

This is structurally a port in form. I am NOT raising `mathematica_transliteration`, consistent
with the established pass-2 treatment of the directly-analogous numerical-fixed-point stage 156
(and stages 150/158/159), because:

1. The shared skeleton is **forced by the physics** — the stage's content is "iterate the
   nonlinear operator `Σ ↦ e^{-Φ}/∫e^{-Φ}` to convergence." There is no second algorithmic
   pathway to a numerical fixed point; both engines must implement the same Picard iteration.
   There is no answer-baked substitution (no hand-typed `*_expected` closed form fed back as the
   target) that would betray a transliteration — every reported moment emerges from the converged
   `sig`, not from a copied literal. This is the precise distinction from stage 154/157, which
   WERE flagged because they re-typed hand-derived symbolic `*_expected` expressions.
2. The `.wl` carries substantive independent content the `.py` lacks: the overflow-stability shift
   `phShift = ph - Min[ph]` before exponentiation (line 64), which cancels under `normalize` and
   so is a genuinely different numerical implementation of the same map.
3. The two engines land at the same fixed point with `|Δ| ≲ 2.2e-16` on every moment — the
   correct agreement level for two independent discretizations of the same integral equation, and
   meaningful redundancy (a bug in either grid/weight/kernel/iteration would surface as
   disagreement).

Borderline-port, recorded, but treated as legitimate dual-engine cross-confirmation rather than a
defect — the same disposition pass-2 reached for stage 156.

## Engine cross-check

Both engines converge in 16 iterations to residual 6.084e-14. Reported moments:

| Quantity | SymPy | Mathematica | \|Δ\| |
|---|---|---|---|
| g_fp | 0.693352419668063 | 0.6933524196680632 | ≲1.1e-16 |
| S_fp | 0.6216013167514007 | 0.621601316751401 | ≲2.2e-16 |
| R_fp | 0.2827139049082381 | 0.2827139049082381 | 0 |
| Pi_fp | 1.4885734438300713 | 1.488573443830071 | ≲2.2e-16 |
| ΔR_direct | 0.03271390490823811 | 0.03271390490823811 | ≲1.7e-17 |
| ΔR_transport | 0.032713904908237662 | 0.032713904908237605 | ≲5.7e-17 |

The Mathematica transcript also reports all five `expectApprox` diffs at ≤2.22e-16 with PASS. No
`engine_disagreement`.

## Verdict justification

`clean`. The math holds under attack: both engines independently converge to the same fixed-point
moments at 1e-16 agreement; both anchor the moment quartet `g_fp, S_fp, R_fp, Pi_fp` to the
paper-stated values; both confirm the transport-law identity. R_fp = 0.2827139… matches the card's
bottom-line `\(\mathcal R_{\rm fp}\approx0.282714\)` and the notes' precise 0.2827139049082381, and
the appendix block independently carries `g_fp`/`S_fp`/`R_fp`. Attacks tried that failed:
(i) is the transport-law check tautological? — no; it requires `g_* = r − √(1+r²)/2`, which I
verified the pinned constants satisfy (|Δ|≈9e-16), so it ties together `g_star`, `rF1`, and the
converged solve. (ii) Is `rF1=1.77799353547498` stale / does the canonical Family-1 radius use
the bad `168π²`? — no; it matches the canonical closed form `√(4107−100π²)/(10π)` to 2e-15 (same
literal carried consistently across stages 154–159), no `168π²`/`168%` artifact present.
(iii) Symbol-assumption / parity trap? — the iteration is purely numerical on `[0,1]` with explicit
weights; no `simplify`-under-assumptions and no symmetric-domain parity claim to break.
(iv) Is the `.wl` a transliteration that would copy a SymPy bug? — the shared skeleton is forced by
the integral equation, contains no answer-baked substitution, and the `.wl` differs in the overflow
guard; both engines independently discretize and agree at FP precision. I read the card, notes, and
appendix block; the script's verified claims match all of them. The first-pass two findings are
already remediated. Outputs are fresh (sympy .txt 11:30 > .py 10:02; mathematica .txt 11:31 >
.wl 23:11 prior day). One harmless observation, not raised as a finding: `S_star = 0.658075937605429`
(`.py` line 27) is dead/unused scaffolding (never referenced; absent from the `.wl`); it corrupts no
check.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit (source: script text +
committed `.txt` transcripts; nothing executed):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `g_fp ≈ 0.693352419668063` | py:89 / wl:94; sympy.txt:7, math.txt:8 | notes §1 line 41; appendix line 947 (`0.693352419668063`) | MATCH |
| `S_fp ≈ 0.6216013167514007` | py:90 / wl:95; sympy.txt:8, math.txt:9 | notes §1 line 44 (`0.6216013167514007`); appendix line 949 (`0.621601316751401`, rounded) | MATCH |
| `R_fp ≈ 0.2827139049082381` | py:91 / wl:96; sympy.txt:9, math.txt:10 | card line 16 (`0.282714`); notes §2 line 71; appendix line 951 (`0.282713904908238`) | MATCH |
| `Pi_fp ≈ 1.4885734438300713` | py:92 / wl:97; sympy.txt:10, math.txt:11 | notes §2 line 94 (`1.4885734438300713`) | MATCH |
| `δR_fp = R_fp − 1/4 ≈ 0.0327139049082381` | py:102–103 / wl:88,98; sympy.txt:11, math.txt:11 | notes §2 line 83 (`0.0327139049082381`) | MATCH |

Input constants carried as literals (reconciled against upstream/appendix, not stage-155
deliverables but checked for the stale-radius/value traps):
- `rF1 = 1.77799353547498` (py:25, wl:28) → notes (stage 154 line 41) and stages 154–159; canonical
  closed form `√(4107−100π²)/(10π)` (stage 157 py:65). MATCH.
- `g_star = 0.758035078944663` (py:28, wl:29) → appendix line 571. MATCH.
- `Pi_star = 1.50882951349316` (py:29, wl:31) → appendix line 663. MATCH.
- `Sigma0_star = 1.80594111095636` (py:29, wl:31) → notes Goal line 9; appendix line 773. MATCH.

INTERNAL items (verification scaffolding / dead code, no prose expected, no finding): iteration
count `16`; `max residual` `6.08e-14`; `δR from exact transport law ≈ 0.032713904908237605`
(a verification residual that matches the direct value to ~5.7e-17 — the transport-law check, not a
reported deliverable); pass/fail flags and `expectApprox` diffs; tolerances `1e-12`/`1e-10`/`1e-13`;
`S_star = 0.658075937605429` (`.py` line 27, unused dead code).

reconciliation: complete; 5 deliverable values + 4 input constants checked, 0 misaligned

## Self-test notes

Variable-independence trap: the new SymPy asserts (94–97) are plain `abs(value − constant) < tol`
on values from the converged solve; no `diff(EXPR, VAR)` anywhere, so no identically-zero-derivative
trap. Parity/symmetry trap: the iteration is numerical on `[0,1]` with explicit trapezoidal weights —
no unbounded-domain symmetric integral, no vanishing-by-parity claim. Trivial-case / anchoring trap:
I confirmed the transport-law identity is non-vacuous by checking the pinned constants satisfy
`g_* = r − √(1+r²)/2` (|Δ|≈9e-16) and that `rF1` matches `√(4107−100π²)/(10π)` (|Δ|≈2e-15), ruling
out both a tautological check and the stale-`168π²` radius. No directive written (zero findings).
