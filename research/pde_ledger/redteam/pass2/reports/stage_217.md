---
unit_id: 217
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00-06:00
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
  notes_stage_files: [moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction.md]
  paper_appendix: present
---

# Audit unit 217 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_217.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.txt`

## What the paper claims

Stage 217 reduces the unique support-cardinality-5 mixed simplex's four-parameter interior optimization to a finite algebraic candidate set. The card's `\stagefield{Output}` reads: "Finite support-five interior candidate set, preferred budget \(324\), fallback budget \(1500\), and support-five interior improvement/non-improvement theorem." The `\stagefield{Derivation ledger}` states the derivation "obtains degree pattern \((3,3,3,3,2)\), proves the preferred \(162\)-candidate bound per envelope, and gives a square-root-free fallback bound." Distinct deliverables (from the notes §1–§8): (1) the unique positive spherical five-simplex with five codimension-one quadruple faces; (2) the exact stationary-numerator identity ∂Φ=0 ⇔ N=0; (3) the lifted polynomial system with degree pattern (3,3,3,3,2) and lifted Bézout bound 162 = 3⁴·2; (4) the projected square-root-free elimination with degrees (5,5,5,6) and bound 750; (5) two special reductions (diagonal-isotropic gradient-optimal ray, fivefold-symmetric equal-mix barycenter); (6) the local improvement/non-improvement theorem. The appendix (part06) carries the budget arithmetic: 162 per envelope → 324 across {lo,hi}; fallback 750 per envelope → 1500; and the Stage-218-level totals 1140+324=1464 (preferred) and 1140+1500=2640 (fallback).

## What the script claims to verify

The SymPy script verifies: (I) exactly five codimension-one faces of the 5-simplex; (II) the four exact stationary-numerator identities `2√Δ·S^(3/2)·∂_xΦ − (2M_x√Δ + L_x) = 0` for x∈{r,s,t,u}; (III) the lifted polynomials F_r,F_s,F_t,F_u,F_Δ have total degrees (3,3,3,3,2) and product 162; (IV) the projected polynomials C_rs,C_rt,C_ru,S_r have degrees (5,5,5,6) and product 750; (V) the diagonal-isotropic reduction L_x = 2·K_lin·M_x and that the gradient-optimal ratios annihilate both M_x and L_x; (VI) the fivefold-symmetric equal-mix barycenter (r=s=t=u=1) annihilates both M_x and L_x. The Mathematica script verifies the same six deliverables (sections M1–M6), with the face/degree facts in M2/M3/M4 and the stationary-numerator identity in M1.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Five codim-1 quadruple faces of Δ₅⁺ | py I (len(faces)==5); wl M2 (count 5, each card. 4) | match |
| Stationary-numerator identity ∂Φ=0 ⇔ N=0 | py II (4 expect_zero); wl M1 (4 expectZero) | match |
| Degree pattern (3,3,3,3,2), lifted Bézout 162 | py III (deg tuple + bezout_lift==162); wl M3 (tuple + product==3^4·2) | match |
| Projected degrees (5,5,5,6), bound 750 | py IV (tuple + bezout_proj==750); wl M4 | match |
| Diagonal-isotropic gradient-optimal ray exact | py V; wl M5 | match |
| Fivefold-symmetric equal-mix barycenter exact | py VI; wl M6 | match |
| Preferred budget 324 = 2×162 (across envelopes) | not script-emitted; prose-level arithmetic over the per-envelope 162 | match (downstream arithmetic) |
| Fallback budget 1500 = 2×750 | not script-emitted; prose-level arithmetic over the per-envelope 750 | match (downstream arithmetic) |
| Stage-218 totals 1140+324=1464, 1140+1500=2640 | not script-emitted (these belong to Stage 218) | n/a (out-of-stage) |
| Improvement / non-improvement theorem (interval comparison) | not script-tested (it is a logical implication on imported brackets β₅) | partial — symbolic, no numeric data exists yet |

`paper_alignment: aligned`. The one finding is a `stale_output` on the SymPy `.txt`, which is not a paper-misalignment.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 54-55 | `if len(faces) != 5: raise` | deliverable 1 (face count) | yes |
| A2 | sympy | 92-107 | 4× `expect_zero(2√Δ S^(3/2)∂Φ − N)` | deliverable 2 (stationary identity) | yes (non-trivial; involves √Δ and S^(3/2)) |
| A3 | sympy | 133-134 | `if (deg...) != (3,3,3,3,2): raise` | deliverable 3 (degree pattern) | yes (degrees computed from constructed polys) |
| A4 | sympy | 138-139 | `if bezout_lift != 162: raise` | deliverable 3 (Bézout 162) | yes (product of computed degrees, not literal) |
| A5 | sympy | 162-168 | `(5,5,5,6)` + `bezout_proj != 750` | deliverable 4 (projected 750) | yes |
| A6 | sympy | 211-223 | 12× `expect_zero` diag reduction | deliverable 5 | yes |
| A7 | sympy | 256-263 | 8× `expect_zero` equal-mix | deliverable 6 | yes |
| M1 | mathematica | 117-125 | 4× `expectZero[scaledDerivNum − stationaryNum]` | deliverable 2 | yes |
| M2 | mathematica | 134-135 | `expectTrue` count 5 + card. 4 | deliverable 1 | yes |
| M3 | mathematica | 147-148 | `expectTrue` tuple {3,3,3,3,2} + product==3^4·2 | deliverable 3 | yes |
| M4 | mathematica | 162-163 | `expectTrue` tuple {5,5,5,6} + product==5·5·5·6 | deliverable 4 | yes |
| M5 | mathematica | 189-211 | 12× `expectZero` diag reduction | deliverable 5 | yes |
| M6 | mathematica | 247-253 | 8× `expectZero` equal-mix | deliverable 6 | yes |

All assertions are non-tautological and anchored. The stationary-numerator identities (A2/M1) are the load-bearing checks: they verify that the symbolic derivative `∂_xΦ` (computed by the CAS from the actual Φ = (K_lin+√Δ)/√S) reduces to the claimed boxed form `2M_x√Δ + L_x` up to the `S^(3/2)` Jacobian factor. They CAN fail if the boxed M_x/L_x definitions in the notes are wrong, so they genuinely exercise the paper claim.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.txt` (mtime 2026-05-11, epoch 1778525386)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py` (mtime 2026-06-03, epoch 1780523951)

**What's wrong:**
The committed SymPy output `.txt` predates the current script by ~23 days. The captured transcript carries a stale stage label — line 11 prints `STAGE 200 — UNIQUE FIVE-COORDINATE SIMPLEX INTERIOR OPTIMIZER` and the header date is `2026-05-11` — whereas the current `.py` banner (line 35) prints `STAGE 217 — UNIQUE FIVE-COORDINATE SIMPLEX INTERIOR OPTIMIZER`. The numeric results in the transcript (162, 750, all `= 0` residuals) are still correct and consistent with the current script's checks, so this is an informational freshness flag, not a content disagreement. The Mathematica output (mtime 2026-06-02) is newer than its `.wl` (also 2026-06-02) and is fresh; its banner already reads `STAGE 217`.

**Why this matters:**
The committed transcript does not reflect the current banner text, and the stale `STAGE 200` label could be mistaken for a numbering drift in the active script (it is not — the live `.py` already says 217). A fresh capture removes the ambiguity.

**Required change:**
Re-run the SymPy script and re-capture its output to refresh the `.txt` so the banner reads `STAGE 217` and the header date/mtime are current. No script-logic edit is needed.

**Verification:**
After re-run, the `.txt` header date is current and line ~11 reads `STAGE 217 — UNIQUE FIVE-COORDINATE SIMPLEX INTERIOR OPTIMIZER`; all residuals remain `= 0`, Bézout bounds remain 162 and 750, `EXIT_CODE: 0`.

## Independent-derivation check (Mathematica)

Verdict: **INDEPENDENT** (borderline-leaning-independent on the degree/Bézout sections, but the load-bearing M1 identity is genuinely independent).

Single discriminating operation: the construction of the envelope-transport numerators `L_x`. The SymPy script hardcodes the Jacobian term as `2*r*Delta`:
```
Lr = sp.simplify(S * sp.diff(Delta, r) - 2 * r * Delta)        # py:85
```
The Mathematica script instead writes the same quantity through the symbolic derivative `D[S, #]` rather than substituting the closed form `2r`:
```
envelopeNumerators = FullSimplify[S D[delta, #] - D[S, #] delta, ...] & /@ ratioVars   # wl:100-103
```
That is a derive-vs-derive divergence (Mathematica lets the CAS differentiate `S` rather than positing `∂_r S = 2r`), not a line-for-line port.

The load-bearing stationary-numerator identity (M1) is also reached by a different route. SymPy scales by `S^(3/2)` and uses `sp.simplify`:
```
2 * sp.sqrt(Delta) * S**sp.Rational(3,2) * sp.diff(Phi, r) - (2*Mr*sp.sqrt(Delta) + Lr)   # py:94
```
Mathematica forms `2 sqrtDelta S^(3/2) D[phi,#]`, then applies `Cancel[Together[...]]` and prints the explicit cleared symbolic numerator (the long `M1 cleared derivative numerators = {...}` block in the output) before differencing against the stationary numerators with a `ConditionalExpression`-stripping `cleanExpr`. Different simplification machinery, different intermediate representation surfaced.

The M3/M4 degree/Bézout checks ARE the same operation in both engines (compute `total_degree` of the constructed polynomials, multiply). This is the borderline part — but computing a polynomial's total degree is a canonical, engine-agnostic fact (not a posited value), so corroborating it across two CAS implementations is legitimate cross-checking, not transliteration. The .wl additionally prints `3^4*2` and `5*5*5*6` as separate literal evaluations and checks the computed product equals them, which is an independent re-statement rather than an echo of the .py's single `!= 162` guard. Overall the .wl derives the same results through its own algebra; it is not a port.

## Engine cross-check

Both engines agree fully. SymPy output: all four stationary-numerator identities `= 0`; `(deg_Fr..deg_FDelta) = (3,3,3,3,2)`, `Lifted Bezout bound = 162`; `(deg_Crs..deg_Sr) = (5,5,5,6)`, projected bound `750`; all 12 diag-reduction and 8 equal-mix residuals `= 0`; `EXIT_CODE: 0`. Mathematica output: M1 four identities `= 0` (PASS); M3 `{3,3,3,3,2}`, product `162` == `3^4*2`; M4 `{5,5,5,6}`, product `750` == `5*5*5*6`; M5 (12) and M6 (8) all `= 0` (PASS); "All Stage 217 Mathematica audit checks passed." No residual, sign, or factor disagreement.

## Verdict justification

`findings: 1`. The mathematics holds up under attack. The load-bearing 162 (= 3⁴·2) is genuinely DERIVED in both engines — the per-factor degrees (3,3,3,3,2) are computed from the actually-constructed lifted polynomials via `total_degree`, then multiplied; neither engine hardcodes 162 as a posited answer (the literal `162` appears only as the value the computed product is asserted to equal). The four stationary-numerator identities are non-tautological (they confront the CAS-computed `∂_xΦ` against the boxed M_x/L_x forms) and both engines clear them to zero by different routes. The two special reductions (diagonal-isotropic and fivefold-symmetric) substitute concrete Hessian-envelope structures and confirm the gradient-optimal ray and equal-mix barycenter are exact stationary points. The 324/1500/1464/2640 budget figures are prose-level arithmetic combinations (324=2×162, 1500=2×750, and the Stage-218 totals 1140+324, 1140+1500) consistent across card/appendix/notes. The only finding is a low-severity `stale_output` on the SymPy transcript (older mtime + stale `STAGE 200` banner text); its numbers still match the current script. Attacks tried that failed: (a) is 162 hardcoded? — no, it's a product of computed degrees; (b) is the .wl a port? — no, the L_x construction and M1 route differ; (c) any surviving 179/230 anywhere? — none (see RE-CONFIRMATION); (d) do the engines disagree? — no. Paper, notes, and appendix all consistently state 162/324, matching what the scripts derive.

## 217 RE-CONFIRMATION (published value_mismatch carry-over)

The pass-1 user-resolved correction holds. I independently re-confirmed:

(a) **Scripts emit/derive 162** — YES. SymPy `bezout_lift = deg_Fr·deg_Fs·deg_Ft·deg_Fu·deg_FDelta` (py:136) with `if bezout_lift != 162: raise` (py:138), output line 42 `Lifted Bezout bound = 162`. Mathematica `liftedProduct = Times @@ liftedDegrees` (wl:142) with `expectTrue["M3 lifted product equals 3^4*2", liftedProduct == 3^4*2]` (wl:148), output line 33 `M3 lifted Bezout product = 162`, line 34 `M3 literal 3^4*2 evaluates to 162`. The 162 is DERIVED (product of computed degrees), not hardcoded.
  - **324** — NOT script-emitted; it is `2×162` (across {lo,hi} envelopes). Stated in card line 15 (`preferred budget \(324\)`), appendix line 1205 (`2\times162=324`), notes lines 620. Arithmetically forced from the script-derived 162.
  - **1464** — NOT script-emitted (it is a Stage-218 total). Appendix line 1261 `1140+324=1464`; appendix row line 67 (Stage 218) `preferred budget \(1464\)`. Arithmetically forced: 1140 (imported support-≤4 ledger) + 324 = 1464.
  - **2640** — NOT script-emitted (Stage-218 fallback total). Appendix line 1266 `1140+1500=2640`. Arithmetically forced: 1140 + 1500 (= 2×750) = 2640.

(b) **Card + appendix + notes all say 162, no surviving 179/230** — YES. `grep` for `179`/`230` across stage_217.tex, the notes .md, both scripts, both outputs, and part06 appendix returned ZERO matches. 162 present in card (line 13), appendix (lines 65, 1200), notes (lines 31, 409, 506, 616, 634).

**Arithmetic relation forcing 162:** the lifted stationary system has degree pattern (3,3,3,3,2); by Bézout the candidate bound is the product `3·3·3·3·2 = 3⁴·2 = 162` per envelope. Then `324 = 2·162` (lo+hi envelopes); `1464 = 1140 + 324` and `2640 = 1140 + 1500` (= 1140 + 2·750) are the Stage-218 grand totals over the imported 1140-candidate support-≤4 ledger. 162 is the unique value consistent with the (3,3,3,3,2)→324→1464 chain, exactly as the pass-1 user resolution recorded; the wrong pass-1 `179/230` would break `2×=324` and `1140+=1464`.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation (the scripts emit symbolic identities + integer degree/Bézout facts; the only standalone numbers they pin are 162 and 750):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| degree pattern (3,3,3,3,2) | py:133 assert, out L37-41; wl:147, out L32 | card line 13; appendix line 1195; notes lines 30, 404 | MATCH |
| lifted Bézout bound `162` | py:136-138, out L42; wl:142-148, out L33-34 | card line 13; appendix lines 65, 1200; notes lines 31, 409, 506, 616, 634 | MATCH |
| projected degrees (5,5,5,6) | py:162 assert, out L47-50; wl:162, out L43 | appendix line 1207 (implicit via 750); notes line 502 (`5·5·5·6`) | MATCH |
| projected one-chart bound `750` | py:165-167, out L51; wl:157-163, out L44-45 | appendix line 1207; notes lines 32, 502, 506 | MATCH |
| preferred budget `324` (=2×162) | not emitted (prose arithmetic) | card line 15; appendix line 1205; notes line 620 | MATCH (derived from emitted 162) |
| fallback budget `1500` (=2×750) | not emitted (prose arithmetic) | card line 15; appendix line 1207 | MATCH (derived from emitted 750) |
| Stage-218 preferred total `1464` (=1140+324) | not emitted (Stage-218 prose) | appendix lines 67, 1261 | MATCH (out-of-stage; arithmetically consistent) |
| Stage-218 fallback total `2640` (=1140+1500) | not emitted (Stage-218 prose) | appendix line 1266 | MATCH (out-of-stage; arithmetically consistent) |
| 5 codim-1 faces, each cardinality 4 | py:54 assert, out L19-24; wl:134-135, out L24-27 | notes §1 lines 79-91; appendix line 1215 (strata 5+10+10+5=30) | MATCH |
| diagonal-isotropic gradient-optimal ray exact (12 residuals 0) | py:211-223; wl:189-211 | notes §6.1 | MATCH (symbolic) |
| equal-mix barycenter r=s=t=u=1 exact (8 residuals 0) | py:256-263; wl:247-253 | notes §6.2 lines 559-567 | MATCH (symbolic) |

INTERNAL scaffolding (no finding): `S`, `Klin`, `Delta`/`delta`, `Phi`/`phi`, `Mr/Ms/Mt/Mu`, `Lr/Ls/Lt/Lu`, `Fr/Fs/Ft/Fu/FDelta`, `Crs/Crt/Cru/Sr`, `diag_subs`/`diagRules`, `sym_subs`/`symRules`, `grad_ratios`, `bary_subs`, per-residual `= 0` outputs, `EXIT_CODE`/PASS flags, the `M1 cleared derivative numerators` printed intermediate.

reconciliation: complete; 11 deliverable values checked, 0 misaligned.

## Self-test notes

Checked: (1) variable-independence — every `sp.diff(Phi, var)` and `D[phi, var]` uses a `var` that Φ genuinely depends on (r,s,t,u all appear in S, Klin, and Δ), so no derivative is identically zero and the `expect_zero` checks are non-trivial. (2) No unbounded integrals here, so parity is N/A. (3) Trivial-case: at the gradient ratios r=k_c/k_L etc. under diag-isotropic curvature, M_r = S·k_c − r·K_lin with r=k_c/k_L collapses correctly (output confirms 0); at r=s=t=u=1 under fivefold symmetry the M_x are antisymmetric pairs that cancel (output confirms 0). (4) No missing-script finding (both engines present). (5) Paper round-trip: the only fix (stale_output re-run) introduces no new constant — the refreshed transcript will still print 162/750 and zero residuals, matching the card/appendix/notes.
