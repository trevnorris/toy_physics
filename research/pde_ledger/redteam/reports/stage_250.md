---
unit_id: 250
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-03T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 250 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_250.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (row at line 98; main-text edge/speed equations at lines 286-303)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 250 is a "timescale compiler" that turns the Stage-248 relaxed dynamic branch into a one-sided survival inequality "crossing outruns collapse." The card and notes enumerate distinct deliverables: (1) the exact event-chain transit integral `T_traj(E;r_a→r_b)=sqrt(m_s/2)∫dr/sqrt(E−V(r))`; (2) the characteristic crossing-time compiler `t_cross(E)=lambda_eff·sqrt(m_s/(2(E−Vpeak)))`, monotone decreasing in E (notes §1.2, `dt_cross/dE<0`); (3) the collapse-time compiler `t_collapse=sqrt(mu_eta/(g_UV·chi_peak))`; (4) the survival ratio `S(E)=t_cross/t_collapse` and the exact lower-edge `E_edge=Vpeak+(m_s/mu_eta)(g_UV·chi_peak·lambda_eff²/2)`, with the **one-sided / unique-edge structure** `S(E)<1 ⟺ E>E_edge` resting on `S(E)` being monotone decreasing (notes §3.1 "at most one lower survival edge", §3.2, appendix line 303 "this is a one-sided lower edge", card line 67 "t_collapse energy-independent and t_cross(E) decreases monotonically"); (5) the heavy-throat cancellation `∂E_edge/∂m_s=0` when `mu_eta=alpha·m_s`; (6) the speed-space identity `v_safe,min² = v_crit,new² + lambda_eff²·g_UV·chi_peak/mu_eta` (main-text eq `goldilocks-speed`, line 296); (7) the Session-III benchmark numerics (§7-§8) and the trigger-width/steepness sensitivity (`E_edge−Vpeak ∝ lambda_eff²·chi_peak`). The `\stagefield{Verification}` line states "SymPy audit: …; Mathematica audit: none yet."

## What the script claims to verify

The docstring lists 7 items matching deliverables 1-7. Sections 1-3 print the transit/crossing/collapse forms (no asserts). Section 4 derives `E_edge` via `sp.solve(Eq(S**2,1),E)[0]` and asserts it equals the closed form, plus `S_inf==0`. Section 5 substitutes `mu_eta→alpha·m_s` and asserts `dE_edge/dm_s==0` (the cancellation theorem). Section 6 asserts the speed-form crossing time and the speed-space identity `v_safe² == v_crit² + lambda_eff²·g_UV·chi_peak/mu_eta`. Section 7 prints sensitivity derivatives. Sections 8-9 evaluate the Session-III benchmark and raw-width numerics and assert them against literals matching the notes, plus `Smax_num<1.0`, `transit_max_num<tcollapse_num`, `Eedge_raw>Eedge_num`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| 1. Exact transit integral | §1 prints `Ttraj` (no assert) | partial (print-only; definitional, acceptable) |
| 2. Crossing-time compiler + monotonicity (§1.2) | §2 prints `tcross`, prints `dtcross_dE` (no assert) | partial (form printed; monotonicity printed not asserted) |
| 3. Collapse-time compiler | §3 prints `tcollapse` (no assert) | partial (definitional, acceptable) |
| 4a. Lower-edge `E_edge` | line 92 `assert E_edge − closedform == 0` (solved) | match |
| 4b. One-sided/unique edge ⟺ global monotone `dS/dE<0` (§3.1, §3.2, appendix L303) | only `S_inf==0` (line 93) + single-point `Smax_num<1.0` (line 190); `dS_dE` printed not asserted | **mismatch/missing** — see F2 |
| 5. Heavy-throat cancellation `∂E_edge/∂m_s=0` | line 110 `assert dEedge_dm == 0` | match |
| 6. Speed-space identity | lines 132-133 asserts | match |
| 7a. Benchmark numerics | lines 183-191, 208-210 asserts vs notes values | match |
| 7b. Sensitivity `E_edge−Vpeak ∝ lambda_eff²·chi_peak` | §7 prints `dE_dlam/dE_dchi/dE_dmu` (no assert); `Eedge_raw>Eedge_num` (line 210) | partial (printed; raw-width inequality asserted) |
| Mathematica second engine | absent | **missing** — see F1 |

`paper_alignment: aligned` — every load-bearing numeric literal in the script (`9.430664762121758`, `5.322659467573074`, `0.07469790710905001`, `1.259480370925967`, `6.171635157122822`, `27.532730948999998`, `chi_raw=50.74399964`, `Emax=80.93332737`) matches the notes §7-§8 boxed benchmark values and the appendix. No `paper_misalignment`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 92 | `simplify(E_edge − closedform)==0` (E_edge from solve) | claim 4a (lower edge) | yes |
| A2 | sympy | 93 | `S_inf==0` | claim 4b (partial: limit only) | partial |
| A3 | sympy | 110 | `dEedge_dm==0` after `mu_eta→alpha·m_s` | claim 5 (cancellation) | yes |
| A4 | sympy | 132 | `simplify(tcross_v − lam_eff/sqrt(v0²−vcrit²))==0` | claim 6 (speed crossing time) | yes |
| A5 | sympy | 133 | `simplify(v_safe_sq − v_safe_expected)==0` | claim 6 (speed-space identity) | yes |
| A6 | sympy | 183-189 | `abs(num − literal)<tol` (7 asserts) | claim 7a benchmark | yes (paper-anchored regression) |
| A7 | sympy | 190 | `Smax_num < 1.0` (single point) | claim 4b (window safe at one E) | partial |
| A8 | sympy | 191 | `transit_max_num < tcollapse_num` | claim 7 (dynamic cross-check) | yes |
| A9 | sympy | 208-209 | `abs(num − literal)<tol` (raw width) | claim 7b sensitivity | yes |
| A10 | sympy | 210 | `Eedge_raw > Eedge_num` | claim 7b sensitivity direction | yes |
| — | mathematica | — | none | all | missing |

17 asserts total; comfortably exceeds the docstring's 7 stated claims by count. The gap is qualitative, not count-based.

## Findings

### F1 — missing_verification_script (missing_mathematica)

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/` (no `moving_throat_pde_stage250_*_mathematica_audit.wl` exists; MANIFEST `mathematica.path: null`)

**What's wrong:**
The unit is `is_checkpoint: false` but `is_status_only_candidate: false`, so the dual-engine rule applies: a `.wl` is required wherever Mathematica CAN independently verify the stage. Stage 250 is entirely algebraic/symbolic — solving `S(E)²=1` for the edge, simplifying the survival ratio, taking derivatives for the cancellation and sensitivity theorems, verifying the speed-space identity, and evaluating closed-form numerics. All of these are squarely within native Mathematica primitives (`Solve`, `Reduce`, `Simplify`/`FullSimplify`, `D`, `Resolve[ForAll]`, `Limit`). The upstream Stage 248 (`...stage248_dynamic_event_chain..._mathematica_audit.wl`) and downstream Stage 253 both have `.wl` audits, confirming this family of stages is Mathematica-verifiable. There is no genuine single-engine impossibility here.

**Why this matters:**
A single-engine stage cannot catch transliteration-class or engine-specific simplification errors, and the dual-engine policy is a user-level requirement for any stage Mathematica can verify. The most load-bearing claim — the one-sided window via global monotonicity — is exactly the kind of `ForAll` statement Mathematica's `Resolve`/`Reduce` proves natively and SymPy proves awkwardly.

**Required change:**
Add a NEW independent-route Mathematica audit (not a transliteration of the `.py`) at `mathematica/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_mathematica_audit.wl`. Use a different decomposition: prove the window/edge globally with `Reduce`/`Resolve[ForAll]` rather than evaluating at sample points. See the directive's Claim manifest.

**Verification:**
`redteam exec-mathematica 250` runs and exits 0; the new `.wl` exists with `expectZero`/`expectTrue`-style checks covering M1-M7 below.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.py:83-93` (monotonicity printed not asserted) and `:190` (window strictness sampled at one point)

**What's wrong:**
The paper's headline deliverable is the **one-sided Goldilocks window with a unique lower edge**: notes §3.1 "Because S(E) is monotone decreasing in E, there is at most one lower survival edge"; §3.2 items 1-2; appendix line 303 "Under the declared characteristic closure this is a one-sided lower edge"; card line 67 "t_collapse is energy-independent and t_cross(E) decreases monotonically." This entire structure rests on the **global strict monotonicity** `dS/dE < 0` for all `E > Vpeak` (and equivalently `dt_cross/dE < 0`, notes §1.2). The script computes `dS_dE` (line 83) and `dtcross_dE` (line 57) and merely **prints** them (output lines 36, 22); it never asserts their sign. The only ratio-side asserts are `S_inf==0` (line 93, a limit, not a monotonicity statement) and `Smax_num < 1.0` (line 190) — a strict inequality checked at a **single sampled energy** `E_max,scan`, not globally. So the load-bearing "window is a half-line / edge is unique" claim is exercised only at one point plus a limit, which does not establish monotonicity or uniqueness.

Note also the relevant symbol-domain context: `E` is declared only `real=True` (line 37), with no `E > Vpeak`; a sign assertion on `dS_dE` requires SymPy to know `E - Vpeak > 0`.

**Why this matters:**
A non-monotone `S(E)` could dip below 1, rise back above 1, and create a closed safe island or a finite upper edge — exactly the structure the paper explicitly rules out "within this closure" (card line 67, notes §3.2). Without asserting `dS/dE < 0` globally, the script's checks would still pass even if monotonicity failed, so the most distinctive physics claim of the stage is unverified.

**Required change:**
Add a global strict-monotonicity assertion in section 4 of the `.py`. Because the sign needs `E > Vpeak`, introduce the gap variable `x = E - Vpeak` with `x` positive, rewrite `dS/dE` in terms of `x`, and assert it is strictly negative (e.g. assert the numerator/`x`-dependence forces a negative sign for all positive `x` with positive params). Equivalently, assert `sp.simplify(dS_dE * (E - Vpeak)**sp.Rational(3,2))` is a negative constant in the remaining positive parameters, OR (cleaner) assert that `S` is strictly decreasing by checking the sign of `dS_dE` under positivity assumptions. Also add `assert dtcross_dE` is negative under the same domain (covers notes §1.2). See directive F2.

**Verification:**
New asserts appear in section 4 (and section 2) of the `.py`; the script still exits 0; a `dS/dE < 0` strict-monotonicity check is present and non-vacuous (it must reference `E`/`x` so it is not a vacuous independence test).

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration is moot for now. The directive's prescribed `.wl` MUST use an independent decomposition (global `Reduce`/`Resolve[ForAll]` for the window and monotonicity; `Solve` for the edge) rather than echoing the `.py`'s `subs`/`simplify` choreography, to satisfy the second-engine independence policy.

## Engine cross-check

Not applicable — only one engine present. (Becomes applicable once F1's `.wl` lands; the orchestrator's verifier will compare residuals/edge formula/benchmark numerics.)

## Verdict justification

The SymPy script is paper-aligned on every numeric and the closed-form edge: I attacked the `E_edge` solve (single root after squaring `S²=1`; `[0]` is safe), the heavy-throat `dEedge_dm==0` (genuinely non-vacuous because `mu_eta→alpha·m_s` makes `m_s` cancel, so the derivative tests a real cancellation, not a missing variable), the `subs(E,E_launch)` in section 6 (a no-op since `E_edge` is E-free, but the resulting algebraic identity `2(E_edge−V0)/m_s = v_crit² + lambda_eff²·g_UV·chi_peak/mu_eta` is still independently constructed and correct), and the benchmark literals (they match notes §7-§8 exactly, so they are paper-anchored regression anchors, not self-confirming hardcodes). Those hold up. Two real gaps remain: the stage is single-engine when Mathematica can independently verify it (F1), and the headline one-sided/unique-edge claim is exercised only via a limit plus a single sample point, never via the global strict monotonicity `dS/dE<0` it actually rests on (F2). Hence `verdict: findings`, no stop-cold (neither fix changes a downstream-quoted constant; the edge formula and benchmarks are unchanged).

## Self-test notes

- **Variable-independence trap:** Checked the `dEedge_dm==0` assert (line 110) — NOT vacuous: `E_edge` depends on `m_s` before `mu_eta→alpha·m_s`, and the substitution is precisely what cancels `m_s`, so the zero derivative is the physics, not a missing variable. The F2 monotonicity fix I prescribe references `E`/`x` explicitly (the integrand depends on `E−Vpeak`), so it is not a vacuous-derivative trap.
- **Strict-inequality / window trap:** This is the core of F2 — the load-bearing strict window/uniqueness is currently tested non-strictly-globally (single point `Smax_num<1.0` + limit), so I prescribe a global `dS/dE<0` assertion (SymPy) and `Reduce`/`Resolve[ForAll]` (Mathematica).
- **Trivial-case pre-check:** For `dS/dE<0`, with all params positive and `x=E−Vpeak>0`, `dS_dE = −(positive)/x^{3/2} < 0` — confirmed negative, so the proposed assert is satisfiable and meaningful, not trivially true.
- **Output freshness:** output mtime (2026-05-11 12:52) is newer than script mtime (11:58); not stale.
