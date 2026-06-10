---
unit_id: 250
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 250 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_250.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (row 250 at line 98; main goldilocks-edge eqs at lines 288/296)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_mathematica_audit.txt`

## What the paper claims

Stage 250 is a "timescale compiler" that turns the Stage-248 relaxed dynamic branch into a survival inequality "crossing outruns collapse." The card and notes (Section 11 enumerates 8 deliverables) state the stage proves: (1) the exact event-chain transit integral `T_traj = sqrt(m_s/2) ∫ dr/sqrt(E-V)`; (2) the characteristic crossing-time compiler `t_cross = λ_eff sqrt(m_s/(2(E-Vpeak)))`; (3) the unstable-leg growth rate `Γ_coll = sqrt(g_UV χ_peak/μ_η)` and `t_collapse = sqrt(μ_η/(g_UV χ_peak))`; (4) the survival ratio `S(E)` and the lower-edge `E_edge = Vpeak + (m_s/μ_η)(g_UV χ_peak λ_eff²/2)` (card eq:app-part08-stage250-survival-ratio + main goldilocks-edge); (5) heavy-throat cancellation `∂E_edge/∂m_s = 0` when `μ_η=α m_s`; (6) the speed-space identity `v_safe,min² = v_crit,new² + λ_eff² g_UV χ_peak/μ_η`; (7) the full Session-III benchmark numerics; (8) the dynamic transit-time cross-check `transit_max < t_collapse`. The one-sidedness claim (card line 67) rests on `t_collapse` energy-independent AND `t_cross(E)` monotone decreasing. The card `\stagefield{Verification}` line states "Mathematica audit: none yet."

## What the script claims to verify

The SymPy script derives all 8 deliverables symbolically (sections 1–9) plus a sensitivity ledger, and pins the Session-III benchmark + raw-width numerics with `assert abs(...) < tol`. The Mathematica script (M1–M7) re-derives the crossing form and — crucially — establishes the load-bearing global claims by quantifier elimination: M1 `Resolve[ForAll[..., D[t_cross,En] < 0]]` (global monotonicity), M3 `m3Unique = Resolve[ForAll[..., S²==1 ⟹ En==Eedge]]` (unique edge), and M4 `m4Window = Resolve[ForAll[..., Equivalent[S<1, En>Eedge]]]` (global one-sided window). It also derives the edge independently via `Solve`/`Reduce` and checks the same benchmark numerics.

## Paper ↔ script cross-check

| Deliverable | Script check | Status |
|---|---|---|
| T_traj integral | py §1 prints exact form | match |
| t_cross compiler | py §2, wl M1 `expectZero` | match |
| t_collapse / Γ_coll | py §3, wl M2 | match |
| S(E), E_edge | py §4 assert, wl M3 Solve+Reduce | match |
| Monotonicity → one lower edge | wl M1/M3/M4 `Resolve[ForAll]` (GLOBAL) | match |
| Heavy-throat cancellation | py §5 `dEedge_dm==0`, wl M5 | match |
| Speed-space identity | py §6 assert, wl M6 | match |
| Benchmark numerics | py §8, wl M7 | match |
| Raw-width sensitivity | py §9, wl M7 | match |
| transit_max < t_collapse | py §8 `assert transit_max_num < tcollapse_num` | match |
| Card "Mathematica audit: none yet" | a full passing `.wl` exists | mismatch (F1) |

Dominant pattern is match; one paper-text mismatch → `partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 67 | `dtcross_dE_gap == -√2 λ √m/(4 x^{3/2})` | claim 2 | yes |
| A2 | sympy | 70 | `(dtcross_dE_gap < 0) == True` | monotonicity (point in gap-var, but symbolic over x>0) | yes |
| A3 | sympy | 103 | `E_edge - (Vpeak + λ²g χ m/(2μ)) == 0` | claim 4 | yes |
| A4 | sympy | 112 | `(dS_dE_gap < 0) == True` | monotonicity | yes |
| A5 | sympy | 113 | `S_inf == 0` | high-E limit | yes |
| A6 | sympy | 130 | `dEedge_dm == 0` | claim 5 | yes |
| A7 | sympy | 152/153 | `tcross_v` & `v_safe_sq` identities | claim 6 | yes |
| A8 | sympy | 203–211 | benchmark + transit `assert abs<tol` | claim 7/8 | yes |
| A9 | sympy | 228–230 | raw-width numerics | claim 7 | yes |
| M1 | math | 74–93 | `Resolve[ForAll[..., D<0]]` ×2 | monotonicity (GLOBAL) | yes |
| M2 | math | 101 | `expectZero[Γ - √(gχ/μ)]` | claim 3 | yes |
| M3 | math | 115–129 | `Resolve[ForAll[S²==1⟹En==Eedge]]` + Solve | claim 4 (GLOBAL unique) | yes |
| M4 | math | 135–145 | `Resolve[ForAll[Equiv[S<1,En>Eedge]]]` | one-sided window (GLOBAL) | yes |
| M5 | math | 151–155 | `expectZero[D[EedgeHeavy,ms]]` | claim 5 | yes |
| M6 | math | 161–167 | `expectZero[speedResidual]` | claim 6 | yes |
| M7 | math | 209–216 | `expectApprox` benchmark | claim 7 | yes |

## Findings

### F1 — paper_misalignment

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_250.tex:4`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_mathematica_audit.wl` (full file, passing)

**What's wrong:**
The stage card's `\stagefield{Verification}` line (stage_250.tex:4) reads: "SymPy audit: ... Mathematica audit: none yet." But a complete, independent, passing Mathematica audit `.wl` exists for this stage (M1–M7 all `PASS` in the saved output). The card understates the verification coverage.

**Why this matters:**
This is a documentation/coverage statement, not a math error. It does not affect any derived result. Subtype `paper_missing_script_claim` (script does more than the card records). Direction of resolution is the user's call (update the card to cite the `.wl`).

**Required change:**
None for Codex — this is a paper-side text update routed to the user. See directive `## Resolve before fix_loop`.

**Verification:**
After user resolution, card line 4 should cite the `.wl` path instead of "Mathematica audit: none yet."

## Independent-derivation check (Mathematica)

The `.wl` is genuinely independent, NOT a transliteration. Evidence:
1. The `.py` solves the edge via `sp.solve(sp.Eq(S**2, 1), E)[0]` and proves monotonicity only by `(dS_dE_gap < 0) == True` (a sign test on the symbolic gap-variable expression). The `.wl` instead establishes the SAME claims by quantifier elimination it never borrows: `m1Monotone = Resolve[ForAll[{En,Vpeak,ms,lambdaEff}, Implies[..., D[...,En] < 0]], Reals]` (M1, lines 74–82), `edgeReduce = Reduce[Sratio^2==1 && En>Vpeak && ..., En, Reals]` plus `m3Unique = Resolve[ForAll[..., S²==1 ⟹ En==Eedge]]` (M3, lines 109–124), and `m4Window = Resolve[ForAll[..., Equivalent[Sratio<1, En>Eedge]]]` (M4, lines 135–145).
2. M4's global window-equivalence certificate has NO counterpart in the `.py` at all — the SymPy script never proves the `S<1 ⟺ E>E_edge` biconditional globally; it only asserts the edge formula and pointwise sign of the derivative. The `.wl` proves more, via different machinery (ForAll over the whole positive domain).
3. The `.wl` derives `t_collapse` and the edge from scratch with `FullSimplify`/`Solve`/`Reduce` rather than echoing the SymPy `1/Gamma_coll` choreography.

This is an independent re-derivation; no `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree exactly. Shared symbolic results: `S(E) = √2 √χ √g λ √m /(2 √μ √(E-Vpeak))` (py output line 26 ≡ wl M3 `(lambdaEff Sqrt[(chiPeak gUV ms)/(En muEta - muEta Vpeak)])/Sqrt[2]`); `E_edge = Vpeak + χ g λ² m/(2μ)` (py line 28 ≡ wl M3 Reduce line); `dt_cross/dE < 0` (py asserts, wl proves globally). Benchmark numerics agree across both within tol: `t_collapse≈9.43066476`, `E_safe,min≈5.32265943`, `v_safe,min≈0.07469791`, `ratio≈1.25948037`, `t_collapse_raw≈6.17163516`, `E_safe,min_raw≈27.53273095`, `S(E_edge)=1`. No disagreement.

## Verdict justification

The math holds. Attacks tried that failed: (a) checked the heavy-throat cancellation `∂E_edge/∂m_s=0` is non-tautological — it is, because `E_edge` depends on `m_s` through `m_s/μ_η` and the cancellation only occurs after the substitution `μ_η=α m_s` (genuine, not by construction); (b) checked the monotonicity / one-sided window claim — the load-bearing GLOBAL claim is established by `Resolve[ForAll]` (M1, M3, M4) over the full positive domain, NOT by a single sample point, so the pass-1 fix is present and intact; (c) checked variable-independence of every derivative (all `D[expr,var]` act on exprs that genuinely depend on the variable — no identically-zero traps); (d) checked the benchmark asserts are anchored to externally-imported Session-III values, not hardcoded self-confirmations of in-script algebra (they re-derive `E_edge`, `t_collapse`, `v_safe` from the closed forms and compare to the published constants). The single finding is a low-severity paper-text mismatch (card says "Mathematica audit: none yet" but a full passing `.wl` exists), routed to the user. Verdict: `findings` (1, paper-side, user-gated).

## Value Reconciliation (pass-2 augmentation)

Every script-emitted deliverable value was located in the notes (Section 7, the benchmark carrier). The terse `.tex` card legitimately omits the numerics; they all live correctly in the `.md`, which is a MATCH per the guards.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| v_crit,p = 0.05930851 | py out L61 / py L203 | notes:431 | MATCH |
| t_collapse = 9.43066476 | py out L62 / wl M7 | notes:443 | MATCH |
| E_safe,min = 5.32265943 | py out L63 / wl M7 | notes:457, 490 | MATCH |
| v_safe,min = 0.07469791 | py out L64 / wl M7 | notes:469, 496 | MATCH |
| v_safe/v_crit = 1.25948037 | py out L65 / wl M7 | notes:476 | MATCH |
| v_safe,max(scan) = 0.29654256 | py out L66 / py L208 | notes:496 | MATCH |
| E_max,scan = 80.93332737 | py L179 | notes:490 | MATCH |
| t_collapse_raw = 6.17163516 | py out L75 / wl M7 | notes:533 | MATCH |
| E_safe,min_raw = 27.53273095 | py out L76 / wl M7 | notes:543 | MATCH |
| χ_peak = 21.73204372 | py L176 / wl L177 | notes:417 | MATCH |
| χ_raw = 50.74399964 | py L218 / wl L179 | notes:525 | MATCH |
| E_edge = Vpeak + m_s g_UV χ λ²/(2μ) (symbolic) | py L103 / wl M3 | notes:228–235; card eq L46–49 & appendix:289 | MATCH |
| v_safe² = v_crit² + λ² g χ/μ (symbolic) | py L153 / wl M6 | notes:286–292; card eq L286–293; appendix:297 | MATCH |
| Mathematica audit existence | full passing `.wl` | card L4 "none yet" | MISMATCH → F1 |

INTERNAL (no finding expected in prose): residual/diff check values, tolerances, `transit_min=0.204`/`transit_max=4.054` (sampled-trajectory bounds, present in notes:504–507 anyway), `dE_edge/d{λ,χ,μ}` sensitivity derivatives, `S(E_max,scan)` intermediate, pass/fail flags, `v_safe,min(raw)=0.17254490` (intermediate, py-only print).

reconciliation: complete; 13 deliverable values checked, 0 numeric misaligned (1 paper-text coverage mismatch → F1, the `.wl`-existence statement).

## Self-test notes

Checked the variable-independence trap: every `D[expr,var]`/`sp.diff` acts on an expression that genuinely depends on the variable (e.g. `EedgeHeavy` depends on `ms` only through the residual `Vpeak` term after `μ=α m_s`, so `D=0` is the real physics result, and the M5 second `expectZero` independently pins the heavy-edge value so the zero is not a vacuous-derivative artifact). Checked the single-sample-point trap: the global monotonicity / one-sided-window claims are proven by `Resolve[ForAll]` over the whole positive domain (M1/M3/M4), not at one numeric point — pass-1 fix confirmed present. No new `paper_misalignment` introduced (the F1 fix is a paper-text update routed to the user, not a Codex edit).
