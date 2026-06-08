---
unit_id: 144
batch: IV.5
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-07T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage144_unique_regular_canonical_branch.md]
  paper_appendix: present
---

# Audit unit 144 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_144.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage144_unique_regular_canonical_branch.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows/eqs at lines 533-573, 655-775, 787-790, 1178)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.txt`

## What the paper claims

The stage card (`stage_144.tex:7,13-17`) states Stage 144 is the "coupled mouth fixed point and gain selection" ledger step whose audit target is the quoted verification output, and boxes the bottom line: "Lower compensated branch is the unique finite-bias, finite-traction branch in the explicit positive exponential closure." The notes make this concrete with four deliverables: (1) the upper Family-1 compensated branch `g_+^F1 ≈ 2.79795199200529 > 1` is impossible by positivity (`0 ≤ g_c ≤ 1`); (2) the lower branch `g_-^F1 ≈ 0.758035078944663` lies strictly in `(2/π, 1) ≈ (0.63662, 1)`, so the monotone map `g_Π` reaches it at a unique finite `Π_* ≈ 1.50882951349316` (notes:32,64); (3) at `Π_*` the traction/amplitude are finite and moderate, `Σ₀(Π_*) ≈ 1.80594111095636` and `T̂_m(Π_*) ≈ 0.901484054174205` (notes:75-79), reached before the self-matched derivative point `Π_match ≈ 1.90848600654854`, `T̂_m(Π_match) ≈ 1.01132972803599` (notes:84-86); and (4) the branch-selection theorem that the lower branch is uniquely selected as the regular finite positive-source branch (notes:94-107, appendix Theorem at lines 787-790). The appendix supplies the underlying closed forms: `r_F1 = √(4107−100π²)/(10π) ≈ 1.77799` (appendix:560-563), `g_±(r) = r ± ½√(1+r²)` (appendix:542,570), `g_Π = 2Π(2Πe^Π+π)/((4Π²+π²)(e^Π−1))` (appendix:658), the S-kernel at κ=π/2 (appendix:692-694), `R_q = (g_Π−r_F1)²/(1+r_F1²)` (appendix:759), `Σ₀(Π)=Π/(1−R_q S_q)` (appendix:766), and `Σ₀ = (20/9) T̂_m²` (appendix:751).

## What the script claims to verify

Both scripts construct the same physical premises symbolically — `r` (=r_F1), `g_±`, `g_Π`, the S-kernel `S_q`, `R_q`, `Σ₀`, and `T̂_m = √((9/20)Σ₀)` — and then (a) assert the branch inequalities `g_+^F1 > 1` and `2/π < g_-^F1 < 1` (sympy:53-56; wl:72-73); (b) locate `Π_*` (where `g_Π = g_-^F1`) and `Π_match` (where `g_Π = π/4`) and assert each numeric target to a `1e-12` tolerance against the notes/appendix values, plus the ordering `0 < Π_* < Π_match` (sympy:32-67; wl:53-88); and (c) assert `Σ₀(Π_*)`, `T̂_m(Π_*)`, and `T̂_m(Π_match)` against the notes targets. The Mathematica script additionally asserts that the cleared-denominator residual genuinely vanishes at its independently-found `Π_*` (wl:79-82) and anchors `Π_*` to Stage 131's independently-owned value (wl:83-84). The verdict applies to: the two branch inequalities, the uniqueness/ordering of the two fixed points, and the moderate finite values of the canonical-point quantities.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Upper branch impossible: `g_+^F1 ≈ 2.79795 > 1` | sympy:53-54 `gplus_N > 1`; wl:72 `N[gPlus,30] > 1` | match |
| Lower branch in `(2/π, 1)`: `g_-^F1 ≈ 0.758035` | sympy:55-56 `2/pi < gminus < 1`; wl:73 | match |
| Unique finite `Π_* ≈ 1.50882951349316` | sympy:32,58-59 nsolve + drift assert; wl:53-54,79-84 cleared-denom FindRoot + residual + owner anchor | match |
| `Σ₀(Π_*) ≈ 1.80594111095636` | sympy:40,62-63; wl:67,86 | match |
| `T̂_m(Π_*) ≈ 0.901484054174205` | sympy:33,60-61; wl:55,85 | match |
| `Π_match ≈ 1.90848600654854` (reached after Π_*) | sympy:35,47-48,64-65; wl:57,65,87 | match |
| `T̂_m(Π_match) ≈ 1.01132972803599` | sympy:36,66-67; wl:58,88 | match |
| Branch-selection theorem (narrative) | sympy:70-79 / wl:90-99 ledger Prints (narrative restatement, not an assertion) | match (narrative) |

`paper_alignment: aligned` — every numeric/symbolic deliverable in the card+notes+appendix has a corresponding non-tautological script assertion, and every script constant matches the paper-stated value.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 47-48 | `Pi_star>0 and Pi_match>Pi_star` (raise) | ordering (deliv. 2/3) | yes |
| A2 | sympy | 53-54 | `gplus_N > 1` (raise) | upper-branch impossible | yes |
| A3 | sympy | 55-56 | `2/pi < gminus_N < 1` (raise) | lower-branch reachability | yes |
| A4 | sympy | 58-59 | `\|Pi_star − 1.50882951349316\| < 1e-12` | Π_* value | yes |
| A5 | sympy | 60-61 | `\|That_star − 0.901484054174205\| < 1e-12` | T̂_m(Π_*) | yes |
| A6 | sympy | 62-63 | `\|Sigma0_star − 1.80594111095636\| < 1e-12` | Σ₀(Π_*) | yes |
| A7 | sympy | 64-65 | `\|Pi_match − 1.90848600654854\| < 1e-12` | Π_match | yes |
| A8 | sympy | 66-67 | `\|That_match − 1.01132972803599\| < 1e-12` | T̂_m(Π_match) | yes |
| A9 | mathematica | 72 | `N[gPlus,30] > 1` else Exit[1] | upper-branch impossible | yes |
| A10 | mathematica | 73 | `2/Pi < gMinus < 1` else Exit[1] | lower-branch reachability | yes |
| A11 | mathematica | 79-82 | `Chop[residual] === 0 \|\| Abs<1e-25` | Π_* is a true root (independent route) | yes |
| A12 | mathematica | 83-84 | `\|piStar − 1.5088295134931555830…\| < 1e-12` (Stage-131 owner) | Π_* value | yes |
| A13 | mathematica | 85-88 | four `Abs[… − target] < 1e-12` | T̂_m(Π_*), Σ₀(Π_*), Π_match, T̂_m(Π_match) | yes |
| A14 | mathematica | 65 | `piStar>0 && piMatch>piStar` else Exit[1] | ordering | yes |

No assertion is tautological: each `Π_*`/`Π_match` value is found by an independent root solve and then compared to a target that is NOT how it was constructed (the SymPy nsolve solves `gPi==gminus`/`gPi==π/4`, not `Π==1.5088…`; the Mathematica FindRoot solves the cleared-denominator residual, not the literal). `Σ₀`, `T̂_m` are derived from the located root through the full kernel chain, so their drift asserts can fail if the kernel formulas drift.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an **INDEPENDENT** derivation for the load-bearing `Π_*` step, not a transliteration. Three corresponding sections:

1. **Π_* root route (the discriminating step).** SymPy: `Pi_star = sp.N(sp.nsolve(sp.Eq(gPi, gminus), 1.5), 30)` (sympy:32) — solves the *rational* equation `gPi == gminus` from a single scalar seed. Mathematica: it does NOT call `FindRoot[gPi == gMinus, …]`; instead it defines a cleared-denominator polynomial-in-(p, e^p) residual `gThresholdResidual[p_] := 2*p*(2*p*Exp[p] + Pi) - gMinus*(4*p^2 + Pi^2)*(Exp[p] - 1)` (wl:52) and roots it with a **bracketing seed pair** `FindRoot[…, {piM, 1.4, 1.6}]` (wl:53-54), exploiting monotonicity. This is a genuinely different formulation (multiply through by the strictly-positive denominator `(4p²+π²)(e^p−1)` and root the polynomial form), explicitly commented as such (wl:46-51), with precedent in stages 131-F3 and 142.
2. **Residual self-check + cross-stage anchor (.wl-only).** The `.wl` adds an assertion absent from SymPy: `piStarResidual = Chop[N[gThresholdResidual[piStar],30], 10^-30]` must be `=== 0` or `< 1e-25` (wl:79-82), then anchors `piStar` to Stage 131's independently-owned value `1.50882951349315558300555075595` (wl:83-84), which matches Stage 131's own SymPy output line 2. This is independent corroboration, not an echo of Stage 144's SymPy.
3. **Shared physics premises (legitimately common).** Both engines encode the same `r`, `gPi`, `Sq`, `Rq`, `Sigma0`, `That` (sympy:18-25 vs wl:32-39). This is the stage's physical setup that BOTH engines must encode from the appendix — it is the premise, not transliterated derivation algebra; the discriminating derivation (the Π_* solve) differs.

**Sign verification (the §131-F3/144 trap):** the cleared-denominator residual on wl:52 uses `(Exp[p] - 1)` = `(e^p − 1)`, the **correct** sign. A `(1 − e^p)` flip would put a ~6366-magnitude nonzero residual at the true root and trip the `Abs < 1e-25` guard at wl:80. The only `(1 - Exp…)` occurrence is `(1 - Exp[-piM])` in the `sQ` denominator (wl:36), which is the legitimate `(1 − e^{−Π})` factor from the S-kernel (appendix:694), not the residual. The saved `.wl` output line 20 (`PASS: Pi_* solves cleared-denominator residual (independent root)`) confirms the residual vanishes. The stale `168π²` → `100π²` watch is clean: `grep 168` returns nothing in either script, and `r` uses `4107 - 100*Pi^2` (sympy:18, wl:32), matching appendix:562.

## Engine cross-check

Both engines agree to all printed digits on every deliverable:

| value | SymPy (txt) | Mathematica (txt) |
|---|---|---|
| g_-^F1 | 0.758035078944662826919680890414 (l.9) | 0.7580350789446628269196808904141… (l.9) |
| g_+^F1 | 2.79795199200529341011158893417 (l.10) | 2.7979519920052934101115889341720… (l.10) |
| Π_* | 1.50882951349315552747043511772 (l.12) | 1.5088295134931555830055507559542… (l.12) |
| T̂_m(Π_*) | 0.901484054174204022702401688674 (l.13) | 0.9014840541742040389645127111412… (l.13) |
| Σ₀(Π_*) | 1.80594111095635380721796724713 (l.14) | 1.8059411109563538723736729091995… (l.16) |
| Π_match | 1.90848600654854538838378630317 (l.16) | 1.9084860065485453947610306052087… (l.14) |
| T̂_m(Π_match) | 1.01132972803599475860454058210 (l.17) | 1.0113297280359947602555983007368… (l.15) |
| Σ₀(Π_match) | 2.27286181957635360635379332349 (l.15) | 2.2728618195763536137749655615305… (l.17) |

Agreement is to ~16-22 significant figures (limited by SymPy's nsolve precision-30 vs Mathematica's WorkingPrecision-80 / N[...,30]); both clear the asserted `1e-12` tolerance comfortably. No `engine_disagreement`.

## Verdict justification

**clean.** I read the stage card, the notes, and the relevant appendix equations (533-573, 655-775, 787-790) before the scripts, and the script's claim matches the paper's claim deliverable-for-deliverable. Attacks tried and failed: (1) constant-drift — checked `r = √(4107−100π²)/(10π)` against appendix:562, no stale `168π²`; (2) the §131-F3/144 sign trap — the cleared-denominator residual uses `(e^p − 1)` (correct), and the residual-zero assertion would catch a flip; (3) transliteration — the Mathematica Π_* route is an independently-formulated cleared-denominator bracketed FindRoot (plus a residual self-check and a cross-stage anchor SymPy lacks), not a port of SymPy's `nsolve(gPi==gminus)`; (4) tautology — every numeric target is compared against a value reached by an independent solve through the full kernel chain, not against its own construction; (5) freshness — both `.txt` outputs are newer than their scripts and the engines agree to ~16+ figures.

## Self-test notes

I checked: (1) Variable independence — no `diff`/`D` derivatives in this stage; the load-bearing operations are root-solves of `g_Π(Π)=target`, whose residuals genuinely depend on Π, and `Σ₀`/`T̂_m` depend on the located Π through the kernel, so the drift asserts are non-trivial. (2) Trivial-case — the Π_* drift assert (sympy:58, wl:84) compares an nsolve/FindRoot output to a literal that is NOT the solve target, so it can fail; confirmed non-tautological. (3) Sign — re-derived the cleared denominator `(4Π²+π²)(e^Π−1)` and confirmed wl:52 carries `(Exp[p]-1)` with the correct sign, distinct from the legitimate `(1-Exp[-piM])` S-kernel factor at wl:36. No directive is written (zero findings).

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 8 deliverable values checked, 0 misaligned.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| g_-^F1 ≈ 0.758035078944663 | py:19,28 / wl:33,42; sympy txt:9, wl txt:9 | notes:23, appendix:571 | MATCH |
| g_+^F1 ≈ 2.79795199200529 | py:20,29 / wl:34,43; sympy txt:10, wl txt:10 | notes:24, appendix:573 | MATCH |
| 2/π ≈ 0.636619772367581 | py:30 / wl:44; sympy txt:11, wl txt:11 | notes:51 (2/π ≈ 0.636619772367581) | MATCH |
| Π_* ≈ 1.50882951349316 | py:32,38 / wl:53,60; sympy txt:12, wl txt:12 | notes:64, appendix:663 | MATCH |
| Σ₀(Π_*) ≈ 1.80594111095636 | py:40,42 / wl:67,69; sympy txt:14, wl txt:16 | notes:78, appendix:773 | MATCH |
| T̂_m(Π_*) ≈ 0.901484054174205 | py:33,39 / wl:55,61; sympy txt:13, wl txt:13 | notes:79, appendix:775 | MATCH |
| Π_match ≈ 1.90848600654854 | py:35,44 / wl:57,62; sympy txt:16, wl txt:14 | notes:85 | MATCH |
| T̂_m(Π_match) ≈ 1.01132972803599 | py:36,45 / wl:58,63; sympy txt:17, wl txt:15 | notes:86 | MATCH |

INTERNAL (printed but genuine comparison/scaffolding, not stated deliverables — no finding):
- `Σ₀(Π_match) ≈ 2.27286181957635` (py:41,43 / wl:68,70; sympy txt:15, wl txt:17) — printed as a comparison-only intermediate at the self-matched derivative point; the notes report `Π_match` and `T̂_m(Π_match)` for that "for comparison" point (§3) but never `Σ₀` there, so it is not a stage deliverable. Absent from both docs but excluded by the augmentation's deliverable-only guard.
- `r` (=r_F1) — computed internally (py:18, wl:32) but NOT printed as a labeled result by either script, so it is not in the emitted-value set; its numeric value `≈1.77799353547498` is nonetheless correctly carried at appendix:563.

All 8 emitted deliverable values reconcile against the notes and/or appendix. No MISMATCH and no MISSING-DELIVERABLE; the augmentation adds no findings.
