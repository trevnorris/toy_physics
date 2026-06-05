---
unit_id: 080
batch: III.4
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
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage080_family1_zeta_thresholds.md]
  paper_appendix: present
---

# Audit unit 080 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_080.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage080_family1_zeta_thresholds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 138: "080 & Demand thresholds in \(\zeta\) & \StatusExactClosure{} & Family--1 wall-depth windows converted to support-ratio demand windows.")
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.txt`

## What the paper claims

Stage 080 (`\stagefield{Output}{Demand thresholds \eqref{eq:app-stage080-zeta-chi}--\eqref{eq:app-stage080-zeta-J}.}`) converts the Stage-078 Family-1 wall-depth windows (expressed in the transport-bias variable `Pe_req`) into explicit quadrupole-demand windows in `zeta_req`, using the Stage-079 demand map `zeta_F1(Pe) = A_F1 Omega_Pe^2`. The bottom-line deliverables are four boxed numeric thresholds at `lambda_mu = 1`: `zeta_suff^(chi) ≈ 2.46622291347846`, `zeta_fail^(chi) ≈ 2.46752913273870` (natural shell-weighted branch), and `zeta_suff^(J) ≈ 2.44257571477179`, `zeta_fail^(J) ≈ 2.46752736855058` (conservative lower envelope). The notes add two further claims: each threshold is `zeta_F1(Pe_threshold)` with the four Stage-078 Pe constants (`96.5285247264386`, `11220.5441626259`, `22.0062226330754`, `2558.01892349205` `lambda_mu^2`); and all four thresholds saturate the hard Family-1 ceiling `zeta_max^(F1) = A_F1 pi^2/4 ≈ 2.46752922945601` as `lambda_mu → ∞`. The success threshold sits ≈0.00131 below the ceiling and the failure threshold is essentially saturated.

## What the script claims to verify

The sympy script builds `zeta(Pe) = A_F1 * Omega(Pe)^2` symbolically (with `A_F1` from the `kappa_F1 = 12321/5`, `eta_F1 = 37` Family-1 constants and `Omega` the Stage-079 closed form), substitutes the four Stage-078 Pe floats, prints the four `zeta_*(1)` values, then (a) numerically re-checks each substituted value against a frozen 25-digit `_expected` literal to tolerance `1e-14`, and (b) takes `lambda_mu → ∞` limits of all four thresholds and asserts each saturates to the independently-derived `zeta_max = limit(zeta, Pe, oo)` to tolerance `1e-10`. The Mathematica script mirrors the four-threshold construction, re-checks each against an independently re-evaluated target (`zetaTarget*`) to `1e-14`, checks `lambda_mu → ∞` saturation to `zetaMax = aF1*Pi^2/4` (closed form, not a Pe-limit) to `1e-14`, and adds four ordering assertions (`suff < fail < max` for both branches; `(J) <= (chi)` for both suff and fail).

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| `zeta_suff^(chi)(1) ≈ 2.46622291347846` | sympy line 41/54 print + line 78-83 numeric check; .wl line 50/68 + line 78 | match |
| `zeta_fail^(chi)(1) ≈ 2.46752913273870` | sympy line 42/55 + 79-83; .wl line 51/69 + 79 | match |
| `zeta_suff^(J)(1) ≈ 2.44257571477179` | sympy line 43/56 + 80-83; .wl line 52/70 + 80 | match |
| `zeta_fail^(J)(1) ≈ 2.46752736855058` | sympy line 44/57 + 81-83; .wl line 53/71 + 81 | match |
| thresholds = `zeta_F1(Pe_thresh)` with the 4 Stage-078 Pe constants | sympy lines 36-44 substitute the 4 Pe floats; .wl lines 45-53 | match |
| `zeta_max^(F1) = A_F1 pi^2/4 ≈ 2.46752922945601` | sympy line 33 `limit(zeta,Pe,oo)`; .wl line 43 `aF1*Pi^2/4` | match (cross-validated by two routes) |
| all 4 thresholds → `zeta_max` as `lambda_mu→∞` | sympy lines 86-95; .wl lines 83-96 | match |

Every paper deliverable maps to a script check. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 82-83 | `if diff > 1e-14: raise` (4×, vs frozen `_expected`) | the 4 boxed `zeta_*(1)` values | partial (regression-anchor; `_expected` is a frozen copy of the formula output, so confirms formula-vs-frozen-literal, not an independent route) |
| A2 | sympy | 94-95 | `if abs(lim - zeta_max) > 1e-10: raise` (4×) | `lambda_mu→∞` saturation | yes (`zeta_max` derived via `limit` over Pe, non-tautological) |
| A3 | math | 78-81 | `expectApprox[..., zetaTarget*, 1e-14]` (4×) | the 4 boxed `zeta_*(1)` values | partial (same formula re-evaluated; regression-anchor) |
| A4 | math | 93-96 | `expectApprox[lim*, zetaMax, 1e-14]` (4×) | `lambda_mu→∞` saturation | yes (`zetaMax` = closed form `aF1*Pi^2/4`; cross-route to sympy's Pe-limit) |
| A5 | math | 97-100 | `expectTrue[suff < fail < max; (J)<=(chi)]` (4×) | ordering/structure of the window | yes (would fail under any threshold mis-mapping) |

## Findings

### F1 — (label-only stale cross-stage references; not one of the ten math categories)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py:5` ("Stage-61 Family-1 Pe_req thresholds")
- `.../moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py:6` ("Stage-62 Family-1 demand map")
- `.../moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py:27` ("# Family-1 constants from Stage 62.")
- `.../moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py:35` ("# Stage-61 explicit Pe thresholds.")
- `.../moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py:65` ("# ...the four Stage-61 numerical Pe thresholds...")

**What's wrong:**
The sympy docstring and four inline comments cite upstream sources as "Stage-61" (Pe thresholds) and "Stage-62"/"Stage 62" (the `A_F1` constants and demand map). The paper card states the real upstream sources are **Stage 078** (`\stagefield{Inputs}{... the Stage~078 Peclet windows.}`) and **Stage 079** (the demand map `zeta_F1 = A_F1 Omega_Pe^2`, where `kappa_F1 = 12321/5`, `eta_F1 = 37` are defined — confirmed in `notes/stages/moving_throat_pde_stage079_family1_quadrupole_pe_map.md:45,47`). These "61/62" labels are stale pre-renumber references from the known +17 EM-extension realignment drift (61→078, 62→079); the canonical card/notes are ground truth. They are unambiguous because the paper card and notes name the upstream stages explicitly, and the four Pe constants are byte-identical to Stage 078's. No math is affected — the constants and formulas used are correct.

**Why this matters:**
Stale cross-stage labels in committed scripts mislead anyone tracing provenance and perpetuate the numbering drift the project is actively eliminating. Per the in-loop Reading-2 policy, a verdict:findings stage gets its unambiguous self/cross labels fixed.

**Required change:**
Label-only edits, no math change:
- line 5: "Stage-61" → "Stage-078"
- line 6: "Stage-62" → "Stage-079"
- line 27: "Stage 62" → "Stage 079"
- line 35: "Stage-61" → "Stage-078"
- line 65: "Stage-61" → "Stage-078"
(The banner on line 21 already correctly reads "STAGE 080"; the `.wl` carries no stage-number labels and needs no change.)

**Verification:**
`grep -n "Stage.6[12]\|Stage 6[12]\|Stage-6[12]" scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py` returns nothing; the script still runs and exits 0 (comment-only edit cannot change behavior).

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration of the `.py`. Decisive evidence of independent derivation:
1. **`zeta_max` via a different route.** SymPy (`.py:33`) derives the ceiling by a genuine symbolic limit `zeta_max = sp.limit(zeta, Pe, sp.oo)`. Mathematica (`.wl:43`) uses the closed form `zetaMax = N[aF1*Pi^2/4, 50]` (i.e., it assumes the analytic value `Omega(∞) = π/2`). The two routes agree to 13 digits (`2.46752922945601…`), cross-validating `Omega(∞) = π/2` rather than echoing one algebra.
2. **Different saturation harness.** SymPy checks `lambda_mu→∞` limits and compares to its Pe-derived `zeta_max`; Mathematica checks `lambda_mu→∞` limits and compares to its closed-form `zetaMax` — so the two saturation checks anchor to ceilings obtained by independent means.
3. **Extra structural assertions only in `.wl`.** Lines 97-100 add four ordering checks (`suff < fail < max`, `(J) <= (chi)`) absent from the sympy script — a real second-engine contribution, not a port.
4. Root-finder differs (`FindRoot` `WorkingPrecision->80` vs `nsolve tol=1e-34`). No shared variable choreography beyond the genuinely-shared physical formula `A_F1 Omega^2`, which both must use because it is the stage's definition.

## Engine cross-check

| quantity | sympy output | mathematica output | agree? |
|---|---|---|---|
| `zeta_max^(F1)` | 2.4675292294560112350 | 2.46752922945601223332… | yes (~13 digits; sympy `1e-10`-tol limit is the looser route) |
| `zeta_suff^(chi)(1)` | 2.4662229134784638979 | 2.46622291347846457779… | yes |
| `zeta_fail^(chi)(1)` | 2.4675291327387028754 | 2.46752913273870334015… | yes |
| `zeta_suff^(J)(1)` | 2.4425757147717912819 | 2.44257571477179109710… | yes |
| `zeta_fail^(J)(1)` | 2.4675273685505776147 | 2.46752736855057822496… | yes |
| 4× `lambda_mu→∞` limits | 2.46752922945600… (≈zeta_max) | 2.46752922945601223… (≈zetaMax) | yes |

Both transcripts end cleanly (`Stage 080 Mathematica audit passed.` / `FINAL LEDGER`), all PASS. The `Limit::alimv` warnings in the `.wl` output are benign (the `lambdaMu>0` assumption is ignored for the limit variable; the limits still evaluate correctly because each threshold is a monotone function of `lambda_mu` and the closed-form result is unambiguous).

## Verdict justification

The scripts faithfully and exactly verify every Stage-080 paper deliverable: the four boxed `zeta_*(1)` values, their construction as `zeta_F1` of the four Stage-078 Pe constants, the `zeta_max = A_F1 π²/4` ceiling, and the `lambda_mu→∞` saturation — all cross-checked by two engines that derive `zeta_max` by genuinely different routes. I attacked: (1) tautology — the numeric `_expected`/`zetaTarget` checks are weak (formula re-evaluated against a frozen copy of its own output) but the load-bearing physics is in the non-tautological limit-saturation and ordering checks, where `zeta_max` is independently derived, so I did not raise a `tautological_check` or `insufficient_verification` finding; (2) `Omega(∞)=π/2` — confirmed analytically and matched across the limit and closed-form routes; (3) constant provenance — `kappa_F1=12321/5`, `eta_F1=37`, and all four Pe constants match Stage 079/078 exactly; (4) value reconciliation — all five emitted deliverable values match the `.tex` and `.md` to stated precision. The single finding is label-only: the sympy docstring/comments cite stale pre-renumber "Stage-61/62" instead of the canonical Stage-078/079. Verdict is `findings` solely on that label fix; the math is clean and `paper_alignment` is `aligned`.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 5 deliverable values checked, 0 misaligned.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `zeta_suff^(chi)(1) ≈ 2.46622291347846` | py:54 + sympy.txt:10; wl:73 + math.txt:10 | `stage_080.tex:17`; `...stage080...md:31,93` | MATCH |
| `zeta_fail^(chi)(1) ≈ 2.46752913273870` | py:55 + sympy.txt:11; wl:74 + math.txt:11 | `stage_080.tex:19`; `...stage080...md:95` | MATCH |
| `zeta_suff^(J)(1) ≈ 2.44257571477179` | py:56 + sympy.txt:12; wl:75 + math.txt:12 | `stage_080.tex:25`; `...stage080...md:99` | MATCH |
| `zeta_fail^(J)(1) ≈ 2.46752736855058` | py:57 + sympy.txt:13; wl:76 + math.txt:13 | `stage_080.tex:27`; `...stage080...md:101` | MATCH |
| `zeta_max^(F1) ≈ 2.46752922945601` | py:46 + sympy.txt:5; wl:62 + math.txt:5 | `notes ...md:35,105`; carried from `stage_079.tex:46` (`\eqref{eq:app-stage079-ceiling}`, referenced by `stage_080.tex:29`) | MATCH |

Note on `zeta_max`: the Stage-080 card does not re-box `zeta_max` (it references the Stage-079 ceiling via `\eqref{eq:app-stage079-ceiling}` on line 29), so its absence from a boxed equation in `stage_080.tex` is a legitimate terse-card omission, not MISSING — the value is present in the Stage-080 notes (`md:35,105`) and in the upstream Stage-079 card it cites. The four Stage-078 Pe input constants (`96.5285247264386`, etc.) are inputs, not Stage-080 deliverables, and reconcile to `stage_078.tex:17-27` and the Stage-078 notes; they are printed as labeled map-arguments, not new results.

INTERNAL (scaffolding, no prose expected): `y_F1`/`yRoot` (root of `y tan y = 37`), `A_F1`/`aF1` and `aF1Indep`/`A_F1_recomputed`, `Omega`/`omegaFn`/`omegaIndep`/`_omega_explicit` intermediate expressions, the per-check `diff` residuals (all ~1e-16 or 0), the four `_expected`/`zetaTarget*` frozen reference literals, tolerances `1e-14`/`1e-10`, PASS/FAIL flags, and the `Limit::alimv` warnings.

## Self-test notes

I checked: (1) variable independence — the two `Limit`/`limit` calls are over `lambda_mu`/`Pe` and each integrand/expression genuinely depends on that variable, so no identically-zero-derivative trap (there are no `diff` derivative assertions here). (2) Symmetry/parity — n/a, no symmetric-domain integrals. (3) Trivial-case pre-check — `Omega(∞)=π/2` verified by leading-order asymptotics, giving `zeta_max = A_F1 π²/4`, matching both engines; the `lambda_mu→∞` limits correctly saturate because every threshold ∝ `lambda_mu²` drives `Pe→∞`. (4) Paths — n/a, no missing-script finding; both engines present. (5) Paper round-trip — the single label-only fix changes comments only and introduces no constant, so it cannot create a new `paper_misalignment`.
