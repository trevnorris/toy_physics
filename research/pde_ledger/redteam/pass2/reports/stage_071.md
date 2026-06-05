---
unit_id: 071
batch: III.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage071_tanh_wall_branch.md]
  paper_appendix: present
---

# Audit unit 071 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_071.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage071_tanh_wall_branch.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows at lines 120, 260, 339 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.txt`

## What the paper claims

Stage 071 picks the canonical smooth finite-thickness tanh wall `f(xi) = (1+tanh xi)/2` (so `chi_phi = f'(xi) = (1/2) sech^2 xi`) and records the resulting branch controls. The card's boxed moments (eq:app-stage071-moments) are `I_f = 1/3`, `I_g = 4/15`, `I_g/I_f = 4/5`. It then imposes the first natural local Robin mouth closure `K_m = T_X/ell` (eq:app-stage071-mouth), which gives `eta = K_m L/T_X = L/ell = Lambda_ell`. The stage `\stagefield{Output}` is the boxed canonical branch controls (eq:app-stage071-control-laws): `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2`, `eta = Lambda_ell`, `W_wall = Upsilon_w Lambda_ell^2`, with the parent ratios `chi_s = m c_{s,w} L/hbar`, `Lambda_ell = L/ell`, `Upsilon_w = 4 rho_w^2 V0^2/(hbar^2 c_{s,w}^2)`. The notes add the explicit branch coefficients (§2): `T_X = pi a^2 ell hbar^2/(3 m rho_w)`, `K_X = 4 pi a^2 (5 m^2 c_{s,w}^2 ell^2 + hbar^2)/(15 ell m rho_w)`, `J_1 = 1/(3 H_w) = rho_w/(3 m c_{s,w}^2)` with `H_w = m c_{s,w}^2/rho_w`. The appendix row (line 120) summarizes it as "Exact I_f=1/3, I_g=4/15, local mouth closure, and branch controls"; the numeric Family-1 closure at line 339 (Lambda_ell=eta=37, chi_s=37/2, kappa=12321/5) is a *downstream* application (the card's "Downstream use" assigns the numeric fixing to Stages 073-074), not a 071 deliverable.

## What the script claims to verify

Both scripts construct `f`, `fp`, `fpp` symbolically and compute the two shell moments by genuine integration over `(-inf, inf)`, cross-checked against an independent `t = tanh xi` substitution form. They assert `I_f = 1/3`, `I_g = 4/15`, `I_g/I_f = 4/5`, and the two direct-vs-substitution agreements. They then build `T_X`, `K_X`, `J_1`, `W_wall` from those moments, derive `K_m = T_X/ell` and `eta = L/ell`, and reduce `kappa = K_X L^2/T_X` and `W_wall` to the parent-ratio control laws. SymPy asserts the moments, `K_m`, `eta`, and the reduced `kappa`/`W_wall` laws (the latter via `.subs` pattern replacement of `m c_sw L/hbar -> chi_s`, `L/ell -> Lambda_ell`). Mathematica asserts the same moments and additionally pins the closed forms of `T_X`, `K_X`, `J_1` and verifies the `kappa`/`W_wall` reductions as algebraic identities in the original symbols (the stronger route). Bottom line tested: the canonical tanh-wall moments and the three control laws.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| `I_f = 1/3` | sympy:48 `If - 1/3 == 0`; wl:50 | match |
| `I_g = 4/15` | sympy:50 `Ig - 4/15 == 0`; wl:52 | match |
| `I_g/I_f = 4/5` | sympy:52; wl:54 | match |
| `K_m = T_X/ell` ⟹ `eta = L/ell = Lambda_ell` | sympy:65-71 (`Km`, `eta`, both asserted); wl:77-83 | match |
| `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2` | sympy:81-82 (via `.subs`); wl:85-91 (identity in original symbols) | match |
| `W_wall = Upsilon_w Lambda_ell^2` | sympy:84-89; wl:87,94 | match |
| `T_X = pi a^2 ell hbar^2/(3 m rho_w)` (notes §2) | sympy:60 print only; wl:70 asserted | match |
| `K_X = 4 pi a^2 (5 m^2 c_sw^2 ell^2 + hbar^2)/(15 ell m rho_w)` (notes §2) | sympy:61 print only; wl:71-74 asserted | match |
| `J_1 = rho_w/(3 m c_sw^2)` (notes §2) | sympy:62 print only; wl:75 asserted | match |

`paper_alignment: aligned`. Every paper-side deliverable has a faithful script-side check; the load-bearing control laws are asserted by both engines, and the §2 branch coefficients are at minimum pinned by Mathematica.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 48 | `If - 1/3 == 0` (If from real integral, line 36) | I_f | yes |
| A2 | sympy | 49 | `If - If_sub == 0` (independent t-sub route) | I_f cross-check | yes |
| A3 | sympy | 50 | `Ig - 4/15 == 0` (Ig from real integral, line 39) | I_g | yes |
| A4 | sympy | 51 | `Ig - Ig_sub == 0` | I_g cross-check | yes |
| A5 | sympy | 52 | `Ig/If - 4/5 == 0` | I_g/I_f | yes |
| A6 | sympy | 70 | `Km - pi a^2 hbar^2/(3 m rho_w) == 0` | K_m closure | yes |
| A7 | sympy | 71 | `(Km_expected L/Tx) - L/ell == 0` | eta = L/ell | yes |
| A8 | sympy | 81 | `kappa_red - (4 chi_s^2 + 4/5 Lambda_ell^2) == 0` | kappa law | yes (substitution-dependent) |
| A9 | sympy | 89 | `W_red - Upsilon_w Lambda_ell^2 == 0` | W_wall law | yes (substitution-dependent) |
| B1 | wl | 50-54 | five moment identities (Integrate-based) | I_f, I_g, ratio | yes |
| B2 | wl | 70 | `tx - pi a^2 ell hbar^2/(3 m rhoW) == 0` | T_X | yes |
| B3 | wl | 71-74 | `kx - 4 pi a^2(5 m^2 cSw^2 ell^2 + hbar^2)/(15 ell m rhoW) == 0` | K_X | yes |
| B4 | wl | 75 | `j1 - rhoW/(3 m cSw^2) == 0` | J_1 | yes |
| B5 | wl | 82-83 | K_m / eta closure | K_m, eta | yes |
| B6 | wl | 91 | `kappa - kappaExpected == 0` (original-symbol identity) | kappa law | yes |
| B7 | wl | 94 | `wwall - wExpected == 0` (original-symbol identity) | W_wall law | yes |

No tautological rows: each moment is produced by an actual `integrate`/`Integrate`, not defined-then-reasserted. The kappa/W_wall reductions are real: SymPy's `.subs` is the weaker mechanism but the output (lines 25, 27) confirms the substitution actually fired; Mathematica confirms the same reduction as an algebraic identity in original symbols (the independent, stronger route).

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.txt:3` and `:31`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.txt:3` and `:46`

**What's wrong:**
Both committed `.txt` outputs have mtime `May 22 20:11`, predating both scripts (`Jun 3 15:59`). The captured banners do not match what the current scripts print. SymPy output line 3 says `STAGE 54 — CANONICAL TANH-WALL BRANCH` and line 31 `STAGE 54 THEOREM LEDGER`, but the current SymPy script prints `STAGE 71 ...` (script lines 24, 91). Mathematica output line 3 says `STAGE 054 ...` and line 46 `STAGE 054 THEOREM LEDGER`, but the current `.wl` prints `STAGE 071 ...` (script lines 26, 96). All numeric/symbolic *results* in the saved outputs agree with the current scripts (I_f=1/3, I_g=4/15, I_g/I_f=4/5, T_X, K_X, J_1, W_wall, K_m, eta, the two reductions — sympy output lines 7-28, wl output lines 8-43), so the math is unaffected; only the stage-label banner is stale.

**Why this matters:**
The transcript is the auditable record cited by the card's `\stagefield{Verification}` line. A reader following the citation sees "STAGE 54/054", which contradicts the stage 071 card and looks like a wrong-stage transcript.

**Required change:**
Re-run both scripts and overwrite the saved `.txt` outputs so the banners read STAGE 71 / STAGE 071 and the transcripts reflect the current source.

**Verification:**
After refresh, sympy output line 3 reads `STAGE 71 — CANONICAL TANH-WALL BRANCH`, line 31 `STAGE 71 THEOREM LEDGER`; mathematica output line 3 reads `STAGE 071 ...`, the ledger banner `STAGE 071 THEOREM LEDGER`. Result lines unchanged.

### F2 — stale numbering self-labels (SymPy docstring)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py:3` and `:5`

**What's wrong:**
The SymPy module docstring carries stale stage-54 self-labels that disagree with the on-disk filename (`...stage071...`) and the card. Line 3: `moving_throat_pde_stage54_tanh_wall_branch_sympy_audit.py` (wrong embedded filename). Line 5: `SymPy audit for Stage 54:`. The runtime banner (`banner("STAGE 71 ...")`, lines 24, 91) and the actual filename are 071, so these are leftover pre-renumber labels (the known +17 / stale-stage drift). The Mathematica script has no analogous stale self-label (its banners are already 071).

**Why this matters:**
Self-label/filename drift inside the canonical script; cosmetic but it is exactly the in-loop self-label class the second-pass policy fixes on any verdict:findings stage. Material change: none (comment/docstring only).

**Required change:**
In `scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py`, edit line 3 `moving_throat_pde_stage54_tanh_wall_branch_sympy_audit.py` → `moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py`, and line 5 `SymPy audit for Stage 54:` → `SymPy audit for Stage 071:`. Optionally also relabel the two banner strings (lines 24, 91) from `STAGE 71` to `STAGE 071` to match the `.wl` zero-padded convention; do not change any math.

**Verification:**
Docstring filename and "Stage" label read 071; script still exits 0 with all `expect_zero` checks passing; refreshed output banner reflects the change.

## Independent-derivation check (Mathematica)

The `.wl` shares the same physical choreography as the `.py` (same `f`, `fp`, `fpp`, same two moment routes, same `T_X/K_X/J_1/W_wall` build), so the *structure* is parallel. But it is not a blind line-by-line port: (1) Mathematica computes the moments with its own `Integrate` (wl:38-41), not by echoing SymPy's numbers; (2) it asserts the closed forms of `T_X`, `K_X`, `J_1` that SymPy only prints (wl:70-75); (3) it verifies the `kappa`/`W_wall` reductions as algebraic identities in the *original* symbols (`kappa - kappaExpected == 0`, wl:85-91), whereas SymPy reduces via fragile `.subs` pattern replacement (py:75-82). The reduction mechanism differs materially. I judge this `independent` (transliteration-leaning at the moment-build level, but with a genuinely distinct and stronger reduction/closed-form layer) — it does not rise to a `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree on every shared quantity (saved outputs side by side):

| quantity | sympy output | mathematica output |
|---|---|---|
| f'(xi) | `1/(2*cosh(xi)**2)` (L5) | `Sech[xi]^2/2` (L5) — identical |
| f''(xi) | `-sinh(xi)/cosh(xi)**3` (L6) | `-(Sech[xi]^2*Tanh[xi])` (L6) — identical |
| I_f | `1/3` (L7) | `1/3` (L7) |
| I_g | `4/15` (L9) | `4/15` (L9) |
| T_X | `pi*a**2*ell*hbar**2/(3*m*rho_w)` (L16) | `(a^2*ell*hbar^2*Pi)/(3*m*rhoW)` (L21) |
| K_X | `4*pi*a**2*(5*c_sw**2*ell**2*m**2 + hbar**2)/(15*ell*m*rho_w)` (L17) | matching (L22) |
| J_1 | `rho_w/(3*c_sw**2*m)` (L18) | `rhoW/(3*cSw^2*m)` (L23) |
| kappa reduced | `4*Lambda_ell**2/5 + 4*chi_s**2` (L25) | `4*chiS^2 + (4*lambdaEll^2)/5` (L38) |
| W_wall reduced | `Lambda_ell**2*Upsilon_w` (L27) | `lambdaEll^2*upsilonW` (L41) |

No disagreement. (Both saved outputs are stale only in their banner label, per F1, not in any result.)

## Verdict justification

The scripts faithfully verify the paper's stage-071 claims. Attacks tried that failed: (a) tautology — the moments are produced by real integrals and asserted against independent literals AND an independent substitution route, not defined-then-reasserted; (b) the `.subs`-based SymPy kappa/W_wall reduction could have silently failed to fire (leaving the assertion trivially true on an un-substituted expression), but the saved output (lines 25, 27) shows the substitution actually produced the reduced form, and Mathematica independently confirms the reduction as an original-symbol identity; (c) symbol domains are physical (`positive=True`/`>0` on all scale parameters; `xi` real), and the improper integrals converge under those assumptions; (d) the appendix Family-1 numbers (Lambda_ell=eta=37, chi_s=37/2, kappa=12321/5) check out against the control law (4·(37/2)² + (4/5)·37² = 1369 + 5476/5 = 12321/5) and belong to the downstream numeric fixing, not to 071's symbolic deliverables. The two findings are label-only (stale `.txt` banners and a stale SymPy docstring self-label), material_change false; they set the verdict to `findings` under the Reading-2 in-loop self-label policy but carry no math impact. I confirm I read the card, the notes, and the appendix rows, and the script's verified claim matches the paper.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit against the `.tex` card and `.md` notes:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `f'(xi) = (1/2) sech^2 xi` | py:42 / wl:43 (out py L5, wl L5) | `.tex:18`, `.md:47` | MATCH |
| `f''(xi) = -sech^2 xi tanh xi` | py:43 / wl:44 (out py L6, wl L6) | `.md:49` | MATCH |
| `I_f = 1/3` | py:48 / wl:50 (out py L7, wl L7) | `.tex:24` (boxed), `.md:53` | MATCH |
| `I_g = 4/15` | py:50 / wl:52 (out py L9, wl L9) | `.tex:25` (boxed), `.md:55` | MATCH |
| `I_g/I_f = 4/5` | py:52 / wl:54 (out py L10/15, wl L19) | `.tex:27` (boxed), `.md:59` | MATCH |
| `T_X = pi a^2 ell hbar^2/(3 m rho_w)` | py:60 / wl:70 (out py L16, wl L21) | `.md:69` | MATCH |
| `K_X = 4 pi a^2(5 m^2 c_sw^2 ell^2 + hbar^2)/(15 ell m rho_w)` | py:61 / wl:71-74 (out py L17, wl L22) | `.md:71` | MATCH |
| `J_1 = rho_w/(3 m c_sw^2)` ( =1/(3 H_w) ) | py:62 / wl:75 (out py L18, wl L23) | `.md:73-75` | MATCH |
| `W_wall = 4 rho_w^2 V0^2 L^2/(hbar^2 c_sw^2 ell^2)` | py:58 / wl:63 (out py L19, wl L24) | `.md:83` | MATCH |
| `K_m = pi a^2 hbar^2/(3 m rho_w)` (= T_X/ell) | py:67 / wl:79 (out py L20, wl L31) | `.tex:33` (`K_m=T_X/ell`), `.md:93` | MATCH |
| `eta = L/ell = Lambda_ell` | py:68 / wl:80 (out py L21, wl L32) | `.tex:34-35` (boxed), `.md:99,123` | MATCH |
| `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2` | py:81 / wl:90-91 (out py L25, wl L38) | `.tex:50` (boxed), `.md:121` | MATCH |
| `W_wall = Upsilon_w Lambda_ell^2` | py:88-89 / wl:93-94 (out py L27, wl L41) | `.tex:53` (boxed), `.md:125` | MATCH |

INTERNAL items (scaffolding, no prose expectation): `If_sub`, `Ig_sub` (independent substitution cross-check routes); `Hw` (intermediate `= m c_sw^2/rho_w`, though it does appear in `.md:75`); `kappaExpected`/`wExpected`/`Km_expected`/`kappa_red`/`W_red` (comparison targets); pass/fail flags and `expectZero`/`expect_zero` residuals.

Note: appendix `.tex:339` numeric Family-1 closure (`Lambda_ell=eta=37`, `chi_s=37/2`, `kappa=12321/5`) is a downstream numeric application (card "Downstream use" → Stages 073-074), not a 071 script deliverable; it is self-consistent with the 071 control law and is not in scope for this reconciliation.

reconciliation: complete; 13 values checked, 0 misaligned.

## Self-test notes

Variable-independence trap: no `diff`/`D` against a non-present symbol — `fp`, `fpp` are derivatives w.r.t. `xi`, on which `f` genuinely depends; no identically-zero derivative dressed as a nonzero check. Parity/symmetry trap: the two improper integrals are of `fp^2` (even·even, sech^4) and `fpp^2` (even, sech^4 tanh^2) over a symmetric domain — both legitimately nonzero, matching the asserted 1/3 and 4/15; the substitution routes integrate `(1-t^2)/4` and `t^2(1-t^2)` over `[-1,1]`, both even, both nonzero. Trivial-case trap: each `expect_zero`/`expectZero` target is a real residual against an external literal or original-symbol identity, confirmed nonzero-before-cancellation by the saved outputs. Both findings are label-only with material_change false; no paper_misalignment and no value mismatch were found.
