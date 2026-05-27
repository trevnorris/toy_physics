---
unit_id: 071
batch: III.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-26T00:00:00Z
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
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage071_tanh_wall_branch.md
  paper_appendix: present
---

# Audit unit 071 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_071.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage071_tanh_wall_branch.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (line 260 `\input{stages/stage_071}` — no per-stage prose row beyond the include)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.txt`

## What the paper claims

Stage 071 chooses the canonical tanh wall `f(xi) = (1 + tanh xi)/2` and the natural local Robin mouth closure `K_m = T_X/ell`, and records the resulting branch controls. The boxed shell-moment ledger (paper eq. `app-stage071-moments`) is `I_f = 1/3`, `I_g = 4/15`, `I_g/I_f = 4/5`. With definitions `chi_s = m c_{s,w} L/hbar`, `Lambda_ell = L/ell`, `Upsilon_w = 4 rho_w^2 V0^2/(hbar^2 c_{s,w}^2)`, the stage's `\stagefield{Output}` is the canonical branch-control laws (boxed eq. `app-stage071-control-laws`): `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2`, `eta = Lambda_ell`, `W_wall = Upsilon_w Lambda_ell^2`. The notes additionally pin the closed forms `T_X = pi a^2 ell hbar^2/(3 m rho_w)`, `K_X = 4 pi a^2 (5 m^2 c_{s,w}^2 ell^2 + hbar^2)/(15 ell m rho_w)`, and `J_1 = 1/(3 H_w) = rho_w/(3 m c_{s,w}^2)` with `H_w = m c_{s,w}^2/rho_w` as derivation-support results.

## What the script claims to verify

Both engines independently (i) compute `f'`, `f''` from the tanh profile and evaluate the moment integrals on the full real line; (ii) cross-check those integrals against a `t = tanh xi` substitution form on `[-1, 1]`; (iii) assert `I_f = 1/3`, `I_g = 4/15`, `I_g/I_f = 4/5`; (iv) build `T_X`, `K_X`, `J_1`, `W_wall` from the integral results, set `K_m = T_X/ell`, and assert that the closed-form `K_m = pi a^2 hbar^2 / (3 m rho_w)` matches and that `Km_expected * L / Tx = L/ell`; and (v) assert that `kappa = K_X L^2/T_X` reduces to `4 chi_s^2 + (4/5) Lambda_ell^2` and that `W_wall` reduces to `Upsilon_w Lambda_ell^2`. The Mathematica script additionally asserts the closed forms of `T_X`, `K_X`, and `J_1` against the notes' formulas.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `I_f = 1/3` | sympy L48 `expect_zero("I_f - 1/3", ...)`; .wl L50 `expectZero["I_f - 1/3", ...]` | match |
| `I_g = 4/15` | sympy L50; .wl L52 | match |
| `I_g/I_f = 4/5` | sympy L52; .wl L54 | match |
| `eta = Lambda_ell` (i.e. `L/ell`) | sympy L71 (closed-form K_m route); .wl L83 (closed-form K_m route) | match |
| `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2` | sympy L81–82 (via pattern subs); .wl L91 (direct identity in physical vars) | match |
| `W_wall = Upsilon_w Lambda_ell^2` | sympy L89; .wl L94 | match |
| `T_X`, `K_X`, `J_1` closed forms (notes-level) | sympy: built but not asserted in closed form; .wl L70/L73/L75 explicit `expectZero` | match (the .wl pins all three; the .py builds them and uses them downstream, which is equivalent verification through `kappa` and `eta`) |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 48 | `expect_zero("I_f - 1/3", If - 1/3)` | `I_f = 1/3` | yes |
| A2 | sympy | 49 | `expect_zero("I_f direct - substitution", If - If_sub)` | cross-check of integration form | yes (non-tautological: two independent integrals — `int_{-oo,oo} f'^2 dxi` vs. `int_{-1,1} (1-t^2)/4 dt`) |
| A3 | sympy | 50 | `expect_zero("I_g - 4/15", Ig - 4/15)` | `I_g = 4/15` | yes |
| A4 | sympy | 51 | `expect_zero("I_g direct - substitution", Ig - Ig_sub)` | cross-check | yes |
| A5 | sympy | 52 | `expect_zero("I_g/I_f - 4/5", ...)` | `I_g/I_f = 4/5` | yes |
| A6 | sympy | 70 | `expect_zero("K_m - pi a^2 hbar^2 / (3 m rho_w)", Km - Km_expected)` | closed-form `K_m` (cross-checks closed-form against symbolic `T_X/ell`) | yes |
| A7 | sympy | 71 | `expect_zero("eta - L/ell (from closed-form K_m)", (Km_expected*L/Tx) - L/ell)` | `eta = Lambda_ell`, non-tautologically pinned via closed-form `K_m` and symbolic `Tx` | yes |
| A8 | sympy | 81–82 | `expect_zero("kappa_reduced - [4 chi_s^2 + 4/5 Lambda_ell^2]", ...)` | `kappa` control law | yes |
| A9 | sympy | 89 | `expect_zero("W_wall_reduced - Upsilon_w Lambda_ell^2", ...)` | `W_wall` control law | yes |
| B1 | mathematica | 50 | `expectZero["I_f - 1/3", ...]` | `I_f = 1/3` | yes |
| B2 | mathematica | 51 | `expectZero["I_f direct - substitution", ...]` | cross-check | yes |
| B3 | mathematica | 52 | `expectZero["I_g - 4/15", ...]` | `I_g = 4/15` | yes |
| B4 | mathematica | 53 | `expectZero["I_g direct - substitution", ...]` | cross-check | yes |
| B5 | mathematica | 54 | `expectZero["I_g/I_f - 4/5", ...]` | `I_g/I_f = 4/5` | yes |
| B6 | mathematica | 70 | `expectZero["T_X exact formula", ...]` | closed-form `T_X` (notes) | yes |
| B7 | mathematica | 71–74 | `expectZero["K_X exact formula", ...]` | closed-form `K_X` (notes) | yes |
| B8 | mathematica | 75 | `expectZero["J_1 exact formula", ...]` | closed-form `J_1` (notes) | yes |
| B9 | mathematica | 82 | `expectZero["K_m - pi a^2 hbar^2 / (3 m rhoW)", ...]` | closed-form `K_m` | yes |
| B10 | mathematica | 83 | `expectZero["eta - L/ell (from closed-form K_m)", ...]` | `eta = Lambda_ell` (non-tautological route) | yes |
| B11 | mathematica | 91 | `expectZero["kappa reduced law", kappa - kappaExpected]` | `kappa` control law (direct identity, no pattern-subs) | yes |
| B12 | mathematica | 94 | `expectZero["W_wall reduced law", wwall - wExpected]` | `W_wall` control law | yes |

All assertions trace to specific paper- or notes-side deliverables. The previous prior-tautological `eta - L/ell` check (which was a definitional rewrite of `K_m = T_X/ell`) has been replaced in both engines by the closed-form-`K_m` route shown at A6/A7 and B9/B10, which is non-tautological because it pins `K_m` to the explicit `pi a^2 hbar^2/(3 m rho_w)` and divides by the symbolic `T_X` built from the integrated `I_f`; the assertion fails under any factor or sign error in `T_X`.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` script is not a transliteration of the `.py`. The two engines share the natural calculation path (compute `f`, `f'`, `f''`; integrate; assemble physical quantities), but their verification idioms differ in load-bearing ways:

- The .py reduces `kappa` by SymPy `subs({m*c_sw*L/hbar: chi_s, L/ell: Lambda_ell})` and then asserts equality with `4 chi_s^2 + (4/5) Lambda_ell^2`. The .wl avoids the pattern-subs entirely and writes `kappaExpected = FullSimplify[4*(m*cSw*L/hbar)^2 + (4/5)*(L/ell)^2, ...]` in the same physical variables, asserting `kappa - kappaExpected == 0`. The SymPy path probes whether `subs` pattern-matches; the Mathematica path probes whether the closed-form combination matches.
- The .wl additionally asserts closed forms for `T_X` (L70), `K_X` (L71–74), and `J_1` (L75) — the .py builds these but does not pin them to their notes-level closed forms (it only asserts the downstream `kappa` and `eta` they feed). This is an honest second-engine difference, not a transliteration.

Quoted parallel sections:

- sympy 32–34: `f = sp.Rational(1, 2) * (1 + sp.tanh(xi))` / `fp = sp.simplify(sp.diff(f, xi))` / `fpp = sp.simplify(sp.diff(fp, xi))`
- mathematica 34–36: `f = (1 + Tanh[xi])/2;` / `fp = FullSimplify[D[f, xi], ...];` / `fpp = FullSimplify[D[fp, xi], ...];`

The shared structure here is inevitable given the profile is `f = (1 + tanh xi)/2`; downstream the engines diverge in verification approach (subs-based reduction vs. direct closed-form identity). Independent-derivation: confirmed; not a `mathematica_transliteration` finding.

## Engine cross-check

Both transcripts show every check passing with residual `0` in the same theorem ledger. Side-by-side highlights:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `f'(xi)` | `1/(2*cosh(xi)**2)` | `Sech[xi]^2/2` |
| `f''(xi)` | `-sinh(xi)/cosh(xi)**3` | `-(Sech[xi]^2*Tanh[xi])` |
| `I_f` | `1/3` | `1/3` |
| `I_g` | `4/15` | `4/15` |
| `T_X` | `pi*a**2*ell*hbar**2/(3*m*rho_w)` | `(a^2*ell*hbar^2*Pi)/(3*m*rhoW)` |
| `K_X` | `4*pi*a**2*(5*c_sw**2*ell**2*m**2 + hbar**2)/(15*ell*m*rho_w)` | `(4*a^2*(hbar^2 + 5*cSw^2*ell^2*m^2)*Pi)/(15*ell*m*rhoW)` |
| `J_1` | `rho_w/(3*c_sw**2*m)` | `rhoW/(3*cSw^2*m)` |
| `W_wall` | `4*L**2*V0**2*rho_w**2/(c_sw**2*ell**2*hbar**2)` | `(4*L^2*rhoW^2*V0^2)/(cSw^2*ell^2*hbar^2)` |
| `kappa` (raw) | `4*L**2*c_sw**2*m**2/hbar**2 + 4*L**2/(5*ell**2)` | `(4*L^2*(ell^(-2) + (5*cSw^2*m^2)/hbar^2))/5` |
| `eta` | `L/ell` | `L/ell` |

Engines agree at the level claimed. No `engine_disagreement` finding. `engines_agree: true`.

Output freshness: both `.txt` mtimes (`2026-05-22 20:11`) post-date their script mtimes (`2026-05-22 20:10`) by ~50 s. `outputs_fresh: true`.

## Verdict justification

`clean`. I attacked the audit on five fronts and each held:

1. **Tautology check on `eta`** (the original finding from the prior audit pass): the current scripts no longer write `expect_zero("eta - L/ell", eta - L/ell)` where `eta = Km*L/Tx` and `Km = Tx/ell` (which would be a definitional rewrite). They write `expect_zero("eta - L/ell (from closed-form K_m)", (Km_expected*L/Tx) - L/ell)` — i.e. they pin `Km_expected` to the closed form `pi a^2 hbar^2/(3 m rho_w)` and divide by the symbolic `Tx` (which is itself built from the integrated `I_f`). This fails if `I_f != 1/3`, if a factor of `ell` is dropped from `Tx`, or if `Km` is mis-defined. The earlier finding's required change has been applied correctly and the new assertion has a genuine failure mode.
2. **`kappa` reduction by pattern-subs**: SymPy's `.subs` of `m*c_sw*L/hbar -> chi_s` into `m^2 c_sw^2 L^2/hbar^2` is fragile, but the transcript confirms it worked, and the assertion `kappa_red - (4 chi_s^2 + 4/5 Lambda_ell^2) == 0` non-trivially exercises the result. The Mathematica engine corroborates by an entirely different idiom (direct identity in physical variables, no pattern-subs), eliminating subs-fragility as a hidden risk.
3. **Moment integrals**: Hand-checked `I_f = (1/4) int sech^4 xi dxi = (1/4) int_{-1}^{1} (1-t^2) dt = (1/4)(4/3) = 1/3` and `I_g = int sech^4 xi tanh^2 xi dxi = int_{-1}^{1} t^2(1-t^2) dt = 4/15`. The substitution-form checks reproduce these.
4. **Closed-form `K_X`**: Computed `4 pi a^2 ell (1/3) Hw + pi a^2 (4/15) hbar^2/(m rho_w ell) = (4 pi a^2/(15 ell m rho_w))(5 m^2 c_{s,w}^2 ell^2 + hbar^2)`. Matches the notes form pinned by the Mathematica assertion at L71–74.
5. **Paper alignment**: every paper-side deliverable (boxed moments and the three control laws) maps to a script-side assertion in at least one engine; the Mathematica engine pins additional notes-level closed forms.

A minor cosmetic note (not a finding): the SymPy docstring header (`stage54`) and both engines' banners (`STAGE 54` / `STAGE 054`) reflect the unit's prior numbering before the part-III reorder. Filenames and content are correct for stage 071; the discrepancy is metadata-only and does not fit any of the ten finding categories.

## Self-test notes

I checked: (1) no `sp.diff(expr, var)` where `var` does not appear in `expr` — all derivatives are wrt `xi` against tanh-expressions, so derivatives are nonzero by construction; (2) the integrands `f'^2 = sech^4 xi /4` and `f''^2 = sech^4 xi * tanh^2 xi` are both even in `xi`, so the symmetric `(-oo, oo)` integrals genuinely sample both halves; (3) the substitution-form integrands `(1 - t^2)/4` and `t^2(1 - t^2)` are also even, so the `[-1, 1]` evaluation is non-degenerate; (4) the closed-form `K_m` check `Km - Km_expected = 0` fails under any factor or sign error in `Tx` (substituting the explicit `Tx = pi a^2 ell hbar^2/(3 m rho_w)` gives `Km = pi a^2 hbar^2/(3 m rho_w) = Km_expected`, matching only for the correct `I_f = 1/3`); (5) the follow-up `(Km_expected * L / Tx) - L/ell` collapses to `L/ell - L/ell = 0` only when `Tx` carries the correct `pi a^2 ell hbar^2/(3 m rho_w)` factor structure. No traps tripped.
