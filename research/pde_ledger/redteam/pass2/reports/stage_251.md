---
unit_id: 251
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-10T00:00:00Z
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 251 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_251.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.md`
- part appendix: cited a wrong project root (mis-rooted) → actual: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_mathematica_audit.txt`

## What the paper claims

Stage 251 replaces the phenomenological Session-IV envelope law `gamma_tot*Vdot` by the first microscopic odd export kernel seen by the active `V` leg. Deliverables: (1) exact cubic coefficient `Gamma_{3,0}=gamma_1 eta_0^2 Omega_{U,0}^4/Delta_0^2`, projected `Gamma_3=Pi_V0^2 Gamma_{3,0}`; (2) exact quintic coefficient `Gamma_{5,-}=a^5/(27 c_s^5) P_{0,-}`, `P_{0,-}=beta_0 s_-/lambda_-`, projected `Gamma_5=Pi_V-^2 Gamma_{5,-}`; (3) the super-Ohmic kernel `K_exp(s)=Gamma_3 s^3+Gamma_5 s^5`; (4) the Schott power identity giving `P_exp=Gamma_3 Vddot^2+Gamma_5 Vdddot^2>=0`; (5) the characteristic polynomial `F(s)=Gamma_5 s^5+Gamma_3 s^3+mu_eta s^2-kappa_V` with exactly one positive root (no finite microscopic `gamma_crit`); (6) the small-kernel slowdown `s_+ = s_0 - (Gamma_3 s_0^2+Gamma_5 s_0^4)/(2 mu_eta)`; (7) the exact event-safe half-plane (`\stagefield`-equivalent boxed): `Gamma3hat + s_c^2 Gamma5hat >= (s_0^2-s_c^2)/s_c^3`; (8) the Session-IV benchmark `Gamma3hat + 0.3013336471 Gamma5hat >= 289.61004918`, with `G5hat_safe ~ 961.09429528`. Appendix table row 100 and eqs `app-part08-main-kexp` / `app-part08-main-safe-half-plane` carry the same forms.

## What the script claims to verify

The SymPy docstring enumerates exactly the seven items above. The assertions verify: the omega^2 coefficient of `N0` equals `eta_0^2 Omega_U0^4/Delta_0^2` (line 59); the projected kernel coefficients (85-86); homogeneity degrees of `Gamma_5`, `P0_minus` under symbol rescaling (89-92); the Schott power balance reduces exactly to `-G3 Vddot^2 - G5 Vdddot^2` (120); `F'(s)` form (156); the safe-`G3` solution and its mass-normalized half-plane form (157-158); the weak-root shift matches the hand form (159); the benchmark numerics (198-203); and that the actual safe-surface polynomials have exactly one positive real root equal to `s_c` (224-227). The `.wl` mirrors these targets with native machinery.

## Paper ↔ script cross-check

| paper deliverable | script check | status |
|---|---|---|
| Gamma_{3,0} cubic coeff | py:59 / wl:90-91 | match |
| Gamma_3 = Pi_V0^2 Gamma_{3,0} | py:72 / wl:92 | match |
| Gamma_{5,-}, P_{0,-} | py:69-71 + homogeneity 89-92 / wl:96-107 | match |
| Gamma_5 = Pi_V-^2 Gamma_{5,-} | py:73,89 / wl:98,106 | match |
| K_exp(s)=Gamma_3 s^3+Gamma_5 s^5 | py:82-86 / wl:99,108-121 | match |
| Schott P_exp>=0 | py:114-120 / wl:135-151 | match |
| F(s), one positive root | py:130-156 / wl:155-186 | match |
| small-kernel slowdown | py:140-159 / wl:190-201 | match |
| safe half-plane | py:132-158,166 / wl:205-209 | match |
| Session-IV benchmark (289.61, 0.3013, 961.09) | py:176-227 / wl:211-248 | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 59 | `simplify(N0_coeff2 - eta0^2 OmU0^4/Delta0^2)==0` | deliverable 1 | yes (coeff from Series, not pre-set) |
| A2 | sympy | 85-86 | kernel coeff == Gamma5/Gamma3 | 3 | yes |
| A3 | sympy | 89-92 | rescale-homogeneity of Gamma5,P0_minus | 2 | yes (symbols present) |
| A4 | sympy | 120 | Schott balance == -G3 Vddot^2-G5 Vdddot^2 | 4 | yes (built from F_odd, not asserted form) |
| A5 | sympy | 156 | `Fprime - (5G5 s^4+3G3 s^2+2mu s)==0` | 5 | yes |
| A6 | sympy | 157-158 | safe-G3 solve vs closed form | 7 | yes (solve of char poly) |
| A7 | sympy | 159 | weak-root eps1_sol == root_shift | 6 | yes (series-solve vs hand form) |
| A8 | sympy | 198-203 | benchmark numerics | 8 | yes |
| A9 | sympy | 224-227 | NSolve roots == sc, count==1 | 5,8 | yes (independent root solve) |
| M1 | math | 90-92 | expectZero coeff/Gamma30/Gamma3 | 1 | yes |
| M2 | math | 104-109 | homogeneity + kernel coeffs | 2,3 | yes |
| M3 | math | 118-125 | K_exp/Sigma coeff presence+absence | 3 | yes |
| M4 | math | 141-151 | Schott residual + ForAll positivity | 4 | yes (Resolve quantifier) |
| M5 | math | 181-186 | F' form + diff-quotient + ForAll monotonicity | 5 | yes (independent certificates) |
| M6 | math | 201 | weak root-shift | 6 | yes |
| M7 | math | 209-248 | safe half-plane + benchmark + NSolve roots | 7,8 | yes |

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration. It shares the physics skeleton (M1–M7 mirror sections 1–8) but uses genuinely Mathematica-native, independent machinery in the load-bearing places:
- M4/M5 prove the positivity and strict-monotonicity facts with `Resolve[ForAll[...], Reals]` quantifier elimination (wl:144-151, 162-177); the SymPy side merely asserts the algebraic form of `F'(s)` (py:156). These are different proof routes for the same claim.
- M5 adds a **difference-quotient certificate** `F(x)-F(y)-(x-y)*diffQuot == 0` (wl:158-161,185) establishing strict increase via a factored Newton-quotient that has NO SymPy counterpart.
- M7 confirms the safe-surface roots with `NSolve[..., Reals, WorkingPrecision->50]` (wl:228-248) at 50-digit precision, vs SymPy's `nroots` at machine precision (py:213-227) — independent numeric solvers, agreeing to ~1e-49.
The shared structure is the physical derivation order, not echoed algebra. Independent.

## Engine cross-check

Final symbolic and numeric outputs agree:
- `Gamma_{3,0} = eta_0^2 OmU0^4 gamma_1/(OmU0^2 OmW0^2 - r0^2)^2` — both (py.txt:11, wl prints coeff at wl.txt:9-11, residual 0).
- `Gamma_{5,-} = a^5 beta_0 s_-/(27 c_s^5 lambda_-)`, `Gamma_5 = Pi_V-^2 (...)` — both (py.txt:18-20, wl.txt:22-23).
- Schott `power - dSchott/dt = -G3 Vddot^2 - G5 Vdddot^2` — both (py.txt:32, wl.txt:58 residual 0).
- Benchmark: `s_c=0.5489386551`, `s_c^2=0.3013336471`, `rhs=289.61004918`, `G5hat=961.09429528`, positive root = `0.5489386551` — both engines (py.txt:56-69, wl.txt:94-113). SymPy machine-precision vs Mathematica 50-digit, agree to displayed precision. No disagreement.

## Verdict justification

Clean. I read the card, notes, and appendix first and confirmed all ten paper deliverables map to anchored, non-tautological script checks in both engines. Attacks tried and failed: (a) the pass-1 tautology/round-trip concern — A6 (`safe_eq` solved from the characteristic polynomial vs an independently hand-written closed form) and A7 (`eps1_sol` from a Series-solve vs `root_shift` hand form) are genuine cross-derivations, not define-then-assert; A9/M7 independently NSolve the actual safe-surface polynomials and confirm the root is `s_c`, which would fail if the closed forms were wrong — this decisively retires the round-trip risk. (b) Variable-independence trap — the rescaling-homogeneity subs (py:89-92) all test symbols that genuinely appear, so no identically-zero/trivial-pass; the Schott `D[v,{t,n}]` derivatives all act on `v(t)` which depends on `t`. (c) Symbol assumptions — `positive=True`/`Reals && >0` are physically justified (frequencies, projection amplitudes, growth rates, mass/stiffness). (d) Symmetry — no unbounded symmetric-domain integrals. Engines agree and the `.wl` is independent. Outputs are fresh (scripts mtime 13:13/13:15, outputs 13:16).

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Gamma_{3,0}=gamma_1 eta_0^2 OmU0^4/Delta_0^2 | py:49 / py.txt:11; wl.txt:9-11 | notes:155-160 (boxed); tex:30-31 | MATCH |
| Gamma_3=Pi_V0^2 Gamma_{3,0} | py:72 / py.txt:19 | notes:162-166 (boxed); tex:32 | MATCH |
| P_{0,-}=beta_0 s_-/lambda_- | py:69 / py.txt:16 | notes:181-183; tex:40 | MATCH |
| Gamma_{5,-}=a^5/(27 c_s^5) P_{0,-} | py:71 / py.txt:18 | notes:178-184 (boxed); tex:37-38 | MATCH |
| Gamma_5=Pi_V-^2 Gamma_{5,-} | py:73 / py.txt:20 | notes:188-192 (boxed); tex:42 | MATCH |
| K_exp(s)=Gamma_3 s^3+Gamma_5 s^5 | py:82 / py.txt:34 (wl) | notes:245-247 (boxed); tex:eq-laplace-d/appendix-kexp | MATCH |
| P_exp=Gamma_3 Vddot^2+Gamma_5 Vdddot^2 | py:120 / py.txt:32 | notes:295-300 (boxed); tex:74-75 | MATCH |
| F(s)=Gamma_5 s^5+Gamma_3 s^3+mu_eta s^2-kappa_V | py:130 / py.txt:37 | notes:318-320 (boxed); tex:80-81 | MATCH |
| s_+ = s_0-(Gamma_3 s_0^2+Gamma_5 s_0^4)/(2 mu_eta) | py:137,144 / py.txt:44-45 | notes:366-372 (boxed); tex:85-88 | MATCH |
| safe half-plane Gamma3hat+s_c^2 Gamma5hat>=(s_0^2-s_c^2)/s_c^3 | py:166 / py.txt:50 | notes:413-417 (boxed); tex:appendix-safe-half-plane | MATCH |
| s_c = 0.5489386551062235 | py:180,199 / py.txt:56 | notes:459 | MATCH |
| s_c^2 = 0.3013336471 | py:200 / py.txt:57 | notes:466; tex:95 | MATCH |
| rhs = 289.61004918 | py:201 / py.txt:59,61 | notes:468,480; tex:95 | MATCH |
| G5hat_safe = 961.09429528 | py:202 / py.txt:62 | notes:482 | MATCH |
| s0_from_t ~ 6.94311 (=gamma_crit to 1e-6) | py:181,198 / py.txt:57 | notes:461-463 (explicitly to displayed precision) | MATCH |
| quotient 961.09/289.61 ~ 3.318 | implied (not asserted) | notes:487-489 | MATCH (prose only) |

INTERNAL (scaffolding, no finding): `N0_series` / `n0Series` (intermediate expansion), `Sigma_exp(omega)` (the i-prefixed kernel, a restatement of K_exp — appears notes:199-205/tex:eq-main-kexp anyway), `F'(s)` form, `Schott storage S_odd`, all `expectZero`/`expectTrue`/`PASS` flags, tolerances, NSolve root lists, `diffQuot` difference-quotient certificate.

One emitted value, `dimensionless rhs = 158.978150899687` (py:187,203, py.txt:63; wl computes weightNum etc. but does NOT assert it), is absent from both `.tex` and `.md`. Assessment: this is `(s_0^2-s_c^2)/s_c^2 = s_0^2/s_c^2 - 1`, a secondary bookkeeping normalization, NOT one of the stage's five stated main results nor a boxed/`\stagefield{Output}` quantity — the notes deliberately frame the safe surface via the `/s_c^3` rhs (which IS present). Per the augmentation guard ("only deliverables absent from BOTH count"), this is INTERNAL, not a MISSING-DELIVERABLE. No finding.

reconciliation: complete; 16 deliverable values checked, 0 misaligned

## Self-test notes

Variable-independence trap: checked every rescaling-substitution (py:89-92) and every `D[v,{t,n}]`/`sp.diff(V,t,n)` — all act on symbols/functions that genuinely appear, so no identically-zero derivative giving a trivial pass. Symmetry/parity: no unbounded symmetric-domain integrals in this stage. Trivial-case: the safe-surface NSolve checks (py:224-227, wl:245-248) substitute the actual benchmark numbers into the real polynomials and recover `s_c` to ~1e-49, confirming the closed forms are not self-referential. No directive written (zero findings).
