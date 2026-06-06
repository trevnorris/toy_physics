---
unit_id: 092
batch: IV.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage092_dynamic_geometry_obstruction.md]
  paper_appendix: present
---

# Audit unit 092 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_092.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage092_dynamic_geometry_obstruction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows at lines 1175, 1184-1186, 1216-1218 reference this stage / MTDC-T8.1)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage092_dynamic_geometry_obstruction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage092_dynamic_geometry_obstruction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.txt`

## What the paper claims

The card's boxed deliverable (stage_092.tex:15-17): "If geometry carries dynamic even moments, the pole fraction becomes \((1+\epsilon_4)/[4(1+\epsilon_2)^2]\)." The Derivation-ledger field names the forced conservative carrier and the obstruction variables \((\epsilon_2,\epsilon_4)\). The Checks field demands three things: (i) static limit \(\epsilon_2=\epsilon_4=0\) returns \(c_{\rm pole}=1/4\); (ii) \(l=0\)/\(l=2\) orthogonality before applying the firewall; (iii) any support/source success carries the minimal-module hypothesis. The notes give the full derivation: low-frequency coefficients `K0=K_g0+K_pole`, `K2=K_g2+K_pole/Omega_Q^2`, `K4=K_g4+K_pole/Omega_Q^4`; the minimal branch identity `K0 K4 = 4 K2^2`; the obstruction formula solving for `K_g0`; the dimensionless variables `eps_2=Omega_Q^2 K_g2/K_pole`, `eps_4=Omega_Q^4 K_g4/K_pole`; `K0 = 4 K_pole (1+eps_2)^2/(1+eps_4)`; `c_pole = K_pole/K0 = (1+eps_4)/[4(1+eps_2)^2]`; `c_geom = 1 - c_pole`; the static-limit values `c_pole=1/4`, `c_geom=3/4`; and the first-order expansion `c_pole = 1/4[1 + eps_4 - 2 eps_2 + O(eps^2)]`. Deliverables: (D1) the obstruction formula for c_pole; (D2) the static-limit recovery 1/4 (and 3/4 for c_geom); (D3) the obstruction formula for K_g0 (geometry contact = 3 K_pole iff static); (D4) the first-order sensitivity expansion.

## What the script claims to verify

The SymPy script (docstring lines 2-7) builds `K_Q^cons(omega) = K_g0 + K_g2 w^2 + K_g4 w^4 + K_pole/(1 - w^2/Omega_Q^2)`, series-expands to O(w^6), reads off K0/K2/K4, forms the branch obstruction `K0 K4 - 4 K2^2`, solves it for `K_g0`, then checks: the static limit of the K_g0 solution equals `3 K_pole`; the dimensionless `c_pole` equals `(1+eps_4)/(4(1+eps_2)^2)`; and the three first-order Taylor coefficients (1/4, -1/2, 1/4). The Mathematica script takes an independent route: it works directly in eps variables from the start, computes `K_0` from the branch identity as `4 K_2^2/K_4` (no omega-series, no Solve), checks it against the notes' closed form `4 K_pole (1+eps2)^2/(1+eps4)`, reads off c_pole, checks the static limits of both c_pole (=1/4) and c_geom (=3/4), and verifies the same three first-order coefficients.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| D1: c_pole = (1+eps4)/[4(1+eps2)^2] | py:67 `expect_zero(cpole_dimless - cpole_expected)`; wl:63 `expectZero[cPole - cPoleExpected]` | match |
| D2: static limit c_pole=1/4 (c_geom=3/4) | wl:66 `c_pole=1/4`, wl:70 `c_geom=3/4`; py:55 static-limit on K_g0 implies 1/4 | match |
| D3: K_g0 obstruction; contact = 3 K_pole iff static | py:55 `expect_zero(Kg0_sol.subs(Kg2:0,Kg4:0) - 3*Kp)` | match |
| D4: first-order expansion 1/4[1+eps4-2eps2] | py:80-82 and wl:79-81 three coefficient checks | match |
| Check (ii): l=0/l=2 orthogonality | none | see note below |
| Check (iii): success carries minimal-module hypothesis | none (prose firewall, not algebraic) | n/a |

Note on Check (ii): the card's `l=0`/`l=2` orthogonality is the upstream isotropic-decoupling hypothesis (the grouped-P_2-as-single-isotropic-pole premise carried in via `\stagefield{Inputs}` from Part III / Stage 74). It is an imported assumption of this stage, not a deliverable this stage derives — the card frames the stage as testing the obstruction formula GIVEN that decoupling. Treating it as a missing script-side check would penalize 092 for not re-deriving an upstream result it explicitly imports; it is correctly out of scope here. Check (iii) is a prose claim-firewall directive, not an algebraic identity. paper_alignment = aligned.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 55 | `expect_zero(Kg0_sol{Kg2,Kg4->0} - 3 Kp)` | D3 | yes |
| A2 | sympy | 67 | `expect_zero(cpole_dimless - (1+eps4)/(4(1+eps2)^2))` | D1 | yes |
| A3 | sympy | 80 | `expect_zero(eps^0 coeff - 1/4)` | D4 | yes |
| A4 | sympy | 81 | `expect_zero(eps2 coeff - (-1/2))` | D4 | yes |
| A5 | sympy | 82 | `expect_zero(eps4 coeff - 1/4)` | D4 | yes |
| A6 | mathematica | 57 | `expectZero[k0FromBranch - 4 kPole(1+eps2)^2/(1+eps4)]` | D1 (K0 closed form) | yes |
| A7 | mathematica | 63 | `expectZero[cPole - (1+eps4)/(4(1+eps2)^2)]` | D1 | yes |
| A8 | mathematica | 66 | `expectZero[cPole{eps2,eps4->0} - 1/4]` | D2 | yes |
| A9 | mathematica | 70 | `expectZero[cGeom{eps2,eps4->0} - 3/4]` | D2 | yes |
| A10 | mathematica | 79-81 | three first-order coefficient checks | D4 | yes |

No assertion is orphaned; every check traces to a stated deliverable. No tautological rows: each subtracts an independently-formed expression from a derived quantity.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration of the `.py`. The two engines take materially different derivation routes:

- SymPy (py:35-50) constructs the full frequency-dependent module `K(omega)`, performs `sp.series(K, omega, 0, 6)`, extracts K0/K2/K4 from the series coefficients, then `sp.solve(Eq(K0*K4-4*K2^2, 0), Kg0)` to get the obstruction solution, and only afterward substitutes the dimensionless variables into `c_pole = Kp/(Kg0_sol+Kp)`.
- Mathematica (wl:38-57) deliberately avoids both the omega-series and `Solve[]` (comment wl:38-39, wl:51: "substitute eps2, eps4 directly without first solving for kg0" / "Compute this directly without Solve[]"). It writes K_2, K_4 directly in eps form, then forms `k0FromBranch = 4*k2Eps^2/k4Eps` algebraically and checks it against the notes' closed form. So Mathematica verifies the branch identity → K0 closed form in one algebraic step, whereas SymPy reconstructs everything from the original omega-expansion and a symbolic solve. These are independent corroborations of the same `c_pole`, not echoed algebra. The only structural overlap (the first-order coefficient block) is a trivial Taylor extraction both can run independently.

## Engine cross-check

Both engines agree at the level they claim. SymPy output line 14: `c_pole in (eps2, eps4) variables = (eps_4 + 1)/(4*(eps_2 + 1)**2)`. Mathematica output line 10: `c_pole = (1 + eps4)/(4*(1 + eps2)^2)`. Identical. c_geom: SymPy line 16 `(4*eps_2**2 + 8*eps_2 - eps_4 + 3)/(4*(eps_2+1)**2)` equals Mathematica line 15 `(3 + 8*eps2 + 4*eps2^2 - eps4)/(4*(1+eps2)^2)` term-by-term. First-order tail: both report dropped tail `-eps2*eps4/2` and the same linear part `1/4 - eps2/2 + eps4/4`. All `expect_zero`/`expectZero` residuals are 0 in both transcripts; the Mathematica transcript shows `PASS` for all seven gated checks and exits 0.

## Verdict justification

Clean. I read the card, the single notes file, and the part-04 appendix rows before opening the scripts, and the script bottom-line assertions verify exactly the card's boxed claim `c_pole=(1+eps4)/[4(1+eps2)^2]` plus all notes deliverables (K_g0 obstruction with the 3 K_pole static contact, the static limits 1/4 and 3/4, and the first-order expansion). Attacks I tried and that failed: (1) tested whether the static-limit assertion py:55 is tautological — it is not, because `Kg0_sol` is the nontrivial output of `solve(K0 K4 - 4 K2^2)` and reducing it to `3 Kp` at Kg2=Kg4=0 genuinely exercises the branch algebra; (2) checked whether the `eps` substitution merely defines the answer into existence — it does not, because the substitution `Kg2->eps2 Kp/OmegaQ^2`, `Kg4->eps4 Kp/OmegaQ^4` follows the notes' definitions and the resulting `c_pole` is independently compared against a separately-formed `cpole_expected`; (3) checked the symbol domains — `Kp, OmegaQ` positive (justified: a physical pole strength and pole frequency), the geometry moments and eps real and unrestricted (correct, contamination may be either sign), no over-strong assumption hides a branch; (4) checked the Mathematica route for transliteration and found it independent (no series, no Solve). Outputs are fresh (both .txt mtimes 14:28/14:29 newer than the 11:13 script mtimes) and their content matches the current scripts. Both engines present and agreeing; checkpoint=False but the higher-bar criteria (both engines, substantive non-tautological assertions, exact paper alignment) are met.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `c_pole = (1+eps4)/[4(1+eps2)^2]` | py out:14; wl out:10 | stage_092.tex:16 (boxed); notes:19,95 | MATCH |
| `K0 = K_g0 + K_pole` | py out:7 | notes:51 | MATCH |
| `K2 = K_g2 + K_pole/Omega_Q^2` | py out:8 | notes:52 | MATCH |
| `K4 = K_g4 + K_pole/Omega_Q^4` | py out:9 | notes:53 | MATCH |
| `K_g0` obstruction solution | py out:11 | notes:72-74 (solved form) | MATCH (equivalent form) |
| static contact `K_g0 = 3 K_pole` | py out:12 (residual 0) | notes:68 | MATCH |
| `K0 = 4 K_pole (1+eps2)^2/(1+eps4)` | wl out:7 | notes:90 | MATCH |
| static limit `c_pole = 1/4` | wl out:13 (residual 0) | stage_092.tex:22; notes:107 | MATCH |
| `c_geom = (3+8eps2+4eps2^2-eps4)/[4(1+eps2)^2]` | py out:16; wl out:15 | notes:99 (`c_geom=1-c_pole`, form) | MATCH (= 1 - c_pole) |
| static limit `c_geom = 3/4` | wl out:16 (residual 0) | notes:109 | MATCH |
| first-order expansion `1/4[1+eps4-2eps2]` | py out:18; wl out:19 | notes:127 | MATCH |
| eps2 = Omega_Q^2 K_g2/K_pole (defn) | py:62; wl:44 | notes:84 | MATCH |
| eps4 = Omega_Q^4 K_g4/K_pole (defn) | py:62; wl:46 | notes:86 | MATCH |

INTERNAL scaffolding (no finding): `K_Q^cons(omega)` printed form, `Series` printed form, `Branch obstruction` numerator/denominator printed form, `K_2/K_4 in eps variables` intermediate prints, "Linear part" and "Dropped higher-order tail" prints, all `expect_zero`/`expectZero` residual flags and `PASS` markers.

reconciliation: complete; 13 deliverable values checked, 0 misaligned.

## Self-test notes

I checked (1) variable independence — no `sp.diff`/`D[]` derivatives are taken in either script, so the "derivative-of-independent-variable → identically zero" trap does not apply; the series/coefficient extraction operates on expressions that genuinely depend on eps2/eps4. (2) Trivial-case pre-check — substituting eps2=eps4=0 into c_pole gives 1/4 and into c_geom gives 3/4 (matching the static-limit asserts), and the K_g0 static substitution gives 3 K_pole as asserted; all `expect_zero` residuals are genuinely the difference of two independently-formed expressions, none collapses to 0 by construction. (3) Symmetry/parity — no unbounded-domain integrals here. No directive written: zero findings, verdict clean.
