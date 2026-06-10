---
unit_id: 233
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 233 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_233.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (row 78 + §"Rigid-mouth orbit-lock compiler" lines 883-942)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: "Rigid-mouth compiler: mouth tracking alone is not enough; the actual branch must control the transfer-shape/nontracking packet that carries $\Xi_1$." The `Derivation ledger` states the proof is that rigid mouth tracking $\delta\ln R_{\rm tr}=0$ implies $\Xi_1=\delta\ln\mathfrak N_*$, that the static ceiling is read as an internal transfer/turbulence gate, and the remaining finite packet is identified. The notes §9 ("SymPy-backed status") enumerate eight executable deliverables: (1) the Stage-188 observable compiler $\Theta_1,\Xi_1,\mathcal R_1$; (2) track-locked specialization $\delta\ln R_{\rm tr}=0\Rightarrow\Theta_1=0,\Xi_1=\delta\ln\mathfrak N_*$; (3) prefactor identity $\Xi_{\rm load}=N_{01}/N_0-D_{01}/D_0=P_1/P_0$; (4) operator-rigid $D_{01}=0\Rightarrow\Xi_{\rm load}=N_{01}/N_0$; (5) transported ceiling $|\epsilon\Xi_1|\le P_{\rm crit}\hat m_0^2/(\Delta_{\rm norm}+T_{\rm quad})-1$; (6) calibrated $\Delta_{\rm norm}=0$ form; (7) equivalent $\bar P_0$ form $|\epsilon\Xi_1|\le P_{\rm crit}/\bar P_0-1$; (8) numerical recovery of the two carried Stage-224 budgets `0.367930328492646` and `0.737619063660757` from `0.002069792318062885`. Status is "Exact within reduced branch language; turbulence/choke are interpretive labels."

## What the script claims to verify

The SymPy script verifies exactly notes §9's eight deliverables: it builds the compiler symbolically, substitutes $\delta\ln R_{\rm tr}=0$ and asserts the collapsed forms (asserts $\Theta_1=0$, $\Xi_1-\delta\ln\mathfrak N_*=0$, $\mathcal R_1$ collapse), constructs $P_0,P_1$ from independent $N_0,N_{01},D_0,D_{01}$ symbols and asserts $P_1/P_0-\Xi_{\rm load}=0$ (non-tautological — $P_1$ and $\Xi_{\rm load}$ are written via different algebra), specializes $D_{01}=0$, proves the ceiling↔$\bar P_0$ equivalence and the $\Delta_{\rm norm}=0$ reduction, and numerically recovers both Stage-224 budgets via `Pcrit/Pbar - 1` against the carried literals to tol 1e-15. The Mathematica script verifies the same eight via an independent route (Coefficient extraction).

## Paper ↔ script cross-check

| Paper/notes deliverable | Script check | Status |
|---|---|---|
| (1) Stage-188 compiler forms | .py L27-30 build, .wl M1 Coefficient extraction | match |
| (2) track-lock $\Theta_1=0,\Xi_1=\delta\ln\mathfrak N_*$ | .py L45-47 asserts; .wl M2 | match |
| (3) prefactor $\Xi_{\rm load}=P_1/P_0$ | .py L64; .wl M3 | match |
| (4) operator-rigid $D_{01}=0$ | .py L75; .wl M4 | match |
| (5) transported ceiling | .py L87,92; .wl M5 | match |
| (6) calibrated $\Delta_{\rm norm}=0$ | .py L94-95; .wl M6 | match |
| (7) $\bar P_0$ equivalent form | .py L88-92; .wl M7 | match |
| (8) numeric budgets 0.3679.., 0.7376.. | .py L126-127; .wl M8 | match |

Dominant pattern: all match → `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 45 | `assert_zero(Theta1_rm)` | claim 2 | yes |
| A2 | sympy | 46 | `assert_zero(Xi1_rm - dln_Nstar)` | claim 2 | yes |
| A3 | sympy | 47 | `assert_zero(R1_rm + c_eta*dln_eps_eta + dln_Nstar)` | claim 2 | yes |
| A4 | sympy | 64 | `assert_zero(P1/P0 - Xi_load)` | claim 3 | yes |
| A5 | sympy | 75 | `assert_zero(Xi_load_or - N01/N0)` | claim 4 | yes |
| A6 | sympy | 92 | `assert_zero(gate_rhs - (Pcrit/Pbar_expr - 1))` | claim 5/7 | yes |
| A7 | sympy | 95 | `assert_zero(gate_rhs_cal - (Pcrit*mhat0**2/Tquad - 1))` | claim 6 | yes |
| A8 | sympy | 126 | `assert_close(recovered_robust, robust_budget)` | claim 8 | yes |
| A9 | sympy | 127 | `assert_close(recovered_nonempty, nonempty_budget)` | claim 8 | yes |
| M1 | mathematica | 100-101 | `expectZero` Coefficient residuals | claim 1 | yes |
| M2 | mathematica | 102-105 | `expectZero` track-lock from coeff data | claim 2 | yes |
| M3 | mathematica | 114 | `expectZero` quotientResidual | claim 3 | yes |
| M4 | mathematica | 115 | `expectZero` operator-rigid | claim 4 | yes |
| M5-M7 | mathematica | 124,127,129 | `expectZero` ceiling/Pbar forms | claim 5/6/7 | yes |
| M8 | mathematica | 142-144 | `expectZero`+`expectSmall` budgets | claim 8 | yes |

## Findings

None. All script-side checks are non-tautological, all eight paper/notes deliverables are exercised in both engines, no symbol-assumption error, no variable-independence self-test trap.

## Independent-derivation check (Mathematica)

The `.wl` is genuinely independent, not a transliteration. Where the `.py` proves the track lock by direct substitution `Theta1.subs({dln_Rtr:0})` and asserting on the result (L40-47), the `.wl` instead extracts polynomial structure: `thetaCoeff = Coefficient[theta1, dlnRtr]` and `thetaConstant = Coefficient[theta1, dlnRtr, 0]` (L94-97), then asserts `thetaCoeff - 1 == 0`, `xiCoeff + bstar == 0`, and reconstructs the track-locked value from the constant term `xiConstant - dlnNstar == 0` (M2, L101-105). That is a different derivation route to the same claim (coefficient extraction vs. specialization). M3 uses `Together[p1/p0 - xiLoad]` vs the .py's `simplify`. The budget literals in M8 are the same independent carried inputs (legitimately shared physical constants from Stage 224), not algebra echoed from the .py. Verdict: independent.

## Engine cross-check

Both engines agree. SymPy output: all symbolic residuals 0, `Recovered robust budget = 0.367930328492646003...` and `Recovered nonempty budget = 0.737619063660756683...` matching the carried targets to tol 1e-15. Mathematica output: every `residual = 0` for M1-M7 and M8-symbolic; M8 numeric residuals `-900068804371/206979231806288500000000000000` (≈ -4.3e-18) and `-134296606540789/4.14e29` (≈ -3.2e-16), both below `budgetTol = 1e-15`, both PASS. Same identities, same numbers, exact rational arithmetic in Mathematica corroborating SymPy's float recovery.

## Verdict justification

Clean. Attacks tried that failed: (a) the prefactor identity A4/M3 is NOT tautological — `P1 = (N01*D0 - N0*D01)/D0**2` and `Xi_load = N01/N0 - D01/D0` are constructed via independent expressions, and `P1/P0 - Xi_load` simplifying to 0 is a real algebraic fact, not construction-guaranteed. (b) No `sp.diff`/`D[]` anywhere, so the variable-independence self-test trap cannot apply. (c) Symbol domains are sound: $N_0,D_0$ declared `nonzero` (needed for the quotient), the gate symbols `positive` consistent with the physical setup (pressures, masses). (d) The track-lock asserts are substantive: A2 would fail if $\Xi_1$ carried any residual $B_*\,\delta\ln R_{\rm tr}$ term after the substitution. (e) Budget recovery is anchored to genuinely independent carried Pcrit literals, not back-solved to be trivially true. Paper card, notes §9, and appendix row 78 all agree with what both scripts verify. I confirm I read the paper card, notes, and appendix before the scripts, and the script claims match the paper claim.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| $\Xi_1=\delta\ln\mathfrak N_* - B_*\delta\ln R_{\rm tr}$ | .py L28 / .wl L90 / out L5 | notes L87 / card "Derivation ledger" | MATCH |
| track lock $\Xi_1\to\delta\ln\mathfrak N_*$ | .py L46 / out L10 | notes L126, L218; card L13 | MATCH |
| $\Xi_{\rm load}=N_{01}/N_0-D_{01}/D_0=P_1/P_0$ | .py L62-64 / .wl L111-112 / out L16-17 | notes L150-155 | MATCH |
| operator-rigid $\Xi_{\rm load}\to N_{01}/N_0$ | .py L75 / out L20 | notes L163-167 | MATCH |
| gate RHS $P_{\rm crit}\hat m_0^2/(\Delta_{\rm norm}+T_{\rm quad})-1$ | .py L87 / .wl L120 / out L23 | notes L196-201 | MATCH |
| $\bar P_0$-form gate $P_{\rm crit}/\bar P_0-1$ | .py L90 / out L24 | notes L398-403 | MATCH |
| $\bar P_0$ compat point `0.002069792318062885` | .py L115 / .wl L135 / out L33 | notes L250, L412 | MATCH |
| robust budget `0.367930328492646` | .py L116 / .wl L138 / out L34 | notes L254, L268, L286 | MATCH |
| nonempty budget `0.737619063660757` | .py L117 / .wl L139 / out L35 | notes L257, L274, L292 | MATCH |

Internal (scaffolding / intermediate carried inputs, no prose finding): implied `Pcrit_robust = 0.0028313316855593175`, `Pcrit_nonempty = 0.0035965105896846573` (back-computed Stage-224 ceiling inputs, sourced via comment to stage224 ceilings dict — intermediate, not a Stage-233 deliverable), `c_eta`, `P0`, `P1` symbolic intermediates, `budgetTol`, residual values.

reconciliation: complete; 9 deliverable values checked, 0 misaligned

## Self-test notes

Checked: (1) Variable-independence trap — N/A, no derivatives in either script. (2) Symmetry/parity — N/A, no integrals. (3) Trivial-case — substituting `dln_Rtr=0` genuinely zeroes A1 and reduces A2's residual to 0 only because the $B_*$ term drops; A4 verified by hand: $(N_{01}D_0-N_0D_{01})/D_0^2 \div (N_0/D_0) = (N_{01}D_0-N_0D_{01})/(N_0 D_0) = N_{01}/N_0 - D_{01}/D_0 = \Xi_{\rm load}$, residual 0, non-trivial. (4)/(5) N/A (no directive). No directive written — zero findings.

## Stale self-label note (deferred numbering band, informational — NOT a math finding)

The `.py` comments/prints carry forward-renumbered source attributions that disagree with the notes' canonical sourcing and with the `.wl`:
- `.py` L20 "Exact Stage 239 branch-observable compiler" and L32 print "Stage 239 observable compiler" — notes L82 attribute this compiler to **Stage 188**.
- `.py` L113 "carried Stage 241 budgets" and L129 print "Stage 240 / Stage 241 carried numbers" — notes L146/L248/L404 attribute the budgets and compatibility point to **Stage 224** (compat point Stage 223).
The `.wl` is canonical (L133-134 comment correctly cites the Stage-224 ceilings dict; M8 banner "Carried budget recovery from independent critical pressures" carries no stale forward-label). These `.py` self-labels (188→239, 224→240/241) are the known deferred script/output-band numbering work (see redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md), not stage-233 math defects; left as-is, noted here per the deferred-numbering policy.
