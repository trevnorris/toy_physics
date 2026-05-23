---
unit_id: 066
batch: III.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 066 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage066_wall_figure_of_merit_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage066_wall_figure_of_merit_mathematica_audit.txt`

## What the script claims to verify

The unit packages the Stage-48 explicit thin-wall thresholds into a single dimensionless wall figure of merit `W_wall = 4*pi*a^2*L^2*J1*V0^2/(T_X*ell)` and verifies that comparison of `W_wall` against `W_fail = Pe_req/Delta_inf` and `W_suff = Pe_req/Delta_0` is equivalent to comparing the explicit `V0^2` thresholds. The scripts substitute the Stage-48 `V0_fail^2` and `V0_suff^2` expressions back into `W_wall` and confirm the residual versus the dimensionless thresholds is identically zero. They additionally verify that under the constant-compressibility identification `J1 = I_f/H_w` the same chain applies to `W_H = 4*pi*a^2*L^2*I_f*V0^2/(H_w*T_X*ell)`. The Mathematica side adds monotonicity sign assertions on the partial derivatives of `W_wall` and a bridge identity `W_wall = kappa * G_eq^(tw)` under `kx -> tx*kappa/len^2`.

## Assertion inventory

| #  | Script       | Line  | Form                                                                              | Anchored to claim? |
|----|--------------|-------|-----------------------------------------------------------------------------------|--------------------|
| A1 | sympy        | 63    | `expect_zero("W_wall(V0_fail)^2 - W_fail", W_wall.subs(V0**2, V0_fail_sq) - W_fail)` | yes                |
| A2 | sympy        | 64    | `expect_zero("W_wall(V0_suff)^2 - W_suff", W_wall.subs(V0**2, V0_suff_sq) - W_suff)` | yes                |
| A3 | sympy        | 79    | `expect_zero("J1 = I_f/H_w reduction", W_wall.subs(J1, If/Hw) - W_H)`              | yes                |
| A4 | sympy        | 80    | `expect_zero("W_H(V0_fail)^2 - W_fail", ...)`                                      | yes                |
| A5 | sympy        | 81    | `expect_zero("W_H(V0_suff)^2 - W_suff", ...)`                                      | yes                |
| A6 | mathematica  | 48    | `expectZero["W_wall - kappa G_eq^(tw)", wWall - (gEqTw /. kx -> tx*kappa/len^2)*kappa]` | yes (bridge)       |
| A7 | mathematica  | 53    | `expectZero["W_wall(V0_fail)^2 - W_fail", ...]`                                    | yes                |
| A8 | mathematica  | 54    | `expectZero["W_wall(V0_suff)^2 - W_suff", ...]`                                    | yes                |
| A9 | mathematica  | 73    | `expectTrue["dW/d(V0^2) > 0", dV0Sq > 0]`                                          | yes                |
| A10| mathematica  | 74    | `expectTrue["dW/da > 0", dA > 0]`                                                  | yes                |
| A11| mathematica  | 75    | `expectTrue["dW/dL > 0", dLen > 0]`                                                | yes                |
| A12| mathematica  | 76    | `expectTrue["dW/dell < 0", dEll < 0]`                                              | yes                |
| A13| mathematica  | 77    | `expectTrue["dW/dJ1 > 0", dJ1 > 0]`                                                | yes                |
| A14| mathematica  | 78    | `expectTrue["dW/dT_X < 0", dTx < 0]`                                               | yes                |
| A15| mathematica  | 85    | `expectZero["J1 = I_f/H_w reduction", (wWall /. j1 -> ifMom/hw) - wH]`             | yes                |
| A16| mathematica  | 86-89 | `expectZero["W_H(V0_fail)^2 - W_fail", ...]`                                       | yes                |
| A17| mathematica  | 90-93 | `expectZero["W_H(V0_suff)^2 - W_suff", ...]`                                       | yes                |

## Findings

No findings.

Attempts to break the assertions:

1. **Tautology probe on A1/A2**: `V0_fail_sq` is defined at sympy lines 60-61 as the closed Stage-48 expression `T_X*ell*Pe_req/(4*pi*a^2*L^2*J1*Delta_inf)`. The check substitutes this into `W_wall = 4*pi*a^2*L^2*J1*V0^2/(T_X*ell)` and tests whether the residual against `W_fail = Pe_req/Delta_inf` is zero. If the prefactor of `W_wall` had been mis-set (e.g., extra factor of 2, missing pi, swapped numerator/denominator), the residual would be a nonzero rational in the symbols rather than identically zero. The check therefore exercises the packaging factor non-trivially.
2. **Bridge identity A6**: `wWall - kappa*(gEqTw /. kx -> tx*kappa/len^2)` reduces correctly only when the prefactor of `gEqTw = 4*pi*a^2*j1*v0^2/(kx*ell)` and the rule `kx -> tx*kappa/len^2` combine to cancel against `wWall`. Any sign flip in the kx-substitution or off-by-len-power factor would leave a non-zero residual.
3. **J1 = I_f/H_w reduction**: not tautological — `W_wall` is built with `J1` symbolic and `W_H` with `If/Hw` symbolic; the substitution rule must align with the algebraic relation that the underlying physical identification claims. A wrong identification (e.g., `J1 = If*Hw`) would fail.
4. **Monotonicity sign checks (A9-A14)**: `expectTrue` is evaluated under `$Assumptions` that mark every variable positive. The Mathematica engine cannot return `True` unless the derivative's sign is structurally determined by those assumptions, so the checks are non-trivial.
5. **Symbol-assumption probe**: all symbols are physical positive reals (lengths, areas, times, dimensionless thresholds, squared velocities); the `positive=True` declarations in SymPy and `$Assumptions` in Mathematica are consistent with the physical setup and do not silence branch ambiguities for these rational forms. No hidden domain pitfalls.
6. **Stale output**: both saved outputs are timestamped after their corresponding scripts (sympy script Apr 1, output May 11; mathematica script May 11 11:56, output May 11 12:55). Outputs fresh.

## Independent-derivation check (Mathematica)

The Mathematica script is not a pure transliteration. It shares the same physical premises (the closed-form Stage-48 thresholds and the dimensionless rewrite) but introduces:

- a separate `gEqTw` expression and a bridge identity `W_wall = kappa * G_eq^(tw)` (sympy has no analog),
- monotonicity sign checks via `expectTrue` on each partial derivative (sympy only prints, does not assert sign),
- additional partial-derivative checks (`dW/dJ1`, `dW/dT_X`) that sympy does not even print.

Corresponding SymPy / Mathematica sections that overlap (e.g., the V0_fail/V0_suff substitutions and the J1 reduction) implement the same algebraic identities, but the Mathematica side does not echo SymPy's choreography — it factors `W_wall` and the thresholds independently with `FullSimplify[..., Assumptions -> $Assumptions]` rather than `sp.simplify` of an expanded form, and includes assertions absent from SymPy.

## Engine cross-check

Final outputs (from the saved transcripts) agree on every overlapping quantity:

- `W_wall`: sympy `4*pi*J1*L^2*V0^2*a^2/(T_X*ell)` vs mathematica `(4*a^2*j1*len^2*Pi*v0^2)/(ell*tx)` — identical up to renaming `L -> len`, `J1 -> j1`, `T_X -> tx`.
- `W_fail`: sympy `Pe_req/Delta_inf` vs mathematica `peReq/deltaInf` — identical.
- `W_suff`: sympy `Pe_req/Delta_0` vs mathematica `peReq/delta0` — identical.
- `W_H`: sympy `4*pi*I_f*L^2*V0^2*a^2/(H_w*T_X*ell)` vs mathematica `(4*a^2*ifMom*len^2*Pi*v0^2)/(ell*hw*tx)` — identical.
- All five overlapping residuals (`W_wall(V0_*)^2 - W_*`, `J1 = I_f/H_w reduction`, `W_H(V0_*)^2 - W_*`) report 0 in both engines.
- Partial derivatives `dW/d(V0^2)`, `dW/da`, `dW/dL`, `dW/dell` agree in sign and magnitude between engines.

## Verdict justification

The unit's central claim — that comparing `W_wall` against `W_fail = Pe_req/Delta_inf` and `W_suff = Pe_req/Delta_0` is equivalent to comparing the explicit Stage-48 `V0^2` thresholds — is verified non-tautologically in both engines by substituting the threshold formulas into `W_wall` and confirming the residual is identically zero. The constant-compressibility extension via `J1 = I_f/H_w` is checked independently. The Mathematica side additionally asserts the monotonicity signs the SymPy docstring claims ("manifestly monotone"), so the joint coverage is complete. Symbol domains are consistent with the physical positive-real setup; no branch ambiguities are silenced. Attempts to find tautological structure, hidden hardcoded prefactors, or branch errors all failed.

## Self-test notes

Walked through traps: (1) variable independence — every `D[..., var]` or `sp.diff(..., var)` targets a variable that appears in the expression, so none of the derivatives are trivially zero; (2) no unbounded-domain integrals appear, parity is not at issue; (3) trivial-case pre-check on `W_wall.subs(V0**2, V0_fail_sq) - W_fail`: substituting `Pe_req=1, Delta_inf=1` gives `4*pi*a^2*L^2*J1*(T_X*ell/(4*pi*a^2*L^2*J1))/(T_X*ell) - 1 = 1 - 1 = 0`, confirming the residual collapses correctly. No script edits needed; no directive written.
