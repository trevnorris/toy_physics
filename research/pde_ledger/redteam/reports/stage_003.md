---
unit_id: 003
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
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
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage003_bdg_coupling.md
  paper_appendix: present
---

# Audit unit 003 red-team report (v2 paper-grounded re-audit)

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_003.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage003_bdg_coupling.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 28 for stage 003: "Minimal BdG--wall coupling — StatusReduced / StatusExactClosure — Stable-mode Schur complement, low-frequency wall moments, pole shift formula, and conservative grouped-isotropy diagnostic.")
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage003_bdg_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage003_bdg_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage003_bdg_mathematica_audit.txt`

mtime check: SymPy script 2026-05-21 00:54, output 2026-05-21 11:26 (fresh). Mathematica script 2026-05-21 01:03, output 2026-05-21 11:50 (fresh).

## What the paper claims

The paper card (`stage_003.tex`) `\stagefield{Output}` states verbatim:

> "Stage 003 outputs the scalar Schur complement \eqref{eq:app-stage003-scalar-kernel}, the grouped kernels \eqref{eq:app-stage003-p2-kernel}, the low-frequency coefficients \eqref{eq:app-stage003-d0d2d4}, the exact pole shift \eqref{eq:app-stage003-pole-shift}, and the first microscopic grouped-isotropy criterion."

Concrete deliverables (from boxed equations in the card body and the notes file):

1. **Scalar Schur kernel** (axisymmetric `(a, L)` wall + scalar BdG modes):
   `D_0^eff(omega) = K_0 - omega^2 M_0 - C (Omega_0^2 - omega^2 I)^{-1} C^T`, with `C` a 2xN coupling matrix and `Omega_0^2 = diag(varpi_{0,alpha}^2)`.
2. **Low-frequency scalar moments**: `K_0^eff = K_0 - C Omega_0^{-2} C^T`, `M_0^eff = M_0 + C Omega_0^{-4} C^T`, `N_0^eff = C Omega_0^{-6} C^T`.
3. **Grouped real P_2 kernel** (per channel A in {20, 21, 22}): `D_A^eff(omega) = K_A - M_A omega^2 - sum_alpha g_{A,alpha}^2 / (varpi_{A,alpha}^2 - omega^2)`.
4. **Grouped low-frequency coefficients**: `d_{0A}^eff = K_A - sum_alpha g_{A,alpha}^2/varpi_{A,alpha}^2`, `d_{2A}^eff = -(M_A + sum_alpha g_{A,alpha}^2/varpi_{A,alpha}^4)`, `d_{4A}^eff = -sum_alpha g_{A,alpha}^2/varpi_{A,alpha}^6`.
5. **Exact one-mode pole formula**: `omega_pm^2 = (Omega_eta^2 + varpi^2 +- sqrt((Omega_eta^2 - varpi^2)^2 + 4 g^2/M))/2`.
6. **Perturbative pole shifts** (for `varpi^2 > Omega_eta^2`): `delta Omega_eta^2 = -g^2/M / (varpi^2 - Omega_eta^2) + O(g^4)`, `delta varpi^2 = +g^2/M / (varpi^2 - Omega_eta^2) + O(g^4)`.
7. **Grouped isotropy criterion**: with `a_x = (2 x_20 - x_21 - x_22)/10` and `b_x = (x_21 - x_22)/2`, if all channelwise inputs (`K_A`, `M_A`, `g_{A,alpha}`, `varpi_{A,alpha}`) are channel-independent then `a_x = b_x = 0`.

Notes section 9 enumerates an additional script-backed deliverable not boxed in the .tex but listed as a verification target:

8. **Harmonic selection rule** on an isotropic reference throat: the angular overlap is diagonal in (l, m), so axisymmetric (l=0) wall motion couples only to l=0 matter modes, and grouped real l=2 wall motion couples channelwise inside l=2.

The paper card body (paragraph after eq. `app-stage003-p2-kernel`) explicitly states the audit verifies **finite-mode witness cases** and that the displayed `sum_alpha` formulas then follow by linear superposition. This is a paper-side green light for the script's witness-case approach to (3) and (4).

The card also lists four `\stagefield{Checks}` items: (a) series expansion gives the right K_0^eff / d_{0A}^eff / etc., (b) static stiffness softens (subtract `g^2/varpi^2`), (c) pole shift has the avoided-crossing sign, (d) isotropy follows from lane equality of microscopic inputs.

## What the script claims to verify

Per the SymPy docstring (lines 7-18) and the final ledger print block (lines 401-412): the unit verifies the axisymmetric `(a, L)` wall + scalar BdG Euler-Lagrange equations and the exact effective wall kernel `D_eff(omega) = K - omega^2 M - C (Omega_m^2 - omega^2 I)^{-1} C^T`; the low-frequency renormalizations `K_eff`, `M_eff`, `N_eff`; the exact two-pole spectrum for one wall + one BdG mode and its perturbative shift; channelwise grouped-P_2 self-energies; the trace/`a2`/`b2` anisotropy invariants; preservation of `a2 = b2 = 0` on the isotropic matter-coupled branch; and the harmonic selection rule enforcing the l=0 / grouped l=2 block structure. The Mathematica script (lines 40-264) mirrors the same scope.

The script's witness instantiation is `N_modes = 2` for the axisymmetric scalar sector and `N_modes = 1` for each of the three grouped P_2 channels.

## Paper <-> script cross-check

| # | Paper deliverable | Script-side check (sympy file:line / mathematica file:line) | Status |
|---|---|---|---|
| 1 | Scalar Schur kernel `D_0^eff` | sympy 145-156 (derived by EL elimination, compared to Schur formula); 162-174 (vs hand-typed manual form). mathematica 100, 108-117 (LinearSolve form, EL-elimination derived, compared) | match (2-mode witness; paper-card-approved) |
| 2 | K_0^eff, M_0^eff, N_0^eff | sympy 176-200 (series expansion compared to explicit closed forms). mathematica 132-145 (same series check) | match (2-mode witness) |
| 3 | Grouped P_2 kernels D_A^eff | sympy 285-292 (single-mode form for A in {20, 21, 22}). mathematica 197-201 | match (1-mode witness per channel; paper-card-approved) |
| 4 | d_{0A}, d_{2A}, d_{4A} | sympy 296-318 (coefficient extraction). mathematica 201-217 (dP2 series, coefficient extraction, projection labels) | match (1-mode witness) |
| 5 | Exact pole formula omega_pm^2 | sympy 234-250 (sp.solve of derived dispersion vs closed form). mathematica 154-178 (derived dispersion vs closed form, direct root substitution into dispersion, plus Vieta sum/product) | match |
| 6 | Perturbative pole shifts | sympy 252-269 (series expansion to O(g^2) of closed-form roots, compared to `+/- g^2/M/delta`). mathematica 180-188 | match |
| 7 | Grouped isotropy a_x = b_x = 0 | sympy 320-342 (definitions match paper exactly, isotropic substitution drives both to 0). mathematica 213-230 | match |
| 8 | Harmonic selection rule (notes section 9) | sympy 352-386 (real-Y spherical harmonics, six cross-integral checks plus norms). mathematica 232-260 (full 4x4 overlap matrix compared to identity) | match |

Every paper-side deliverable is covered by a non-tautological script-side check. The `Inputs` field of the stage card (Stage 001 coupling, Stage 002 wall coords, stable BdG mode assumption with `varpi > 0`, harmonic selection on isotropic throat) is honored by the script: `wa, wb, w20, w21, w22 > 0` declared positive, wall coords carry the Stage 002 names `(qa, qL)`, the `varpi^2 > Omega_eta^2` branch is enforced by `delta > 0` in the perturbative substitution.

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 113-128 | EL equation residuals = 0 for qa, qL, xa, xb (vs hand-typed expected form, after sign convention) | Lagrangian (1) | yes |
| A2 | sympy | 156 | `derived D0_eff vs Deff = 0` (eliminate Xa, Xb from EL, compare to Schur formula) | claim 1 | yes |
| A3 | sympy | 174 | `D0_eff - manual form = 0` (Schur formula vs explicit 2x2 closed form) | claim 1 | yes |
| A4 | sympy | 200 | `series match = 0` (low-freq series of D_eff vs Keff/Meff/Neff target) | claim 2 | yes |
| A5 | sympy | 237 | `derived dispersion vs (K - M w2)(varpi2 - w2) - g^2 = 0` (EL elimination vs hand dispersion) | claim 5 | yes |
| A6 | sympy | 249-250 | `roots - closed_form = 0` (sp.solve roots vs paper boxed formula) | claim 5 | yes |
| A7 | sympy | 264-265 | wall/matter-like shift O(g^2) series match | claim 6 | yes |
| A8 | sympy | 339-342 | isotropic a2, b2, D20-D21, D21-D22 = 0 | claim 7 | yes |
| A9 | sympy | 374-383 | individual Y_lm cross-integrals and norms vanish/are unity | claim 8 | yes (six checks plus four norms) |
| B1 | mathematica | 85-94 | EL equation residuals = 0 for qa, qL, xa, xb | Lagrangian (1) | yes |
| B2 | mathematica | 117 | `derived D0_eff vs Deff = 0` (EL elimination vs Schur formula via LinearSolve) | claim 1 | yes |
| B3 | mathematica | 130 | `D0_eff - manual form = 0` | claim 1 | yes |
| B4 | mathematica | 145 | `series match = 0` | claim 2 | yes |
| B5 | mathematica | 168-169 | `derived dispersion vs ... = 0` | claim 5 | yes |
| B6 | mathematica | 175-176 | each root substituted into dispersion gives 0 | claim 5 | yes (genuinely different from sympy's sp.solve path) |
| B7 | mathematica | 177-178 | Vieta sum/product of roots | claim 5 | yes (independent check) |
| B8 | mathematica | 187-188 | wall/matter-like O(g^2) shift | claim 6 | yes |
| B9 | mathematica | 227-230 | isotropic a2, b2, D20-D21, D21-D22 = 0 | claim 7 | yes |
| B10 | mathematica | 259-260 | 4x4 overlap matrix = identity | claim 8 | yes (structurally distinct from sympy's per-pair approach) |

All rows in the "Anchored?" column are `yes`. No tautological or hardcoded-result rows. The 2-mode (axisymmetric) and 1-mode (grouped P_2) witness instantiations are explicitly authorized by the paper card.

## Independent-derivation check (Mathematica)

The Mathematica script's Section I derives EL equations via a hand-coded `timeD[expr]` operator (lines 71-74) rather than calling SymPy's `euler_equations` analogue (Mathematica has `VariationalMethods\`EulerEquations` but the script does not use it). The Schur elimination uses `LinearSolve[oMat - omega^2 IdentityMatrix[2], Transpose[cMat]]` rather than SymPy's `(Omat - omega^2 I).inv() * Cmat.T`. These are genuinely independent algorithms producing the same symbolic result.

Section II derives the dispersion from a `lOne` Lagrangian via Mathematica's standard `D[D[L, D[q,t]], t] - D[L, q]` (independent of SymPy's `sp.euler_equations`). The closed-form roots are then verified by direct substitution into the dispersion **and** by Vieta sum/product identities — neither of which appears in SymPy. This is the most strongly independent section.

Section III re-uses the same `D_A` formula and the same `a2`, `b2` linear combinations as SymPy. The script defines `T0 = (1/Sqrt[5]) DiagonalMatrix[{1, Sqrt[2], Sqrt[2]}]`, `Ta = (1/Sqrt[10]) DiagonalMatrix[{2, -1, -1}]`, `Tb = (1/Sqrt[2]) DiagonalMatrix[{0, 1, -1}]` on lines 204-206 — labeled as "projections onto representation-theoretic basis" — but the projection matrices are never actually applied to anything. `d2Bar` is computed via `Tr[d2coeffMat]/5` (weighted-trace path, distinct from SymPy's scalar expression) but `a2` and `b2` on lines 214-217 are still computed as the same direct scalar combinations as SymPy. This is the v1 `mathematica_transliteration` (F4) finding partially addressed (Section IV consolidated to a single 4x4 overlap matrix; Section III restructured for `d2Bar` only) and accepted as `material_change: false` at v1 verification time. Not re-raised here: the paper-side claim (`a_x = (2 x_20 - x_21 - x_22)/10` and `b_x = (x_21 - x_22)/2`) is the exact form the script computes, so this is a second-engine-policy concern rather than a paper alignment defect, and the v1 audit cycle already adjudicated it.

Section IV computes the full 4x4 overlap matrix `Table[Integrate[...] ...]` and asserts equality with `IdentityMatrix[4]`. This is structurally distinct from SymPy's pick-each-pair approach.

## Engine cross-check

Both engines exit 0 and pass every assertion in their output transcripts.

Section I, D_eff:
- SymPy prints `D0_eff(omega) =` with entries of the form `Kaa - Maa*omega^2 + c1a^2/(omega^2 - wa^2) + c1b^2/(omega^2 - wb^2)`.
- Mathematica prints `D0_eff(omega) = {{kaa - maa*omega^2 + c1a^2/(omega^2 - wa^2) + c1b^2/(omega^2 - wb^2), ...}, ...}`.
- Algebraically identical (modulo `1/(w^2 - omega^2) = -1/(omega^2 - w^2)` sign convention; the manual form check then closes the loop in each engine).

Section II, roots:
- SymPy: `Omega_eta2/2 + varpi2/2 -+ sqrt(M*Omega_eta2^2 - 2 M Omega_eta2 varpi2 + M varpi2^2 + 4 g^2)/(2 sqrt(M))`.
- Mathematica: `(om2 -+ Sqrt[4 g^2/m + (om2 - varpi2)^2] + varpi2)/2`.
- These are the same expression once `sqrt(M*X + 4 g^2)/sqrt(M) = sqrt(X + 4 g^2/M)`. Both pass the `expect_zero` against the paper's closed form.

Section II, perturbative shifts:
- SymPy: `Omega_eta2 - eps^2 g^2/(M delta)` and `Omega_eta2 + delta + eps^2 g^2/(M delta)`.
- Mathematica: `om2 - eps^2 g^2/(delta m)` and `delta + om2 + eps^2 g^2/(delta m)`.
- Identical.

Section III, anomalies:
- SymPy: `a2 = -M20/5 + M21/10 + M22/10 - g20^2/(5 w20^4) + g21^2/(10 w21^4) + g22^2/(10 w22^4)`, `b2 = -M21/2 + M22/2 - g21^2/(2 w21^4) + g22^2/(2 w22^4)`.
- Mathematica: `a2 = (-2 m20 + m21 + m22 - 2 g20^2/w20^4 + g21^2/w21^4 + g22^2/w22^4)/10`, `b2 = (-m21 + m22 - g21^2/w21^4 + g22^2/w22^4)/2`.
- Identical (modulo notation).

Section IV: all per-pair integrals zero in SymPy; full 4x4 overlap matrix is identity in Mathematica. Equivalent statements.

`engines_agree: true`. No `engine_disagreement`.

## Verdict justification

`clean` (paper_alignment: aligned, findings_count: 0).

Every paper-side deliverable — the scalar Schur kernel, the K/M/N moments, the grouped P_2 kernels, the low-frequency d_{0A}/d_{2A}/d_{4A} coefficients, the exact pole formula, the perturbative shifts, the isotropy criterion, and the harmonic selection rule — has a non-tautological script-side check in both engines, with assertions anchored to the named claim. The witness-case instantiation (2 modes for the scalar sector, 1 mode per grouped P_2 channel) is explicitly authorized by the paper card body ("the executable audit closes finite-mode witness cases explicitly; the displayed sum_alpha formulas then follow by linear superposition because each stable BdG mode contributes additively to the Schur complement"). The two checks `Inputs` carries forward — Stage 001 confinement coupling sign and Stage 002 wall coord names — are honored by the script's symbol choices.

Adversarial attacks attempted:

- **Sign convention on the Schur subtraction**: paper has `K_0 - C Omega_0^{-2} C^T` (subtract). Script's `Keff = Kaa - c1a^2/wa^2 - c1b^2/wb^2` (subtract). Match. A flipped sign in the script would fail the series-match assertion (line 200 of `.py`, line 145 of `.wl`). Both engines pass.
- **Sign of pole shift**: paper has `delta Omega_eta^2 = -g^2/(M (varpi^2 - Omega_eta^2))` (negative for `varpi^2 > Omega_eta^2`). Script's wall-like shift is `Omega_eta2 - eps^2 g^2/(M delta)` with `delta > 0`. The corresponding `expect_zero` would fail under a sign flip. Both engines pass.
- **Witness-case multi-mode additivity**: paper allows witness-case audit because additivity over alpha follows by linear superposition. Section I (N = 2 modes) directly verifies the additivity for the scalar sector. Section III (N = 1 per channel) does not separately re-verify additivity for the grouped P_2 channels, but the algebraic structure is identical (single oscillator per channel coupled to a wall coordinate via the same Schur form), so the additivity established in Section I propagates. Not a finding.
- **Coupling matrix convention**: paper's `C` is 2xN (wall rows by mode columns). Script's `Cmat = [[c1a, c1b], [c2a, c2b]]` is 2x2 — first wall coord couples to modes a, b; second wall coord couples to modes a, b. The `Cmat * (Omat - omega^2 I)^{-1} * Cmat.T` then yields the per-element form the paper's closed-form `D_0^eff` expects. Convention matches.
- **Y_lm normalization conventions**: SymPy uses the real-Y_lm forms `Y00 = 1/(2 sqrt(pi))`, `Y20 = sqrt(5)/(4 sqrt(pi)) (3 cos^2 - 1)`, `Y21c = sqrt(15)/(2 sqrt(pi)) sin cos cos(phi)`, `Y22c = sqrt(15)/(4 sqrt(pi)) sin^2 cos(2 phi)`. The orthonormality checks (norms all = 1, cross-overlaps all = 0) confirm these are the standard real STF basis — the paper's selection-rule notes (section 9 of the .md) require exactly this. Match.

The Mathematica multi-line `lRed = ...` continuation that captured only kinetic terms in the v1 pre-fix code is now correctly compensated by the `lRed = lRed + (...)` block (lines 62-67) that adds back the missing potential and coupling terms; the resulting EL equations on lines 80-83 then match the paper's reduced Lagrangian exactly, confirmed by the four PASS lines in the Mathematica output transcript. This was the v1 F2 fix and has been verified.

Confirmed reading order: paper card (`stage_003.tex`) -> notes (`moving_throat_pde_stage003_bdg_coupling.md`) -> part appendix (`stage_appendix_part01.tex`, row 28) -> SymPy script -> Mathematica script -> both saved outputs. Paper claim and script claim match.

## Self-test notes

Walked through the failure modes flagged by the v2 prompt's self-test list:

1. **Variable independence in derivatives**: every `sp.diff(EXPR, VAR)` and `D[expr, var]` in both scripts is on a variable that actually appears in `EXPR` (verified by inspection — e.g., `sp.diff(EL_qa_red, Qa)` after substitution `qa -> Qa exp(-i omega t)` correctly has `Qa` as a dependency).
2. **Symmetry/parity in integrals**: the spherical harmonic overlap integrals span `(theta, phi) in [0, pi] x [0, 2 pi]`. Cross-overlaps vanish by m-orthogonality (phi integral) or l-orthogonality (theta integral); norms are positive even integrands. Both engines obtain the correct values.
3. **Trivial-case pre-check**: for the isotropic substitution `K20 -> K2, K21 -> K2, ...`, the residual `a2` reduces to `-3 M2/10 + 3 M2/10 - (3 g2^2)/(10 w2^4) + (3 g2^2)/(10 w2^4) = 0` and similarly for `b2 = 0`. Matches the assertion.
4. **Multi-line `lRed` continuation**: re-checked Mathematica's parsing rule (newline terminates expression when syntactically complete) — line 55 ends with `1/2 mLL D[qL, t]^2` (a complete subexpression), so the `\n` terminates `lRed`'s assignment with only kinetic terms in the initial assignment. The `lRed = lRed + (...)` extension on lines 62-67 then correctly adds back the missing potential and coupling, restoring the full Lagrangian before the EL operator is applied. Verified consistency with the v1-iter1 failure transcript at log line 5327, which shows the exact residual `-(kaa qa) - kaL qL + c1a xa + c1b xb` that confirms the pre-fix `lRed` had only kinetic terms.
5. **Paper round-trip**: re-read the boxed equations in `stage_003.tex` after auditing the script. The script's symbolic forms for `D_0^eff`, `K_eff/M_eff/N_eff`, `D_A^eff`, `d_0/d_2/d_4`, `omega_pm^2`, and `delta Omega_eta^2 / delta varpi^2` all match the boxed forms verbatim (modulo trivial notation: `varpi^2` vs `varpi2`, `Omega_eta^2` vs `Omega_eta2 / Om2`).
