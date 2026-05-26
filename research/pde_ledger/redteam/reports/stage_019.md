---
unit_id: 019
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
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
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 019 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_019.tex`
- notes: (none)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row at line 60; `\input{stages/stage_019}` at line 117)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_mathematica_audit.txt`

## What the paper claims

The stage card (`paper/stages/stage_019.tex`) states the stage rewrites the isotropic one-pole and outgoing-normalization targets in parent-wall variables. The four exported equations are:
(i) parent-complete denominator moments `D_0=K_\Sigma-B_0-Z_0`, `D_2=-(M_\Sigma+B_2+Z_2)`, `D_4=-(B_4+Z_4)` (eq. `stage019-parent-bundle`);
(ii) one-pole closure `K_\Sigma = B_0+Z_0 + 3(M_\Sigma+B_2+Z_2)^2/(B_4+Z_4)` (eq. `stage019-one-pole-parent`);
(iii) outgoing normalization `K_\Sigma = B_0+Z_0 + N_0/P_{0,\rm target}` with `P_{0,\rm target} = 54 G c_s^5 / (5 a^5 c^5 \widehat m_0^{\,2})` (eq. `stage019-normalization-parent`);
(iv) compatibility identity `N_0/P_{0,\rm target} = 3(M_\Sigma+B_2+Z_2)^2/(B_4+Z_4)` (eq. `stage019-isotropic-compatibility`).
The Output paragraph reads verbatim: "Stage 019 exports (parent-bundle)--(isotropic-compatibility)." The card also notes the audit "checks the response-sign criterion selecting the positive one-pole branch and the constant-prefactor conditions." The appendix row (line 60) corroborates the status as `\StatusExactClosure{}`.

## What the script claims to verify

The SymPy and Mathematica scripts both: (a) construct D0, D2, D4 verbatim; (b) build the bundle pole series in `x` and check the one-pole numerator equivalence `u_4 - 4 u_2^2 \propto D_0(B_4+Z_4) - 3(M_\Sigma+B_2+Z_2)^2`; (c) solve for `K_\Sigma` from the one-pole condition and from the normalization condition and confirm each matches the paper's closed forms; (d) derive a compatibility expression from the difference of the two `K_\Sigma` forms; (e) solve `P_2 = 0` and `P_4 = 0` for closed-form constant-prefactor expressions `N_2^{\rm closed}`, `N_4^{\rm closed}`, verify the Jacobian determinant, factorize `P_2` and `P_4`, and include mutation guards against epsilon-perturbed closures; (f) factor the one-pole numerator into M-roots, verify Vieta identities, and evaluate `u_2` on both roots; (g) numerically sample three `(D_0, B_4+Z_4)` triples to verify the positive-vs-negative branch sign convention; (h) evaluate a concrete Gaussian wall profile to produce closed numeric values for `M_\Sigma` and `K_\Sigma` as a sanity instantiation. Both engines report STATUS: PASS with zero residuals on all assertions.

## Paper to script cross-check

| Paper-side deliverable | Script coverage | Verdict |
|---|---|---|
| (i) D0/D2/D4 definitions, eq. `stage019-parent-bundle` | sympy lines 31-33; mathematica lines 46-48 (verbatim) | match |
| (ii) one-pole `K_\Sigma`, eq. `stage019-one-pole-parent` | sympy line 51 (`K_from_one_pole` check); mathematica lines 93-98 (M2) | match |
| (iii) normalization `K_\Sigma` with `P_{0,target}`, eq. `stage019-normalization-parent` | sympy lines 47, 52 (`K_from_norm`); mathematica lines 50, 101-105 (M3) | match |
| (iv) compatibility identity, eq. `stage019-isotropic-compatibility` | sympy lines 53-56 derives `3(M+B2+Z2)^2/(B4+Z4) - N0/P0_target` as the algebraic difference of (ii) and (iii); mathematica does not have a separate compatibility check, but (ii)-(iii) imply it | match |
| Response-sign criterion (audit-only) | sympy lines 103-173 (symbolic + 3 numeric samples); mathematica lines 156-179 (M10, M11) | match |
| Constant-prefactor conditions (audit-only) | sympy lines 58-101 (N2/N4 closed-form solves, Jacobian, P2/P4 factorizations, mutation guards); mathematica lines 108-153 (M4-M9) | match |

`paper_alignment: aligned` — every paper-side deliverable has a corresponding non-tautological script-side check, and the constants used by the script (`P_{0,target}` coefficients 54, 5, 5, 5, 2) match the paper card verbatim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 39 | `simplify(one_pole_defect - one_pole_numerator/D0**2) == 0` | (i)+(ii) preamble: M-quadratic structure of `u_4-4u_2^2` | yes |
| A2 | sympy | 51 | `K_from_one_pole - (B0+Z0+3(MSigma+B2+Z2)^2/(B4+Z4)) == 0` | (ii) | yes |
| A3 | sympy | 52 | `K_from_norm - (B0+Z0+N0/P0_target) == 0` | (iii) | yes |
| A4 | sympy | 56 | `compatibility - (3(MSigma+B2+Z2)^2/(B4+Z4) - N0/P0_target) == 0` | (iv) | yes (derives the compatibility relation) |
| A5 | sympy | 85 | `N2_const - N2_const_closed == 0` | constant-prefactor (audit-only) | yes |
| A6 | sympy | 86 | `N4_const - N4_const_closed == 0` | constant-prefactor (audit-only) | yes |
| A7 | sympy | 87 | `const_prefactor_matrix.det() - D0**3 == 0` | constant-prefactor (audit-only) | yes |
| A8 | sympy | 88 | `P2 - (N2 - N2_const_closed)/D0 == 0` | constant-prefactor (audit-only) | yes |
| A9 | sympy | 89 | `P4.subs(N2, N2_const_closed) - (N4 - N4_const_closed)/D0 == 0` | constant-prefactor (audit-only) | yes |
| A10 | sympy | 90-91 | solve results match closed forms | constant-prefactor (audit-only) | yes |
| A11 | sympy | 92-99 | mutation guards (epsilon shift detected) | constant-prefactor (audit-only) | yes |
| A12 | sympy | 101 | `N4_const.subs(KSigma,K_one_pole) - N4_md_one_pole.subs(KSigma,K_one_pole) == 0` | md-equivalent form (audit-only) | yes |
| A13 | sympy | 107-114 | M-root factorization, Vieta sum, Vieta product | response-sign criterion | yes |
| A14 | sympy | 118-120 | `u2_root_positive - root_gap/D0 == 0`; negative-root counterpart; discrimination | response-sign criterion | yes |
| A15 | sympy | 121-173 | 3 numeric sample positivity guards | response-sign criterion | yes |
| A16 | sympy | 184-185 | Gaussian wall integrals `MSigma=sqrt(pi)`, `KSigma=3*sqrt(pi)/2` | not in paper card (concrete example) | n/a (informational) |
| A17 | mathematica | 84 | `P0series - P0 == 0` (Series consistency) | series infrastructure | yes |
| A18 | mathematica | 86-90 | M1 one-pole numerator | (ii) preamble | yes |
| A19 | mathematica | 93-98 | M2 KSigma from one-pole | (ii) | yes |
| A20 | mathematica | 101-105 | M3 KSigma from normalization | (iii) | yes |
| A21 | mathematica | 108-119 | M4/M5 N2/N4 closures from solve | constant-prefactor | yes |
| A22 | mathematica | 122-131 | M6 Jacobian det = D0^3 | constant-prefactor | yes |
| A23 | mathematica | 134-144 | M7/M8 factorizations | constant-prefactor | yes |
| A24 | mathematica | 146-154 | M9 mutation guards | constant-prefactor | yes |
| A25 | mathematica | 156-168 | M10 M-root + Vieta identities | response-sign criterion | yes |
| A26 | mathematica | 171-178 | M11 u2 on M roots | response-sign criterion | yes |
| A27 | mathematica | 187-197 | M12 Gaussian wall integrals | not in paper card (concrete example) | n/a (informational) |

No tautological rows; every assertion exercises a non-trivial algebraic claim. The two `n/a` rows (A16, A27) implement a concrete-instantiation sanity check that is not load-bearing — they confirm one specific Gaussian profile produces specific numeric values, which does not break paper alignment because the paper card does not assert any specific wall profile.

## Findings

### F1 — symbol_assumption_error

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py:25-29`

**What's wrong:**
The SymPy script declares every symbol with only `nonzero=True`:

```
KSigma, MSigma = sp.symbols('KSigma MSigma', nonzero=True)
B0, B2, B4, Z0, Z2, Z4 = sp.symbols('B0 B2 B4 Z0 Z2 Z4', nonzero=True)
N0, N2, N4 = sp.symbols('N0 N2 N4', nonzero=True)
mhat0, G, cs, a, c = sp.symbols('mhat0 G cs a c', nonzero=True)
eps = sp.symbols('eps', nonzero=True)
```

The physical setup is real (wall constants, projector moments, target-surface parameters), and several quantities (`D_0`, `B_4+Z_4`, the dimensionful prefactors `G, c_s, a, c, \widehat m_0`) are physically positive — the numeric sample loop at line 165 explicitly enforces `D0_value > 0 and tail_value > 0`. The Mathematica mirror at `.wl:62-69` correctly declares `Element[{...}, Reals]` plus the same nonzero set. Without `real=True`, SymPy carries the `sqrt` arguments in their unsigned form: the saved output line 23 reports `u2_on_positive_root = -sqrt(3)*sqrt(-(B4+Z4)*(B0-KSigma+Z0))/(3*B0-3*KSigma+3*Z0)`, which is algebraically equivalent to `+sqrt(D0(B4+Z4)/3)/D0` only when `D0>0` and `B4+Z4>0` — exactly the assumption the script then has to re-establish numerically.

**Why this matters:**
Every current assertion passes because the residuals algebraically cancel without needing positivity. The latent risk is that any future edit that introduces a `simplify`, `radsimp`, `together`, or `sqrt`-manipulation step that *does* depend on sign would silently produce a wrong simplification, and the assertion could pass for the wrong reason. The Mathematica/SymPy declaration asymmetry also undermines the second-engine independence guarantee: the two engines verify the same identities under different domain assumptions, so a sign-sensitive disagreement would not surface as `engine_disagreement`.

**Required change:**
Add `real=True` to all physical-symbol declarations and `positive=True` to the symbols the physical setup makes positive (the dimensionful prefactors `G, c_s, a, c, \widehat m_0`). Concretely, replace lines 25-29 with:

```python
KSigma, MSigma = sp.symbols('KSigma MSigma', real=True, nonzero=True)
B0, B2, B4, Z0, Z2, Z4 = sp.symbols('B0 B2 B4 Z0 Z2 Z4', real=True, nonzero=True)
N0, N2, N4 = sp.symbols('N0 N2 N4', real=True, nonzero=True)
mhat0, G, cs, a, c = sp.symbols('mhat0 G cs a c', positive=True)
eps = sp.symbols('eps', real=True, nonzero=True)
```

Leave KSigma/MSigma/B*/Z* without `positive=True` because the M-root analysis genuinely sweeps both signs of `MSigma` (the positive and negative one-pole branches), and the script must not pre-commit `D_0` positivity into the algebra — that is what the response-sign criterion is testing. `positive=True` on the dimensionful prefactors is safe and unambiguously physical.

**Verification:**
After the edit, re-run `redteam exec-sympy 019`. All existing assertions must still pass. The saved-output line currently reading `u2_on_positive_root = -sqrt(3)*sqrt(-(B4+Z4)*(B0-KSigma+Z0))/(3*B0-3*KSigma+3*Z0)` is permitted to change form (e.g., the sign cosmetics may collapse), but the line `positive-root numeric u2 = 0.6324555320336759` must remain identical, as must every `... = PASS` line and `STATUS: PASS`.

## Independent-derivation check (Mathematica)

The two engines derive `u_2` and `u_4` by genuinely different routes:
- SymPy hardcodes the closed forms (sympy lines 35-36): `u2 = -D2/D0`, `u4 = (D2**2 - D0*D4)/D0**2`.
- Mathematica builds the pole series and reads off coefficients (mathematica lines 71-75): `den = D0 + D2*x + D4*x^2`, `poleSeries = Normal[Series[1/den, {x, 0, 2}]]`, then `u2 = FullSimplify[Coefficient[normalizedPoleSeries, x, 1], ...]` and `u4 = FullSimplify[Coefficient[normalizedPoleSeries, x, 2], ...]`.

Similarly for the bundle moments: SymPy writes `P2 = sp.factor(sp.together((D0*N2 - 2*D2*N0)/D0**2))` (closed-form, sympy lines 42-45), whereas Mathematica builds `bundleSeries = Normal[Series[D0(N0 + N2*x + N4*x^2)/den^2, {x, 0, 2}]]` and extracts `P2 = FullSimplify[Coefficient[bundleSeries, x, 1], ...]` (mathematica lines 77-82). The Mathematica derivation route then independently confirms via M7 that `P_2 = (N_2 - N_2^{\rm closed})/D_0`.

The M-root and Vieta assertions are stated symbolically the same way in both engines, but that is the natural structure of a quadratic factorization, not a transliteration tell. The two scripts are independent re-derivations, not a line-by-line port. No `mathematica_transliteration` finding.

## Engine cross-check

Both engine outputs end with `STATUS: PASS` after zero-residual reports on the entire assertion family. Side-by-side:

| Claim | SymPy outcome | Mathematica outcome |
|---|---|---|
| one-pole numerator equivalence | "one-pole numerator equivalence" passes (no failure raised in output) | `M1 one-pole numerator residual = 0` |
| K from one-pole | `K_from_one_pole = (B0*B4 + B0*Z4 + 3*B2**2 + 6*B2*MSigma + 6*B2*Z2 + B4*Z0 + 3*MSigma**2 + 6*MSigma*Z2 + Z0*Z4 + 3*Z2**2)/(B4 + Z4)` (output line 11) | `M2 one-pole KSigma residual = 0` |
| K from normalization | `K_from_norm = B0 + Z0 + 5*N0*a**5*c**5*mhat0**2/(54*G*cs**5)` (output line 12) | `M3 normalization KSigma residual = 0` |
| N2/N4 closures | matches `N2_const_closed`, `N4_const_closed` (output lines 18-19) | `M4`, `M5` pass |
| M-root, Vieta | passes | `M10` triple passes |
| u2 on roots | passes | `M11` pair passes |
| Gaussian wall example | `MSigma=sqrt(pi)`, `KSigma=3*sqrt(pi)/2` (output line 33) | `M12 Gaussian wall inertia` and `... stiffness` both pass |
| mutation guards | "constant-prefactor mutation guards = PASS" (output line 22) | `M9 mutated N2 closure guard residual = -(eps/(B0 - KSigma + Z0))`, nonzero as required |

Engines agree. The Mathematica M9 residual `-(eps/(B0-KSigma+Z0))` is the expected nonzero perturbation, consistent with the SymPy mutation guard claim. No `engine_disagreement` finding.

## Verdict justification

The paper card lists four exported equations plus two audit-only checks (response-sign and constant-prefactor). The SymPy and Mathematica scripts both exercise all six items with non-tautological assertions, independent derivation paths for the pole series, and a Mathematica-side explicit `Reals` declaration that the SymPy side under-uses. Attacks tried and failed: (a) checking whether the compatibility assertion at sympy:56 is tautological — it is not, because it depends on the prior solve-based assertions at lines 51-52 that produce K_from_one_pole and K_from_norm independently; (b) checking the one-pole numerator equivalence at sympy:39 for circularity — it is a real structural check that `u_4-4u_2^2` is the discriminant of `den` and matches the M-quadratic; (c) verifying numeric constants `54`, `5`, `5`, `5`, `2` in `P_{0,target}` against paper eq. `stage019-normalization-parent` — exact match; (d) attempting to find a paper deliverable without script coverage — none; (e) attempting to find a script assertion without paper basis — A16/A27 (Gaussian wall integrals) are concrete instantiations explicitly framed as "Concrete wall-integral example" and "M12 Gaussian wall ... integral", not load-bearing, so not a `paper_missing_script_claim`. The one finding raised is a SymPy assumption hygiene issue worth fixing for future-proofing and for second-engine symmetry, but it does not invalidate any current assertion.

## Self-test notes

I checked: (1) variable independence — every `sp.diff` in the Jacobian at lines 78-79 takes derivatives of `D0**2 * P2` and `D0**3 * P4` w.r.t. `N2` and `N4`, both of which genuinely appear in those expressions per the construction at lines 42-45 (so the determinant `D0**3` is a non-trivial structural result, not a hidden 0=0). (2) Symmetry/parity of the Gaussian integrals: `mu_eta * beta^2 = exp(-w^2)` is even and integrates to `sqrt(pi)`; `T_w*(d beta/dw)^2 = w^2*exp(-w^2)` is even and integrates to `sqrt(pi)/2`; `(K_eta + 6*T_omega)*beta^2 = 1*exp(-w^2)` integrates to `sqrt(pi)`; total `KSigma = sqrt(pi)/2 + sqrt(pi) = 3*sqrt(pi)/2` matches the asserted value. (3) Trivial-case substitution for the `M_root_positive` branch using the baseline sample `(B0,Z0,KSigma,B2,Z2,B4,Z4)=(1,2,13,3,4,5,7)`: `D0 = 10`, `B4+Z4 = 12`, `root_gap = sqrt(10*12/3) = sqrt(40)`, `u2_root_positive = sqrt(40)/10 ≈ 0.632` — matches the output value `0.6324555320336759`. The required change for F1 was paper-round-tripped: adding `real=True`/`positive=True` does not change any constant carried into the assertions and does not introduce any new symbol the paper does not also use.
