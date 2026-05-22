---
unit_id: 038
batch: III.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 038 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage038_dimensionless_continuum_placement_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage038_dimensionless_continuum_placement_mathematica_audit.txt`

## What the script claims to verify

Per the SymPy docstring, the stage "compress[es] the Stage-20 continuum formulas into the exact dimensionless kernel map, verif[ies] the product relation, and factor[s] the one-way parameter tendencies." Concretely:
1. Substituting the dimensionless ratios (`eps_eta = c_etaU^2/(K_U Keta_eff)`, `eps_W = sigma c_UW^2/(K_U KWeff)`, `Z_W = c_etaW^2/(Keta_eff KWeff)`, `rho = c_UW c_etaU/(K_U c_etaW)`, `delta0 = pi^2 T_w/(L^2 Keta_eff)`, and `Lambda = 27 pi^2 c_s^5 G KWeff/(20 a^5 c_light^5 mu_W)`) reduces the four operators `delta`, `M_mix`, `R_target`, `beta0` to closed forms in those ratios.
2. The product `R_target * M_mix = 8 Lambda (1 - eps_W)/pi^2`.
3. The nine first-order partial derivatives of the closed forms factor into specific manifestly-signed forms ("one-way parameter tendencies").

## Assertion inventory

| #   | Script      | Line  | Form                                                                                                          | Anchored to claim? |
| --- | ----------- | ----- | ------------------------------------------------------------------------------------------------------------- | ------------------ |
| A1  | sympy       | 109   | `expect_zero("delta - delta0/(1-eps_eta)", ...)`                                                              | yes                |
| A2  | sympy       | 110-113 | `expect_zero("M_mix - 8 Z_W (1+rho)^2/[pi^2 (1-eps_eta)(1-eps_W)]", ...)`                                   | yes                |
| A3  | sympy       | 114-117 | `expect_zero("R_target - Lambda (1-eps_eta)(1-eps_W)^2/[Z_W (1+rho)^2]", ...)`                                | yes                |
| A4  | sympy       | 118-121 | `expect_zero("beta0 - (mu_W/mu_eta)(Keta_eff/KWeff) Z_W (1+rho)^2/(1-eps_W)^2", ...)`                         | yes                |
| A5  | sympy       | 135   | `expect_zero("R_target M_mix - 8 Lambda (1-eps_W)/pi^2", ...)`                                                | yes                |
| A6  | sympy       | 141-143 | `expect_zero("8 Lambda (1-eps_W)/pi^2 - NQ * KWeff(1-eps_W)/mu_W", ...)`                                    | partial (consistency of Lambda's definition) |
| A7-A15 | sympy    | 167-175 | nine `expect_zero` calls for d{delta,M,R}/d{eps_eta,eps_W,Z_W,rho} factorizations                           | partial (factorization only; sign claim unchecked) |
| B1  | mathematica | 106   | `expectZero["delta - delta0/(1-epsEta)", ...]`                                                                | yes                |
| B2  | mathematica | 107-110 | `expectZero["M_mix - 8 Z_W (1+rho)^2/[Pi^2 (1-epsEta)(1-epsW)]", ...]`                                      | yes                |
| B3  | mathematica | 111-114 | `expectZero["R_target - lambda (1-epsEta)(1-epsW)^2/[Z_W (1+rho)^2]", ...]`                                  | yes                |
| B4  | mathematica | 115-118 | `expectZero["beta0 - (muW/muEta)(KetaEff/KWEff) Z_W (1+rho)^2/(1-epsW)^2", ...]`                              | yes                |
| B5  | mathematica | 128   | `expectZero["R_target M_mix - 8 lambda (1-epsW)/Pi^2", ...]`                                                  | yes                |
| B6  | mathematica | 129-137 | `expectZero["8 lambda (1-epsW)/Pi^2 - NQ KWeff(1-epsW)/muW", ...]`                                          | partial            |
| B7-B15 | mathematica | 157-165 | nine `expectZero` calls for derivative factorizations                                                      | partial            |

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage038_dimensionless_continuum_placement_mathematica_audit.wl:1-175`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py:1-185`

**What's wrong:**
The `.wl` script is structurally a line-by-line port of the `.py` script: same four banner sections in the same order, same `sigma = 88/(9 pi^2)`, same `kEtaEff = kEta + 6 tOmega`, same `kWEff = kW + Pi^2 tW/(4 ell^2)`, same `A = (kU*kEtaEff - cEtaU^2)/(muEta*kU)`, same builds of `delta`, `mMix`, `beta0`, `rTarget`, and the same 15 assertion targets in the same order with the same exact phrasing (banner strings like `"d M / d epsW factorization"` match `"d M / d eps_W factorization"` modulo underscore). The closing commentary block ("On the natural nonvanishing transfer branch 1+rho > 0") is reproduced verbatim.

Concrete corresponding sections:

SymPy (lines 50-58):
```
sigma = sp.Rational(88, 9) / sp.pi**2
Keta_eff = sp.simplify(K_eta + 6 * T_Omega)
KWeff = sp.simplify(K_W + sp.pi**2 * T_W / (4 * L**2))
A = sp.simplify((K_U * Keta_eff - c_etaU**2) / (mu_eta * K_U))
delta = sp.simplify(sp.pi**2 * T_w * K_U / (L**2 * (K_U * Keta_eff - c_etaU**2)))
M_mix = sp.simplify(8 * (K_U * c_etaW + c_UW * c_etaU) ** 2 / (sp.pi**2 * (K_U * Keta_eff - c_etaU**2) * (K_U * KWeff - c_UW**2 * sigma)))
```

Mathematica (lines 43-50) — same formulas, same variable choreography:
```
sigma = 88/(9 Pi^2);
kEtaEff = FullSimplify[kEta + 6 tOmega, ...];
kWEff = FullSimplify[kW + Pi^2 tW/(4 ell^2), ...];
a = FullSimplify[(kU*kEtaEff - cEtaU^2)/(muEta*kU), ...];
delta = FullSimplify[Pi^2 tw*kU/(ell^2 (kU*kEtaEff - cEtaU^2)), ...];
mMix = FullSimplify[8 (kU*cEtaW + cUW*cEtaU)^2/(Pi^2 (kU*kEtaEff - cEtaU^2) (kU*kWEff - cUW^2 sigma)), ...];
```

Furthermore, the Mathematica `applyDimless` (lines 77-99) introduces a rule that bakes the *answer* into the substitution rather than independently deriving it:
```
(cEtaU*cUW + cEtaW*kU)^2 -> zW kEtaEff kWEff kU^2 (1 + rho)^2
```
This is the target factored form; the SymPy script substitutes only the atomic monomials (`c_etaU^2`, `c_UW * c_etaU`, `c_etaW^2`) and lets SymPy expand. The `.wl` script's wholesale-substitution shortcut is the giveaway of a port engineered after the SymPy answer was known, not an independent derivation from the same physical premises.

**Why this matters:**
The two-engine policy demands that SymPy and Mathematica each derive the result independently, so that an agreement constitutes genuine cross-validation. A transliteration provides no information beyond the SymPy result. If the SymPy script had a sign error (or an incorrect substitution rule), the transliterated Mathematica would inherit it.

**Required change:**
Refactor the Mathematica script so it derives the dimensionless-kernel reductions from independent algebra, with at least the following differentiations from the SymPy approach:
1. Do not pre-introduce the squared-numerator-as-product rule `(cEtaU*cUW + cEtaW*kU)^2 -> zW kEtaEff kWEff kU^2 (1+rho)^2`. Substitute only the atomic monomials (`cEtaU^2`, `cUW*cEtaU`, `cEtaW^2`, `tw`, `gNewton`) and let `FullSimplify`/`Together` assemble the result.
2. Express the dimensionless ratios as *symbol introductions* via `Solve` or explicit definitions from their inverse forms (e.g., `cEtaU^2 = epsEta kU kEtaEff` solved as `epsEta = cEtaU^2/(kU kEtaEff)`), and then verify the dimensionless closed forms by substituting back.
3. Drop the cloned section banners ("STAGE 21" / "Stage-20 continuum formulas" / "Exact derivative factors") if you want to keep a parallel structure; otherwise restructure the Mathematica script's order of operations (e.g., do the derivative factorization first using the closed forms, then verify the closed forms reduce from the continuum operators — the reverse of the SymPy order).

The minimum substantive change is item (1): remove the answer-baked substitution rule for `(cEtaU*cUW + cEtaW*kU)^2`, and rely on the atomic substitutions for `cUW*cEtaU`, `cEtaW^2`, `cEtaU^2`, `cUW^2`. After that change, the same `expectZero` checks must still pass (else the original SymPy result was unsupported and a separate finding applies).

**Verification:**
After the edit, line containing `(cEtaU*cUW + cEtaW*kU)^2 -> zW kEtaEff kWEff kU^2 (1 + rho)^2` (and its `^(-2)` partner) should be removed from `applyDimless`. The four `expectZero` checks in section 2 must still print `PASS:` in the regenerated `.txt`. If they fail, the rule was carrying load — meaning the atomic substitutions alone are not sufficient and the script's claim relied on the pre-baked answer.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py:153-180`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage038_dimensionless_continuum_placement_mathematica_audit.wl:143-170`

**What's wrong:**
Section 4 of both scripts claims to "factor the one-way parameter tendencies" and prints commentary:
```
print("On the natural nonvanishing transfer branch 1+rho > 0:")
print("  delta increases with eps_eta")
print("  M_mix increases with eps_eta, eps_W, Z_W, rho")
print("  R_target decreases with eps_eta, eps_W, Z_W, rho")
```
But no assertion verifies these sign claims. The nine `expect_zero` calls (sympy 167-175) only check that SymPy's symbolic derivative of the manually-typed closed form equals a manually-typed factored expression. For example, line 167 asserts `sp.diff(delta0/(1-eps_eta), eps_eta) - delta0/(1-eps_eta)**2 == 0`. That is a routine calculus identity for SymPy; it does not address the physical "one-way" claim.

Moreover, the symbol declarations (sympy line 80: `eps_eta, eps_W, rho, Z_W, delta0, Lambda` declared `positive=True, real=True`) do not impose `eps_eta < 1` or `eps_W < 1`. The factored derivative forms have `(1-eps_eta)` and `(1-eps_W)` in denominators; the sign claim in the commentary requires both to be in `(0, 1)`. Without that constraint, e.g. `dM_deps_eta = 8 Z_W (1+rho)^2/[pi^2 (1-eps_eta)^2 (1-eps_W)]` is not manifestly positive: if `eps_W > 1`, it is negative.

**Why this matters:**
The script's docstring and commentary present the sign of each derivative as a stage result. As written, the assertions verify factorization identities (which can fail only via author typo) but do not verify the stated tendencies. A reader (or downstream user) consuming the stage's PASS status as evidence for the tendency claim is misled.

**Required change:**
After the existing factorization assertions in both scripts, add explicit sign assertions under the stated transfer-branch assumption (`0 < eps_eta < 1`, `0 < eps_W < 1`, `1 + rho > 0`, `Z_W > 0`, `Lambda > 0`, `delta0 > 0`). Two acceptable approaches:

Approach A (preferred): tighten the symbol declarations to include `eps_eta < 1` and `eps_W < 1`, then assert each derivative's sign via `sp.refine` / `sp.simplify` on `sp.Abs` of a normalized form. Concretely, for each derivative `D`, factor it as `D = sign_factor * positive_template`, and assert `expect_zero` for the difference `D / positive_template - expected_sign`, where `positive_template` is a manifestly positive product (e.g. `Z_W*(1+rho)**2/((1-eps_eta)**2*(1-eps_W))`) and `expected_sign` is `+1` or `-1`.

Approach B (lighter): add nine sign-marker assertions, one per derivative, of the form
```python
expect_zero("sign(d delta / d eps_eta)", sp.sign(d_delta_deps_eta) - 1)
```
under a SymPy `Q.positive` refinement context that imposes `1 - eps_eta > 0`, `1 - eps_W > 0`, etc. (Note: `sp.sign` will not auto-resolve without those assumptions; use `sp.refine(sp.sign(expr), Q.positive(1-eps_eta) & ... )`.)

Pick whichever the script's author tooling supports cleanly. For Mathematica, the parallel is `Simplify[Sign[...], 0 < epsEta < 1 && 0 < epsW < 1 && rho > -1 && zW > 0 && lambda > 0 && delta0 > 0] == 1` (or `-1` for the R derivatives), asserted via `expectZero` against the expected sign.

**Verification:**
After the patch, sympy script lines ~167-175 should still be present; immediately after them, nine new `expect_zero` lines (or one consolidated block) should appear asserting the sign of each derivative under the branch assumptions. The output `.txt` should print nine new lines (e.g., `sign(d delta / d eps_eta) - 1 = 0`, `sign(d R / d eps_eta) + 1 = 0`). Exit code remains 0. The Mathematica script gets the parallel addition with the analogous output.

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration of the SymPy script — see F1 for details. The two scripts share section structure, banner text, formula choreography, variable definition order, and the exact assertion set. The Mathematica `applyDimless` further bakes the squared-numerator answer into a substitution rule (`(cEtaU*cUW + cEtaW*kU)^2 -> zW kEtaEff kWEff kU^2 (1 + rho)^2`), which is the target factored form. This is not an independent derivation; it is algebra engineered backward from the known SymPy result.

## Engine cross-check

Both engines produce the same final results across all 15 corresponding assertions:

SymPy output (lines 26-33, 38-40, 45-53):
- `delta - delta0/(1-eps_eta) = 0`
- `M_mix - 8 Z_W (1+rho)^2 / [pi^2 (1-eps_eta)(1-eps_W)] = 0`
- `R_target - Lambda (1-eps_eta)(1-eps_W)^2 / [Z_W (1+rho)^2] = 0`
- `beta0 - (mu_W/mu_eta)(Keta_eff/KWeff) Z_W (1+rho)^2/(1-eps_W)^2 = 0`
- `R_target M_mix - 8 Lambda (1-eps_W)/pi^2 = 0`
- nine derivative-factorization zeros

Mathematica output (lines 26-33, 42-45, 51-68): same five primary checks PASS, same nine derivative-factorization checks PASS.

Numerical/symbolic agreement at the residual level: both engines return exactly `0` on every check after `simplify`/`FullSimplify`. No engine disagreement.

The agreement, however, is weakened by the transliteration (F1): the agreement is at the level of the same algorithm running in two CAS, not of two independent derivations converging.

## Verdict justification

The five primary substitution checks (A1-A5/B1-B5) are non-tautological and pass cleanly under inspection: the algebra of substituting `c_etaU^2`, `c_UW*c_etaU`, `c_etaW^2`, `T_w`, and `G` into `delta`, `M_mix`, `R_target`, `beta0` produces the stated closed forms in the dimensionless ratios. The product relation `R_target * M_mix = 8 Lambda (1-eps_W)/pi^2` follows from multiplying the two closed forms, and the `Lambda`-definition check (A6/B6) is a numerical-prefactor consistency check (54/5 = 8*27/20, modulo pi^2 absorbed into the definition) that holds.

Two issues prevent a `clean` verdict:
1. The Mathematica script is a transliteration of the SymPy script, weakening the cross-engine guarantee (F1).
2. Section 4's "one-way parameter tendencies" claim is asserted only as a factorization identity (basic calculus), not as a sign claim under the stated transfer-branch assumption; the commentary about which parameters cause increases vs. decreases is unverified (F2).

Attacks that failed: I tried to break the dimensionless reductions by checking whether the substitution `c_UW*c_etaU -> rho * K_U * c_etaW` mishandles signs (it does not, because every appearance traces to `(K_U c_etaW + c_UW c_etaU)^2`, which is sign-insensitive once squared). I checked whether the `Lambda` definition's numerical prefactor `20/27` cancels correctly against `NQ = 54/5` and the `8/pi^2` factor (it does: `54/5 * 20/27 = 8`). I checked whether the derivative factorizations have sign typos (they do not; each matches elementary calculus). I checked whether `eps_eta -> 1` causes a removable singularity issue (it does not within the scope of the assertions, which never substitute that value). The math is correct; the verification framework is what needs hardening.

Verdict: `findings` with two findings, neither stop-cold. Outputs are fresh (sympy output mtime 2026-05-11 12:41 vs. script mtime 2026-04-01 12:39; mathematica output mtime 2026-05-11 12:48 vs. script mtime 2026-05-11 11:56 — both `.txt` files are newer than their respective scripts).

## Self-test notes

I checked four traps before finalizing the directive:
1. **Variable independence:** The proposed sign-of-derivative additions (F2) use `sp.sign(expr)` or `Simplify[Sign[...]]` on expressions that genuinely depend on the named symbols, so there is no zero-derivative trap.
2. **Symmetry/parity:** No integrals in this stage; not applicable.
3. **Trivial-case substitution:** For F2, substituting `eps_eta=1/2, eps_W=1/2, Z_W=1, rho=0, Lambda=delta0=1` into the nine derivative forms gives finite, nonzero values with the claimed signs (e.g., `dM_deps_eta = 8*1*1/(pi^2 * 1/4 * 1/2) = 64/pi^2 > 0`; `dR_deps_eta = -1*1/4/1 < 0`), confirming the sign assertions are exercising the right physics.
4. **Path specifications:** F1 targets `.wl` file at `/var/projects/toy_physics/research/pde_ledger/mathematica/...`. F2 targets both `.py` at `/var/projects/toy_physics/research/pde_ledger/scripts/...` and `.wl` at `/var/projects/toy_physics/research/pde_ledger/mathematica/...`. Both directories named explicitly in the directive.
