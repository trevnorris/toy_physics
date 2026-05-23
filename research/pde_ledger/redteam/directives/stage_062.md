---
unit_id: 062
batch: III.3
created_at: 2026-05-22T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-22T19:39:41-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 062

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:46-81`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:42-68`

**Issue:**
Three `expect_zero` assertions in each engine (`G_micro - G_expected`, the coherence-factor decomposition, and `Xi_micro - Xi_expected`) consist of defining a quantity and then comparing it to a hand-typed algebraic re-write of itself. They cannot fail under any choice of symbol values, so they verify nothing about the parent-action projection. The Mathematica file mirrors the same pattern.

**Required change:**

Replace the three tautological blocks in BOTH engines with an explicit parent-action derivation. The substantive content the unit is supposed to verify is the chain "parent quadratic action + heavy-channel elimination → effective phi gain `G_micro = rho_star g_phi^2 Osp^2 / (m cs_star_sq KX Nss)`". Implement it as a real derivation.

For the SymPy file (`scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py`), edit lines 46-81 to:

1. **Set up a parent quadratic action** symbolically. Introduce two channel variables `sigma` and `phi` and write:
   ```python
   sigma, phi = sp.symbols("sigma phi", real=True)
   S_parent = sp.Rational(1, 2) * Theta_sigma * sigma**2 + Lambda_phi * sigma * phi + sp.Rational(1, 2) * KX * phi**2
   ```
   (Use the existing definitions `Theta_sigma = (m * cs_star_sq / rho_star) * Nss` and `Lambda_phi = g_phi * Osp` that were already declared at lines 51-52; keep those declarations.)

2. **Eliminate sigma via the stationary condition** and substitute back:
   ```python
   sigma_star = sp.solve(sp.diff(S_parent, sigma), sigma)[0]
   S_eff_phi = sp.expand(S_parent.subs(sigma, sigma_star))
   gain_from_action = -2 * sp.simplify(S_eff_phi.coeff(phi, 2))   # sign per parent-action convention
   ```
   Note: The sign in `gain_from_action` is `-2 * coefficient of phi^2` because eliminating sigma at its minimum subtracts `Lambda_phi^2/(2 Theta_sigma) * phi^2` from the parent, leaving `(KX/2 - Lambda_phi^2/(2 Theta_sigma)) phi^2`. The microscopic *gain* (susceptibility to phi when sigma is integrated out) is conventionally `Lambda_phi^2/(KX * Theta_sigma)`. Use whichever sign convention matches the script's previous `G_micro` formula — verify by substituting the existing closed-form values for `Theta_sigma`, `Lambda_phi`, `KX` and checking the result reduces to `rho_star g_phi^2 Osp^2 / (m cs_star_sq KX Nss)`. Adjust the sign/prefactor in `gain_from_action` to match the script's stated convention.

3. **Assert that the action-derived gain matches the closed-form `G_micro`:**
   ```python
   G_micro_closed = rho_star * g_phi**2 * Osp**2 / (m * cs_star_sq * KX * Nss)
   expect_zero("G_micro from parent action vs closed form",
               gain_from_action - G_micro_closed)
   ```
   This is a non-tautological check: `gain_from_action` comes from `sp.solve` and substitution, while `G_micro_closed` is the claimed formula. A wrong sign on `Lambda_phi`, a wrong factor in `Theta_sigma`, or a wrong exponent on `Osp` would now make the assertion fail.

4. **Drop the coherence-factor check.** `C2 := Osp^2 / (Nss Npp)` is a definition, not a theorem. Replace the existing block at lines 62-71 with a single print line documenting the definition, no `expect_zero`:
   ```python
   print("Coherence factor (definition):  C_(sigma phi)^2 := Osp^2 / (Nss Npp)")
   ```

5. **Replace the Xi_micro tautology with a dimensional-substitution check.** `kappa = KX L^2 / TX` is presented as a definition; instead, introduce `kappa` as an independent symbol up front and assert that the *only* substitution making `Xi_micro` reduce to the `KX, L, TX, Osp, ...` form is `kappa = KX L^2 / TX`:
   ```python
   Xi_micro = sp.simplify(kappa * gain_from_action)
   Xi_target = rho_star * g_phi**2 * Osp**2 * L**2 / (m * cs_star_sq * TX * Nss)
   kappa_solved = sp.solve(Xi_micro - Xi_target, kappa)
   assert kappa_solved == [KX * L**2 / TX], f"Unexpected kappa solution: {kappa_solved}"
   print("kappa solved from Xi_micro = Xi_target:", kappa_solved[0])
   ```
   This is substantive: if `Xi_target` were wrong, `kappa_solved` would not equal `KX L^2/TX`.

Keep the docstring (lines 2-11) and the n=5 EOS check (lines 33-44) for now; F2 addresses those.

For the Mathematica file (`mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl`), edit lines 42-68 to perform the **same physical derivation** but via Mathematica's `Solve` and `Coefficient` primitives, NOT a line-by-line copy. Concretely:

1. Build the parent action `sParent = (1/2)*thetaSigma*sigma^2 + lambdaPhi*sigma*phi + (1/2)*kX*phi^2`.
2. Use `sigmaStar = sigma /. First[Solve[D[sParent, sigma] == 0, sigma]]`.
3. Compute `sEff = Expand[sParent /. sigma -> sigmaStar]`.
4. Extract `gainFromAction = -2*Coefficient[sEff, phi, 2]` (adjust sign to match SymPy's convention so both engines compute the same gain).
5. Assert `expectZero["gMicro from parent action vs closed form", gainFromAction - rhoStar*gPhi^2*oSP^2/(csStarSq*kX*m*nSS)]`.
6. Drop the coherence-factor `expectZero` (line 63); replace with a `Print` statement documenting the definition only.
7. For `xiMicro`, solve `Solve[xiMicro == xiTarget, kappa]` and assert the solution equals `kX*ell^2/tX` via an `If[...] === True` test, not by re-substituting and comparing.

The Mathematica file must use `Solve` and `Coefficient` primitives that are NOT used in the SymPy file's matching block; if SymPy uses `sp.solve` and `.coeff`, Mathematica's `Solve` and `Coefficient` count as the corresponding-but-distinct calls. The key is the algebraic *path* must differ — both arrive at the same gain via independent symbolic chains. (See F3 for the broader transliteration concern.)

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 062` and `redteam exec-mathematica 062`. Both scripts must exit 0, both must contain the new `expect_zero` / `expectZero` call comparing the action-derived gain to the closed-form, and the SymPy file must contain a `kappa_solved` block (or equivalent) replacing the old Xi tautology.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl`
- summary: Replaced tautological gain, coherence, and Xi checks with parent-action stationary elimination, a definition-only coherence print, and a solved kappa dimensional check.
- deviation: none

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:33-44`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:28-40`

**Issue:**
The "n=5 EOS identity" check `h'(rho) - m c_s^2 / rho == 0` is a degree-4-monomial algebraic identity that holds for any `h(rho) = A rho^4` when `c_s^2` is defined as `(4/m) h` (i.e., as `(1/m) d(K rho^5)/drho`). It is not specific to n=5 and does not exercise the polytropic index. A wrong polytropic exponent that propagates consistently through `U` and `c_s^2` would not be detected.

**Required change:**

For the SymPy file, modify lines 33-44 so that the polytropic index appears as a separate symbol whose mismatch can break the assertion. Specifically:

1. Introduce a symbol `n_poly` (the polytropic index) and rewrite `U` as `K * rho**n_poly / n_poly`:
   ```python
   n_poly = sp.symbols("n_poly", positive=True, integer=True)
   U_general = K * rho**n_poly / n_poly
   h_general = sp.diff(U_general, rho)
   hprime_general = sp.diff(h_general, rho)
   ```
2. Define `c_s^2` via the standard thermodynamic chain from a pressure `p(rho)`. For a barotropic monomial EOS, `p = (n_poly - 1) K rho**n_poly / n_poly` (or the convention this unit uses — pick one and state it inline above the code), giving `c_s^2 = dp/d(m rho) = (n_poly - 1) K rho**(n_poly - 1) / m`:
   ```python
   p_general = (n_poly - 1) * K * rho**n_poly / n_poly
   cs_sq_general = sp.diff(p_general, rho) / m   # c_s^2 = dp/d(m*rho) with d(m*rho)/drho = m
   ```
3. Verify the identity `h'(rho) = m c_s^2 / rho` in general:
   ```python
   expect_zero("h'(rho) = m c_s^2 / rho (general polytrope)",
               hprime_general - m * cs_sq_general / rho)
   ```
   This residual evaluates `(n_poly-1) K rho^(n_poly-2) - (n_poly-1) K rho^(n_poly-2)` symbolically, which simplifies to zero ONLY when the pressure-EOS relation `p = (n-1)U` is consistent with `h = U'`. Verify by mental substitution: `h = K rho^(n-1)`, `h' = (n-1) K rho^(n-2)`; `c_s^2 = (n-1) K rho^(n-1)/m`, `m c_s^2/rho = (n-1) K rho^(n-2)`. Equal. Good.
4. Then specialize to n=5 and reproduce the original check:
   ```python
   subs_n5 = {n_poly: 5}
   print("Specializing to n=5:")
   print("  U(rho) =", sp.simplify(U_general.subs(subs_n5)))
   print("  h(rho) =", sp.simplify(h_general.subs(subs_n5)))
   print("  h'(rho) =", sp.simplify(hprime_general.subs(subs_n5)))
   print("  c_s^2(rho) =", sp.simplify(cs_sq_general.subs(subs_n5)))
   expect_zero("n=5 specialization of h' = m c_s^2/rho",
               (hprime_general - m * cs_sq_general / rho).subs(subs_n5))
   ```

Note: the original script had `U = K rho^5/4`, i.e., used `K/4` not `K/5`. That's a *different* normalization than the standard `K rho^n / n`. Inspect the existing line 35 carefully:

```python
U = K * rho**5 / 4
```

The `/4` (not `/5`) means the original `U` is NOT the standard polytropic form. Codex should preserve the existing normalization choice: rewrite as `U_general = K * rho**n_poly / (n_poly - 1)` and recompute the pressure consistently — verify the algebra so the residual still vanishes for general `n_poly`. Mentally: with `U = K rho^n / (n-1)`, `h = U' = n K rho^(n-1)/(n-1)`. For n=5, `h = 5 K rho^4 / 4`, matching the original. Then `h' = n(n-1) K rho^(n-2)/(n-1) = n K rho^(n-2)`. For c_s^2, the original used `c_s^2 = (1/m) d(K rho^n)/drho = n K rho^(n-1)/m`, so `m c_s^2/rho = n K rho^(n-2) = h'`. So the identity holds for the existing normalization with `U = K rho^n/(n-1)` and `c_s^2 = (1/m) d(K rho^n)/drho`, for general n. Codex should use that form.

The substantive content is now: the identity holds because the script's `c_s^2` and `h'` are linked by the same polytropic structure. The check becomes non-trivial when an inconsistent index is fed to one definition but not the other — Codex should add a final sanity check:
```python
# Sanity: identity breaks if c_s^2 uses a different exponent
cs_sq_wrong = sp.diff(K * rho**(n_poly + 1), rho) / m
residual_wrong = (hprime_general - m * cs_sq_wrong / rho).subs(subs_n5)
assert sp.simplify(residual_wrong) != 0, "Inconsistency check failed to detect wrong exponent"
print("Inconsistency probe (n+1 in c_s^2 only):", sp.simplify(residual_wrong))
```
This confirms the assertion machinery would catch an exponent mismatch.

For the Mathematica file, perform the analogous edit at lines 28-40: introduce `nPoly`, define `uGeneral = capitalK*rho^nPoly/(nPoly - 1)`, etc., verify the general identity via `expectZero`, then specialize to `nPoly -> 5`, and add the corresponding inconsistency probe using `If[FullSimplify[residualWrong /. nPoly -> 5] === 0, fail[...], Print["inconsistency probe nonzero: ", ...]]`.

**Verification command:**
The verifier confirms the SymPy file contains the symbol `n_poly` declared as `sp.symbols(...)`, that the EOS identity is asserted for general `n_poly` (not pre-substituted), and that an inconsistency probe (`residual_wrong` or equivalently named) appears and is asserted nonzero. The Mathematica file must have the analogous structure with `nPoly`. Both scripts must exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl`
- summary: Generalized the EOS identity to a symbolic polytropic index, preserved the stage normalization at n=5, and added wrong-exponent inconsistency probes.
- deviation: none

## F3 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:1-73`

**Issue:**
The Mathematica file mirrors the SymPy file's exact algebraic choreography (same variable order, same intermediate steps, same `expectZero` sequence, same comparisons against hand-typed twins). It does not derive the result independently. This violates the second-engine policy and makes both engines blind to the same potential errors in lockstep.

**Required change:**

After F1 and F2 have been applied, the Mathematica file must derive its gain via a *different algebraic path* than the SymPy file. Two options; choose one and apply consistently:

**Option A — Symbolic-extremization path:** Use `Solve` and `Coefficient` (as outlined in F1) to extract `gainFromAction`. The SymPy file uses `sp.solve` + `.coeff(phi, 2)`. The Mathematica file uses `Solve` + `Coefficient`. These are corresponding-but-distinct primitives. To further differentiate, the Mathematica file should additionally verify the result by computing the gain a second way: via `Series` expansion around `phi = 0` of the integrated-out action, reading the quadratic coefficient. So the Mathematica file performs the derivation *twice* and checks both routes match, where SymPy only does it once. Concretely:

```
sEff = Expand[sParent /. sigma -> sigmaStar];
gainFromAction = -2*Coefficient[sEff, phi, 2];
gainFromSeries = -2*SeriesCoefficient[Series[sEff, {phi, 0, 2}], 2];
expectZero["Mathematica two-route consistency", gainFromAction - gainFromSeries];
expectZero["gMicro from parent action vs closed form",
           gainFromAction - rhoStar*gPhi^2*oSP^2/(csStarSq*kX*m*nSS)];
```

**Option B — Numerical cross-check path:** In addition to the symbolic derivation (which can mirror SymPy's), add a numerical sweep where Mathematica randomly samples positive rational values for all parameters and verifies the closed-form `G_micro` agrees with the action-derived gain to high precision:

```
nTrials = 20;
SeedRandom[42];
trials = Table[
  With[{vals = {rhoStar -> RandomReal[{0.1, 10}, WorkingPrecision -> 50], ...}},
    Chop[N[(gainFromAction - rhoStar*gPhi^2*oSP^2/(csStarSq*kX*m*nSS)) /. vals, 30]]
  ],
  {nTrials}
];
If[Max[Abs[trials]] > 10^-20, fail["Numerical sweep failed", trials], pass["20 random trials agree"]];
```

(SymPy can stay purely symbolic; Mathematica adds numerical evidence.)

Pick Option A if Codex wants pure symbolic-only verification; Option B if numerical confirmation is acceptable. Either way, the .wl file must contain at least one call (`SeriesCoefficient`, `Series`, `RandomReal`, `N`, `Chop`, or numerical equivalent) that is NOT present in the .py file.

Apply this change AFTER F1 lands. If F1 introduces the parent-action derivation in both files, F3 ensures the two derivations diverge structurally.

**Verification command:**
After Codex applies, the verifier diffs the structural list of operations in `.py` vs. `.wl`. The `.wl` must contain at least one operation type not present in the `.py` (e.g., `Series`, `SeriesCoefficient`, `RandomReal`, `Chop`, `Table` over numerical trials). The Mathematica script must exit 0 and the new check(s) must appear in its output.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl`
- summary: Added a Mathematica-only series-coefficient route that cross-checks the action-derived gain before comparing it to the closed form.
- deviation: none
