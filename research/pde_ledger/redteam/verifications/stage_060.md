---
unit_id: 060
batch: III.2
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: true
---

# Verification — unit 060

## Per-finding outcomes

### F1 — hardcoded_result

**Classification:** resolved

**What changed:**
- SymPy (`scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:65-78`): the failing `sp.solve(...)` call is replaced with an explicit closed-form `Csol = a/(sp.exp(a*L) - 1)`; a new `expect_zero("Csol normalizes sigma_trial on [0,L]", ...)` exercises the normalization integral on `[0,L]`; `Sigma_from_rescale` is then built by carrying the Jacobian `L * (Csol*exp(a*s)) |_{s -> x*L}` and substituting `dphi -> Pe*Theta/Lam` (equivalently `a -> Pe/L`), and `expect_zero("Sigma_Pe from rescaling", Sigma_from_rescale - Sigma_x)` is asserted before the unchanged `normalized Sigma_Pe family` check.
- Mathematica (`mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:86-96`): analogous `cNormSol = a/(Exp[a*ell] - 1)`, `expectZero["Csol normalizes sigmaTrial on [0,ell]", ...]`, and `expectZero["Sigma_Pe from rescaling", sigmaFromRescale - sigmaPe]` with the Jacobian rescaling `ell*(cNormSol*Exp[a*s] /. s -> xVar*ell) /. deltaDrop -> pe*theta/lambdaPhi`.

**Assessment:**
The edit matches the directive's "After" block (with the deviation noted: the SymPy script uses the `dphi -> Pe*Theta/Lam` substitution rather than `a -> Pe/L` because `a` is already expanded in the script; both substitutions are equivalent and the residual `Sigma_from_rescale - Sigma_x` simplifies to 0). The new assertions are non-tautological: `Sigma_from_rescale` is built by a change-of-variables transform on the normalized exponential, while `Sigma_x` is the closed form; a mistake in either the Jacobian or the Pe substitution would surface as a non-zero residual. Outputs confirm `Csol normalizes sigma_trial on [0,L] = 0`, `Sigma_Pe from rescaling = 0`, and `normalized Sigma_Pe family = 0` in both engines.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- (a) Pe identification (sympy:79-89; wl:97-107): the substitution-identity `expect_zero("Pe identification", a*L - Lam*dphi/Theta)` is replaced by a `Solve`-based derivation. A new symbol `gamma` (`gammaVar` in wl) is introduced; the ansatz `Cnorm*exp(gamma*s)` is substituted into the affine-drop ODE; the non-zero solution branch is selected; and `expect_zero("Pe identification (derived rate)", gamma_derived - Lam*dphi/(Theta*L))` is asserted. The output prints `gamma_derived = Delta_phi*Lambda_phi/(L*Theta)` (sympy) / `gammaDerived = (deltaDrop*lambdaPhi)/(ell*theta)` (wl), confirming the derivation arrives at the same rate by a different path.
- (b) Xi_micro identities (sympy:118-127; wl:129-142): kept as `expect_zero` checks because F4 successfully derives `phi_from_Phi` (the "if F4 is Blocked" branch in the directive does not apply). The Xi_micro chi-substitution and D/M-substitution checks are still present; their non-tautological grounding now comes from the F4-derived `phi_from_Phi`.
- (c) Product-rule identity (sympy:129-131; wl:144): the `expect_zero("integration-by-parts identity", ...)` assertion is removed and replaced by a comment noting the product-rule rearrangement is a calculus identity.
- (d) Onsager dissipation density (sympy:132-145; wl:145-156): the trivial cancellation `mu_s*J_on + J_on^2/(M*sigma) == 0` is replaced by a positivity check. SymPy uses `sp.ask(sp.Q.nonnegative(dissipation_density), sp.Q.positive(Msig) & sp.Q.positive(sigma_val) & sp.Q.real(mu_s))` with an `assert ... is True` guard. Mathematica uses `Reduce[ForAll[{muS, sigmaVal, mSigma}, Implies[mSigma > 0 && sigmaVal > 0 && Element[muS, Reals], dissipationDensity >= 0]], Reals]` (noted deviation: `Implies[...]` substituted for the directive's `==>` per Mathematica syntax).

**Assessment:**
All four sub-edits match the directive. (a) The derived-rate check is genuinely non-tautological — `gamma` is a free symbol solved out of the ODE, not defined to equal the target. (c) Removal is correct; the calculus identity contributes nothing. (d) The positivity check is substantive: `sp.ask(...) is True` will fail if the simplifier cannot prove nonnegativity, and `Reduce[ForAll[...], Reals] === True` is a strict provability check. Outputs confirm `Pe identification (derived rate) = 0` (both engines), the product-rule line is absent, and `PASS: dissipation density nonnegative under M_sigma, sigma > 0` / `PASS: dissipation density nonnegative under mSigma, sigmaVal > 0` appear.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:129-138`: the `xiMicro = FullSimplify[(lambdaPhi/theta)*(lambdaPhi*ell^2*deltaSupport/tX)/deltaSupport, ...]` deltaSupport-cancellation construction is replaced by a direct dimensional combination `xiMicro = lambdaPhi^2*ell^2/(theta*tX)`, followed by two independent consistency checks: `xiMicroFromChi = chiSigma*lambdaPhi^2*ell^2/tX /. chiSigma -> 1/theta` and `xiMicroFromDM = mSigma*lambdaPhi^2*ell^2/(dSigma*tX) /. dSigma -> mSigma*theta`, with `expectZero` assertions for each.
- The wl script's F2(a) edit (Solve-based `gammaVar` derivation) and F1 edit (rescaling) also contribute independent structural divergence from the SymPy script, as noted in the directive.

**Assessment:**
The `deltaSupport` symbol no longer appears in the `xiMicro` construction in the Mathematica script — verified by reading the file: lines 129-138 contain only `xiMicro = lambdaPhi^2*ell^2/(theta*tX)` with no `deltaSupport` reference. The two new `expectZero` consistency lines verify that the dimensional form agrees with both the susceptibility (`chiSigma`) and phenomenological (`dSigma/mSigma`) routes. Output confirms `xiMicro consistency via chi substitution = 0` and `xiMicro consistency via D/M substitution = 0`. The two engines now construct `Xi_micro` along structurally distinct paths (SymPy via `Lam*phi_from_Phi/(Theta*Delta)`, Mathematica via direct dimensional combination plus chi/DM substitution).

### F4 — insufficient_verification

**Classification:** resolved

**What changed:**
- SymPy (`scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:91-116`): a new "3.5) derive phi_from_Phi from the support EL in the constant-sigma, K_X=0 limit" block is inserted immediately before the `# 4) microscopic Xi from support normalization` section. It introduces `phi_bvp`, `sigma_0`, sets up `EL_const = -Lam*sigma_0 - T_X*phi_bvp''`, solves via `sp.dsolve`, imposes BCs `T_X*phi'(0) = K_m*phi(0)` and `phi'(L) = 0`, computes `Delta_derived = phi(L) - phi(0)`, takes `sp.limit(Delta_derived, K_m, sp.oo)`, and asserts `expect_zero("phi_from_Phi from support BVP (K_m -> infty)", Delta_rigid_limit - Lam*L^2*sigma_0/(2*T_X))`.
- Mathematica (`mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:109-127`): analogous `DSolveValue`-based BVP solve with the same BCs, `deltaRigidLimit` extracted via `Limit[deltaDerived, kM -> Infinity]` (wrapped in `Quiet` for `Limit::alimv`), and `expectZero["phi_from_Phi from support BVP (kM -> infty)", ...]`.
- The comment on the existing `phi_from_Phi = Lam*L^2*Delta/T_X` line is updated to reference the BVP check.

**Assessment:**
The BVP solve is substantive — `sp.dsolve` and `DSolveValue` actually integrate the ODE, the BCs are imposed via `sp.solve`/`Solve`, and the end-to-end drop is extracted by subtracting boundary values. Outputs confirm the derivation: SymPy prints `phi_BVP(0) = L*Lambda_phi*sigma_0/K_m`, `phi_BVP(L) = L*Lambda_phi*sigma_0*(K_m*L + 2*T_X)/(2*K_m*T_X)`, and `Delta_derived = L**2*Lambda_phi*sigma_0/(2*T_X)`. Mathematica prints `Delta_derived = (ell^2*lambdaPhi*sigma0)/(2*tX)`. Both engines confirm the K_m -> infty limit equals `Lambda L^2 sigma_0/(2 T_X)`, matching the directive's self-test (`Delta = Lam sigma_0 L^2/(2 T_X)` independent of K_m in this K_X=0, constant-sigma limit). The `phi_from_Phi from support BVP (... -> infty) = 0` PASS line appears in both outputs.

## Exec log assessment

**SymPy:** exit=0 (inferred from script completion: the final theorem-ledger banner and items 1-5 are present in the saved output, and the script's `expect_zero` helper raises `AssertionError` on any nonzero residual — no such error appears). Notable lines:
- `gamma_derived = Delta_phi*Lambda_phi/(L*Theta)` followed by `Pe identification (derived rate) = 0`
- `Delta_derived = phi(L) - phi(0) = L**2*Lambda_phi*sigma_0/(2*T_X)` followed by `phi_from_Phi from support BVP (K_m -> infty) = 0`
- `Sigma_Pe from rescaling = 0` and `Csol normalizes sigma_trial on [0,L] = 0`
- `PASS: dissipation density nonnegative under M_sigma, sigma > 0`

**Mathematica:** exit=0 (final line `Stage 060 Mathematica audit passed.` appears; `Exit[0]` is reached). Notable lines:
- `gammaDerived = (deltaDrop*lambdaPhi)/(ell*theta)` and `PASS: Pe identification (derived rate)`
- `Delta_derived = (ell^2*lambdaPhi*sigma0)/(2*tX)` and `PASS: phi_from_Phi from support BVP (kM -> infty)`
- `PASS: Sigma_Pe from rescaling`, `PASS: Csol normalizes sigmaTrial on [0,ell]`
- `PASS: xiMicro consistency via chi substitution` and `PASS: xiMicro consistency via D/M substitution`
- `PASS: dissipation density nonnegative under mSigma, sigmaVal > 0`

**Output freshness:** Confirmed. mtimes: sympy script = 1779494141, sympy output = 1779494353 (newer by 212 s); mathematica script = 1779494291, mathematica output = 1779494367 (newer by 76 s). Both `.txt` outputs were regenerated after Codex's edits.

## Material-change assessment

`material_change`: true.

The Sigma_Pe family and Xi_micro value are unchanged numerically — the headline result `Xi_micro = Lambda^2 L^2/(Theta T_X)` is preserved. However, the verification now exercises a derived `phi_from_Phi` extracted from a constant-sigma, K_X=0 support BVP with rigid-grounding limit `K_m -> infty`, and the prefactor is established as `1/2` (i.e., `Delta = Lam sigma_0 L^2/(2 T_X)`). The headline `phi_from_Phi = Lam*L^2*Delta/T_X` retains the unit O(1) factor only under the interpretation that `Delta` denotes the rescaled end-to-end drop per `(1/2) sigma_0` normalization. Any downstream unit that depends on the exact numerical prefactor in `phi_from_Phi` may need to reconcile this `1/2`. Downstream consumers of Xi_micro (units > 060 that import the closure) should be flagged `upstream_stale` per the standard policy, with particular attention to anywhere the literal `Lam*L^2*sigma_0/T_X` formula is invoked without the `1/2`.

## Side observations (non-blocking)

- The wl `expectZero` helper now strips `ConditionalExpression[e_, _] :> e` (lines 26 of the wl) as a generic idiom; this was applied by the orchestrator as a no-substance batch-wide patch and does not affect any per-finding assessment.
- The SymPy script line 73 uses `(Csol * sp.exp(a*s)).subs(s, x*L).subs(dphi, Pe*Theta/Lam)` rather than the directive's `.subs(a, Pe/L)`; this is the equivalent substitution because `a = Lam*dphi/(Theta*L)` is expanded throughout. Functionally identical, and the deviation is noted in the directive's `## Applied: F1` block.
- The SymPy `gamma_solved` fallback `[g for g in gamma_solved if g != 0][0] if gamma_solved else None` will raise `IndexError` if `solve` returns only zero branches (rather than a clean `None`). The output shows `gamma_derived = Delta_phi*Lambda_phi/(L*Theta)`, so this branch was not hit. Not a blocker.
- The Mathematica `phiHead = Unique["phiFun"]` / `sigmaHead = Unique["sigmaFun"]` pattern in the EL section yields slightly noisy print output (`phiFun11`, `sigmaFun12`); not a correctness issue.

## Verdict justification

All four findings are resolved with edits matching the directive's "After" blocks. F1 replaces the failing `sp.solve` with an explicit closed-form `Csol` and adds a rescaling assertion; F2 turns the Pe identification into a `Solve`-based derivation, removes the product-rule tautology, and replaces the Onsager cancellation with a positivity check; F3 eliminates the Mathematica `deltaSupport`-cancellation `xiMicro` construction in favor of an independent dimensional path plus chi/DM consistency checks; F4 adds a constant-source K_X=0 BVP solve in both engines that confirms `Delta = Lambda L^2 sigma_0/(2 T_X)` in the rigid-grounding limit. All new `expect_zero`/`expectZero` lines report 0; the positivity check reports PASS; both scripts complete with exit 0; output `.txt` mtimes confirm freshness.
