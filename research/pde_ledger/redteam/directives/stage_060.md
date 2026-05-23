---
unit_id: 060
batch: III.2
created_at: 2026-05-22T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-22T17:57:04-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 060

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — hardcoded_result

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:66-74`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:80-81`

**Issue:**
`Sigma_Pe(x) = Pe * exp(Pe*x)/(exp(Pe)-1)` is hand-typed as the normalized exponential family in both engines without derivation. In sympy, the intervening `sp.solve` returns `Piecewise((1/L, Eq(Delta_phi*Lambda_phi, 0)), (nan, True))` for `Csol` (visible in the saved .txt line 22) — i.e. solve picked the degenerate branch and returned nan for the generic case. The downstream integral check (sympy line 73-74, wl line 81) merely verifies the hand-typed form normalizes to 1, which it always does. The actual derivation step (normalizing `Cnorm*exp(a*s)` over `[0,L]` then rescaling to `s = x*L, Pe = a*L`) is not performed.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py`:

Replace the block at lines 65-74 (the failing `sp.solve` and the hand-typed `Sigma_x`) with an explicit symbolic derivation:

Before (lines 65-74):
```
# normalization on [0,L]
Csol = sp.simplify(sp.solve(sp.Eq(sp.integrate(sigma_trial, (s, 0, L)), 1), Cnorm)[0])
print("Normalization constant C =", Csol)
Pe = sp.symbols('Pe', positive=True, real=True)
Csol_Pe = sp.simplify(Csol.subs(a, Pe/L))
x = sp.symbols('x', real=True)
Sigma_x = sp.simplify(Pe * sp.exp(Pe*x) / (sp.exp(Pe) - 1))
print("Sigma_Pe(x) =", Sigma_x)
expect_zero("normalized Sigma_Pe family",
            sp.simplify(sp.integrate(Sigma_x, (x, 0, 1)) - 1))
```

After:
```
# normalization on [0,L] — compute C explicitly (sympy's solve degenerates on a=0)
Csol = a / (sp.exp(a*L) - 1)
print("Normalization constant C =", Csol)
expect_zero("Csol normalizes sigma_trial on [0,L]",
            sp.simplify(sp.integrate(Csol*sp.exp(a*s), (s, 0, L)) - 1))
# rescale to x = s/L with Pe = a*L; carry the Jacobian via Sigma_x dx = sigma(s) ds
Pe = sp.symbols('Pe', positive=True, real=True)
x = sp.symbols('x', real=True)
Sigma_from_rescale = sp.simplify(L * (Csol * sp.exp(a*s)).subs(s, x*L).subs(a, Pe/L))
Sigma_x = sp.simplify(Pe * sp.exp(Pe*x) / (sp.exp(Pe) - 1))
print("Sigma_Pe(x) =", Sigma_x)
expect_zero("Sigma_Pe from rescaling", Sigma_from_rescale - Sigma_x)
expect_zero("normalized Sigma_Pe family",
            sp.simplify(sp.integrate(Sigma_x, (x, 0, 1)) - 1))
```

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl`:

Insert immediately before the existing `sigmaPe = ...` definition at line 80 the analogous derivation. Replace the block at lines 80-81:

Before (lines 80-81):
```
sigmaPe = FullSimplify[pe*Exp[pe*s]/(Exp[pe] - 1), Assumptions -> pe > 0];
expectZero["normalized Sigma_Pe family", Integrate[sigmaPe, {s, 0, 1}] - 1];
```

After:
```
cNormSol = a/(Exp[a*ell] - 1);
expectZero["Csol normalizes sigmaTrial on [0,ell]",
  FullSimplify[Integrate[cNormSol*Exp[a*s], {s, 0, ell}] - 1, Assumptions -> $Assumptions]];
xVar = Symbol["xVar"];
sigmaFromRescale = FullSimplify[
  ell*(cNormSol*Exp[a*s] /. s -> xVar*ell) /. a -> pe/ell,
  Assumptions -> $Assumptions && pe > 0];
sigmaPe = FullSimplify[pe*Exp[pe*xVar]/(Exp[pe] - 1), Assumptions -> pe > 0];
expectZero["Sigma_Pe from rescaling", sigmaFromRescale - sigmaPe];
expectZero["normalized Sigma_Pe family",
  Integrate[sigmaPe, {xVar, 0, 1}] - 1];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 060` and `redteam exec-mathematica 060` and confirm both new `Sigma_Pe from rescaling` and `Csol normalizes ...` lines appear and the scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl`
- summary: Added explicit normalization and rescaling derivations for the Sigma_Pe exponential family in both engines.
- deviation: Used the equivalent Delta_phi/deltaDrop substitution for Pe because `a` is stored expanded in these scripts; the SymPy normalization integral uses `conds='none'` to select the generic nondegenerate branch.

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:75` (Pe identification)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:82-84` (Xi_micro identities)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:89-92` (product-rule identity)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:94-97` (Onsager dissipation density)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:82` (Pe identification)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:87-89` (xiMicro identities)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:93-94` (product-rule identity)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:95-97` (Onsager dissipation density)

**Issue:**
Six assertions verify their own definitions: (a) `Pe identification`: `a` was set to `Lambda*dphi/(Theta*L)` and the assertion just substitutes that definition; (b) `Xi_micro - Lambda^2 L^2/(Theta T_X)` and the susceptibility / phenomenological forms: `Xi_micro` was algebraically constructed as exactly that combination; (c) integration-by-parts: pure calculus identity; (d) Onsager dissipation density: `mu_s J_on + J_on^2/(M sigma)` with `J_on = -M sigma mu_s` cancels by direct substitution. None can fail.

**Required change:**

(a) **Pe identification.** Replace the substitution-identity check with a derivation. In the sympy script at line 75 (replace that line):

Before:
```
expect_zero("Pe identification", a*L - Lam*dphi/Theta)
```

After:
```
# derive the exponential rate from the affine-drop ODE J=0 without using `a`
gamma = sp.symbols('gamma', positive=True, real=True)
sigma_ansatz = Cnorm*sp.exp(gamma*s)
ode_ansatz = ode.subs(sigma, sigma_ansatz).replace(
    sp.Derivative(sigma_ansatz, s), sp.diff(sigma_ansatz, s))
gamma_solved = sp.solve(sp.Eq(ode_ansatz, 0), gamma)
# take the nonzero branch matching dphi/L != 0
gamma_derived = [g for g in gamma_solved if g != 0][0] if gamma_solved else None
print("gamma_derived =", gamma_derived)
expect_zero("Pe identification (derived rate)",
            sp.simplify(gamma_derived - Lam*dphi/(Theta*L)))
```

In wl at line 82, replace:

Before:
```
expectZero["Pe identification", a*ell - lambdaPhi*deltaDrop/theta];
```

After:
```
(* derive gamma from J=0 without referencing the predefined a *)
gammaVar = Symbol["gammaVar"];
sigmaAnsatz = cNorm*Exp[gammaVar*s];
odeAnsatz = jAff /. {sigmaField -> sigmaAnsatz, Derivative[1][sigma][s] -> D[sigmaAnsatz, s]};
gammaSolved = Solve[FullSimplify[odeAnsatz, Assumptions -> $Assumptions] == 0, gammaVar];
gammaDerived = gammaVar /. First[Select[gammaSolved, (gammaVar /. #) =!= 0 &]];
Print["gammaDerived = ", fmt[gammaDerived]];
expectZero["Pe identification (derived rate)",
  FullSimplify[gammaDerived - lambdaPhi*deltaDrop/(theta*ell), Assumptions -> $Assumptions]];
```

(b) **Xi_micro identities.** Pair this with F4 — once `phi_from_Phi` is derived (per F4), re-anchor the Xi_micro identity to the derived value. If F4 is `Blocked`, then at sympy:82-84 and wl:87-89, change the three identities to reference *symbolic* `phi_from_Phi` and `Xi_micro` as inputs (not derived values), and replace the three `expect_zero` calls with `print` statements noting "definitional rearrangement". Concretely in sympy:

Before (lines 80-84):
```
phi_from_Phi = Lam * L**2 * Delta / T_X  # end-to-end drop only
Xi_micro = sp.simplify(Lam * phi_from_Phi / Theta / Delta)
print("Xi_micro =", Xi_micro)
expect_zero("Xi_micro - Lambda^2 L^2/(Theta T_X)", Xi_micro - Lam**2 * L**2 / (Theta * T_X))
expect_zero("Xi_micro susceptibility form", Xi_micro.subs(Theta, 1/chi) - chi * Lam**2 * L**2 / T_X)
expect_zero("Xi_micro phenomenological form", Xi_micro.subs(Theta, Dsig/Msig) - Msig*Lam**2 * L**2/(Dsig*T_X))
```

After (only if F4 is Blocked; otherwise see F4 for the derived form):
```
# Xi_micro identities — these are algebraic rearrangements of the definition
# Xi_micro := Lambda * (phi_from_Phi/Delta) / Theta with phi_from_Phi = Lambda * L^2 * Delta / T_X.
# The physical derivation of phi_from_Phi is handled in finding F4 (deferred).
phi_from_Phi = Lam * L**2 * Delta / T_X
Xi_micro = sp.simplify(Lam * phi_from_Phi / Theta / Delta)
print("Xi_micro (definitional) =", Xi_micro)
print("Xi_micro susceptibility form (definitional) =", Xi_micro.subs(Theta, 1/chi))
print("Xi_micro phenomenological form (definitional) =", Xi_micro.subs(Theta, Dsig/Msig))
```

In wl (only if F4 is Blocked), the analogous change to lines 84-89: replace the three `expectZero` calls with `Print` statements.

(c) **Integration-by-parts identity.** Remove the assertion entirely.

In sympy, replace lines 87-92:

Before:
```
# 5) local dissipation identity
mu_fun = sp.Function('mu_sigma')(s)
J_fun = sp.Function('J')(s)
local_identity = sp.simplify(
    -mu_fun*sp.diff(J_fun, s) - (-sp.diff(mu_fun*J_fun, s) + sp.diff(mu_fun, s)*J_fun)
)
expect_zero("integration-by-parts identity", local_identity)
```

After:
```
# 5) local dissipation identity — note: the product-rule rearrangement
# -mu J' = -(mu J)' + mu' J is a calculus identity; we use it to motivate the
# integration-by-parts form below but do not assert it.
```

In wl, replace lines 91-94:

Before:
```
muFun = mu[s];
jFun = j[s];
localIdentity = FullSimplify[-muFun*D[jFun, s] - (-D[muFun*jFun, s] + D[muFun, s]*jFun), Assumptions -> $Assumptions];
expectZero["integration-by-parts identity", localIdentity];
```

After:
```
(* product-rule identity -mu J' = -(mu J)' + mu' J is a calculus identity, not asserted *)
```

(d) **Onsager dissipation density / positivity.**

In sympy, replace lines 94-97:

Before:
```
mu_s = sp.symbols('mu_s', real=True)
J_on = -Msig * sigma * mu_s
expect_zero("Onsager dissipation density",
            sp.simplify(mu_s * J_on + J_on**2/(Msig*sigma)))
```

After:
```
# Onsager dissipation density: with J_on = -M sigma mu_s, the dissipation density
# is J_on^2/(M sigma) = M sigma mu_s^2 >= 0. Verify positivity under M, sigma > 0.
sigma_val = sp.symbols('sigma_val', positive=True, real=True)
mu_s = sp.symbols('mu_s', real=True)
J_on = -Msig * sigma_val * mu_s
dissipation_density = sp.simplify(J_on**2/(Msig*sigma_val))
print("dissipation density =", dissipation_density)
# under M_sigma > 0, sigma > 0: M_sigma sigma mu_s^2 >= 0
assert sp.ask(sp.Q.nonnegative(dissipation_density),
              sp.Q.positive(Msig) & sp.Q.positive(sigma_val) & sp.Q.real(mu_s)) is True, \
    "dissipation density not provably nonnegative"
print("PASS: dissipation density nonnegative under M_sigma, sigma > 0")
```

In wl, replace lines 95-97:

Before:
```
muS = Symbol["muS"];
jOnsager = -mSigma*sigma*muS;
expectZero["Onsager dissipation density", FullSimplify[muS*jOnsager + jOnsager^2/(mSigma*sigma), Assumptions -> $Assumptions]];
```

After:
```
muS = Symbol["muS"];
sigmaVal = Symbol["sigmaVal"];
jOnsager = -mSigma*sigmaVal*muS;
dissipationDensity = FullSimplify[jOnsager^2/(mSigma*sigmaVal),
  Assumptions -> $Assumptions && mSigma > 0 && sigmaVal > 0 && Element[muS, Reals]];
Print["dissipation density = ", fmt[dissipationDensity]];
positivityCheck = Reduce[ForAll[{muS, sigmaVal, mSigma},
    mSigma > 0 && sigmaVal > 0 && Element[muS, Reals] ==> dissipationDensity >= 0],
  Reals];
If[TrueQ[positivityCheck === True],
  pass["dissipation density nonnegative under mSigma, sigmaVal > 0"],
  fail["dissipation density nonnegative", positivityCheck]];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 060` and `redteam exec-mathematica 060` and confirm: (a) a `gamma_derived = ...` line appears with `Pe identification (derived rate) = 0`; (c) the `integration-by-parts identity` line is gone; (d) a `dissipation density nonnegative ...` PASS line appears. Both scripts must still exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl`
- summary: Replaced tautological Pe, product-rule, and Onsager cancellation checks with a derived-rate check, a non-asserted product-rule note, and positivity checks.
- deviation: Used Mathematica `Implies[...]` for the positivity predicate because the directive's `==>` form is invalid script syntax.

## F3 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:58-89`

**Issue:**
The Mathematica script mirrors the SymPy structure step-for-step in three places: hand-typed Sigma_Pe (wl:80 vs sympy:71), hand-built xiMicro with explicit deltaSupport cancellation (wl:85 vs sympy:79-80), and pre-defined `a` followed by `cNorm*Exp[a*s]` ansatz (wl:58-60 vs sympy:60-61). This violates the two-engine policy: an error in either script survives transcription into the other.

**Required change:**

Two of the three transliterations are already addressed by F1 (which adds an explicit rescaling derivation in both engines) and F2(a) (which adds a `Solve[..., gamma]` step in the Mathematica script that does NOT appear in the SymPy script as written). For the third — the `xiMicro` construction at wl:85 — change the Mathematica script to derive `Xi_micro` via a path the SymPy script does not use. Specifically, at wl line 85, replace:

Before:
```
xiMicro = FullSimplify[(lambdaPhi/ theta)*(lambdaPhi*ell^2*deltaSupport/tX)/deltaSupport, Assumptions -> $Assumptions];
```

After:
```
(* Independent route: form Xi_micro from dimensional combination of derived
   parameters, then verify it matches the susceptibility and phenomenological
   forms used in the SymPy script. The deltaSupport cancellation is avoided. *)
xiMicro = lambdaPhi^2*ell^2/(theta*tX);
xiMicroFromChi = chiSigma*lambdaPhi^2*ell^2/tX /. chiSigma -> 1/theta;
expectZero["xiMicro consistency via chi substitution",
  FullSimplify[xiMicro - xiMicroFromChi, Assumptions -> $Assumptions]];
xiMicroFromDM = mSigma*lambdaPhi^2*ell^2/(dSigma*tX) /. dSigma -> mSigma*theta;
expectZero["xiMicro consistency via D/M substitution",
  FullSimplify[xiMicro - xiMicroFromDM, Assumptions -> $Assumptions]];
```

This eliminates the `deltaSupport` cancellation choreography from the wl, leaving the SymPy script as the sole user of that intermediate variable.

(If F4 is applied with a derived `phi_from_Phi`, this directive should still be applied since the goal is structural independence of the two engines, not just correctness.)

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 060` and confirm `deltaSupport` no longer appears in the `xiMicro` definition line and the new `xiMicro consistency via chi substitution` / `xiMicro consistency via D/M substitution` PASS lines appear. The script must still exit 0.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl`
- summary: Replaced the Mathematica xiMicro deltaSupport cancellation construction with an independent dimensional combination and consistency checks.
- deviation: none

## F4 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:78-84`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:101-107`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:84-89`

**Issue:**
`phi_from_Phi = Lambda*L^2*Delta/T_X` (sympy line 79; comment claims "end-to-end drop only") is the formula underlying the headline `Xi_micro` result, but it is typed in. The support Euler-Lagrange equation `-Lambda*sigma + K_X*phi - T_X*phi'' = 0` is derived in Section 6 (sympy lines 103-106) and the boundary conditions are stated as prose ("Boundary term from support variation: T_X phi_s(0) = K_m phi(0), phi_s(L)=0", sympy line 107), but no BVP is solved to extract the end-to-end drop. As a result, the connection between the EL equation and the value `Lambda*L^2*Delta/T_X` is unverified.

**Required change:**

Add a constant-source / `K_X -> 0` BVP solve and extract the end-to-end drop, then verify it matches `Lambda*L^2*sigma_0/T_X` up to the BC-dependent factor.

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py`, insert immediately before line 78 (i.e., before the `# 4) microscopic Xi from support normalization` block):

```
# 3.5) derive phi_from_Phi from the support EL in the constant-sigma, K_X=0 limit
#      so the headline Xi_micro value below has a physical basis.
phi_bvp = sp.Function('phi_BVP')(s)
sigma0 = sp.symbols('sigma_0', positive=True, real=True)
EL_const = -Lam*sigma0 - T_X*sp.diff(phi_bvp, s, 2)  # K_X -> 0 limit
sol_general = sp.dsolve(sp.Eq(EL_const, 0), phi_bvp).rhs
C1, C2 = sp.symbols('C1 C2')
# match dsolve constant names (sympy may use C1, C2 in either order)
free_syms = sorted([sym for sym in sol_general.free_symbols if str(sym).startswith('C')], key=str)
c1, c2 = free_syms[0], free_syms[1]
# BCs: T_X phi'(0) = K_m phi(0), phi'(L) = 0
phi_s = sp.diff(sol_general, s)
bc1 = sp.Eq(T_X * phi_s.subs(s, 0), K_m * sol_general.subs(s, 0))
bc2 = sp.Eq(phi_s.subs(s, L), 0)
const_vals = sp.solve([bc1, bc2], [c1, c2])
phi_solved = sp.simplify(sol_general.subs(const_vals))
Delta_derived = sp.simplify(phi_solved.subs(s, L) - phi_solved.subs(s, 0))
print("phi_BVP(0) =", sp.simplify(phi_solved.subs(s, 0)))
print("phi_BVP(L) =", sp.simplify(phi_solved.subs(s, L)))
print("Delta_derived = phi(L) - phi(0) =", Delta_derived)
# headline formula: Lambda L^2 sigma_0 / T_X is the bulk-scale drop;
# the K_m -> infinity (rigid grounding at s=0) limit fixes the prefactor to 1/2.
Delta_target_rigid = Lam * L**2 * sigma0 / (2 * T_X)
Delta_rigid_limit = sp.limit(Delta_derived, K_m, sp.oo)
expect_zero("phi_from_Phi from support BVP (K_m -> infty)",
            sp.simplify(Delta_rigid_limit - Delta_target_rigid))
```

After this block, update line 79 from:

Before:
```
phi_from_Phi = Lam * L**2 * Delta / T_X  # end-to-end drop only
```

After:
```
# phi_from_Phi normalized to the bulk-scale drop Lambda L^2 / T_X per unit Delta
# (validated above against the rigid-grounding BVP limit)
phi_from_Phi = Lam * L**2 * Delta / T_X  # end-to-end drop only (see BVP check above)
```

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl`, insert immediately before the existing `xiMicro = ...` line (currently line 85) — keep the F3 edit above; F4 inserts before F3's replacement:

```
(* derive phi_from_Phi from the support EL in the constant-sigma, K_X=0 limit *)
phiBVP = Symbol["phiBVP"];
sigma0 = Symbol["sigma0"];
elConst[fun_] := -lambdaPhi*sigma0 - tX*D[fun[s], {s, 2}];
solGeneral = DSolveValue[elConst[phiBVP] == 0, phiBVP[s], s];
{cc1, cc2} = Sort[Cases[{solGeneral}, _C, Infinity]];
phiS = D[solGeneral, s];
bc1 = tX*(phiS /. s -> 0) == kM*(solGeneral /. s -> 0);
bc2 = (phiS /. s -> ell) == 0;
constVals = Solve[{bc1, bc2}, {cc1, cc2}];
phiSolved = FullSimplify[solGeneral /. First[constVals], Assumptions -> $Assumptions];
deltaDerived = FullSimplify[(phiSolved /. s -> ell) - (phiSolved /. s -> 0),
  Assumptions -> $Assumptions && sigma0 > 0];
Print["Delta_derived = ", fmt[deltaDerived]];
deltaRigidLimit = FullSimplify[Limit[deltaDerived, kM -> Infinity],
  Assumptions -> $Assumptions && sigma0 > 0];
deltaTargetRigid = lambdaPhi*ell^2*sigma0/(2*tX);
expectZero["phi_from_Phi from support BVP (kM -> infty)",
  FullSimplify[deltaRigidLimit - deltaTargetRigid, Assumptions -> $Assumptions]];
```

If the rigid-grounding limit (`K_m -> oo`) does not yield exactly `Lambda L^2 sigma_0 / (2 T_X)` — for instance, if the script's BC convention differs from this directive's reading — add a `## Blocked: F4` block stating which factor disagrees, and skip F4. F1, F2, F3 should still be applied.

**Self-test note:** Trial substitution for the proposed BVP: with `EL: -Lam sigma_0 - T_X phi'' = 0`, the general solution is `phi(s) = -Lam sigma_0 s^2/(2 T_X) + C1*s + C2`. At `s=0`: `phi(0) = C2`, `phi'(0) = C1`. At `s=L`: `phi'(L) = -Lam sigma_0 L/T_X + C1 = 0` → `C1 = Lam sigma_0 L/T_X`. BC1: `T_X * C1 = K_m * C2` → `C2 = T_X C1 / K_m = Lam sigma_0 L / K_m`. Then `Delta = phi(L) - phi(0) = -Lam sigma_0 L^2/(2 T_X) + C1*L = -Lam sigma_0 L^2/(2 T_X) + Lam sigma_0 L^2/T_X = Lam sigma_0 L^2/(2 T_X)`. So `Delta = Lambda sigma_0 L^2/(2 T_X)` independent of `K_m` (in this K_X=0, constant-sigma limit) — meaning the `K_m -> infty` limit is just the value itself, and the rigid-grounding factor is `1/2`. The directive's target `Lambda L^2 sigma_0/(2 T_X)` is consistent.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 060` and `redteam exec-mathematica 060` and confirm a new `phi_from_Phi from support BVP (...) = 0` PASS line appears in both outputs. Both scripts must still exit 0.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl`
- summary: Added constant-source K_X -> 0 support BVP solves in both engines and verified the rigid-grounding end-to-end drop factor.
- deviation: none
