---
unit_id: 060
batch: III.2
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 060 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.txt`

## What the script claims to verify

The script claims to derive an entropic source microclosure: (i) the chemical potential mu_sigma from a positive-density free energy is `Theta log(sigma/sigma_*) - Lambda phi`; (ii) the Onsager current `J = -D_sigma sigma_s + M_sigma Lambda sigma phi_s` with the Einstein relation `D_sigma = M_sigma Theta`; (iii) on an affine support drop, the zero-flux branch is the exponential family `Sigma_Pe(x) = Pe e^{Pe x}/(e^{Pe}-1)` with `Pe = Lambda Delta_phi/Theta`; (iv) the microscopic gain is `Xi_micro = Lambda^2 L^2/(Theta T_X)` in three equivalent forms; (v) the closure is dissipative (`dF/dt <= 0`); (vi) the support Euler-Lagrange equation is `-Lambda sigma + K_X phi - T_X phi'' = 0`. The verdict is therefore on these six chained claims.

## Assertion inventory

| #  | Script | Line | Form | Anchored to claim? |
|----|--------|------|------|--------------------|
| A1 | sympy  | 42-43 | `simplify(mu_expr - (Theta log(sigma/sigma_*) - Lambda phi)) == 0` | yes (verifies sympy's differentiation of f) |
| A2 | sympy  | 48-49 | `J_expr - (-M Theta sigma' + M Lam sigma phi') == 0` | partial (just expansion of chain rule) |
| A3 | sympy  | 51-52 | `J_expr.subs(M*Theta, D) - (-D sigma' + M Lam sigma phi') == 0` | no (pure substitution) |
| A4 | sympy  | 62-63 | exponential trial satisfies J=0 on affine branch | yes (substantive) |
| A5 | sympy  | 73-74 | `integrate(Pe e^{Pe x}/(e^{Pe}-1), 0..1) - 1 == 0` | partial (Sigma_x hardcoded; see F1) |
| A6 | sympy  | 75 | `a*L - Lam*dphi/Theta == 0` | no (a defined as Lam*dphi/(Theta*L)) |
| A7 | sympy  | 82 | `Xi_micro - Lam^2 L^2/(Theta T_X) == 0` | no (Xi_micro is `Lam*phi_from_Phi/(Theta*Delta)` with `phi_from_Phi = Lam L^2 Delta/T_X` — algebraic identity) |
| A8 | sympy  | 83 | `Xi.subs(Theta,1/chi) - chi Lam^2 L^2/T_X == 0` | no (substitution identity) |
| A9 | sympy  | 84 | `Xi.subs(Theta,D/M) - M Lam^2 L^2/(D T_X) == 0` | no (substitution identity) |
| A10| sympy  | 89-92 | product-rule `-mu J' - (-(mu J)' + mu' J) == 0` | no (algebraic identity, always 0) |
| A11| sympy  | 96-97 | `mu_s J_on + J_on^2/(M sigma) == 0` with `J_on = -M sigma mu_s` | no (direct substitution: `-M sigma mu_s^2 + M sigma mu_s^2`) |
| A12| sympy  | 106 | EL[phi] of f_full equals `-Lam sigma + K_X phi - T_X phi''` | yes (verifies sympy's variational derivative) |
| B1 | math   | 42 | mu identity | yes |
| B2 | math   | 49-51 | Onsager current decomposition | partial |
| B3 | math   | 53-56 | D_sigma substitution | no (pure substitution) |
| B4 | math   | 69-78 | exponential trial annihilates jAff | yes (substantive) |
| B5 | math   | 81 | normalized Sigma_Pe integral | partial (sigmaPe hardcoded; see F1) |
| B6 | math   | 82 | `a*ell - lambdaPhi*deltaDrop/theta == 0` | no (a defined to make this zero) |
| B7 | math   | 87 | `xiMicro - lambdaPhi^2 ell^2/(theta tX) == 0` | no (algebraic identity) |
| B8 | math   | 88 | susceptibility-form substitution | no |
| B9 | math   | 89 | phenomenological-form substitution | no |
| B10| math   | 94 | product-rule identity | no |
| B11| math   | 97 | Onsager dissipation density | no |
| B12| math   | 106 | EL equation form | yes |

## Findings

### F1 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:66-74`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:80-81`

**What's wrong:**
The script attempts to derive the normalized exponential family `Sigma_Pe(x) = Pe e^{Pe x}/(e^{Pe}-1)` from the affine-drop branch. In the SymPy script, line 66 calls

```
Csol = sp.simplify(sp.solve(sp.Eq(sp.integrate(sigma_trial, (s, 0, L)), 1), Cnorm)[0])
```

but the captured output (line 22 of the .txt) shows sympy returned `Piecewise((1/L, Eq(Delta_phi*Lambda_phi, 0)), (nan, True))` — i.e., sympy picked the degenerate `a*Lambda_phi = 0` branch and returned `nan` for the generic case. `Csol_Pe = Csol.subs(a, Pe/L)` then inherits that `nan` and is never asserted on. Line 71 then types in `Sigma_x = Pe*exp(Pe*x)/(exp(Pe)-1)` by hand, and the subsequent integral check (line 73-74) merely verifies the typed-in form integrates to 1. The Mathematica script does the same on lines 80-81: `sigmaPe = pe*Exp[pe*s]/(Exp[pe]-1)` is hardcoded with no derivation from `cNorm*Exp[a*s]` and the normalization `Integrate[sigmaTrial, {s,0,ell}] = 1`. The hardcoded `Sigma_Pe` happens to be correct, but the script never proves it.

**Why this matters:**
A core claim of the theorem ledger ("Under affine support drop, the zero-flux branch is exactly the Stage-39 exponential family") is asserted only by writing the family down and checking it normalizes to 1 — which is true of any hand-typed function of that form. The actual derivation step (normalizing `C*exp(a*s)` and re-scaling to `s = x L`, `Pe = a L`) is unverified. If a future edit changes the `a` definition or the affine-drop sign convention, this script will continue to pass.

**Required change:**
In the SymPy script, replace line 66 with an explicit symbolic solve that excludes the degenerate branch (e.g., compute `Csol = a/(sp.exp(a*L) - 1)` manually and then `expect_zero` that `sp.integrate(Csol*sp.exp(a*s), (s, 0, L)) - 1 == 0`). Then add a new `expect_zero` immediately before line 74 that verifies the change-of-variables `s = x*L, Pe = a*L` carries `C*exp(a*s) * ds` to `Sigma_x dx`:

```
Sigma_from_rescale = sp.simplify(L * Csol.subs(a, Pe/L) * sp.exp((Pe/L)*(x*L)))
expect_zero("Sigma_Pe from rescaling", Sigma_from_rescale - Sigma_x)
```

In the Mathematica script, add an analogous assertion immediately before line 81: solve `Integrate[cNorm*Exp[a*s],{s,0,ell}] == 1` for `cNorm` (manually as `a/(Exp[a*ell]-1)`), verify it normalizes, then `expectZero["sigmaPe from rescaling", FullSimplify[ell * cNormSol /. a -> pe/ell * Exp[(pe/ell)*(s*ell)/ell] /.{...}, ...] - sigmaPe]` (the exact form may need adjustment; the point is the rescaling must be exercised, not the answer hardcoded).

**Verification:**
The new `Sigma_Pe from rescaling` line should appear in the .txt outputs after a re-run, and both scripts must still exit 0. If the rescaling is correct, the residual will simplify to 0; if not, the bug will surface.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:75` (`expect_zero("Pe identification", a*L - Lam*dphi/Theta)`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:80-84` (Xi_micro block)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:89-97` (dissipation identities)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:82` (`a*ell - lambdaPhi*deltaDrop/theta`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:85-89` (xiMicro block)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:93-97` (dissipation identities)

**What's wrong:**
Six assertions are algebraic identities that follow trivially from the surrounding definitions and cannot fail no matter what the physics is.

1. **Pe identification** (sympy:61,75 / wl:58,82). The script defines `a = Lambda*dphi/(Theta*L)` (sympy line 61) and then asserts `a*L - Lambda*dphi/Theta == 0` (line 75). Substituting `a`'s definition gives `Lambda*dphi/Theta - Lambda*dphi/Theta = 0` identically. Same in the wl on lines 58 and 82.

2. **Xi_micro block** (sympy:79-84 / wl:85-89). The script defines

   ```
   phi_from_Phi = Lam * L**2 * Delta / T_X      # sympy line 79
   Xi_micro = sp.simplify(Lam * phi_from_Phi / Theta / Delta)
   ```

   so `Xi_micro = Lam^2 L^2/(Theta T_X)` by construction. The three `expect_zero` lines that follow (Xi_micro definition, susceptibility form `Theta -> 1/chi`, phenomenological form `Theta -> D/M`) are then just substituting symbols into a known algebraic identity. None of these check a physical claim — they check substitution. The Mathematica version on line 85 is even more transparent: `xiMicro = FullSimplify[(lambdaPhi/theta)*(lambdaPhi*ell^2*deltaSupport/tX)/deltaSupport]` literally cancels `deltaSupport` against itself. The "Xi_micro - Lambda^2 L^2/(Theta T_X) == 0" check is unfalsifiable.

3. **Integration-by-parts identity** (sympy:89-92 / wl:93). The expression `-mu J' - (-(mu J)' + mu' J)` simplifies to 0 by the product rule for any choice of mu and J. It is a calculus identity, not a physics claim.

4. **Onsager dissipation density** (sympy:94-97 / wl:95-97). With `J_on = -M*sigma*mu_s` substituted in, `mu_s*J_on + J_on^2/(M*sigma) = -M*sigma*mu_s^2 + M*sigma*mu_s^2 = 0` follows by direct substitution. It cannot fail and does not establish the positivity claim in the docstring ("dF/dt <= 0") — it just shows two expressions for the same quantity cancel.

**Why this matters:**
The theorem ledger lists "Pe identification" and the Xi_micro forms as derived results, but the script only verifies that the author can do algebra. A wrong physical derivation of `a`, `phi_from_Phi`, or `Xi_micro` would never be caught by these assertions because each is checked against its own definition. The dissipation identity claim — `dF/dt <= 0` — requires demonstrating that `J^2/(M*sigma) >= 0`, which is a positivity statement under `M, sigma > 0`, not an algebraic cancellation.

**Required change:**

(a) Pe identification (sympy:75, wl:82): replace with a derivation check. Compute the exponential rate from the affine-drop ODE `Theta*sigma_s = Lambda*sigma*dphi/L` (which gives `sigma_s/sigma = a`) and verify the extracted rate equals `Lambda*dphi/(Theta*L)` *without* first defining `a` that way. Specifically, in sympy: derive `a_derived` by `sp.solve(sp.Eq(sp.diff(sigma_trial, s)/sigma_trial, a_unknown), a_unknown)[0]` substituted into the affine ODE, then `expect_zero("Pe identification", a_derived - Lam*dphi/(Theta*L))`. In wl: the analogous derivation.

(b) Xi_micro block (sympy:79-84, wl:84-89): the value `phi_from_Phi = Lam*L^2*Delta/T_X` is unsupported (see F4). Either (i) derive `phi_from_Phi` from the support EL with the stated BCs (`T_X phi_s(0) = K_m phi(0), phi_s(L) = 0`) under a constant-sigma assumption — e.g., solve `-Lam*sigma_0 + K_X phi - T_X phi'' = 0` and extract the end-to-end drop — and then assert `Xi_micro` from the derivation, OR (ii) remove the Xi_micro identities as substitution-only checks and replace them with a single derivation. The susceptibility and phenomenological forms are then either dropped or re-asserted from the *derived* Xi_micro, not from the typed-in form.

(c) Dissipation positivity (sympy:96-97, wl:97): replace the trivial cancellation with a positivity check. Concretely, drop the `expect_zero` line and add a sign-domain check:

```
assert J_on**2/(Msig*sigma) is nonnegative under (Msig > 0, sigma > 0)
```

In sympy, use `sp.assuming(sp.Q.positive(Msig) & sp.Q.positive(sigma), lambda: sp.ask(sp.Q.nonnegative(J_on**2/(Msig*sigma))))` or substitute concrete positive values and assert numerical nonnegativity. In wl, `Resolve[ForAll[{...}, mSigma > 0 && sigma > 0 ==> jOnsager^2/(mSigma*sigma) >= 0], Reals]` and `If[result =!= True, fail[...]]`.

(d) Integration-by-parts identity (sympy:89-92, wl:93): this is a pure calculus identity and contributes nothing. Remove these assertions or label them explicitly as "calculus identity" comments rather than physics checks.

**Verification:**
After the edits, the .txt outputs should show: (a) a non-trivial residual derivation line for Pe identification (e.g., `a_derived = Lambda_phi*Delta_phi/(L*Theta)`), (b) either a new derivation line for `phi_from_Phi` or removal of the substitution-only Xi_micro lines, (c) a positivity verdict (e.g., `dissipation density positivity: True`) replacing the `Onsager dissipation density = 0` line, (d) the product-rule line removed or relabelled. Both scripts must still exit 0.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:35-106` (entire body)

**What's wrong:**
The Mathematica script is structurally a direct port of the SymPy script rather than an independent re-derivation. Specifically, the same section ordering, the same intermediate variables, the same hand-typed trial functions, and the same hand-typed `phi_from_Phi` are used in both, in the same logical sequence. Three corresponding excerpts:

1. **Identical hand-typed Sigma_Pe with no derivation in either engine.** SymPy line 71: `Sigma_x = sp.simplify(Pe * sp.exp(Pe*x) / (sp.exp(Pe) - 1))`. Mathematica line 80: `sigmaPe = FullSimplify[pe*Exp[pe*s]/(Exp[pe] - 1), Assumptions -> pe > 0]`. Both bypass the failed/missing solve and type the answer directly.

2. **Identical Xi_micro construction with explicit `Delta` cancellation.** SymPy lines 79-80:

   ```
   phi_from_Phi = Lam * L**2 * Delta / T_X
   Xi_micro = sp.simplify(Lam * phi_from_Phi / Theta / Delta)
   ```

   Mathematica line 85: `xiMicro = FullSimplify[(lambdaPhi/ theta)*(lambdaPhi*ell^2*deltaSupport/tX)/deltaSupport, Assumptions -> $Assumptions]`. The exact same algebraic decomposition (the explicit `Delta`/`deltaSupport` cancellation that yields `lambda^2 L^2/(theta T_X)`) is performed, when a truly independent re-derivation in Mathematica would either skip introducing `deltaSupport` at all (since it cancels) or derive `phi_from_Phi` from the support EL.

3. **Same trial-function ansatz with the same definition of `a`.** SymPy line 60: `a = Lam * dphi / (Theta * L); sigma_trial = Cnorm * sp.exp(a*s)`. Mathematica line 58-60: `a = FullSimplify[lambdaPhi*deltaDrop/(theta*ell), Assumptions -> ...]; cNorm = Symbol["Cnorm"]; sigmaTrial = cNorm*Exp[a*s]`. Same definition order, same C * exp(a s) ansatz.

**Why this matters:**
The two-engine policy requires independent derivation paths so a mistake or simplifier artefact in one engine cannot survive in the other by mere transcription. Here, any error in the SymPy script (e.g., the failed `solve` returning `nan`, the unsupported `phi_from_Phi`, the tautological `Pe identification`) is reproduced verbatim in the Mathematica script and would not be caught by cross-engine comparison.

**Required change:**
Restructure the Mathematica script to derive the same claims via a different route. Concrete options (any one is sufficient): (i) in the affine-drop section (wl:58-78), do not pre-define `a`; instead, substitute the ansatz `sigmaTrial = cNorm*Exp[gamma*s]` into the affine-drop ODE and solve `Solve[J_aff == 0, gamma]` for the rate, then verify the solution equals `lambdaPhi*deltaDrop/(theta*ell)`. (ii) in the Xi_micro section (wl:84-89), drop the `deltaSupport`-cancellation form; instead, solve the support EL `-lambdaPhi*sigma0 + kX*phi - tX*phi'' = 0` with the stated BCs in the constant-sigma limit, extract the end-to-end drop `Delta = phi(ell) - phi(0)`, form `Xi = lambdaPhi * Delta/(theta*sigma0*Delta_normalized)` or analogous dimensionless combination, and assert it equals `lambdaPhi^2 ell^2/(theta tX)` in the small-`kX` limit. (iii) at minimum, change at least 2 of the 3 quoted constructions above so the SymPy script's bugs cannot survive transcription.

**Verification:**
After the edit, the Mathematica .txt should show a derivation chain visibly distinct from the SymPy output — at least one of: a `Solve[..., gamma]` result printed, a `DSolve`-derived support profile, or a different intermediate variable for Xi_micro. Both scripts must still exit 0 and their final residuals must still agree.

### F4 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:78-84`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:101-107`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl:84-89`

**What's wrong:**
The docstring/theorem ledger states (item 4): "The microscopic support/source gain is `Xi_micro = Lambda_phi^2 L^2/(Theta T_X) = chi_sigma Lambda_phi^2 L^2/T_X`." The script defines this value via `phi_from_Phi = Lam * L^2 * Delta / T_X` (sympy line 79; called "end-to-end drop only" in the comment), but never derives this `phi_from_Phi`. It is typed in. The support EL equation is derived in Section 6 (sympy:103-106) — `-Lam sigma + K_X phi - T_X phi'' = 0` — and the boundary conditions are *mentioned in a print statement* (line 107: "Boundary term from support variation: T_X phi_s(0) = K_m phi(0), phi_s(L)=0") but never *imposed* in a BVP solve. As a result, the value `Lam L^2 Delta/T_X` for the end-to-end drop is unverified — it is consistent with the EL only in a specific limit (e.g., `K_X = 0`, constant sigma, and a particular boundary regime), but no such limit is identified or checked.

**Why this matters:**
Xi_micro is the headline result of this unit. Its three "alternative forms" all reduce to a hand-typed value with one symbol replaced. If the actual support BVP yields a different `phi_from_Phi` (e.g., one involving `K_X`, the mass term, or sensitive to the boundary stiffness `K_m`), the script will not detect it. Item 6 of the theorem ledger (EL equation) does exist in the script but is not connected to item 4. Worse, the boundary conditions are stated as prose (`print("Boundary term from support variation: ...")`) rather than imposed, so the variational derivation is incomplete on its own terms.

**Required change:**
Add a substantive derivation linking the EL equation (sympy:104-106, wl:104) to the value `phi_from_Phi = Lam L^2 Delta/T_X`. The simplest unambiguous route: in the constant-sigma limit (`sigma_fun = sigma_0`), set `K_X = 0` (or take the relevant limit explicitly), solve the resulting equation `-Lam sigma_0 - T_X phi'' = 0` with BCs `T_X phi_s(0) = K_m phi(0), phi_s(L) = 0`, and verify that the end-to-end drop `phi(L) - phi(0)` reduces to `Lam L^2 sigma_0/T_X` times a O(1) numerical factor that depends only on the dimensionless ratio `K_m*L/T_X`. Then express `Xi_micro` in terms of this derived drop and confirm the headline form is recovered. In sympy:

```
phi_func = sp.Function('phi_BVP')(s)
sigma0 = sp.symbols('sigma_0', positive=True)
EL = -Lam*sigma0 - T_X*sp.diff(phi_func, s, 2)  # K_X=0 limit
sol = sp.dsolve(sp.Eq(EL, 0), phi_func)
# impose BCs, extract Delta = phi(L) - phi(0)
# assert Delta proportional to Lam L^2 sigma0/T_X with the right numeric factor
expect_zero("phi_from_Phi from support BVP",
            Delta_derived - Lam*L**2*sigma0/T_X * <derived_factor>)
```

In wl: the analogous `DSolve`. If introducing this BVP is judged too invasive, the alternative is to *remove* the Xi_micro identities (sympy:80-84, wl:85-89) and the corresponding theorem-ledger items, and replace the docstring claim with one limited to "Xi_micro is *defined* as Lam^2 L^2/(Theta T_X); the support BVP that justifies this is deferred to unit XXX." The current state — assertion without derivation — must not stand.

**Verification:**
After the edit, the .txt outputs should show either (a) a new derivation line `phi_from_Phi from support BVP = 0` corresponding to the solved BVP, or (b) the Xi_micro identity lines removed and a comment in the script's docstring deferring the derivation. Both scripts must still exit 0.

## Independent-derivation check (Mathematica)

The Mathematica script is **not** an independent re-derivation. Three of its key constructions mirror the SymPy script step-for-step: the `Sigma_Pe` family is hand-typed (wl:80) just as it is in sympy (line 71); `xiMicro` is built from the same `lambdaPhi*ell^2*deltaSupport/tX` expression and the same `Delta` cancellation as sympy's `phi_from_Phi` / `Xi_micro` lines (wl:85 versus sympy:79-80); and the affine-drop branch uses the same `Cnorm*Exp[a*s]` ansatz with the same pre-defined `a = lambdaPhi*deltaDrop/(theta*ell)` (wl:58-60 versus sympy:60-61). See F3 for the formal finding.

## Engine cross-check

Both engines produce zero residuals on every assertion (sympy .txt lines 15-34; wl .txt lines 15-41). The numerical / symbolic outputs agree. However, because both scripts assert the *same* expressions in the *same* order (including the same tautological and hardcoded ones), agreement here does not constitute independent verification — see F3.

## Verdict justification

The verdict is `findings`: four findings, all medium severity. What holds up: the chemical-potential derivation (A1/B1), the Onsager current chain-rule expansion (A2/B2), the substantive affine-drop annihilation (A4/B4), and the Euler-Lagrange derivation in Section 6 (A12/B12) are genuine symbolic checks. What does not hold up: the normalized exponential family `Sigma_Pe` is asserted as a hand-typed answer rather than derived from the affine ODE (F1); six assertions across the two scripts are algebraic identities trivially satisfied by their own definitions (F2); the Mathematica script transliterates the SymPy structure rather than re-deriving (F3); and the headline `Xi_micro` value depends on an unjustified end-to-end-drop formula `phi_from_Phi = Lam L^2 Delta/T_X` (F4). Attacks that the script **did** survive: the chain-rule expansion of `J_expr` (A2) does check sympy's algebra; the exponential-ansatz substitution at A4 truly reduces to zero only when `a = Lam*dphi/(Theta*L)`, so that one assertion is real; and the EL derivation at A12 catches the variational derivative correctly. The findings target the unsupported claims, not these.

## Self-test notes

Walked through each proposed directive: (1) **Variable independence** — the `expect_zero` lines I propose for F1 (`Sigma_from_rescale - Sigma_x`) and F2(a) (`a_derived - Lam*dphi/(Theta*L)`) involve `s, L, Pe, x, Lambda, Delta_phi, Theta` only; each variable actually appears in the expression so no spurious "derivative is identically zero" trap. (2) **Symmetry/parity** — no symmetric-domain integrals are introduced; the only integral is `integrate(Sigma_x, x, 0, 1)` over a half-domain, which is the unchanged check. (3) **Trivial-case pre-check** — for F1, substituting `Pe = 0` in `Sigma_from_rescale = L * Pe/L / (exp(Pe) - 1) * exp(Pe*x)` gives the `Pe -> 0` limit `1`, matching `Sigma_x = 0/0` regularized; the proposed check returns zero in the regular case, which is the desired behavior, and the existing degenerate-branch issue is exactly what F1 wants surfaced. For F2(c) positivity, the trivial case `mu_s = 1, M = 1, sigma = 1` gives `J^2/(M*sigma) = 1 > 0`, confirming the assertion is non-trivially nonzero in general. (4) **Path specifications** — no `missing_verification_script` findings here, so target paths are existing files only; both paths confirmed in the file listing.
