---
unit_id: 062
batch: III.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 062 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage062_parent_action_gain_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.txt`

## What the script claims to verify

The docstring lists five checks: (1) an n=5 EOS identity `h'(rho) = m c_s^2 / rho`; (2) parent-channel projection formulas for `Theta_sigma` and `Lambda_phi`; (3) an "exact parent formula" for the microscopic gain `G_micro = (1/Theta_sigma) Lambda_phi^2 / KX`; (4) a coherence-factor decomposition `Osp^2 = C_(sigma phi)^2 Nss Npp`; and (5) the reduction `Xi_micro = kappa G_micro` with `kappa = KX L^2 / TX`. The four substantive `expect_zero` assertions are the bottom line. Each engine prints intermediate symbolic forms and asserts that a residual simplifies to zero.

## Assertion inventory

| #  | Script       | Line | Form                                                                             | Anchored to claim? |
|----|--------------|------|----------------------------------------------------------------------------------|--------------------|
| A1 | sympy        | 44   | `expect_zero("h'(rho) - m c_s^2 / rho", hprime - m*cs_sq/rho)`                  | partial            |
| A2 | sympy        | 61   | `expect_zero("G_micro - expected parent formula", G_micro - G_expected)`         | no                 |
| A3 | sympy        | 68-71| `expect_zero("coherence-factor decomposition", G_expected.subs(...) - G_coh)`    | no                 |
| A4 | sympy        | 78-81| `expect_zero("Xi_micro - parent projected formula", Xi_micro.subs(...) - Xi_expected)` | no             |
| A5 | mathematica  | 40   | `expectZero["h'(rho) - m c_s^2 / rho", hPrime - m*csSq/rho]`                    | partial            |
| A6 | mathematica  | 58   | `expectZero["G_micro - expected parent formula", gMicro - gExpected]`            | no                 |
| A7 | mathematica  | 63   | `expectZero["coherence-factor decomposition", ...gExpected /. oSP^2 -> c2*nSS*nPP - gCoherence]` | no |
| A8 | mathematica  | 68   | `expectZero["Xi_micro - parent projected formula", ...xiMicro /. kappa -> kX*ell^2/tX - xiExpected]` | no |

A1/A5 are flagged "partial" because the identity `h'(rho) = m c_s^2 / rho` is an algebraic consequence of choosing `U = K rho^5/4` together with `c_s^2 := (1/m) d(K rho^5)/drho`: since `c_s^2 = (4/m) h` and `h` is a degree-4 monomial in rho (so `rho * h' = 4 h`), the identity `h' = m c_s^2 / rho = 4h/rho` is forced by the polynomial degree, not by any physical input. It is technically derivable rather than imposed, so I list it as "partial" rather than fully tautological — but it carries essentially no information beyond "h is degree 4 in rho".

A2/A6, A3/A7, A4/A8 are **tautological** in the strict sense (see F1).

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:51-61`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:62-71`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:73-81`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:48-58`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:60-63`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:65-68`

**What's wrong:**
The three "parent-action" assertions (A2/A3/A4 and their mirror A6/A7/A8) consist entirely of defining a quantity and then asserting that it equals its own algebraic re-write. They cannot fail.

A2 (`G_micro - G_expected`), sympy lines 51-61:
```
Theta_sigma = (m * cs_star_sq / rho_star) * Nss
Lambda_phi  = g_phi * Osp
chi_eff     = simplify(1 / Theta_sigma)                # = rho_star / (m cs_star_sq Nss)
G_micro     = simplify(chi_eff * Lambda_phi**2 / KX)   # = rho_star g_phi^2 Osp^2 / (m cs_star_sq KX Nss)
G_expected  = simplify(rho_star * g_phi**2 * Osp**2 / (m * cs_star_sq * KX * Nss))
expect_zero("G_micro - expected parent formula", G_micro - G_expected)
```
`G_micro` is *by construction* `rho_star g_phi^2 Osp^2 / (m cs_star_sq KX Nss)` once you substitute the definitions of `Theta_sigma` and `Lambda_phi`. `G_expected` is then the literal same expression, hand-typed. There is no independent path to the RHS — it's the algebraic simplification of the LHS, retyped. `expect_zero(G_micro - G_expected)` is therefore `0 == 0` for any choice of symbols.

A3 (coherence-factor decomposition), sympy lines 62-71:
```
C2_def = simplify(Osp**2 / (Nss * Npp))
G_coh  = simplify(rho_star * g_phi**2 * Npp * C2 / (m * cs_star_sq * KX))
expect_zero("coherence-factor decomposition",
            G_expected.subs(Osp**2, C2 * Nss * Npp) - G_coh)
```
Applying the substitution `Osp**2 -> C2 * Nss * Npp` to `G_expected = rho_star g_phi^2 Osp^2/(m cs_star_sq KX Nss)` gives `rho_star g_phi^2 (C2 Nss Npp)/(m cs_star_sq KX Nss) = rho_star g_phi^2 C2 Npp/(m cs_star_sq KX)` which is exactly `G_coh` by definition. The check verifies that SymPy can perform a textbook substitution, not that the coherence factorization is a physical claim about the model. Note also that `C2_def` is printed but never asserted to match `C2` anywhere — the "definition" is informational only.

A4 (Xi_micro), sympy lines 73-81:
```
Xi_micro    = simplify(kappa * G_micro)
Xi_expected = simplify(rho_star * g_phi**2 * Osp**2 * L**2 / (m * cs_star_sq * TX * Nss))
expect_zero("Xi_micro - parent projected formula",
            Xi_micro.subs(kappa, KX * L**2 / TX) - Xi_expected)
```
`Xi_micro.subs(kappa, KX L^2/TX) = (KX L^2/TX) * rho_star g_phi^2 Osp^2/(m cs_star_sq KX Nss) = rho_star g_phi^2 Osp^2 L^2/(m cs_star_sq TX Nss)`, which is exactly `Xi_expected` by hand-typing. Tautological.

The Mathematica file replicates each of these constructions identically (A6/A7/A8 at lines 48-58, 60-63, 65-68 respectively), so the issue is mirrored in both engines.

**Why this matters:**
The script's docstring promises "exact parent formula for G_micro" and "Xi_micro = kappa G_micro reduction" as substantive results. As written, the script proves only that SymPy and Mathematica can simplify algebra correctly. If the "parent formula" `G_micro = rho_star g_phi^2 Osp^2/(m cs_star_sq KX Nss)` were wrong (e.g., a missing factor of 2, a sign on `Lambda_phi`, an exponent on `Osp`, or a different combinatorial factor in `Theta_sigma`), the assertions would still pass because the script never tests it against anything external — it tests it against a retyped copy of itself.

The actual physics anchors that should drive these residuals are not derived in the script. For example:
- The form `Theta_sigma = (m cs_star_sq / rho_star) Nss` should be derivable from a stated parent action by taking the appropriate second functional derivative on the sigma channel; the script just writes it down.
- The form `Lambda_phi = g_phi Osp` should be derivable from the projection of the coupling onto the phi support; the script just writes it down.
- The chain `G = chi_eff Lambda^2 / KX` should be derivable from a Gaussian integral over sigma after integrating out the heavy channel; the script just writes it down.
- The relation `kappa = KX L^2 / TX` should follow from dimensional carry of `KX`, `L`, `TX`; the script just substitutes it.

Without one of those derivations in-script, the assertions are not adversarial — they cannot detect a sign or coefficient error in any of the four formulas.

**Required change:**
Replace the three tautological identity checks with substantive ones. Concretely:

1. **Derive `Theta_sigma` and `Lambda_phi` from a parent quadratic action.** Define a parent action `S = (1/2) Theta_sigma sigma^2 + Lambda_phi sigma phi + (KX/2) phi^2` (or whatever the stated parent action of this unit is — read it off the inline comments / docstring; the script must state it explicitly). Then assert by Gaussian elimination of sigma that the effective phi action has gain `G_micro` with the claimed form. The substantive check is then: extremize over sigma, substitute back, read off the coefficient of phi^2/2, and compare to the claimed `G_micro`. That residual can fail if `Theta_sigma` or `Lambda_phi` are wrong by a sign or factor.

2. **Coherence-factor check.** Either (a) drop the assertion (since `Osp^2 = C2 Nss Npp` is a *definition* of `C2`, not a theorem, and asserting it is round-trip nonsense), or (b) replace it with a Cauchy–Schwarz bound `0 <= C2 <= 1` if that is the script's actual claim about the support overlap, or with the operational definition of `C2` in terms of an inner product / integral the script can construct and evaluate.

3. **Xi_micro derivation.** Express `kappa` from its physical definition (e.g., `kappa = KX L^2 / TX` interpreted as a time-scale × length-scale combination), and verify the Xi formula by independently constructing it from the same parent action by switching to time-domain variables (or whichever route the comments claim). The substantive residual must connect to a piece of physics the script does not already assume in identical form.

If the unit really is just "rename variables and substitute", say so in the docstring and replace the assertions with a single sanity check that the symbolic algebra performs as expected — but mark the unit as bookkeeping, not verification.

For the Mathematica file, the same restructuring must follow, and it must derive the parent-action result independently (Cf. F3).

**Verification:**
After the rewrite, the script should contain at least one `expect_zero` whose residual is NOT obtainable by literal substitution of a single variable redefinition. The verifier inspects the file: every `expect_zero` call's LHS-minus-RHS, when written symbolically, must involve a step where (e.g.) `sp.solve`, a functional derivative, or an extremization is invoked — not merely `.subs(symbol, expr)` followed by comparison to a hand-typed twin. Sympy output should still simplify to 0 if the math is right, but the path to 0 must traverse an actual computation.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:32-44`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:28-40`

**What's wrong:**
The single substantive-looking check, `h'(rho) - m c_s^2 / rho == 0` (A1/A5), is presented in the docstring as "n=5 EOS identities". But the relation `h'(rho) = m c_s^2 / rho` is not specific to n=5; it is the trivial algebraic consequence of how the script *defines* `c_s^2` in terms of `K rho^5`. The chain is:

- `U = K rho^5 / 4` (line 35)
- `h = dU/drho = 5 K rho^4 / 4` (line 36)
- `c_s_sq = (1/m) * d(K rho^5)/drho = 5 K rho^4 / m = (4/m) * h` (line 38)

Then `m c_s^2 / rho = 4 h / rho`, and `h'(rho)` for any monomial `h = A rho^n` satisfies `h' = n A rho^(n-1) = n h / rho`. With n=4 (the degree of `h`), we get `h' = 4 h / rho = m c_s^2 / rho` automatically.

So the check is satisfied for any single-monomial polytropic EOS where `c_s^2` is defined via that particular derivative shortcut, not specifically for `n=5`. To call this an "n=5" identity the script should at minimum check a *different* polytropic index (e.g., n=3 or n=7) and observe a coefficient mismatch when the engine still uses the n=5 form, demonstrating that the assertion would catch the wrong polytropic index. As written, it does not.

Furthermore, no check ties `c_s^2` back to the standard thermodynamic definition `c_s^2 = dP/d(rho m)` (or whatever convention this unit uses). `c_s_sq` is just declared via `(1/m) d(K rho^5)/drho` with no derivation that connects this to a pressure or to the parent action's quadratic kernel for `sigma`. Without that anchor, the "n=5 EOS identity" is just an algebra exercise.

**Why this matters:**
The audit's only non-tautological-looking check has no failure mode tied to the physics. A genuine error in the polytropic index, in the prefactor `K/4`, in the chain rule connecting `U`, `p`, `h`, and `c_s^2`, would not be detected.

**Required change:**
Either (a) derive `c_s^2` from a pressure `p(rho)` (e.g., for a barotropic EOS `p = (n-1)/n * U` or whatever convention the unit uses; state it inline) and then assert `h'(rho) = c_s^2 * (dp/dh)^{-1}` or similar non-trivial chain identity, or (b) parametrize the EOS exponent as a symbol `n` and verify that the identity `h' = m c_s^2 / rho` holds with the substituted exponent only when the polytropic relation between `p` and `U` is correctly inserted, then specialize to `n=5`. Either route makes the assertion fail under a wrong exponent or wrong prefactor.

**Verification:**
The amended script should contain at least one parameter substitution (polytropic index, pressure prefactor, or both) whose value is not pre-baked into the EOS definition. Running with the wrong value should make the assertion fail. The verifier inspects the file for that variable substitution and for at least one `subs`/`/.` step that is not a tautology.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:1-73`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:1-83`

**What's wrong:**
The `.wl` file is a structural line-by-line port of the `.py` file. Same variable choreography (`U`, `h`, `h'`, `cs^2`, then `Theta_sigma`, `Lambda_phi`, `chi_eff`, `G_micro`, then coherence decomposition, then `Xi_micro`), same intermediate steps in the same order, same `expect_zero` calls in the same order, with each Python symbol renamed to a Mathematica camelCase identifier:

| sympy (py) | mathematica (wl) |
|------------|------------------|
| `U = K * rho**5 / 4` (l.35) | `uRho = capitalK*rho^5/4` (l.31) |
| `h = sp.diff(U, rho)` (l.36) | `hRho = FullSimplify[D[uRho, rho], ...]` (l.32) |
| `hprime = sp.diff(h, rho)` (l.37) | `hPrime = FullSimplify[D[hRho, rho], ...]` (l.33) |
| `cs_sq = (1/m) * sp.diff(K * rho**5, rho)` (l.38) | `csSq = FullSimplify[(1/m)*D[capitalK*rho^5, rho], ...]` (l.34) |
| `expect_zero("h'(rho) - m c_s^2 / rho", hprime - m * cs_sq / rho)` (l.44) | `expectZero["h'(rho) - m c_s^2 / rho", hPrime - m*csSq/rho]` (l.40) |
| `Theta_sigma = (m * cs_star_sq / rho_star) * Nss` (l.51) | `thetaSigma = (m*csStarSq/rhoStar)*nSS` (l.48) |
| `Lambda_phi = g_phi * Osp` (l.52) | `lambdaPhi = gPhi*oSP` (l.49) |
| `chi_eff = sp.simplify(1 / Theta_sigma)` (l.53) | `chiEff = FullSimplify[1/thetaSigma, ...]` (l.50) |
| `G_micro = sp.simplify(chi_eff * Lambda_phi**2 / KX)` (l.54) | `gMicro = FullSimplify[chiEff*lambdaPhi^2/kX, ...]` (l.51) |
| `Xi_micro.subs(kappa, KX * L**2 / TX)` (l.80) | `xiMicro /. kappa -> kX*ell^2/tX` (l.68) |

There is no place where the Mathematica file performs a step the Python file does not, or vice versa; the two engines proceed through the exact same algebra in the same order. This violates the second-engine policy: both engines must derive the result independently from physical premises, not echo each other's choreography. Note that this also amplifies F1 — both engines being tautological in lockstep makes the duplicated tautology look like agreement.

A truly independent Mathematica derivation would, for example, set up the parent action symbolically, use `Solve`/`Reduce` to extremize over sigma, and read off the phi gain coefficient — taking a path the SymPy file does not take. Or it could verify `G_micro` numerically by sampling `Theta_sigma`, `Lambda_phi`, `KX` at random positive values and checking the formula, providing orthogonal evidence to the SymPy symbolic check.

**Why this matters:**
The second engine exists to catch errors specific to one engine's simplifier or to the author's choice of algebraic path. A transliteration provides neither: any algebraic error common to both files (which is most of them) is invisible.

**Required change:**
After F1's restructuring lands, the Mathematica file must derive the gain via an independent algebraic route from the same parent action. Options include:
- Symbolic extremization: define `S[sigma, phi]` from the parent action, take `D[S, sigma]`, solve for `sigma*[phi]`, substitute back to get `S_eff[phi]`, and read the coefficient of `phi^2/2`. Compare to the claimed `G_micro`.
- Or perform a numerical sanity check: substitute random rational values for `m, cs_star_sq, rho_star, Nss, g_phi, Osp, KX` and verify both the parent extremization route and the closed-form formula give the same number.

Specific edit: at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:48-58`, replace the direct definitions with an extremization computation that produces `gMicro` as an output of `Solve` / `Coefficient`, not as a hand-typed product.

**Verification:**
The new Mathematica file must contain at least one `Solve`/`Reduce`/`Coefficient`/`Series`/`Limit` call that is *not* present in the SymPy file, and must derive `gMicro` (or `xiMicro`) via that call rather than by direct product. The verifier inspects the diff: the .wl file should no longer be the same sequence of definitions as the .py file.

## Independent-derivation check (Mathematica)

Not independent — see F3. The `.wl` is a near-mechanical port of the `.py` with renamed identifiers. The two files do not provide independent confirmation; they re-execute the same algebra.

## Engine cross-check

Both engines produce identical final residuals (all four `expect_zero` residuals print `0` in both outputs), and the printed intermediate expressions match modulo notation:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `h'(rho)` | `5*K*rho**3` | `5*capitalK*rho^3` |
| `c_s^2(rho)` | `5*K*rho**4/m` | `(5*capitalK*rho^4)/m` |
| `Theta_sigma` | `N_ss*cs_star_sq*m/rho_star` | `(csStarSq*m*nSS)/rhoStar` |
| `chi_sigma^(eff)` | `rho_star/(N_ss*cs_star_sq*m)` | `rhoStar/(csStarSq*m*nSS)` |
| `G_micro` | `O_sp**2*g_phi**2*rho_star/(K_X*N_ss*cs_star_sq*m)` | `(gPhi^2*oSP^2*rhoStar)/(csStarSq*kX*m*nSS)` |
| `Xi_micro` | `O_sp**2*g_phi**2*kappa*rho_star/(K_X*N_ss*cs_star_sq*m)` | `(gPhi^2*kappa*oSP^2*rhoStar)/(csStarSq*kX*m*nSS)` |

Numerical agreement is total. This is expected given F3 (they're the same algebra). It does not constitute independent verification.

## Verdict justification

The audit verdict is **findings** (not stop_cold).

What holds up: both scripts run to `exit 0`, both produce identical (renamed) intermediate forms, the outputs are fresh (sympy output mtime 2026-05-11 12:44 ≥ script mtime 2026-04-01 12:39; mathematica output mtime 2026-05-11 12:54 ≥ script mtime 2026-05-11 11:56). The algebra inside each script is internally consistent.

What does not hold up: three of four assertions are textbook tautologies (define `X`, define `X_expected := X`, assert `X - X_expected == 0`); the fourth, framed as an "n=5 EOS identity", is the trivial monomial-degree identity for any degree-4 `h(rho)` and is not specific to the polytrope choice; and the Mathematica file is a line-by-line transliteration of the SymPy file rather than an independent derivation. The audit therefore does not provide adversarial evidence for the unit's claimed "parent-action projection" result.

Attacks attempted and what they found:
- Try to make A2 (`G_micro - G_expected`) fail by perturbing the parent definitions: impossible — `G_expected` is hand-typed as the algebraic simplification of `G_micro`'s definition. Any change to `Theta_sigma` or `Lambda_phi` propagates identically into both LHS and RHS only if the script is rewritten to do so; as written, `G_expected` is a fixed literal copy of the simplified product, so the symbols on both sides are *the same* symbols, not independent quantities.
- Try to make A3 fail by varying the substitution: the substitution `Osp^2 -> C2 Nss Npp` is applied to `G_expected` (one literal) to produce `G_coh` (a second literal). The two are designed to match. No physical claim about the support overlap is being tested.
- Try to make A4 fail by varying `kappa`: `kappa` is substituted to `KX L^2/TX` in `Xi_micro`, which produces an expression that equals the hand-typed `Xi_expected` by inspection. Same structure as A2.
- Try to make A1/A5 fail by varying the polytropic index: the script hardcodes `rho^5`. The identity it asserts holds for any monomial EOS of the appropriate degree, so changing the exponent in *both* `U` and `cs_sq` definitions leaves the check passing. Only an inconsistent change (different exponents in `U` vs. in `cs_sq`'s derivative) would break it, but the script ties the two together via the same `K rho^5` term.

The verdict is `findings` because the math is not wrong — it is genuinely simplifying to zero — but the assertions, as designed, cannot detect any plausible physics error. This is `tautological_check` + `insufficient_verification` + `mathematica_transliteration`. The fixes are well-defined and do not propagate to downstream units (the closed-form `G_micro = rho_star g_phi^2 Osp^2/(m cs_star_sq KX Nss)` is presumably stable; what's needed is a substantive in-script derivation rather than a different formula). Therefore not `CRITICAL_DOWNSTREAM`.

## Self-test notes

- **Tautology check**: I confirmed for each `expect_zero` call that the LHS is constructed by substituting variables that already appear in the RHS — A2, A3, A4 all reduce to `0` by direct substitution, no algebraic identity outside of the script's own definitions is exercised.
- **Monomial-degree trap (A1)**: I confirmed that `h'(rho) = m c_s^2/rho` is the n=4 monomial identity `rho h' = 4 h`, given the script's chain `c_s^2 = (4/m) h`. Not specific to n=5.
- **Path specifications for directive**: F1/F2 target existing files in `scripts/` and `mathematica/`; F3 targets the existing `.wl` file. No `missing_verification_script` finding, so no new file paths to specify.
