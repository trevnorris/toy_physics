---
unit_id: 039
batch: III.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T07:59:49Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 039 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage039_split_u_sector_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage039_split_u_sector_mathematica_audit.txt`

## What the script claims to verify

Per the SymPy docstring and inline comments, the unit asserts that turning on the first axial U-sector stiffness (`deltaU`) (i) splits the flat U-doublet so the direct wall anisotropy ratio becomes an exact shifted form `delta_split = (delta0 + eps_eta*deltaU/(1+deltaU))/(1-eps_eta)`, (ii) renormalizes the mixed U/W blocking ratio to `eps_W_split = eps_W*(1 - 2/11*deltaU/(1+deltaU))`, (iii) yields a non-collinear mixed-loading vector whose direction-splitting invariant `D_dir = kappa0*z1 - kappa1*z0` is exactly `-kappa0*kappa1*g_W*rho0*deltaU/(1+deltaU)`, and (iv) preserves the Stage-21 placement-map product law `M_mix * R_target = 8*Lambda*(1-eps_W_split)/pi^2` after substituting `eps_W -> eps_W_split` and `delta -> delta_split`. Section 22.5 prints (without asserting) the linear-order series of all the above in `deltaU`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 86 | `expect_zero("A0 direct - expected", A0.subs({c_etaU**2: eps_eta*K_U*K_eta_eff}) - A0_expected)` | yes |
| A2 | sympy | 87 | `expect_zero("A1 direct - expected", A1.subs({c_etaU**2: eps_eta*K_U*K_eta_eff}) - A1_expected)` | yes |
| A3 | sympy | 98-101 | `expect_zero("eps_W direct - split formula", eps_W_direct.subs(...) - eps_W_split)` | yes |
| A4 | sympy | 114 | `expect_zero("z1/z0 - (kappa1/kappa0) R_U", z1/z0 - (kappa1/kappa0)*R_U)` | **no (tautological)** |
| A5 | sympy | 119 | `expect_zero("direction-splitting invariant", D_dir - D_dir_expected)` | partial |
| A6 | sympy | 133 | `expect_zero("product law", product - 8*Lambda*(1-eps_W_split)/pi**2)` | **no (tautological by construction)** |
| A7 | mma | 66 | `expectZero["A0 direct - expected", (a0 /. cEtaU^2 -> epsEta*kU*kEtaEff) - a0Expected]` | yes |
| A8 | mma | 67 | `expectZero["A1 direct - expected", (a1 /. cEtaU^2 -> epsEta*kU*kEtaEff) - a1Expected]` | yes |
| A9 | mma | 77 | `expectZero["eps_W direct - split formula", (epsWDirect /. cUW^2 -> ...) - epsWSplit]` | yes |
| A10 | mma | 89 | `expectZero["z1/z0 - (kappa1/kappa0) R_U", z1/z0 - (kappa1/kappa0)*rU]` | **no (tautological)** |
| A11 | mma | 94 | `expectZero["direction-splitting invariant", dDir - dDirExpected]` | partial |
| A12 | mma | 105 | `expectZero["product law", product - 8*lambda*(1-epsWSplit)/Pi^2]` | **no (tautological by construction)** |

The "no/tautological" rows feed F1 and F2 below. The "partial" rows are direct algebraic identities that follow without freedom from the postulated forms, but the postulated `D_dir_expected` formula at least has to come out with the right sign and `deltaU/(1+deltaU)` factor — they get a pass.

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py:107-114`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:82-89`

**What's wrong:**

In Section 22.3 the script constructs

```
z0 = kappa0 * g_W * (1 + rho0)
z1 = kappa1 * g_W * (1 + rho0/(1+deltaU))
R_U = (1 + rho0/(1+deltaU)) / (1 + rho0)
```

then asserts `expect_zero("z1/z0 - (kappa1/kappa0) R_U", z1/z0 - (kappa1/kappa0)*R_U)` (sympy line 114; Mathematica line 89). But `z1/z0` evaluates to `(kappa1/kappa0) * (1 + rho0/(1+deltaU))/(1+rho0) = (kappa1/kappa0)*R_U` by *substitution of the literal definitions*. No physics can make it fail — `R_U` is defined to be exactly the ratio of the rho-factors that appear in `z0` and `z1`. The check just confirms the definition against itself.

**Why this matters:**

The output line `z1/z0 - (kappa1/kappa0) R_U = 0` is currently reported as evidence that the mixed-loading vector satisfies a derived ratio relation, but it is actually evidence of nothing more than `a/b = a/b`. If the verifier ever audits "what is the non-tautological consequence of split-U here," this assertion contributes none.

**Required change:**

Replace the tautological identity with a substantive check that ties `R_U`'s structure to an independently-named physical input. Concretely, derive `z0`, `z1` from the source-vector definition `g_W = c_etaW/sqrt(mu_eta*mu_W)` and the kappa basis, and verify that the *ratio* `z1/z0` equals an expression that does not contain `R_U` as a defined symbol — e.g. assert

```
expect_zero(
    "z1/z0 collinearity gap",
    z1 * (1 + rho0) - z0 * (kappa1/kappa0) * (1 + rho0/(1 + deltaU))
)
```

which, after substituting the constructed `z0`, `z1`, must vanish identically only if the kappa-weighted rho-factor structure is correct. Then *separately* introduce `R_U` as a named shorthand and observe that the gap is `(kappa1/kappa0)*R_U`, without an `expect_zero` on the trivial identity.

Equivalently, drop the existing assertion (sympy:114, mma:89) entirely if it is not replaced — the printout `R_U = ...` already communicates the intended observation.

**Verification:**

The verifier should see either (a) the existing line `expect_zero("z1/z0 - (kappa1/kappa0) R_U", ...)` removed from sympy:114 and mma:89, or (b) replaced with a check whose left-hand side does not equal its right-hand side by direct substitution of named definitions. The new check should pass when the kappa-rho structure is correct and fail under, e.g., perturbing the `(1 + rho0/(1+deltaU))` factor in `z1`.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py:126-133`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:98-105`

**What's wrong:**

In Section 22.4 the script defines

```
M_mix_split   = 8 * Z_W * (1 + rho0)^2 / (pi^2 * (1 - eps_eta) * (1 - eps_W_split))
R_target_split = Lambda * (1 - eps_eta) * (1 - eps_W_split)^2 / (Z_W * (1 + rho0)^2)
```

then asserts `expect_zero("product law", M_mix_split * R_target_split - 8 * Lambda * (1 - eps_W_split) / pi^2)` (sympy line 133; Mathematica line 105). The product literally cancels: `Z_W` cancels, `(1+rho0)^2` cancels, `(1-eps_eta)` cancels, one factor of `(1-eps_W_split)` cancels. The residual is exactly `8*Lambda*(1-eps_W_split)/pi^2` by elementary algebra, irrespective of any physics. The factors `Z_W`, `(1+rho0)^2`, `(1-eps_eta)`, `(1-eps_W_split)` were chosen to appear with the *exact* powers needed for the product to reduce, and the assertion just confirms that the chosen powers indeed cancel.

**Why this matters:**

The script's theorem ledger (sympy:160-161) cites this as "the Stage-21 factorization survives at the scalar placement level," but the assertion currently reduces to `(a*b/c) * (c*b/a) = b^2 * (something) = ...` — pure schoolbook cancellation with no input from the split-U sector. Whatever derivation actually established the placement-map structure happened upstream of this script; A6/A12 verify nothing about Stage 039.

**Required change:**

Replace the identity check with a substitution-based check that links `M_mix_split`, `R_target_split` to their Stage-21 (flat-U) counterparts. Concretely, the claim that the script *says* it is making (sympy docstring item 5: "Stage-21 continuum placement map survives at the scalar level with `eps_W -> eps_W_split` and `delta -> delta_split`") implies there should be flat-U expressions `M_mix_flat` and `R_target_flat` such that

```
M_mix_split   == M_mix_flat   .subs(eps_W, eps_W_split)
R_target_split == R_target_flat.subs(eps_W, eps_W_split)
```

This is what should be asserted. Define `M_mix_flat = 8*Z_W*(1+rho0)^2/(pi^2*(1-eps_eta)*(1-eps_W))` and `R_target_flat = Lambda*(1-eps_eta)*(1-eps_W)^2/(Z_W*(1+rho0)^2)` once, then `expect_zero("M_mix split is M_mix_flat under eps_W->eps_W_split", M_mix_split - M_mix_flat.subs(eps_W, eps_W_split))` and similarly for `R_target`. Keep the product printout but drop the trivial `expect_zero("product law", ...)`.

Equivalent in Mathematica: define `mMixFlat`, `rTargetFlat`, then `expectZero["M_mix split substitution", mMixSplit - (mMixFlat /. epsW -> epsWSplit)]`.

**Verification:**

The verifier should see at sympy:133 (and mma:105) that the `expect_zero("product law", ...)` is either removed or replaced by substitution checks `M_mix_split - M_mix_flat.subs(...)` and `R_target_split - R_target_flat.subs(...)`. The new checks must equal zero only when the substitution structure holds; perturbing the exponent on `(1-eps_W_split)` in either `M_mix_split` or `R_target_split` should make at least one of them fail.

### F3 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:33-126`

**What's wrong:**

The Mathematica script is structurally a line-by-line port of the SymPy script's algebraic choreography, not an independent re-derivation. Side-by-side correspondence:

| SymPy block | Mathematica block | Identical algebra? |
|---|---|---|
| `K_U1 = K_U*(1+deltaU)` (line 71) | `kU1 = kU*(1+deltaU)` (line 55) | yes |
| `A0 = (K_eta_eff - c_etaU**2/K_U)/mu_eta`; `A1 = (...(1+delta0) - c_etaU**2/K_U1)/mu_eta` (72-73) | `a0 = ...`; `a1 = ...` (56-57) | yes |
| `delta_split = (delta0 + eps_eta*deltaU/(1+deltaU))/(1-eps_eta)` (78) | `deltaSplit = (delta0 + epsEta*deltaU/(1+deltaU))/(1-epsEta)` (59) | yes |
| `A0.subs({c_etaU**2: eps_eta*K_U*K_eta_eff})` (86) | `a0 /. cEtaU^2 -> epsEta*kU*kEtaEff` (66) | yes (literal substitution rule, same direction) |
| `S_U = kappa0**2/K_U + kappa1**2/K_U1` (92) | `sU = kappa0^2/kU + kappa1^2/kU1` (71) | yes |
| `eps_W_split = eps_W*(1 - 2/11*deltaU/(1+deltaU))` (94) | `epsWSplit = epsW*(1 - (2/11)*deltaU/(1+deltaU))` (73) | yes (`2/11` postulated identically) |
| `D_dir = kappa0*z1 - kappa1*z0` (116) | `dDir = kappa0*z1 - kappa1*z0` (91) | yes |
| `M_mix_split = 8*Z_W*(1+rho0)**2/(pi**2*(1-eps_eta)*(1-eps_W_split))` (126) | `mMixSplit = 8*zW*(1+rho0)^2/(Pi^2*(1-epsEta)*(1-epsWSplit))` (98) | yes |

Variable names are camelCased and the substitution operator changes from `.subs({a: b})` to `/. a -> b`, but every algebraic step is in the same order with the same intermediate forms. There is no Mathematica-side independent derivation (e.g., constructing `S_U` from a different basis decomposition, or solving the U-sector splitting via `Solve[...]` rather than direct postulation) that could disagree with the SymPy script through any path other than a typo. The second-engine policy requires both engines to derive the result independently from the physical premises so that a mistake in either chain can surface as `engine_disagreement`; this script cannot surface such a disagreement because the chain is shared.

**Why this matters:**

The independent-engine check provides almost no additional confidence. If the SymPy author's choice of postulated forms (e.g., the `2/11` coefficient on `deltaU/(1+deltaU)`, or the `R_U` ratio definition) had a sign or factor error rooted in a misunderstanding of the underlying PDE, the Mathematica script — which copied those postulated forms — would still PASS. This is exactly the failure mode the second-engine policy is meant to catch.

**Required change:**

Restructure the Mathematica script so that the key derived quantities — `epsWSplit`, `deltaSplit`, `dDirExpected`, and at least one of the placement-map factors `mMixSplit` / `rTargetSplit` — are obtained by Mathematica's own algebra from the *direct* expressions, rather than postulated and then checked. Concrete edits:

1. At mma:73, remove the postulated `epsWSplit = FullSimplify[epsW*(1 - (2/11)*deltaU/(1+deltaU)), ...]`. Instead, *derive* it: `epsWSplit = FullSimplify[(epsWDirect /. cUW^2 -> epsW*kU*kWEff/sigma), Assumptions -> $Assumptions]`. Then the assertion at mma:77 should be `expectZero["eps_W split agrees with postulated form", epsWSplit - epsW*(1 - (2/11)*deltaU/(1 + deltaU))]` — i.e., the postulated form is on the right, the Mathematica-derived form is on the left. This way the `2/11` coefficient is computed by Mathematica from `kappa0^2`, `kappa1^2`, `sigma`, not assumed.

2. At mma:59, similarly derive `deltaSplit` from solving `a0Expected*(1 + d) == a1Direct` for `d`, where `a1Direct = (a1 /. cEtaU^2 -> epsEta*kU*kEtaEff)`. Then assert that the postulated closed form matches.

3. At mma:92, derive `dDirExpected` by `FullSimplify[Together[kappa0*z1 - kappa1*z0]]`, no postulated form. Then assert it equals the closed form `-kappa0*kappa1*gW*rho0*deltaU/(1+deltaU)`.

These changes make the Mathematica script independent at exactly the points where the SymPy script postulated answers, so a typo in either side's postulate now produces an `engine_disagreement` instead of mutual silence.

**Verification:**

The verifier should see at mma:59, 73, 92 that the postulated closed forms are no longer the *left-hand side* of the definition but the *right-hand side* of an `expectZero` check, with the left-hand side being obtained by Mathematica's `FullSimplify`/`Solve`/`Together` on the direct expressions. Output should print the Mathematica-derived form before the assertion.

## Independent-derivation check (Mathematica)

The `.wl` script does not derive the claims independently — it transliterates the SymPy script's algebraic choreography step-for-step (see F3 for the side-by-side table). All key derived quantities (`epsWSplit`, `deltaSplit`, `dDirExpected`, the `M_mix`/`R_target` factors) are postulated with the same closed forms as in the SymPy script, then verified against direct constructions that themselves match the SymPy direct constructions. There is no Mathematica-only intermediate step that could disagree with SymPy through an independent algebraic path.

## Engine cross-check

Both engines complete with `Exit 0` and report the same simplified outputs:

| Quantity | SymPy output | Mathematica output | Agree? |
|---|---|---|---|
| `kappa0` | `2*sqrt(2)/pi` | `(2*Sqrt[2])/Pi` | yes |
| `kappa1` | `-4/(3*pi)` | `-4/(3*Pi)` | yes |
| `sigma` | `88/(9*pi**2)` | `88/(9*Pi^2)` | yes |
| `lambda0` | `2/9` | `2/9` | yes |
| `S_U` | `8*(9*deltaU + 11)/(9*pi**2*K_U*(deltaU + 1))` | `(8*(11 + 9*deltaU))/(9*(1+deltaU)*kU*Pi^2)` | yes |
| `eps_W_split` | `eps_W*(9*deltaU+11)/(11*(deltaU+1))` | `epsW - (2*deltaU*epsW)/(11*(1+deltaU))` | yes (algebraically equal) |
| `R_U` | `(deltaU + rho0 + 1)/((deltaU+1)*(rho0+1))` | `(1+deltaU+rho0)/(1+deltaU+rho0+deltaU*rho0)` | yes (`1+deltaU+rho0+deltaU*rho0 = (1+deltaU)(1+rho0)`) |
| `D_dir` | `8*sqrt(2)*c_etaW*deltaU*rho0/(3*pi**2*sqrt(mu_W*mu_eta)*(deltaU+1))` | `(8*Sqrt[2]*cEtaW*deltaU*rho0)/(3*(1+deltaU)*Sqrt[muEta*muW]*Pi^2)` | yes |
| `product` (final) | `8*Lambda*(11*deltaU - eps_W*(9*deltaU+11) + 11)/(11*pi**2*(deltaU+1))` | `(8*(11+11*deltaU - 11*epsW - 9*deltaU*epsW)*lambda)/(11*(1+deltaU)*Pi^2)` | yes |

All six assertions in each script PASS. But — as F3 notes — this agreement is forced by the transliteration, not earned by independent derivation, so the engine-cross-check is weak evidence here.

## Verdict justification

The script is mathematically consistent in the algebra it does perform — the substitution `c_etaU^2 -> eps_eta*K_U*K_eta_eff` correctly yields the postulated `delta_split` and `eps_W_split` closed forms, the direction-splitting invariant comes out with the right sign and `deltaU/(1+deltaU)` factor, and the two engines agree on every simplified expression. However, two of the six assertion pairs (`z1/z0 = (kappa1/kappa0)*R_U` and the product law `M_mix*R_target = 8*Lambda*(1-eps_W_split)/pi^2`) are tautological by construction (F1, F2), and the Mathematica script is a step-by-step transliteration of the SymPy choreography rather than an independent derivation (F3). Verdict is `findings` with three actionable items; no `stop_cold` flag because none of the fixes propagate downstream — they strengthen the verification of *this* unit's existing claims rather than changing any derived constant or sign.

Attacks attempted that failed: (a) trying to break the `2/11` coefficient by recomputing `kappa0^2/sigma`, `kappa1^2/sigma` — both come out exactly to `9/11` and `2/11`, consistent; (b) trying to find a sign error in `D_dir` by direct expansion — sign is correct; (c) trying to find a domain-assumption mismatch (e.g., requiring `eps_eta < 1` for `(1-eps_eta)` in the denominator) — the algebraic identities in the script don't actually require positivity of `(1-eps_eta)`, only that it be nonzero, and the `simplify` calls don't exploit positivity in a way that would hide a branch error; (d) trying to find a parity / domain mismatch in the series expansions (22.5) — these are just `removeO()` prints with no assertion, so they cannot fail in a hidden way.

## Self-test notes

- **Variable independence**: No new `sp.diff` or `D[...]` is proposed in any "Required change"; all proposed checks are algebraic identities of the same symbols already in the script. The substitution-based checks in F2's "Required change" (`M_mix_split - M_mix_flat.subs(eps_W, eps_W_split)`) involve only `eps_W`, `eps_W_split`, both already declared as positive reals — no new symbol introductions, no derivatives, so the variable-independence trap does not apply.
- **Parity / trivial-case**: Setting `deltaU -> 0` in the proposed substitution check `M_mix_split - M_mix_flat.subs(eps_W, eps_W_split)` gives `M_mix_split|_{deltaU=0} = M_mix_flat` (since `eps_W_split|_{deltaU=0} = eps_W`), so the residual collapses to 0 as required. For F1's proposed `z1*(1+rho0) - z0*(kappa1/kappa0)*(1+rho0/(1+deltaU))`: substituting the literal definitions gives `kappa1*g_W*(1+rho0/(1+deltaU))*(1+rho0) - kappa0*g_W*(1+rho0)*(kappa1/kappa0)*(1+rho0/(1+deltaU)) = 0` — correct, this is the *non-tautological* form because the rho-factor structure is no longer encapsulated in a defined symbol `R_U`.
- **Path specifications**: The directive's "Target" lines name `.py` files in `scripts/` and `.wl` files in `mathematica/`, matching the canonical layout.
