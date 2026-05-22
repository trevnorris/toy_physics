---
unit_id: 047
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 047 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.txt`

## What the script claims to verify

The two scripts (filenamed "stage047", docstrings still calling it "Stage 30") claim to verify, for a coherent-kernel map of the moving-throat PDE: (1) the three coherent interference ratios collapse to a common value, rho_0 = sigma_0 = chi_0 = gamma*c_etaU/KU; (2) exact proportionality eps_phi = zeta_def*eps_W and Z_phi = zeta_def*Z_W; (3) the exact total baseline factorization M_tr = M_mix*S(zeta;eps); (4) the support-loaded product law R_target*M_tr = (8 Lambda (1-eps)/pi^2)*S together with the resulting zeta-blindness of the loaded R_target reconstruction; (5) monotonicity dS/dzeta = (1-eps)/(1-zeta*eps)^2 and the boundary value S(0)=1. A numeric "spoiler" probe at the end perturbs M_supp off the exact coefficient to demonstrate the zeta-blindness is non-vacuous.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 50 | `expect_zero("rho_0 - chi_0", rho0 - chi0)` | no (tautological — see F1) |
| A2 | sympy | 51 | `expect_zero("sigma_0 - chi_0", sigma0 - chi0)` | no (tautological — see F1) |
| A3 | sympy | 73 | `expect_zero("eps_phi - zeta_def*eps_W", eps_phi - zeta_def*epsW)` | yes |
| A4 | sympy | 74 | `expect_zero("Z_phi - zeta_def*Z_W", Zphi - zeta_def*ZW)` | yes |
| A5 | sympy | 99 | `expect_zero("M_tr - M_mix*S", Mtr - Mmix*S)` | yes |
| A6 | sympy | 109 | `expect_zero("product law", product_actual - product_expected)` | yes |
| A7 | sympy | 110 | `expect_zero("support-loaded R_target reconstruction", Rtarget_loaded - Rtarget)` | yes |
| A8 | sympy | 111 | `expect_zero("dR_target_loaded/dzeta", sp.diff(Rtarget_loaded, zeta))` | partial (redundant with A5+A6, but non-trivial) |
| A9 | sympy | 112-115 | `expect_zero("R_target_loaded(zeta)-R_target_loaded(0)", ...)` | partial (redundant with A8) |
| A10 | sympy | 116 | `expect_zero("dS/dzeta - (1-eps)/(1-zeta eps)^2", ...)` | yes |
| A11 | sympy | 117 | `expect_zero("S(zeta=0)-1", S.subs(zeta,0)-1)` | yes |
| A12 | sympy | 153-154 | spoiler numeric: `abs(dR_spoiled/dzeta) < 1e-12 -> fail` | yes |
| A13 | mathematica | 43 | `expectZero["rho_0 - chi_0", rho0 - chi0]` | no (tautological — see F1) |
| A14 | mathematica | 44 | `expectZero["sigma_0 - chi_0", sigma0 - chi0]` | no (tautological — see F1) |
| A15 | mathematica | 61 | `expectZero["eps_phi - zeta_def eps_W", epsPhi - zetaDef*epsW]` | yes |
| A16 | mathematica | 62 | `expectZero["Z_phi - zeta_def Z_W", zPhi - zetaDef*zW]` | yes |
| A17 | mathematica | 83 | `expectZero["M_tr - M_mix S", mTr - mMix*sEnhance]` | yes |
| A18 | mathematica | 90 | `expectZero["product law", productActual - productExpected]` | yes |
| A19 | mathematica | 91 | `expectZero["support-loaded R_target reconstruction", loadedRTarget - rTarget]` | yes |
| A20 | mathematica | 92 | `expectZero["dR_target_loaded/dzeta", D[loadedRTarget, zeta]]` | partial |
| A21 | mathematica | 93-96 | `expectZero["R_target_loaded(zeta) - R_target_loaded(0)", ...]` | partial |
| A22 | mathematica | 97 | `expectZero["dS/dzeta - (1-eps)/(1-zeta eps)^2", ...]` | yes |
| A23 | mathematica | 98 | `expectZero["S(zeta=0)-1", ...]` | yes |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py:42-51`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl:36-44`

**What's wrong:**
The "coherent interference identities" rho_0 = sigma_0 = chi_0 are not derived — they are hardcoded with the cancelling factor written in by hand. From the SymPy script:

```
rho0   = sp.simplify((gamma * lamW) * c_etaU / (KU * lamW))     # line 42
sigma0 = sp.simplify(c_etaU * (gamma * lamphi) / (KU * lamphi)) # line 43
chi0   = sp.simplify(gamma * c_etaU / KU)                        # line 44
```

The `lamW/lamW` and `lamphi/lamphi` factors cancel by construction before `simplify` ever runs, so all three expressions are identically `gamma*c_etaU/KU` as constructed. The Mathematica script (lines 36-38) mirrors this exactly:

```
rho0   = FullSimplify[(gamma*lamW)*cEtaU/(kU*lamW), ...];
sigma0 = FullSimplify[cEtaU*(gamma*lamPhi)/(kU*lamPhi), ...];
chi0   = FullSimplify[gamma*cEtaU/kU, ...];
```

Consequently the assertions `expect_zero("rho_0 - chi_0", ...)` and `expect_zero("sigma_0 - chi_0", ...)` (sympy lines 50-51, wl lines 43-44) cannot fail no matter what the upstream physics is: the symbol-cancellation guarantees the equality at construction time. The captured outputs confirm this — both print `rho_0 = c_etaU*gamma/KU` and `sigma_0 = c_etaU*gamma/KU` (identical reduced forms), not the un-cancelled microscopic numerators that a genuine collapse identity would produce as an intermediate.

**Why this matters:**
The unit's docstring lists "Coherent-kernel identities rho_0 = sigma_0 = chi_0" as check #1, but no independent definitions of rho_0 and sigma_0 are constructed from the underlying coherent kernels (e.g., from the W-channel and phi-channel polarisations separately) — they're written directly as chi_0 times a unit-valued factor. The identity claim is a no-op; any putative regression in the coherent-kernel physics that affected rho_0 or sigma_0 via lamW or lamphi would not be detected by these two assertions.

**Required change:**
Construct rho_0 and sigma_0 from independent microscopic definitions that exercise their channel-specific coupling, then assert collapse to chi_0. Replace the trivially cancelling forms in both engines:

- sympy `:42-43`: Define rho_0 as a coherent W-channel polarisation ratio and sigma_0 as a coherent phi-channel polarisation ratio that do NOT pre-cancel their channel coupling. The most defensible minimal change is to introduce independent channel kernels K_W_rho and K_phi_sigma such that rho_0 = (gamma * lamW^2 / (KU * lamW * (lamW/K_W_rho * K_W_rho))) collapses to gamma*c_etaU/KU only after substituting the channel saturation rules (lamW * c_etaU/lamW -> c_etaU). If no such independent definition can be constructed from the unit's stated scope, then mark rho_0 and sigma_0 as auxiliary aliases (`rho0 = chi0`, `sigma0 = chi0`) and remove the two assertions; the misleading "identity check" is worse than no check.
- wl `:36-37`: Mirror the same change independently — derive rho_0 and sigma_0 from a separately stated channel saturation rule (not from the SymPy script's algebra).

**Verification:**
After the fix, the new construction expressions should not have lamW or lamphi appearing as both a numerator and denominator factor in the same product. The captured intermediate `rho_0 = ...` print line should display a non-trivial expression (not the already-reduced `c_etaU*gamma/KU`) before the simplification step, demonstrating that the collapse is performed by simplification rather than by string cancellation. Both engines should exit 0.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl:1-103` (entire file)

**What's wrong:**
The `.wl` script is a line-by-line port of the `.py` script: identical variable order, identical pre-cancelled forms, identical assertion list in the same sequence, no independent derivation of any quantity from physical premises. Three corresponding sections demonstrate this:

(a) Pre-cancelled rho_0/sigma_0 (already cited in F1):
- `.py:42-43`:
  ```
  rho0   = sp.simplify((gamma * lamW) * c_etaU / (KU * lamW))
  sigma0 = sp.simplify(c_etaU * (gamma * lamphi) / (KU * lamphi))
  ```
- `.wl:36-37`:
  ```
  rho0   = FullSimplify[(gamma*lamW)*cEtaU/(kU*lamW), ...];
  sigma0 = FullSimplify[cEtaU*(gamma*lamPhi)/(kU*lamPhi), ...];
  ```
Same numerator/denominator choreography (`(gamma*lamW)*cEtaU/(kU*lamW)`); a Mathematica re-derivation would either expand from a Fourier integral or write the simplified `gamma*cEtaU/kU` directly.

(b) M_supp coefficient identical:
- `.py:88`: `Msupp = sp.simplify(8 * zeta * ZW * (1 + chi0) ** 2 / (sp.pi ** 2 * (1 - eps_eta) * (1 - zeta * eps)))`
- `.wl:73`:  `mSupp = FullSimplify[8*zeta*zW*(1 + chi0)^2/(Pi^2*(1 - epsEta)*(1 - zeta*eps)), ...]`
The pre-factor `8`, the `(1+chi0)^2` structure, and the `(1-zeta*eps)` denominator are byte-identical up to variable renaming. The `.wl` does not derive this from a phi-channel integral; it copies the SymPy expression.

(c) Spoiler/eta_bad probe (.py:121-154) has no .wl counterpart, but the spoiler is the one place a true independent check would help — and it's only present in one engine. The .wl file ends at line 98 with the trivial S(0)-1 check, mirroring `.py:117`.

The order of assertions in both files is identical (A1=A13, A2=A14, … A11=A23), which is the strongest signal of transliteration.

**Why this matters:**
The two-engine policy exists so that a derivation error in one symbolic system gets caught by the other. A transliteration cannot catch any algebraic mistake — both engines will produce the same answer because they're running the same algebra under different syntax. The agreement of outputs (which is what the SymPy and Mathematica `.txt` files show) is therefore not evidence of correctness.

**Required change:**
Re-derive at least one of the load-bearing quantities in the `.wl` script independently of the `.py` script's algebra. The minimum sufficient change: in `mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl`, replace the direct definitions of `mMix` and `mSupp` (currently at lines 72-73) with definitions obtained by independently integrating a coherent-kernel weight, OR by writing M_mix and M_supp from a different but algebraically equivalent factorization (e.g., compute M_supp as `mMix * zeta * (1-eps)/(1-zeta*eps)` from the S-decomposition first, then confirm it equals the original closed form). If no independent derivation is available within the unit's scope, add an explicit per-line comment in the `.wl` file justifying each definition as an independent statement of the coherent kernel result, with a citation to the upstream unit whose first-principles derivation the .wl is restating.

**Verification:**
The verifier confirms by reading the `.wl` file and checking that at least one of `mMix` or `mSupp` has a derivation path (intermediate expressions, comments referencing physical input) not appearing in the `.py` file. The captured output `M_supp = ...` line should match the SymPy output's simplified form (engines should still agree after FullSimplify), but the symbolic construction path leading to it must differ from the `.py` script's path.

## Independent-derivation check (Mathematica)

The `.wl` script is not an independent derivation of the stage-047 claims; it is a structural transliteration of the `.py` script. Both files define the same fifteen intermediate quantities (rho_0, sigma_0, chi_0, eps_eta, eps_W, zeta_def, Z_W, delta0, deltaU, Lambda, eps_phi, Z_phi, eps, R_tr, delta, M_mix, M_supp, S, M_tr, R_target) in the same order with the same closed-form right-hand sides (only the variable name spellings differ, e.g., `KU`/`kU`, `lamphi`/`lamPhi`, `T_w`/`tw`, `Keta_eff`/`kEtaEff`). The assertion list is in the same sequence. The pre-cancelled tautological construction of rho_0 and sigma_0 (see F1) is identical in both files. See F2 for the recommended fix.

## Engine cross-check

Both engines run to `EXIT_CODE: 0` and produce algebraically equivalent simplified forms for every printed intermediate. Concrete comparison:

| quantity | sympy | mathematica |
|---|---|---|
| rho_0 | `c_etaU*gamma/KU` | `(cEtaU*gamma)/kU` |
| zeta_def | `KW_eff*lamphi**2/(Kphi_eff*lamW**2)` | `(kWEff*lamPhi^2)/(kPhiEff*lamW^2)` |
| eps_W | `88*gamma**2*lamW**2/(9*pi**2*KU*KW_eff)` | `(88*gamma^2*lamW^2)/(9*kU*kWEff*Pi^2)` |
| R_target | `G*c_s**5*(KU*Keta_eff - c_etaU**2)*(9*pi**2*KU*KW_eff*(KU*L**2 + pi**2*T_U) - 8*gamma**2*lamW**2*(11*KU*L**2 + 9*pi**2*T_U))**2 / (60*pi**2*KU*a**5*c**5*lamW**2*mu_W*(KU + c_etaU*gamma)**2*(KU*L**2 + pi**2*T_U)**2)` | `(cs^5*gConst*(-cEtaU^2 + kEtaEff*kU)*(ell^2*kU*(-88*gamma^2*lamW^2 + 9*kU*kWEff*Pi^2) + 9*Pi^2*(-8*gamma^2*lamW^2 + kU*kWEff*Pi^2)*tU)^2)/(60*a^5*c^5*kU*(cEtaU*gamma + kU)^2*lamW^2*muW*Pi^2*(ell^2*kU + Pi^2*tU)^2)` |

The two R_target expressions are sign-consistent: `(KU*Keta_eff - c_etaU^2) == (-cEtaU^2 + kEtaEff*kU)`, and the inner quadratic `(9*pi^2*KU*KW_eff*(KU*L^2 + pi^2*T_U) - 8*gamma^2*lamW^2*(11*KU*L^2 + 9*pi^2*T_U))` matches `(ell^2*kU*(-88*gamma^2*lamW^2 + 9*kU*kWEff*Pi^2) + 9*Pi^2*(-8*gamma^2*lamW^2 + kU*kWEff*Pi^2)*tU)` after distributing. All residual zero checks pass identically in both engines. `engines_agree: true` — but see F2: agreement is uninformative because the `.wl` script restates the `.py` script's algebra rather than re-deriving.

## Verdict justification

The two non-tautological substantive checks (M_tr = M_mix*S, support-loaded product law, dS/dzeta closed form, and the spoiler control in sympy) do hold up under attack: the spoiler probe at sympy line 152 returns -44.866, confirming the zeta-blindness of the loaded R_target reconstruction is not a vacuous algebraic identity. The dS/dzeta closed form (line 116) is a non-trivial differentiation check that genuinely exercises the support-enhancement factor's analytic form. However, two findings remain: the rho_0/sigma_0 collapse "identities" are constructed by string cancellation rather than derived (F1: `tautological_check`), and the Mathematica script is a line-by-line port of the SymPy script rather than an independent re-derivation (F2: `mathematica_transliteration`). Verdict: `findings`. Stop-cold: none — fixing either finding is mechanical and does not propagate to downstream units' results (it strengthens but does not change M_tr, R_target, S, or the product law).

## Self-test notes

Variable-independence trap: checked `sp.diff(Rtarget_loaded, zeta)` (line 111). Rtarget_loaded = product_expected/Mtr = 8*Lambda*(1-eps)/(pi^2*Mmix) after cancellation, so it has no zeta dependence by construction once M_tr = M_mix*S holds. The derivative test is therefore equivalent to the M_tr factorization (A5), but it is not a variable-independence trap in the failure-mode sense — Rtarget_loaded does depend on lamW, KU, Keta_eff, KW_eff, etc., and the spoiler probe at line 152-154 confirms zeta-dependence reappears under perturbation, so the check is meaningful as a regression detector. Parity: no integrals on symmetric domains in this unit. Trivial-case substitution: substituted lamphi -> lamW and Kphi_eff -> KW_eff into eps_phi - zeta_def*eps_W; both factors collapse to eps_W with zeta_def=1, residual 0 as expected, so A3/A15 is non-vacuous. Path specifications: directive targets are existing files, no missing-script subtype, so the path-specification trap does not apply.
