---
unit_id: 047
batch: III.1
created_at: 2026-05-26T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-26T01:46:54-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 047

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py:41-74`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl:37-69`

**Issue:**
The `rho_0 = sigma_0 = chi_0` "coherent collapse identities" cannot fail in either engine. In the sympy file the construction is

```python
rho0_num = gamma * lamW * c_etaU        # line 48
rho0_den = KU * lamW                    # line 49
rho0 = sp.simplify(rho0_num / rho0_den) # line 50
...
sigma0_num = c_etaU * gamma * lamphi    # line 53
sigma0_den = KU * lamphi                # line 54
sigma0 = sp.simplify(sigma0_num / sigma0_den)  # line 55
...
chi0 = sp.simplify(gamma * c_etaU / KU) # line 58
```

then the "channel-saturation guard" asserts `simplify(rho0_num*KU - rho0_den*gamma*c_etaU) == 0` and `simplify(sigma0_num*KU - sigma0_den*gamma*c_etaU) == 0`, both of which are `simplify(KU*gamma*lamW*c_etaU - KU*gamma*lamW*c_etaU) == 0` (resp. with `lamphi`) — true by commutativity for arbitrary symbolic content, independent of any physics. The `expect_zero("rho_0 - chi_0", rho0 - chi0)` and `expect_zero("sigma_0 - chi_0", sigma0 - chi0)` calls at lines 73-74 inherit the same tautology because `rho0` and `sigma0` reduce to `chi0` by the `lamW/lamW` (resp. `lamphi/lamphi`) cancellation, which is structural in the script's expression, not derived from any matching condition. The Mathematica file (lines 40-69) mirrors this exactly.

The simplest non-disruptive fix is option (b) from the report — delete the assertions and the no-op guards, and replace them with a comment marking the equality as a notational rename verified upstream at Stage 28. If the unit author later wants to re-establish a substantive collapse identity here, that requires re-derivation from Stage-28 matching conditions and should be scheduled as a separate change.

**Required change:**

Step 1 — SymPy file `scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py`:

Replace lines 41-74 (currently the `# Exact coherent identities.` block plus the two no-op guard `assert` statements plus the banner and the two `expect_zero` calls) with the following block. The block keeps `chi0` as the bare coherent ratio used downstream, removes `rho0`/`sigma0` as separate symbols, prints `chi_0` for trace continuity with the .wl, and inserts a comment explaining where the collapse identity is actually verified.

Before:
```python
# Exact coherent identities.
# rho_0 is the W-channel coherent ratio: cross-correlation strength
# (gamma * lamW * c_etaU) normalised by the W-channel coupling lamW
# composed with the eta-channel propagator KU. The lamW factor enters
# both numerator (through the W-channel polarisation) and denominator
# (through the normalisation), but here we keep them in separate
# subexpressions so the cancellation is a genuine simplification step.
rho0_num = gamma * lamW * c_etaU
rho0_den = KU * lamW
rho0 = sp.simplify(rho0_num / rho0_den)

# sigma_0 is the phi-channel analogue with lamphi in place of lamW.
sigma0_num = c_etaU * gamma * lamphi
sigma0_den = KU * lamphi
sigma0 = sp.simplify(sigma0_num / sigma0_den)

# chi_0 is the bare coherent ratio (no channel coupling on either side).
chi0 = sp.simplify(gamma * c_etaU / KU)

# Sanity: rho_0 and sigma_0 must not be identically chi_0 prior to
# simplification — otherwise the assertions below are tautological.
assert sp.simplify(rho0_num * KU - rho0_den * gamma * c_etaU) == 0, (
    "rho_0 numerator/denominator do not pre-satisfy the channel-saturation rule"
)
assert sp.simplify(sigma0_num * KU - sigma0_den * gamma * c_etaU) == 0, (
    "sigma_0 numerator/denominator do not pre-satisfy the channel-saturation rule"
)

banner("1. Coherent interference ratios")
print("rho_0 =", rho0)
print("sigma_0 =", sigma0)
print("chi_0 =", chi0)
expect_zero("rho_0 - chi_0", rho0 - chi0)
expect_zero("sigma_0 - chi_0", sigma0 - chi0)
```

After:
```python
# Coherent interference ratio.
# On the coherent tracking branch the W-channel and phi-channel
# polarisation amplitudes saturate to the same bare ratio gamma*c_etaU/KU.
# That saturation is established upstream at Stage 28 (matching condition
# for the coherent local D/N kernel); within the scope of Stage 047 the
# equalities rho_0 = sigma_0 = chi_0 are a notational rename rather than
# an independent identity, so we do not assert them here. Any attempt to
# verify them locally reduces to lamW/lamW (resp. lamphi/lamphi) string
# cancellation, which is tautological.
chi0 = sp.simplify(gamma * c_etaU / KU)

banner("1. Coherent interference ratio")
print("chi_0 =", chi0)
```

Step 2 — Mathematica file `mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl`:

Replace lines 37-69 (the `rho0Num/rho0Den/rho0`, `sigma0Num/sigma0Den/sigma0`, `chi0` definitions, the two `If[TrueQ[FullSimplify[...] === 0], Null, ...]` no-op guards, the trace prints, and the two `expectZero` calls) with:

Before:
```mathematica
(* rho_0: W-channel coherent ratio (numerator/denominator kept separate
   so the channel-saturation cancellation is performed by FullSimplify,
   not by string cancellation). *)
rho0Num = gamma*lamW*cEtaU;
rho0Den = kU*lamW;
rho0 = FullSimplify[rho0Num/rho0Den, Assumptions -> $Assumptions];

(* sigma_0: phi-channel analogue. *)
sigma0Num = cEtaU*gamma*lamPhi;
sigma0Den = kU*lamPhi;
sigma0 = FullSimplify[sigma0Num/sigma0Den, Assumptions -> $Assumptions];

chi0 = FullSimplify[gamma*cEtaU/kU, Assumptions -> $Assumptions];

(* Sanity: the channel-saturation rule must hold non-trivially. *)
If[
  TrueQ[FullSimplify[rho0Num*kU - rho0Den*gamma*cEtaU,
    Assumptions -> $Assumptions] === 0],
  Null,
  (Print["FAIL: rho0 channel-saturation rule violated"]; Exit[1])
];
If[
  TrueQ[FullSimplify[sigma0Num*kU - sigma0Den*gamma*cEtaU,
    Assumptions -> $Assumptions] === 0],
  Null,
  (Print["FAIL: sigma0 channel-saturation rule violated"]; Exit[1])
];

Print["rho_0 = ", fmt[rho0]];
Print["sigma_0 = ", fmt[sigma0]];
Print["chi_0 = ", fmt[chi0]];
expectZero["rho_0 - chi_0", rho0 - chi0];
expectZero["sigma_0 - chi_0", sigma0 - chi0];
```

After:
```mathematica
(* Coherent interference ratio.
   On the coherent tracking branch the W-channel and phi-channel
   polarisation amplitudes saturate to the same bare ratio gamma*cEtaU/kU.
   That saturation is established upstream at Stage 28 (matching condition
   for the coherent local D/N kernel); within the scope of Stage 047 the
   equalities rho_0 = sigma_0 = chi_0 are a notational rename, so we do
   not assert them here. Any local verification reduces to lamW/lamW
   (resp. lamPhi/lamPhi) string cancellation, which is tautological. *)
chi0 = FullSimplify[gamma*cEtaU/kU, Assumptions -> $Assumptions];

Print["chi_0 = ", fmt[chi0]];
```

Downstream `chi0` usage in both files is unchanged.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 047` and `redteam exec-mathematica 047` and confirm:
(a) both scripts exit 0;
(b) the captured `.txt` outputs no longer contain `rho_0 =`, `sigma_0 =`, `rho_0 - chi_0 =`, `sigma_0 - chi_0 =`, or `PASS: rho_0 - chi_0` lines — the `chi_0 = ...` print line remains;
(c) every remaining `PASS:` line (eps_phi - zeta_def eps_W, Z_phi - zeta_def Z_W, M_tr - M_mix S, product law, support-loaded R_target reconstruction, dR_target_loaded/dzeta, R_target_loaded(zeta) - R_target_loaded(0), dS/dzeta - (1-eps)/(1-zeta eps)^2, S(zeta=0)-1) is still present.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl`
- summary: Removed the tautological rho_0 and sigma_0 local checks and retained chi_0 as the Stage 047 notational coherent ratio.
- deviation: none

## F2 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl:97-115`

**Issue:**
The .wl script (post the prior F2 patch) derives `mSupp = mMix * supportLoadFactor` where `supportLoadFactor = zeta*(1-eps)/(1-zeta*eps)`, then defines `sEnhance = mTr/mMix` and asserts `sEnhance - sClosedForm == 0` where `sClosedForm = 1 + zeta*(1-eps)/(1-zeta*eps)`. By construction:
- `sEnhance = (mMix + mMix*supportLoadFactor)/mMix = 1 + supportLoadFactor = sClosedForm`,
- so the .wl line currently numbered 115 (`expectZero["S from ratio agrees with closed-form S", sEnhance - sClosedForm]`) is identically true;
- the downstream `expectZero["M_tr - M_mix S", mTr - mMix*sEnhance]` (currently line 124) is also identically true, since `mTr - mMix*(mTr/mMix) = 0` symbolically.

The net effect: the .wl no longer cross-checks the closed-form M_supp formula against the M_mix * S factorization. The two-engine policy on Stage 047's central factorization claim (`M_tr = M_mix * S`) is therefore unsatisfied in the Mathematica engine.

**Required change:**

In `mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl`, replace lines 97-115 (the current `mMix`, support-load-factor derivation of `mSupp`, `mTr`, `sEnhance`, `sClosedForm`, and the `"S from ratio agrees with closed-form S"` assertion) with the following block, which defines `mSupp` from the paper's closed form (notes §4: `M_supp = 8 zeta Z_W (1+chi_0)^2 / (pi^2 (1-eps_eta) (1-zeta eps))`) and `sEnhance` from the paper's closed-form S (Eq. app-stage047-S), so that the downstream `mTr - mMix*sEnhance == 0` assertion (now at the line directly following) is a substantive algebraic identity rather than a tautology.

Before (currently lines 97-115):
```mathematica
mMix = FullSimplify[8*zW*(1 + chi0)^2/(Pi^2*(1 - epsEta)*(1 - eps)), Assumptions -> $Assumptions];

(* Independent derivation of M_supp: the support packet contributes a
   factor zeta*(1-eps)/(1-zeta*eps) on top of M_mix, by the support-
   loading rule of stage 047. We construct M_supp from this rule rather
   than copying the closed form, so this engine has an independent path
   to the result. *)
supportLoadFactor = zeta*(1 - eps)/(1 - zeta*eps);
mSupp = FullSimplify[mMix*supportLoadFactor, Assumptions -> $Assumptions];

(* M_tr from the sum, then S as the ratio M_tr/M_mix. *)
mTr = FullSimplify[mMix + mSupp, Assumptions -> $Assumptions];
sEnhance = FullSimplify[mTr/mMix, Assumptions -> $Assumptions];

(* Cross-check: S as derived from the ratio must equal the closed-form
   1 + zeta*(1-eps)/(1-zeta*eps). If this fails, the support-loading
   rule above is inconsistent with the closed-form S in the paper. *)
sClosedForm = FullSimplify[1 + zeta*(1 - eps)/(1 - zeta*eps), Assumptions -> $Assumptions];
expectZero["S from ratio agrees with closed-form S", sEnhance - sClosedForm];
```

After:
```mathematica
mMix = FullSimplify[8*zW*(1 + chi0)^2/(Pi^2*(1 - epsEta)*(1 - eps)), Assumptions -> $Assumptions];

(* M_supp from the paper's closed form (notes §4): the support lane
   replaces (1-eps) with (1-zeta*eps) in the denominator and acquires a
   prefactor zeta. We write M_supp from the dimensionless ratios directly
   so the factorization M_tr = M_mix * S below is a non-trivial algebraic
   identity rather than a script-built tautology. *)
mSupp = FullSimplify[8*zeta*zW*(1 + chi0)^2/(Pi^2*(1 - epsEta)*(1 - zeta*eps)),
                     Assumptions -> $Assumptions];

(* M_tr as the sum of the two independently defined baselines. *)
mTr = FullSimplify[mMix + mSupp, Assumptions -> $Assumptions];

(* S(zeta;eps) from the closed-form definition (Eq. app-stage047-S). *)
sEnhance = FullSimplify[1 + zeta*(1 - eps)/(1 - zeta*eps),
                        Assumptions -> $Assumptions];
```

The subsequent `Print` lines and assertions are unchanged:
- `Print["M_mix = ", fmt[mMix]];`
- `Print["M_supp = ", fmt[mSupp]];`
- `Print["S(zeta;eps) = ", fmt[sEnhance]];`
- `Print["M_tr = ", fmt[mTr]];`
- `Print["R_target = ", fmt[rTarget]];`
- `expectZero["M_tr - M_mix S", mTr - mMix*sEnhance];`
- (everything from `rTarget = FullSimplify[lambdaScale*...]` line onward, including the product-law and monotonicity checks).

Do not change the `rTarget` definition or any check after `M_tr - M_mix S`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 047` and confirm:
(a) the script exits 0;
(b) the captured `.txt` output's `M_supp = ...` line is algebraically equal (after `FullSimplify`) to the sympy script's `M_supp = ...` line, but its symbolic construction path is now the closed form rather than a multiplication by `supportLoadFactor`;
(c) the captured `.txt` no longer contains `PASS: S from ratio agrees with closed-form S`;
(d) the captured `.txt` still contains `PASS: M_tr - M_mix S`, `PASS: product law`, `PASS: support-loaded R_target reconstruction`, `PASS: dR_target_loaded/dzeta`, `PASS: R_target_loaded(zeta) - R_target_loaded(0)`, `PASS: dS/dzeta - (1-eps)/(1-zeta eps)^2`, `PASS: S(zeta=0)-1`.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl`
- summary: Restored the Mathematica M_supp and S definitions to independent closed forms so the M_tr factorization check is substantive.
- deviation: none
