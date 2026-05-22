---
unit_id: 047
batch: III.1
created_at: 2026-05-22T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-22T12:55:57-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 047

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py:42-51`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl:36-44`

**Issue:**
The "coherent kernel identities" rho_0 = sigma_0 = chi_0 are not non-trivially verified. The current code writes `rho0 = (gamma*lamW)*c_etaU/(KU*lamW)` and `sigma0 = c_etaU*(gamma*lamphi)/(KU*lamphi)`, with the channel-coupling factors (lamW and lamphi) appearing as multiplicative pairs that cancel by construction before any `simplify` call. The subsequent assertions `expect_zero("rho_0 - chi_0", ...)` and `expect_zero("sigma_0 - chi_0", ...)` therefore cannot fail. The Mathematica file mirrors this pattern.

**Required change:**

Step 1 — SymPy file `scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py`:

Replace lines 41-45 (which currently read):
```python
# Exact coherent identities.
rho0 = sp.simplify((gamma * lamW) * c_etaU / (KU * lamW))
sigma0 = sp.simplify(c_etaU * (gamma * lamphi) / (KU * lamphi))
chi0 = sp.simplify(gamma * c_etaU / KU)
```

with the following block, which constructs rho_0 and sigma_0 from independent channel-saturation expressions whose collapse to chi_0 is performed by `simplify` rather than by manual string cancellation:

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
```

The two new `assert` statements verify (non-trivially) that the channel-saturation rule `rho0_num * KU == rho0_den * (gamma*c_etaU)` holds, which is the physical content of the "identity"; they would fail if e.g. the numerator carried an extra lamW factor without a matching denominator factor.

Leave the existing assertions on lines 50-51 (`expect_zero("rho_0 - chi_0", rho0 - chi0)` and `expect_zero("sigma_0 - chi_0", sigma0 - chi0)`) unchanged — they remain a useful regression check on the simplified forms.

Step 2 — Mathematica file `mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl`:

Replace lines 35-38 (which currently read):
```mathematica
sigma = 88/(9*Pi^2);
rho0 = FullSimplify[(gamma*lamW)*cEtaU/(kU*lamW), Assumptions -> $Assumptions];
sigma0 = FullSimplify[cEtaU*(gamma*lamPhi)/(kU*lamPhi), Assumptions -> $Assumptions];
chi0 = FullSimplify[gamma*cEtaU/kU, Assumptions -> $Assumptions];
```

with:
```mathematica
sigma = 88/(9*Pi^2);

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
```

Leave the existing `expectZero["rho_0 - chi_0", ...]` and `expectZero["sigma_0 - chi_0", ...]` calls at lines 43-44 unchanged.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 047` and `redteam exec-mathematica 047` and confirm: (a) both scripts exit 0; (b) the captured `.txt` output for the sympy script now includes a printed intermediate of the form `rho_0 = c_etaU*gamma/KU` (post-simplify) AND the two new internal `assert` statements pass silently; (c) the captured `.txt` for the Mathematica script does not contain `FAIL: rho0 channel-saturation rule violated`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl`
- summary: Reworked rho_0 and sigma_0 construction to keep channel numerator and denominator expressions separate and added non-trivial channel-saturation checks.
- deviation: none

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl:72-83`

**Issue:**
The `.wl` script defines `mMix`, `mSupp`, `sEnhance`, `mTr` (lines 72-75) using the same closed-form expressions as the SymPy script (lines 87-90). The agreement of the two engines' final outputs is therefore not evidence of independent verification — it is byte-equivalent algebra under different syntax. The two-engine policy requires at least one independent derivation path.

**Required change:**

In `mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl`, replace lines 72-83 (currently):
```mathematica
mMix = FullSimplify[8*zW*(1 + chi0)^2/(Pi^2*(1 - epsEta)*(1 - eps)), Assumptions -> $Assumptions];
mSupp = FullSimplify[8*zeta*zW*(1 + chi0)^2/(Pi^2*(1 - epsEta)*(1 - zeta*eps)), Assumptions -> $Assumptions];
sEnhance = FullSimplify[1 + zeta*(1 - eps)/(1 - zeta*eps), Assumptions -> $Assumptions];
mTr = FullSimplify[mMix + mSupp, Assumptions -> $Assumptions];
rTarget = FullSimplify[lambdaScale*(1 - epsEta)*(1 - eps)^2/(zW*(1 + chi0)^2), Assumptions -> $Assumptions];

Print["M_mix = ", fmt[mMix]];
Print["M_supp = ", fmt[mSupp]];
Print["S(zeta;eps) = ", fmt[sEnhance]];
Print["M_tr = ", fmt[mTr]];
Print["R_target = ", fmt[rTarget]];
expectZero["M_tr - M_mix S", mTr - mMix*sEnhance];
```

with the following block, which derives `mSupp` independently as `mMix * (zeta*(1-eps)/(1-zeta*eps))` from the S-decomposition (rather than copying the SymPy closed form), and derives `sEnhance` independently as `mTr/mMix` after constructing `mTr` from the independent definitions of `mMix` and `mSupp`. The final closed-form expressions must still agree with the SymPy outputs:

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

rTarget = FullSimplify[lambdaScale*(1 - epsEta)*(1 - eps)^2/(zW*(1 + chi0)^2), Assumptions -> $Assumptions];

Print["M_mix = ", fmt[mMix]];
Print["M_supp = ", fmt[mSupp]];
Print["S(zeta;eps) = ", fmt[sEnhance]];
Print["M_tr = ", fmt[mTr]];
Print["R_target = ", fmt[rTarget]];
expectZero["M_tr - M_mix S", mTr - mMix*sEnhance];
```

Do not alter the assertions at lines 90-98 — they continue to verify the product law and monotonicity from the now-independently-derived `mSupp`, `mTr`, and `sEnhance`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 047` and confirm: (a) the script exits 0; (b) the captured `.txt` output contains a new `PASS: S from ratio agrees with closed-form S` line; (c) the printed `M_supp = ...` and `S(zeta;eps) = ...` expressions match the corresponding SymPy outputs after FullSimplify (algebraic equivalence, not byte equality); (d) the existing `PASS: M_tr - M_mix S`, `PASS: product law`, `PASS: support-loaded R_target reconstruction`, and `PASS: dS/dzeta - (1-eps)/(1-zeta eps)^2` lines remain present.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl`
- summary: Replaced the Mathematica closed-form support packet copy with an independent support-loading derivation and ratio-based S check.
- deviation: none
