---
unit_id: 047
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage047_coherent_kernel_map.md
  paper_appendix: present
---

# Audit unit 047 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_047.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage047_coherent_kernel_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row line 72; cross-ref line 212)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.txt`

## What the paper claims

Stage 047 introduces the coherent support ratio `zeta` and proves that the support lane only enhances the available loading without moving the selected target ratio. The card's `\stagefield{Output}` is verbatim: "Support enhancement factor \eqref{eq:app-stage047-S} and target independence \eqref{eq:app-stage047-Rtarget}." The two boxed paper equations are `S(zeta;epsilon) = 1 + zeta(1-epsilon)/(1 - epsilon*zeta)` and `R_target = Lambda*(1-eps_eta)*(1-epsilon)^2/(Z_W*(1+chi_0)^2)`, with the factorization `M_tr = M_mix * S(zeta;epsilon)` on Eq. \ref{eq:app-stage047-Mtr}. The notes file supplies the upstream microscopic definitions (sigma = 88/(9*pi^2), eps_eta, eps_W, Z_W, delta0, deltaU, chi_0, zeta, Lambda) and lists additional identities that must hold on the coherent branch: `rho_0 = sigma_0 = chi_0`, `eps_phi = zeta*eps_W`, `Z_phi = zeta*Z_W`, the split mixed blocking `eps = eps_W^(split) = eps_W*(1 - (2/11)*deltaU/(1+deltaU))`, the closed-form `M_mix` and `M_supp`, the product law `R_target*M_tr = (8 Lambda (1-eps)/pi^2)*S`, and `dS/dzeta = (1-eps)/(1-zeta*eps)^2 > 0` (with `S(0;eps) = 1`). The Part III appendix row (line 72) summarizes the deliverable as "Support enhancement factor S(zeta;epsilon) and R_target independent of zeta", consistent with the card.

## What the script claims to verify

Both scripts (sympy `.py` and mathematica `.wl`) — whose internal banners still read "STAGE 30" / "STAGE 030" while the filenames and paper references use "stage 047" — claim to verify: (1) the coherent collapse identities `rho_0 = sigma_0 = chi_0`; (2) the support/mixed proportionalities `eps_phi = zeta_def*eps_W` and `Z_phi = zeta_def*Z_W` with `zeta_def := lamphi^2*KW_eff/(lamW^2*Kphi_eff)`; (3) the exact total baseline factorization `M_tr = M_mix*S(zeta;eps)`; (4) the support-loaded product law `R_target*M_tr = (8 Lambda (1-eps)/pi^2)*S` and the consequent zeta-blindness of the loaded R_target reconstruction; (5) monotonicity `dS/dzeta = (1-eps)/(1-zeta*eps)^2` and the boundary value `S(zeta=0) = 1`. The sympy script also runs a numeric "spoiler" probe at lines 144-177 that perturbs `M_supp` off the exact support-loading coefficient and confirms the loaded-R_target reconstruction's zeta-derivative becomes non-zero. The mathematica file additionally cross-checks `sEnhance = mTr/mMix` against the closed-form `1 + zeta*(1-eps)/(1-zeta*eps)`.

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| `S(zeta;eps) = 1 + zeta(1-eps)/(1-eps*zeta)` (Eq. app-stage047-S) | sympy:139 `dS/dzeta - (1-eps)/(1-zeta eps)^2`, sympy:140 `S(zeta=0)-1`, .wl:138/139 same | match |
| `R_target = Lambda*(1-eps_eta)*(1-eps)^2/(Z_W*(1+chi_0)^2)` (Eq. app-stage047-Rtarget) | sympy:114 defines `Rtarget` from this exact form; sympy:133 reconstructs from product law and confirms agreement; .wl:117 / .wl:132 mirror | match |
| `M_tr = M_mix*S(zeta;eps)` (Eq. app-stage047-Mtr) | sympy:122 `M_tr - M_mix*S`, .wl:124 same | match |
| `rho_0 = sigma_0 = chi_0` (notes §2) | sympy:73-74 `expect_zero("rho_0 - chi_0", ...)`, sympy:74 sigma; .wl:68-69 mirror | partial (assertions exist but cannot fail — see F1) |
| `eps_phi = zeta*eps_W`, `Z_phi = zeta*Z_W` (notes §2) | sympy:96-97; .wl:86-87 | match |
| `M_mix`, `M_supp` closed forms (notes §4) | sympy:110-111; .wl:97/105 (now built via support-load factor rather than copied) | match |
| product law `R_target*M_tr = (8 Lambda (1-eps)/pi^2)*S` (notes §5) | sympy:132 `product law`; .wl:131 same | match |
| `R_target` independence of zeta (notes §5, card Output) | sympy:134 `dR_target_loaded/dzeta = 0` and sympy:135-138 `R_target_loaded(zeta) - R_target_loaded(0)`; .wl:133/134-137 mirror; sympy:144-177 spoiler control | match (with spoiler confirming non-vacuous in sympy only) |
| `dS/dzeta = (1-eps)/(1-zeta*eps)^2 > 0`, `S(0;eps) = 1` (notes §4) | sympy:139, sympy:140; .wl:138, .wl:139 | match |
| `R_tr = (1 + chi_0/(1+deltaU))/(1+chi_0)`, `eps = eps_W*(1 - (2/11)*deltaU/(1+deltaU))`, `delta = ...` (notes §3) | sympy:100-107 print only (no expect_zero); .wl:89-95 print only | extra (printed for trace, no scripted check — paper card does not list these as Outputs for this stage; they are upstream restatements) |
| (none) | sympy:144-177 numeric spoiler probe of zeta-blindness | extra (script-side overhead; serves as a negative control with no direct paper analogue, which is desirable) |

`paper_alignment: aligned` — every paper-side deliverable has a corresponding script-side check, and the script does not assert anything that contradicts the paper. The lone "partial" row (collapse identities) is a script-internal weakness (F1), not a paper-mismatch.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 62-64 | `assert sp.simplify(rho0_num*KU - rho0_den*gamma*c_etaU) == 0` (channel-saturation guard) | notes §2 rho_0 = chi_0 | no (tautological — see F1) |
| A2 | sympy | 65-67 | `assert sp.simplify(sigma0_num*KU - sigma0_den*gamma*c_etaU) == 0` | notes §2 sigma_0 = chi_0 | no (tautological — see F1) |
| A3 | sympy | 73 | `expect_zero("rho_0 - chi_0", rho0 - chi0)` | notes §2 rho_0 = chi_0 | partial (numerator/denominator now formally separated but factors still cancel by commutativity — see F1) |
| A4 | sympy | 74 | `expect_zero("sigma_0 - chi_0", sigma0 - chi0)` | notes §2 sigma_0 = chi_0 | partial (see F1) |
| A5 | sympy | 96 | `expect_zero("eps_phi - zeta_def*eps_W", ...)` | notes §2 eps_phi = zeta*eps_W | yes |
| A6 | sympy | 97 | `expect_zero("Z_phi - zeta_def*Z_W", ...)` | notes §2 Z_phi = zeta*Z_W | yes |
| A7 | sympy | 122 | `expect_zero("M_tr - M_mix*S", ...)` | Eq. app-stage047-Mtr | yes (genuine algebraic factorization of the closed forms) |
| A8 | sympy | 132 | `expect_zero("product law", product_actual - product_expected)` | notes §5 product law | yes |
| A9 | sympy | 133 | `expect_zero("support-loaded R_target reconstruction", Rtarget_loaded - Rtarget)` | card Output / notes §5 | yes |
| A10 | sympy | 134 | `expect_zero("dR_target_loaded/dzeta", ...)` | card Output (zeta-independence) | partial (follows from A7+A8 algebraically; spoiler probe A13 supplies the non-vacuousness control) |
| A11 | sympy | 135-138 | `expect_zero("R_target_loaded(zeta) - R_target_loaded(0)", ...)` | card Output (zeta-independence) | partial (redundant with A10) |
| A12 | sympy | 139 | `expect_zero("dS/dzeta - (1-eps)/(1-zeta eps)^2", ...)` | notes §4 monotonicity closed form | yes |
| A13 | sympy | 140 | `expect_zero("S(zeta=0)-1", ...)` | notes §4 boundary value | yes |
| A14 | sympy | 144-177 | numeric spoiler: perturb `M_supp -> M_supp + eta_bad*zeta*M_mix`, expect non-zero `dR/dzeta` | card Output (non-vacuous zeta-blindness) | yes (genuine negative control) |
| A15 | .wl | 52-57 | `If[TrueQ[FullSimplify[rho0Num*kU - rho0Den*gamma*cEtaU] === 0], Null, ...]` | notes §2 rho_0 = chi_0 | no (tautological — see F1) |
| A16 | .wl | 58-63 | `If[TrueQ[FullSimplify[sigma0Num*kU - sigma0Den*gamma*cEtaU] === 0], Null, ...]` | notes §2 sigma_0 = chi_0 | no (tautological — see F1) |
| A17 | .wl | 68 | `expectZero["rho_0 - chi_0", rho0 - chi0]` | notes §2 | partial (see F1) |
| A18 | .wl | 69 | `expectZero["sigma_0 - chi_0", sigma0 - chi0]` | notes §2 | partial (see F1) |
| A19 | .wl | 86 | `expectZero["eps_phi - zeta_def eps_W", ...]` | notes §2 | yes |
| A20 | .wl | 87 | `expectZero["Z_phi - zeta_def Z_W", ...]` | notes §2 | yes |
| A21 | .wl | 115 | `expectZero["S from ratio agrees with closed-form S", sEnhance - sClosedForm]` | derivation of S from M_tr/M_mix | no (tautological by construction — see F2) |
| A22 | .wl | 124 | `expectZero["M_tr - M_mix S", mTr - mMix*sEnhance]` | Eq. app-stage047-Mtr | partial (by construction `sEnhance := mTr/mMix`, so this is `mTr - mTr == 0` — see F2) |
| A23 | .wl | 131 | `expectZero["product law", productActual - productExpected]` | notes §5 | yes |
| A24 | .wl | 132 | `expectZero["support-loaded R_target reconstruction", loadedRTarget - rTarget]` | card Output / notes §5 | yes |
| A25 | .wl | 133 | `expectZero["dR_target_loaded/dzeta", ...]` | card Output | partial |
| A26 | .wl | 134-137 | `expectZero["R_target_loaded(zeta) - R_target_loaded(0)", ...]` | card Output | partial |
| A27 | .wl | 138 | `expectZero["dS/dzeta - (1-eps)/(1-zeta eps)^2", ...]` | notes §4 | yes |
| A28 | .wl | 139 | `expectZero["S(zeta=0)-1", ...]` | notes §4 | yes |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py:48-74`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl:40-69`

**What's wrong:**
A prior round of red-team fixes (`## Applied: F1` in the directive) reshuffled the rho_0/sigma_0 construction so that numerator and denominator are written as separate Python/Mathematica symbols and the cancellation is nominally performed by `simplify`/`FullSimplify`. The intermediate "channel-saturation guard" assertions were added at sympy:62-67 and .wl:52-63. Inspecting the algebra:

```
sympy:48-50
rho0_num = gamma * lamW * c_etaU
rho0_den = KU * lamW
rho0 = sp.simplify(rho0_num / rho0_den)
...
sympy:62-64
assert sp.simplify(rho0_num * KU - rho0_den * gamma * c_etaU) == 0
```

Substituting:
- `rho0_num * KU = gamma * lamW * c_etaU * KU`
- `rho0_den * gamma * c_etaU = KU * lamW * gamma * c_etaU`

These two expressions are the same product `KU * gamma * lamW * c_etaU` rearranged by commutativity, so their difference is identically zero before `simplify` ever runs. The same is true of the sigma_0 guard (sympy:65-67). The .wl file mirrors this exactly at lines 52-63 with `FullSimplify[rho0Num*kU - rho0Den*gamma*cEtaU] === 0` — again identically zero by commutativity. The downstream assertions at sympy:73-74 and .wl:68-69 (`expect_zero("rho_0 - chi_0", rho0 - chi0)` etc.) therefore inherit the same tautology: `rho0/rho0_den = gamma*lamW*c_etaU/(KU*lamW) = gamma*c_etaU/KU = chi0` because the lamW factors cancel by string structure, not by physics. The captured outputs confirm this — both print `rho_0 = c_etaU*gamma/KU` directly, with no intermediate un-simplified form.

The "physical content" the guard is meant to encode is the channel-saturation rule: that the W-channel polarisation amplitude contributes the same factor in numerator and denominator. But the script does not introduce this as a non-trivial premise; it writes the rule directly into the expression, so the guard's `simplify(...) == 0` is `simplify(KU*gamma*lamW*c_etaU - KU*gamma*lamW*c_etaU) == 0`, which holds for arbitrary symbolic content.

**Why this matters:**
The paper / notes claim `rho_0 = sigma_0 = chi_0` is an exact coherent-kernel identity on the tracking branch. The script's "verification" of this identity is structurally `(x - x) == 0`. If a future refactor accidentally introduces a real asymmetry (e.g., the W-channel polarisation acquires a wall-dressing factor that the phi-channel does not), the existing guards would still pass. The collapse identities are documented as load-bearing in notes §2 ("So the support lane is not independent of the mixed lane at the level of dimensionless couplings"), so trusting an empty assertion is non-trivially dangerous for stages 048-050 which build on this.

**Required change:**
Either:
(a) Define `rho_0` and `sigma_0` from *physically distinct* construction paths — e.g., introduce explicit W-channel and phi-channel polarisation kernels `Pi_W(K_W_eff, lamW)` and `Pi_phi(K_phi_eff, lamphi)` such that the polarisation factor explicitly carries the channel coupling, then state the channel-saturation rule as a separate hypothesis `Pi_W/(lamW * KU) -> gamma*c_etaU/KU under the coherent-branch matching condition` and `simplify` under that substitution. The collapse identity then becomes a non-trivial check on the matching condition. The auditor cannot prescribe the exact algebra without overstepping scope; the unit author should look at the upstream Stage-28 derivation that produced the matching condition.
(b) If no independent definition is available within Stage 047's scope, the honest move is to **delete** the `rho_0 - chi_0` and `sigma_0 - chi_0` assertions (sympy:73-74, .wl:68-69) plus the no-op guards at sympy:62-67 and .wl:52-63, and add a comment stating "rho_0 and sigma_0 are notational aliases for chi_0 on the coherent branch; the equality is a definitional rename, verified upstream at Stage 28, not an independent identity here."

Option (b) is the smaller fix and is acceptable to the auditor; option (a) is more informative but requires re-deriving from Stage-28 material.

**Verification:**
After the fix, the line `expect_zero("rho_0 - chi_0", rho0 - chi0)` either (a) operates on a `rho0` whose construction *does not* contain `lamW/lamW` as a syntactic factor pair, so the residual prior to `simplify` is non-zero; or (b) is absent. The Mathematica analogue should match. Engines should still exit 0.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl:97-124`

**What's wrong:**
A prior round of red-team fixes (`## Applied: F2`) replaced the .wl script's direct closed-form definition of `mSupp` with `mSupp = FullSimplify[mMix * supportLoadFactor]` where `supportLoadFactor = zeta*(1 - eps)/(1 - zeta*eps)` (.wl:104). The script then defines `mTr = mMix + mSupp` (.wl:108) and `sEnhance = mTr/mMix` (.wl:109) and adds the cross-check `expectZero["S from ratio agrees with closed-form S", sEnhance - sClosedForm]` (.wl:115) plus the kept `expectZero["M_tr - M_mix S", mTr - mMix*sEnhance]` (.wl:124).

Trace the algebra:
- `sEnhance = mTr/mMix = (mMix + mMix*supportLoadFactor)/mMix = 1 + supportLoadFactor = 1 + zeta*(1-eps)/(1-zeta*eps)`.
- `sClosedForm = 1 + zeta*(1 - eps)/(1 - zeta*eps)` (.wl:114).

So `sEnhance - sClosedForm` is identically zero by construction; the new "agrees with closed-form S" assertion at .wl:115 cannot fail. Similarly `mTr - mMix*sEnhance = mTr - mMix*(mTr/mMix) = 0` at the symbolic level, regardless of what `mMix` and `mTr` actually are — so the existing `M_tr - M_mix S` assertion at .wl:124 has also become tautological in this engine after the F2 fix (in the sympy script, A7 at line 122 remains substantive because `S` there is defined independently of `Mtr/Mmix`, see sympy:112).

The fix as applied does change the .wl file's surface structure away from a byte-equivalent transliteration of the .py, but it does so by making the .wl's load-bearing assertions algebraically trivial. The .wl no longer verifies the substantive identity (M_mix + M_supp)/M_mix = 1 + zeta(1-eps)/(1-zeta*eps); it asserts (mMix*(1+x))/mMix == 1+x and mMix*(1+x) - mMix*(mTr/mMix) == 0.

**Why this matters:**
The two-engine policy exists to catch algebraic errors. Stage 047's load-bearing claim is that the explicit closed forms of M_mix (8*Z_W*(1+chi_0)^2/(pi^2*(1-eps_eta)*(1-eps))) and M_supp (8*zeta*Z_W*(1+chi_0)^2/(pi^2*(1-eps_eta)*(1-zeta*eps))) sum to M_mix * S where S is the closed-form support-enhancement factor. The sympy script exercises this directly (Msupp defined from its own closed form at sympy:111, then `Mtr = Mmix + Msupp` and `Mtr - Mmix*S == 0` is non-trivial). The .wl script no longer exercises this — it derives M_supp *from* M_mix and the support-load factor, so the agreement is enforced by the script's own construction.

**Required change:**
In `mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl`, replace lines 99-108 to define `mSupp` from the same independent closed form the paper / notes state, then verify the S factorization as a non-trivial identity:

Before (current lines 99-115):
```mathematica
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
(* M_supp from the paper's closed form (notes §4): the support lane
   replaces (1-eps) with (1-zeta*eps) in the denominator and acquires a
   prefactor zeta. We write this from the dimensionless ratios directly
   so the factorization M_tr = M_mix * S below is a non-trivial algebraic
   identity, not a script-built tautology. *)
mSupp = FullSimplify[8*zeta*zW*(1 + chi0)^2/(Pi^2*(1 - epsEta)*(1 - zeta*eps)),
                     Assumptions -> $Assumptions];

(* M_tr as the sum of the two independently defined baselines. *)
mTr = FullSimplify[mMix + mSupp, Assumptions -> $Assumptions];

(* S(zeta;eps) from the closed-form definition (Eq. app-stage047-S). *)
sEnhance = FullSimplify[1 + zeta*(1 - eps)/(1 - zeta*eps),
                        Assumptions -> $Assumptions];
```

Then the existing `expectZero["M_tr - M_mix S", mTr - mMix*sEnhance]` at the line currently numbered 124 is restored to a substantive check: it tests that the sum `mMix + mSupp` equals `mMix * (1 + zeta*(1-eps)/(1-zeta*eps))` as an algebraic identity between independently defined closed forms.

Removing the `sClosedForm` line and the tautological `"S from ratio agrees with closed-form S"` assertion is required as part of this change; that assertion cannot survive after `sEnhance` is defined from the closed form rather than from `mTr/mMix`.

The "transliteration" concern from the prior F2 audit is mitigated separately by the fact that the sympy script's spoiler probe (A14 at sympy:144-177) supplies a non-vacuous numeric attack on the same claim, and by the engine-cross-check on the simplified outputs (see below). A surface-level rearrangement of the algebra to look "different" from the sympy version is *not* a substitute for an actually-failable identity.

**Verification:**
After the fix, the captured `.wl` output should show:
- `M_supp = ...` printed in a form algebraically equal to (but not by-construction identical to) the sympy `M_supp` output;
- `S(zeta;eps) = 1 + zeta*(1-eps)/(1-zeta*eps)` printed (or its FullSimplify normal form);
- `PASS: M_tr - M_mix S` line still present;
- the `PASS: S from ratio agrees with closed-form S` line is **gone** (it is no longer a separate check).
The script exits 0.

## Independent-derivation check (Mathematica)

The .wl script (post the prior F2 patch) is *structurally* different from the .py script — it derives mSupp via the support-load factor and computes sEnhance as mTr/mMix — but as detailed in F2, this surface-level independence comes at the cost of making the load-bearing assertions algebraically tautological in the .wl engine. The independent-derivation requirement is therefore *not* met in a meaningful sense: the .wl no longer encodes a parallel statement of the closed-form M_supp that the .py's `Msupp - 8*zeta*ZW*(1+chi0)^2/(pi^2*(1-eps_eta)*(1-zeta*eps)) == 0` would catch under perturbation. F2 prescribes restoring the closed-form M_supp in the .wl while keeping the assertion that the sum equals M_mix * S — that recovers genuine cross-engine independence on the load-bearing identity. The rho_0/sigma_0 collapse identities are also tautological in both engines (F1).

## Engine cross-check

Both engines run to exit 0 and produce algebraically equivalent simplified forms for every printed intermediate, after accounting for variable-name spelling (`KU` ↔ `kU`, `Keta_eff` ↔ `kEtaEff`, `lamphi` ↔ `lamPhi`, `T_w` ↔ `tw`, etc.) and sign-cleanup (`KU*Keta_eff - c_etaU^2` ↔ `-(cEtaU^2 - kEtaEff*kU)`). Compare the printed `R_target`:

- sympy: `G*c_s**5*(KU*Keta_eff - c_etaU**2)*(9*pi**2*KU*KW_eff*(KU*L**2 + pi**2*T_U) - 8*gamma**2*lamW**2*(11*KU*L**2 + 9*pi**2*T_U))**2/(60*pi**2*KU*a**5*c**5*lamW**2*mu_W*(KU + c_etaU*gamma)**2*(KU*L**2 + pi**2*T_U)**2)`
- .wl: `(cs^5*gConst*(-cEtaU^2 + kEtaEff*kU)*(ell^2*kU*(-88*gamma^2*lamW^2 + 9*kU*kWEff*Pi^2) + 9*Pi^2*(-8*gamma^2*lamW^2 + kU*kWEff*Pi^2)*tU)^2)/(60*a^5*c^5*kU*(cEtaU*gamma + kU)^2*lamW^2*muW*Pi^2*(ell^2*kU + Pi^2*tU)^2)`

After distributing the inner polynomial `ell^2*kU*(-88*gamma^2*lamW^2 + 9*kU*kWEff*Pi^2) + 9*Pi^2*(-8*gamma^2*lamW^2 + kU*kWEff*Pi^2)*tU = 9*pi^2*kU*kWEff*(ell^2*kU + Pi^2*tU) - 8*gamma^2*lamW^2*(11*ell^2*kU + 9*Pi^2*tU)`, the two expressions match. All residual zero-checks return 0 in both engines. `engines_agree: true`. Output transcripts are fresh (`.wl` script mtime 12:55, .wl output 12:57; .py 12:55, .py output 12:57).

## Verdict justification

The paper's three load-bearing claims (Eqs. app-stage047-S, app-stage047-Rtarget, app-stage047-Mtr) are exercised non-tautologically by the sympy script: A7 verifies the M_tr factorization between independent closed forms of M_mix, M_supp, and S; A8 verifies the product law; A12 verifies the closed-form dS/dzeta; and A14 (numeric spoiler) supplies a non-vacuous attack on the zeta-blindness of R_target. The notes' coherent collapse identities `rho_0 = sigma_0 = chi_0` are *not* substantively verified in either engine (F1); the assertions are constructed so the lamW/lamphi factors cancel by string structure rather than by physical content. The .wl script's post-F2-patch construction of `mSupp = mMix * supportLoadFactor` makes the .wl's `M_tr - M_mix*S` assertion tautological in that engine (F2), undermining the two-engine policy on the unit's central factorization claim. Both findings are mechanical to repair within the unit's scope and do not propagate downstream — Stage 048's inversion `S(zeta_req) = required` consumes the closed-form S only, not the .wl's derivation path. Verdict: `findings`. Stop-cold: none. Attacks I tried that did not break the script: (i) substituted lamphi -> lamW and Kphi_eff -> KW_eff in the eps_phi - zeta_def*eps_W identity (A5) — collapses correctly to 0 with zeta_def -> 1; (ii) checked that the spoiler at sympy:144-177 has a non-zero `eta_bad` coefficient sweeping the perturbation off the exact support-loading map and confirmed the printed residual -44.866 is sign-consistent with an over-loaded support packet; (iii) re-derived `R_target * M_tr` by hand from M_mix's closed form and verified `R_target * M_mix = 8 Lambda (1-eps)/pi^2`, then multiplied by `S` to recover the product law (A8) — the algebra holds. The paper card, notes, and appendix row are mutually consistent and the script's claim matches the paper's claim modulo the F1/F2 internal weaknesses.

## Self-test notes

Variable-independence trap: confirmed `Rtarget_loaded` (sympy:127) has no `zeta` symbol after cancellation, so `sp.diff(Rtarget_loaded, zeta) == 0` is by-construction; this is anchored to the paper's zeta-independence claim via the structural argument (zeta cancels iff the product law holds), and the spoiler probe at sympy:144-177 supplies the non-vacuousness. No symmetric-domain integrals in this unit. Trivial-case substitution: substituted lamphi -> lamW into eps_phi - zeta_def*eps_W → both reduce with zeta_def=1, residual 0; A5/A19 non-vacuous. For F1 / F2 directive instructions, confirmed: (F1) deleting the rho_0/sigma_0 assertions does not break any downstream check in either script — A3-A14 do not consume rho_0 or sigma_0; (F2) restoring the closed-form mSupp definition in the .wl keeps `mTr = mMix + mSupp` and the subsequent product-law and monotonicity checks intact, and the directive change instructions reference the correct file paths (`mathematica/` for the `.wl`). Paper round-trip: the F2 fix restores the .wl's M_supp to the *same closed form the notes §4 state*, so no new paper_misalignment is introduced.
