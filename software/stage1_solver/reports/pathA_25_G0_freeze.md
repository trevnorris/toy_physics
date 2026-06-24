SECOND_MEDIUM_DRIFT_AT_FREEZE

# pathA_25 G0 Freeze Artifact: GNLS Polar-Smectic Structure

artifact_status: append-only immutable freeze
directive: `software/stage1_solver/directives/pathA_25_gnls_polar_smectic_consistency_gates.md`
stage: `G0 structure freeze`
date: 2026-06-23
scope: reports-only plus dimensional check; no smectic, spectrum, light-mode, charge, cone-lock, leak, or throat gate solved

Line-1 verdict: `SECOND_MEDIUM_DRIFT_AT_FREEZE`.

Reason: the target-blind finite structure below is admissible, but the honest baseline independent-new-input count is `5`, already above the `>=2` NG5 pressure threshold. The finite sensitivity menu contains additional branch-only inputs. This is a freeze result, not a gate solve.

New G0 combined freeze hash: `f00ee99d465e2e311c68f47fcacf4af0154ca650642271ab66c36d112cb6a290`.

## G0.1 Kept GNLS + T0 Fidelity

The kept GNLS medium is unchanged:

- Bulk coordinates: `X^i=(x,y,z,w)`, `i=1..4`.
- Medium field: `psi=sqrt(rho) exp(i theta)`.
- Existing parent matter action: `research/pde/paper/pde.tex`, equation `parent-Lpsi`.
- Existing EOS: `U(rho)=K rho^5/4`, `P(rho)=K rho^5`, `c_s^2(rho)=5 K rho^4/m`.
- Existing flow: `v_i=(hbar/m) partial_i theta - (q_star/m) A_i`.
- The pathA_25 smectic driver below does not introduce a prescribed `V_conf`, layer profile, wall support, or brane-localized T0 term.

The exact T0 freeze-action bytes from `software/stage1_solver/reports/pathA_24_T0_freeze.md` are embedded unchanged inside the single combined `freeze-action` block in G0.5.

T0 hash recheck command:

```sh
awk '/^```freeze-action$/ {on=1; next} /^```$/ && on {exit} on {print}' software/stage1_solver/reports/pathA_24_T0_freeze.md | sha256sum
```

Recomputed T0 SHA-256:

```text
8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064  -
```

This matches the T0 short hash `8fa41ac51e88`.

## G0.2 Layering Driver Freeze

Target-blind admission criteria used here: O(4) spatial symmetry before spontaneous layering, no preselected layer normal, no prescribed layer profile, finite-dimensional parameterization, no free-form kernel, no added condensate field, locality preferred when otherwise equal, real-matter pedigree used only as non-payoff plausibility, and the strict `k->0` EOS preservation needed for a finite-`k` driver definition.

Forbidden-information statement, binding for every row and branch below:

smectic stability, light mode count, bounded-below, traction-vs-torque, angular-momentum closure, leak, magnetism / bulk-shear-free, charge signs, `c_γ=c_s`/`λγ`, throat radius/stability — none may admit, exclude, rank, tune, or normalize any term.

The older T0 forbidden-information rule also remains active: wall existence, wall lifetime, wall tension, wall-core texture, photon signature, `c_gamma`, bulk shear leakage, Magnus preservation, sector separation, two charge signs, charge universality, Coulomb form, `e`, `alpha`, throat radius, throat stability, dark-energy/cosmology flow behavior, and any T1-T5 outcome label were not used to admit, exclude, tune, rank, or normalize any term.

### Admissibility Audit

| Family | Frozen admissible finite forms | Free parameters and domains | Target-blind audit result | Branch role |
|---|---|---|---|---|
| R: non-local density kernel | `S_R = -1/2 int dt d^4X d^4R delta_rho(X) V_R(|R|) delta_rho(X+R)`, equivalently `tilde V_R(k)=A_R f_R(k/k_R)` with fixed `f_R(x)=(x^4-2x^2) exp(-x^2)`. | `A_R>0`, `[A_R]=M L^6 T^-2`; `k_R>0`, `[k_R]=L^-1`. No free shape function; no free exponent; range `R_R=1/k_R` derived. | Admitted. Finite-dimensional, O(4), finite-range Schwartz kernel by inverse transform, `tilde V_R(0)=0`, negative finite-k annulus. Less local than Family L and has a nonlocal pedigree/implementation cost. | Sensitivity branch `S_R_kernel`, replacing the L baseline. |
| L: Lifshitz/Brazovskii gradient | Baseline energy density `E_L = -1/2 c_L1 (partial_i rho)(partial_i rho) + 1/2 c_L2 (partial_i partial_i rho)^2`; action term `L_L=-E_L`. | `c_L1>0`, `[c_L1]=M L^8 T^-2`; `c_L2>0`, `[c_L2]=M L^10 T^-2`; derived `k_Lstar=sqrt(c_L1/(2 c_L2))`. | Admitted and ranked first. It is local, finite-dimensional, O(4), uses only `rho`, gives a finite-k quadratic feature, and has no support/profile kernel to hide a gate answer. | Pre-committed baseline branch `B0_Lifshitz`. |
| C: polar-density coupling | `E_Cdiv = -lambda_Cdiv delta_rho partial_i P^i`; `E_Cpin = 1/2 chi_Cpin (P^i partial_i rho)(P^j partial_j rho)`. They are admitted only as finite coupled branches with an R or L finite-k regularizer, not as a stand-alone driver. | `lambda_Cdiv` signed, `[lambda_Cdiv]=M L^3 T^-2`; `chi_Cpin` signed, `[chi_Cpin]=M L^8 T^-2`. | Partly admitted. The listed operators are covariant, finite, and use only `rho` and `P^i`. Rejected as stand-alone baseline because C alone supplies no target-blind high-k density regularizer; using it alone would require an uncounted kernel or later rescue. | Sensitivity branches `S_L_plus_Cdiv` and `S_L_plus_Cpin`. |

### Pre-Committed Branch Budget

| Branch | Driver content | Prior/ranking, target-blind | Run budget |
|---|---|---|---|
| `B0_Lifshitz` | Family L only: `c_L1,c_L2>0`; Family R absent; Family C couplings set to zero. | Ranked first by locality, finite-dimensionality, O(4) symmetry, and no kernel-shape function. | Baseline for Gate B4. |
| `S_R_kernel` | Family R kernel only, fixed `f_R`, parameters `A_R,k_R`; Family L absent; Family C zero. | Real-matter finite-range pedigree sensitivity; less local and less minimal. | May be run only as a named sensitivity branch after baseline routing. |
| `S_L_plus_Cdiv` | Family L plus `E_Cdiv`. | Tests the finite admitted divergence coupling without changing the density regularizer. | Named sensitivity only; no silent rescue. |
| `S_L_plus_Cpin` | Family L plus `E_Cpin`. | Tests the finite admitted gradient-alignment coupling and records the NG2 risk. | Named sensitivity only; no silent rescue. |

### Family-C NG2 Note

`E_Cpin` carries a direct NG2 risk. For `chi_Cpin<0`, the energy is lowered by aligning `P^i` with `partial_i rho`, which on a smectic layer is the layer normal; this risks out-of-plane pinning of `P^i` and depletion of in-plane polar order. For `chi_Cpin>0`, the same operator penalizes normal alignment. This risk is recorded only; it did not rank the baseline or exclude the branch. `E_Cdiv` does not directly square `P^i partial_i rho`, so it has no direct out-of-plane-pinning sign channel at the freeze level.

## G0.3 Light-Sector Package Freeze

The light package is `POSTULATED` as a layer-effective constitutive sector tied to the smectic layers that Gate B4 may return. The support is not hand-drawn: `Sigma_n[rho]` denotes the density layers selected by the frozen smectic driver, and `delta_Sigma[rho]` is the corresponding distributional layer measure when a brane term is represented as a 4D bulk density. No Gate B4 result is used here.

Frozen MacCullagh stiffness:

```text
u^a(sigma,t), a=1..3, [u]=L
Omega_u^a := (1/2) epsilon^{abc} partial_b u_c
L_Mac = 1/2 varrho_br[rho] (D_t^Sigma u^a)(D_t^Sigma u^a)
      - 1/2 mu_br (epsilon^{abc} partial_b u_c)(epsilon^{ade} partial_d u_e)
```

`varrho_br[rho]` is the layer mass-density functional `int_layer dn m rho`, dimension `M L^-3`, medium-related and not a free input. `mu_br` is an independent postulated MacCullagh modulus, dimension `M L^-1 T^-2`.

Spin/couple-stress sector:

```text
No independent micro-rotation field is added.
The micro-rotation reservoir reuses P^i from T0.
I_PSigma[rho] := int_layer dn m rho a^2
K_PSigma[rho] := int_layer dn m rho c_s^2(rho) a^2
G_PSigma[rho] := int_layer dn m rho c_s^2(rho)
```

These are surface reductions of the already frozen `P^i` inertia, Frank stiffness, and magnitude stiffness. This is a zero-new-DOF identification; Gate L must still decide whether it actually closes angular momentum and avoids torque-only behavior.

Frozen `P^i <-> u` coupling operator:

```text
P_parallel^a := Pi^a_i P^i on Sigma_n[rho]
delta P_parallel^a := P_parallel^a - Pbar_parallel^a
L_Pu =
    1/2 J_Pu (D_t^Sigma delta P_parallel^a - D_t^Sigma Omega_u^a)^2
  - 1/2 kappa_Pu (delta P_parallel^a - Omega_u^a)^2
```

`J_Pu` is an independent anchoring inertia with dimension `M L^-1`. `kappa_Pu` is an independent anchoring gap/stiffness with dimension `M L^-1 T^-2`. `Pbar_parallel^a` is the local branch background around which the operator is linearized; it is not chosen in G0 and carries no numerical parameter.

C5 stance:

No `phi`-analog and no longitudinal constraint are included in the baseline package. This is the target-blind minimal MacCullagh package with the C5 obstruction deliberately left able to fail. Therefore `FAIL_C5_LONGITUDINAL_ZERO_MODE` is an expected able-to-fail outcome of Gate L(a-iii), not a G0 result.

### Light-Package Provenance and Branch Budget

| Package element | Frozen choice | Provenance | Branch budget |
|---|---|---|---|
| MacCullagh stiffness | Curl-only layer displacement sector with `mu_br` and `varrho_br[rho]`. | `POSTULATED` constitutive package; `varrho_br` is medium-related, `mu_br` independent. | Baseline only. No Cauchy-elastic replacement is admitted in G0; Cauchy appears only as a later negative control if a gate asks for it. |
| Spin/couple reservoir | Reuse T0 `P^i`; no independent micro-rotation field. | `POSTULATED` zero-new-DOF identification, native to the arrows-as-gyrostat picture. | Baseline only. Adding an independent micro-rotation field would require a fresh freeze and DOF/parameter recount. |
| `P^i <-> u` coupling | Rate anchoring plus gap anchoring through `J_Pu` and `kappa_Pu`. | `POSTULATED`; exact operator frozen before Gate L. | Baseline only. No alternate coupling operator is available for silent rescue. |
| C5 handling | No `phi` analog and no longitudinal constraint. | `POSTULATED` minimal-package absence. | Baseline only. A `phi`-completed or constrained branch would be a new structure requiring fresh G0. |

### Field / DOF Inventory

| Field or object | Status | New DOF count | Provenance | Associated constants/functions |
|---|---|---:|---|---|
| `u^a(sigma,t)` | Dynamical layer displacement, 3 raw in-plane components. No G0 gauge removal is asserted. | 3 raw dynamical components | `POSTULATED` light package | `varrho_br[rho]`, `mu_br`, `J_Pu`, `kappa_Pu`, `delta_Sigma[rho]` |
| `phi` analog | Absent. No auxiliary scalar potential and no Lagrange multiplier constraint is frozen. | 0 | Deliberate minimal-package absence | none |
| independent micro-rotation | Absent. Reuses T0 `P^i`; no independent Cosserat rotation field. | 0 beyond T0 `P^i` | `POSTULATED` zero-new-DOF identification | `I_PSigma[rho]`, `K_PSigma[rho]`, `G_PSigma[rho]` |
| `P^i(X,t)` | Already frozen in T0 as a dynamical soft-spin polar vector carried by the medium. | 0 new at G0; inherited T0 field | `INDEPENDENTLY_MOTIVATED_NEW_PARENT_ACTION` from T0 | `m`, `rho`, `a`, `c_s(rho)`; plus `J_Pu`, `kappa_Pu` in the G0 coupling |
| `Sigma_n[rho]`, `delta_Sigma[rho]` | Pure support functional derived from the future smectic profile; not a prescribed layer. | 0 | `POSTULATED` routing convention for layer-effective terms | no free shape function; dimension `[delta_Sigma]=L^-1` |

## G0.4 Complete Parameter Ledger

Independent-new-input count used for the line-1 verdict:

- Active baseline branch `B0_Lifshitz` plus light package: `5`.
- If all named sensitivity branches are activated as a frozen menu, branch-only additional inputs are `4`, for a menu total of `9`.

Since the active baseline already has `5 >= 2`, NG5 pressure is raised immediately as `SECOND_MEDIUM_DRIFT_AT_FREEZE`. This count is not reduced by writing new constants as multiples of `rho,K,m,a,c_s`; no such multiplier has been hidden.

| Parameter/function | Dimension | Role | Classification | Independent-new-input count | Notes |
|---|---:|---|---|---:|---|
| `hbar` | `M L^2 T^-1` | GNLS quantum scale | existing medium-related | 0 | Kept. |
| `m` | `M` | GNLS constituent mass | existing medium-related | 0 | Kept. |
| `q_star` | model charge unit | GNLS gauge coupling in `v_i` | existing medium-related | 0 | Kept. |
| `K` | `M L^18 T^-2` | EOS coefficient | existing medium-related | 0 | Kept. |
| `a` | `L` | existing medium length | existing medium-related | 0 | Kept from T0. |
| `rho`, `rho0` | `L^-4` | density/background density | medium state/background | 0 | `rho0` is a background value, not a new coupling. |
| `c_s^2(rho)=5K rho^4/m` | `L^2 T^-2` | sound speed function | derived medium-related | 0 | Kept. |
| `P^i` normalization `|P|=1` | `1` | T0 polar OP scale | fixed convention | 0 | Byte-identical T0. |
| `c_L1` | `M L^8 T^-2` | Family L negative-gradient driver strength | independent new input | 1 | Active baseline. |
| `c_L2` | `M L^10 T^-2` | Family L stabilizing fourth-gradient/range coefficient | independent new input | 1 | Active baseline. |
| `k_Lstar=sqrt(c_L1/(2 c_L2))` | `L^-1` | Family L roton scale | derived from new inputs | 0 | Not counted separately. |
| `A_R` | `M L^6 T^-2` | Family R Fourier-kernel amplitude | branch-only independent new input | 1 | Sensitivity `S_R_kernel`. |
| `k_R` | `L^-1` | Family R kernel range/scale | branch-only independent new input | 1 | `R_R=1/k_R` derived. |
| `f_R(x)=(x^4-2x^2) exp(-x^2)` | `1` | Family R kernel shape | fixed non-payoff convention | 0 | No shape parameter or free function. |
| `k_Rstar=k_R sqrt(2-sqrt(2))` | `L^-1` | Family R negative-kernel minimum scale | derived from `k_R` | 0 | Not counted separately. |
| `lambda_Cdiv` | `M L^3 T^-2` | Family C divergence coupling | branch-only independent new input | 1 | Sensitivity `S_L_plus_Cdiv`. |
| `chi_Cpin` | `M L^8 T^-2` | Family C gradient-alignment coupling | branch-only independent new input | 1 | Sensitivity `S_L_plus_Cpin`; NG2 risk if negative. |
| `u^a` | `L` | layer displacement field | `POSTULATED` field, not a parameter | 0 in parameter count | Raw DOF counted in G0.3. |
| `varrho_br[rho]=int_layer dn m rho` | `M L^-3` | displacement inertia density | medium-related derived functional | 0 | Depends on the B4 layer profile; no free multiplier. |
| `mu_br` | `M L^-1 T^-2` | MacCullagh rotational modulus | independent new input | 1 | Active baseline. |
| `I_PSigma[rho]=int_layer dn m rho a^2` | `M L^-1` | reused `P^i` rotational inertia | medium-related derived functional | 0 | No independent micro-rotation. |
| `K_PSigma[rho]=int_layer dn m rho c_s^2 a^2` | `M L T^-2` | reused `P^i` surface Frank stiffness | medium-related derived functional | 0 | Gate L decides traction vs torque. |
| `G_PSigma[rho]=int_layer dn m rho c_s^2` | `M L^-1 T^-2` | reused `P^i` surface magnitude/gap scale | medium-related derived functional | 0 | Derived from T0 coefficients. |
| `J_Pu` | `M L^-1` | `P-u` anchoring inertia | independent new input | 1 | Active baseline. |
| `kappa_Pu` | `M L^-1 T^-2` | `P-u` anchoring gap/stiffness | independent new input | 1 | Active baseline. |
| `Pbar_parallel^a` | `1` | branch background for linearized coupling | profile-derived background | 0 | Not selected in G0. |
| `delta_Sigma[rho]` | `L^-1` | layer measure in 4D density representation | derived support functional | 0 | No free support shape. |
| `phi` analog | absent | C5 scalar-potential completion | absent/fixed to none | 0 | C5 remains able to fail. |
| longitudinal constraint multiplier | absent | C5 constraint completion | absent/fixed to none | 0 | C5 remains able to fail. |
| independent micro-rotation field | absent | Cosserat reservoir | absent; reuses `P^i` | 0 | No new micro-rotation DOF. |
| impedance/leak coefficient `epsilon_leak` | absent | leak bound/reservoir impedance | absent/fixed to none | 0 | Leak is not pre-bounded in G0. |

## G0.5 Combined Freeze Block and Hash

This report contains exactly one fence labelled `freeze-action`: the combined block below. The T0 block is embedded inside it byte-for-byte, so the same T0 bytes are protected by both the T0 hash and the new G0 hash.

```freeze-action
PathA_25_G0_combined_action_contract:

Kept GNLS medium:
  Bulk coordinates X^i=(x,y,z,w), i=1..4.
  psi := sqrt(rho) exp(i theta), rho=|psi|^2.
  L_psi =
      (i hbar/2)(psi^* D_t psi - psi D_t psi^*)
    - (hbar^2/(2m))(D_i psi)^*(D_i psi)
    - V_conf(X;Sigma) rho
    - U(rho).
  D_t psi := partial_t psi + (i q_star/hbar) A_0 psi.
  D_i psi := partial_i psi - (i q_star/hbar) A_i psi.
  U(rho) := (K/4) rho^5.
  P_eos(rho) := K rho^5.
  c_s^2(rho) := 5 K rho^4 / m.
  v_i := (hbar/m) partial_i theta - (q_star/m) A_i.
  No new G0 term uses V_conf, a prescribed layer profile, a fixed layer normal,
  or a brane-localized support chosen before Gate B4.

Embedded byte-identical T0 polar-OP freeze-action block:
Field:
  P^i(X,t), i=1..4, a dimensionless polar vector carried by the GNLS medium.

Polarity:
  P and -P are distinct microscopic orientations. There is no director quotient P ~ -P.

Carrier rule:
  The polar field has no standalone vacuum action. Every dynamical OP term is weighted by rho and uses the GNLS material flow v_i. When rho -> 0 the OP energy and inertia vanish with the medium.

Material derivative:
  D_t^v P^i := partial_t P^i + v^j partial_j P^i,
  v_i := (hbar/m) partial_i theta - (q_star/m) A_i.

Local sound-speed function:
  c_s^2(rho) := 5 K rho^4 / m.

Frozen T0 polar-OP Lagrangian density:
  L_pol =
      (1/2) m rho a^2 (D_t^v P^i)(D_t^v P^i)
    - (1/2) m rho c_s^2(rho) a^2 (partial_j P^i)(partial_j P^i)
    - (1/4) m rho c_s^2(rho) (P^i P^i - 1)^2.

Frozen extended parent-action grammar:
  S_T0 = S_GNLS_existing + int dt d^4X L_pol.

Frozen baseline branch:
  O(4)-isotropic soft-spin polar vector with one-constant Frank stiffness,
  single-well magnitude potential |P|=1, no explicit easy axis, no brane-localized term,
  no gauge charge assigned to P except through v_i in the material derivative.

Layering driver family:
  delta_rho := rho - rho0.

  Family L baseline B0_Lifshitz:
    E_L =
        - (1/2) c_L1 (partial_i rho)(partial_i rho)
        + (1/2) c_L2 (partial_i partial_i rho)(partial_j partial_j rho).
    L_L := - E_L.
    c_L1 > 0, c_L2 > 0.
    k_Lstar := sqrt(c_L1/(2 c_L2)).

  Family R sensitivity S_R_kernel:
    S_R =
      - (1/2) int dt d^4X d^4R
        delta_rho(X) V_R(|R|) delta_rho(X+R).
    Equivalently in Fourier variables,
      tilde V_R(k) := A_R f_R(k/k_R),
      f_R(x) := (x^4 - 2 x^2) exp(-x^2).
    A_R > 0, k_R > 0.
    R_R := 1/k_R.
    k_Rstar := k_R sqrt(2 - sqrt(2)).
    No other kernel shape parameter or free-form V_R is admitted.

  Family C sensitivities:
    E_Cdiv := - lambda_Cdiv delta_rho partial_i P^i.
    E_Cpin := (1/2) chi_Cpin (P^i partial_i rho)(P^j partial_j rho).
    L_Cdiv := - E_Cdiv.
    L_Cpin := - E_Cpin.
    Family C terms are admitted only in named branches with Family L or R finite-k regularization.

  Precommitted driver branches:
    B0_Lifshitz: S_driver = int dt d^4X L_L; A_R absent; lambda_Cdiv=0; chi_Cpin=0.
    S_R_kernel: S_driver = S_R; c_L1,c_L2 absent; lambda_Cdiv=0; chi_Cpin=0.
    S_L_plus_Cdiv: S_driver = int dt d^4X (L_L + L_Cdiv).
    S_L_plus_Cpin: S_driver = int dt d^4X (L_L + L_Cpin).

Light-sector package:
  Layer support:
    Sigma_n[rho] denotes the smectic density layers selected by the frozen driver.
    delta_Sigma[rho] is the corresponding layer measure with dimension L^-1.
    No numerical layer profile, layer spacing, layer amplitude, or layer normal is fixed in G0.

  In-plane displacement:
    u^a(sigma,t), a=1..3, is a dynamical layer displacement, [u]=L.
    D_t^Sigma := partial_t + v_parallel^a partial_a.
    Omega_u^a := (1/2) epsilon^{abc} partial_b u_c.

  MacCullagh sector on each Sigma_n:
    L_Mac =
        (1/2) varrho_br[rho] (D_t^Sigma u^a)(D_t^Sigma u^a)
      - (1/2) mu_br (epsilon^{abc} partial_b u_c)(epsilon^{ade} partial_d u_e).
    varrho_br[rho] := int_layer dn m rho.
    mu_br is an independent postulated modulus.

  Spin/couple-stress sector:
    No independent micro-rotation field is added.
    The Cosserat/gyrostat reservoir reuses the T0 polar field P^i.
    P_parallel^a := Pi^a_i P^i on Sigma_n[rho].
    I_PSigma[rho] := int_layer dn m rho a^2.
    K_PSigma[rho] := int_layer dn m rho c_s^2(rho) a^2.
    G_PSigma[rho] := int_layer dn m rho c_s^2(rho).
    These are identified surface reductions of the already-frozen T0 P^i sector,
    not duplicate independent fields.

  P-u coupling:
    delta P_parallel^a := P_parallel^a - Pbar_parallel^a.
    L_Pu =
        (1/2) J_Pu (D_t^Sigma delta P_parallel^a - D_t^Sigma Omega_u^a)^2
      - (1/2) kappa_Pu (delta P_parallel^a - Omega_u^a)^2.
    J_Pu is an independent anchoring inertia.
    kappa_Pu is an independent anchoring gap/stiffness.
    Pbar_parallel^a is a branch background, not a G0-chosen parameter.

  C5 stance:
    No phi-analog and no longitudinal constraint are included.
    A C5 longitudinal zero mode is left as an expected able-to-fail Gate L(a-iii) outcome.

Full action grammar:
  S_G0 =
    S_GNLS_existing
    + int dt d^4X L_pol
    + S_driver
    + sum_n int dt d^3sigma (L_Mac + L_Pu),
  equivalently with layer terms in bulk-density form using delta_Sigma[rho].
```

G0 hash verification command:

```sh
awk '/^```freeze-action$/ {on=1; next} /^```$/ && on {exit} on {print}' software/stage1_solver/reports/pathA_25_G0_freeze.md | sha256sum
```

Frozen G0 SHA-256:

```text
f00ee99d465e2e311c68f47fcacf4af0154ca650642271ab66c36d112cb6a290  -
```

Short G0 hash: `f00ee99d465e`.

## G0.6 Dimensional Check Pointer and k->0 Limit

Dimensional check report: `software/stage1_solver/reports/pathA_25_G0_dimcheck.md`.

Headline result: `G0_DIMCHECK_PASS`; both SymPy and Mathematica exited `0`, reverified the T0 and G0 hashes, agreed on all restored-unit reductions, and verified the strict `k->0` no-EOS-shift limit.

The explicit long-wavelength statement frozen for the dim-check is:

- Family L baseline: the driver quadratic density contribution is `c_L2 k^4 - c_L1 k^2`, so the strict `k->0` contribution to the density compressibility/EOS is zero.
- Family R sensitivity: `tilde V_R(0)=0` and `tilde V_R(k)=-2 A_R k^2/k_R^2 + O(k^4)`, so it also has no `k=0` EOS shift.
- Family C sensitivities enter as `O(k)` or `O(k^2)` couplings and carry no constant density-potential term.
- The light package is layer-effective and does not alter the bulk `U(rho)=K rho^5/4` definition at `k=0`.

Thus the frozen driver is a finite-`k` feature, and the `k->0` `c_s`/EOS remains available for Gate B to test rather than assumed passed.

---

## POST-FREEZE COMMENTARY (orchestrator audit trail — does NOT modify the frozen block or hash)

**Tri-review of G0 (2026-06-23) — verdict accepted; user-gated to proceed to Gate B4.**

1. **Arbiter re-run (orchestrator):** both engines independently re-run — `python3 …_G0_dimcheck_sympy.py` and `math -script
   …_G0_dimcheck.wl` both exit 0 and print `PASS`; both reverify the T0 hash `8fa41ac51e88…` and the G0 hash `f00ee99d465e…`,
   assert the embedded T0 bytes unchanged, and agree on all restored-unit reductions + the strict `k→0` no-EOS-shift limit.
2. **Transliteration-fidelity audit (clean agent): `FIDELITY_CLEAN`.** Every *new* term (driver L/R/C + the 7 light-sector terms +
   constants) is genuinely dimension-reduced from base triples (not hardcoded); hash/byte guards are real and able-to-fail; both
   engines compute the same quantity set. **Two caveats (non-blocking):** (a) the kept GNLS `L_psi` + T0 `L_pol` terms are guarded by
   the hash/byte-identity, **not** re-unit-checked here — so G0.6's "agreed on all restored-unit reductions" is narrower than the
   prose implies (those terms were dim-checked in their own freezes; the auditor hand-verified they remain dimensionally sound);
   (b) the **Family-C `k→0` check is tautological** (`lim λk = 0`, asserted not derived) — but Family C is **sensitivity-only**; the
   **baseline Family-L `k→0` is genuinely able-to-fail** (real limit + derivative-zero at `k*` + 2nd-deriv = `4 c_L1`). **Action:**
   fix the Family-C `k→0` derivation **before any C-sensitivity branch is run**; baseline B4 (Family L) is unaffected.
3. **Adversarial review (clean agent): `ADVERSARIAL_SOUND`.** Independently **reproduced the 5-count** (`c_L1,c_L2,μ_br,J_Pu,κ_Pu`);
   found **no under-counting** (no free knob disguised as `α·a`/`β·c_s`; every `medium-related` functional is a bare `∫_layer` with no
   leading coefficient) and **no padding** (`k*`, `R_R` correctly excluded as derived). Confirmed: target-blind (Family L ranked first
   by locality, not payoff; `f_R` shape + NG2 risk recorded-not-used), complete R/L/C family audit, explicit & coherent C5 stance (no
   φ-analog → C5 left as an expected able-to-fail Gate-L(a-iii) outcome — the harder, more-falsifiable choice), native terms clean.
   Correctly noted the "zero-new-DOF `P^i` reuse" + coefficient-free surface reductions are **conditional on B4/Gate L** (the freeze
   self-flags this; correctly the gate's job, sound at G0).

**Disposition:** `SECOND_MEDIUM_DRIFT_AT_FREEZE` accepted as an honest, expected first-class finding (the "boring" NG5 no-go — 5
calibration inputs, not a requirement contradiction). Per §14: recorded, **not stopping**. User gate (2026-06-23): **proceed to Gate
B4** — the drift becomes the calibration budget, judged later against held-out surplus.
