T0_SHEAR_FROZEN(d9520d3819c3) + SECOND_MEDIUM_DRIFT_AT_FREEZE(11)

# pathA_35 G0 Freeze Artifact: Light-Confining Shear-Surface Brane

artifact_status: append-only immutable freeze
directive: `software/stage1_solver/directives/pathA_35_shear_surface_brane_gates.md`
stage: `G0 structure freeze`
date: 2026-06-26
scope: reports-only plus dual-engine dimensional and flat-brane DOF count checks; no Gate L verdict computed

## Verdict

`T0_SHEAR_FROZEN(d9520d3819c3)`

`SECOND_MEDIUM_DRIFT_AT_FREEZE(11)`

This is a structure freeze, not a gate. The action below specifies the postulated surface and light-sector ingredients, hashes them, dimensionally checks every constructed term in both engines, and reports the flat-brane linear DOF count. It does not assert that any mode is bounded-below, traction-carrying, gauge, gapped in the coupled system, or leak-free.

## G0.1 Kept GNLS and T0 Fidelity

The GNLS medium is kept unchanged:

- `psi=sqrt(rho) exp(i theta)` with `rho=|psi|^2`.
- `rho` is the T0 number density, so the bulk mass density is `m rho`.
- Quantum pressure is retained as `(hbar^2/(8 m rho)) (partial_i rho)(partial_i rho)`.
- `U(rho)=(K/4) rho^5`, `P_eos(rho)=K rho^5`, and `c_s^2(rho)=5 K rho^4/m`.
- `v_i=(hbar/m) partial_i theta - (q_star/m) A_i`.
- Circulation/vortices and the shear-free GNLS bulk are kept.

The T0 polar-OP action is kept byte-for-byte from `software/stage1_solver/reports/pathA_24_T0_freeze.md`, with short hash `8fa41ac51e88`. The exact T0 freeze-action bytes are embedded inside the canonical block below and are rechecked by both scripts.

## G0.2 Postulated Codim-1 Shear Surface

The brane is postulated at `w=0`; the normal `w_hat` and the surface itself are a conceded inherited wall, not derived here. The target-blind confinement profile is the fixed-shape one-width family

`g_ell(w) = exp(-(w/ell_g)^2)/(sqrt(pi) ell_g)`, `ell_g>0`, with `int dw g_ell(w)=1`.

The in-plane displacement `u^a(x,y,z,t)`, `a in {x,y,z}`, is a surface collective/material DOF of the same medium. It is carried by the medium in the material-rotation sense, but tangentially free-slip from bulk mass transport: `dot u^a != v^a`. The out-of-plane component is `u_w`. The brane inertia is the independent surface mass density `rho_br`.

The kept interface traction is

`T_na = T_wa + (T_ww delta_ab - T_ab) partial_b u_w`, with `T_wa = m rho v_w v_a`.

The factor `m` is part of the frozen structure because `rho` is a number density.

## G0.3 Light-Sector Decisions

Baseline `P^i` mass/gap status: `massless` spin waves, inherited from the T0 soft-spin O(4) polar vector around `|P|=1`. This is the target-blind minimal/single-medium baseline because it adds no new `P` mass constant and no hard/slaving constraint. Named alternates, not active in this freeze: `gapped` and `slaved-rigid`, where `slaved-rigid` means `P_parallel = w_hat x (nabla x u)` and carries a `k^4` correction for Gate L to examine.

Gate-L exposure recorded for the baseline only: `FAIL_HIDDEN_PROPAGATING_MODE`, `FAIL_GYROSTAT_NO_CLOSURE`, `FAIL_NOT_BOUNDED_BELOW`, and the linked `FAIL_COUPLE_STRESS_NOGO` chain remain able to fire. This is not a pass/fail statement.

The parity-even `P-u` operator is frozen as

`varpi_a := (w_hat x P_parallel)_a`,

`L_Pu = - lambda_Pu varpi_a Omega_u^a = - lambda_Pu w_hat . (P_parallel x (nabla_parallel x u))`,

with `Omega_u^a := (nabla_parallel x u)^a`. The direct `P . Omega_u` scalar is excluded from the frozen action because it is parity-odd. The repaired operator requires `w_hat`, ties the couple-stress sector to the conceded axis, and re-admits the epsilon-contracted/chiral class that T0 excluded. This is counted as a structural input. No statement is made here about traction closure.

C5 scalar-potential decision: absent in the baseline. No independent scalar `phi`, no `phi=u_w` identification, and no raw `nabla . u=0` projector are frozen. Gate L(a-iii) remains able to test the C5 longitudinal-sector obstruction. Consistency record: the inactive alternate `phi=u_w` would require a massless scalar potential and is in tension with the frozen bare `u_w` gap; the inactive alternate `phi` as a new field would add a new DOF and a new variational term.

The out-of-plane scalar has the bare local gap term

`L_uw = (1/2) rho_br (partial_t u_w)^2 - (1/2) rho_br Omega_w^2 u_w^2`,

with independent scale `Omega_w>0`. This only freezes the term and scale; it does not assert a Gate L(d) coupled-system result.

## G0.4 Parameter Ledger

Field / DOF sub-count, new at G0: `4`.

| item | tag | G0 new DOF count | note |
|---|---|---:|---|
| GNLS `psi,rho,theta` | kept | 0 | unchanged bulk medium |
| T0 `P^i in R^4` | kept | 0 | reused, not duplicated |
| in-plane `u^a` | postulated-ingredient | 3 | surface collective/material components |
| out-of-plane `u_w` | postulated-ingredient | 1 | bending scalar with bare gap term |
| independent micro-rotation | postulated-ingredient | 0 | absent in the baseline; T0 `P^i` is reused |
| C5 scalar `phi` | postulated-ingredient | 0 | absent in the baseline; no scalar-potential analog frozen |

Independent constants sub-count: `4`.

| constant | dimension | tag | count | note |
|---|---:|---|---:|---|
| `rho_br` | `M L^-3` | postulated-ingredient | 1 | surface inertia |
| `mu_R` | `M L^-1 T^-2` | postulated-ingredient | 1 | MacCullagh modulus |
| `lambda_Pu` | `M L^-1 T^-2` | postulated-ingredient | 1 | parity-repaired `P-u` coupling normalization |
| `Omega_w` | `T^-1` | postulated-ingredient | 1 | bare `u_w` gap scale |
| T0 couple-stress inertia/stiffness | kept | 0 | `m rho a^2`, `m rho c_s^2 a^2`, `m rho c_s^2`; no new constants |

Independent functions sub-count: `1`.

| function | tag | count | note |
|---|---|---:|---|
| `g_ell(w)=exp(-(w/ell_g)^2)/(sqrt(pi) ell_g)` | postulated-ingredient | 1 | fixed shape plus one width parameter; no free-form profile/kernel |

Structural constraints / postulates sub-count: `6`.

| structural item | tag | count | note |
|---|---|---:|---|
| imposed `w_hat` axis and `w=0` surface | conceded-wall | 1 | inherited wall, not derived |
| `u^a` as same-medium surface collective but tangentially free-slip from bulk transport | postulated-ingredient | 1 | `dot u^a != v^a` |
| T0 `P^i` reused as Cosserat micro-rotation reservoir | postulated-ingredient | 1 | zero new micro-rotation DOF |
| baseline `P^i` spin-wave status = `massless`; alternates named | postulated-ingredient | 1 | target-blind branch decision |
| `w_hat`-dependent parity-even `P-u` operator | postulated-ingredient / conceded-axis dependent | 1 | re-admits epsilon-contracted/chiral class |
| no C5 `phi` analog and no longitudinal constraint | postulated-ingredient | 1 | leaves C5 exposure to Gate L |

Independent new input count for drift: constants `4` + functions `1` + structural `6` = `11`.

## G0.5 Flat-Brane Linear DOF Count

Dual-engine computed count for the active baseline; reported values are wired to rank outputs, not literals:

| sector | computed source | count |
|---|---|---:|
| in-plane `u^a` | curl rank transverse to `k_a` | 2 |
| in-plane `u^a` | kinetic rank minus curl rank | 1 |
| T0 `P^i` | tangent kinetic/Frank rank | 3 |
| T0 `P^i` | radial Hessian rank | 1 |
| out-of-plane `u_w` | kinetic/gap scalar rank | 1 |
| C5 scalar potential `phi` | absent baseline rank | 0 |

Total computed flat-brane DOF: `8`.

Ranks used by both engines: `u` kinetic rank `3`, curl rank `2`, curl nullity `1`; T0 `P` tangent kinetic/Frank ranks `3/3`, radial kinetic/soft-spin Hessian ranks `1/1`; `u_w` kinetic/gap ranks `1/1`; `phi` rank `0`.

DOF controls now able-to-fail: dropping `u_w` gap term, dropping the T0 soft-spin radial term, or zeroing the `k`-parallel `u` component changes the computed total to `7`, so each FIRED.

The kept bulk remains shear-free GNLS. This count does not classify any component as propagating versus gauge, bounded, leak-free, or coupled-system gapped.

## Canonical Frozen Action Block

Hash target: exact bytes inside the single fenced block labelled `freeze-action` below, excluding the fence lines.

```freeze-action
PathA_35_G0_shear_surface_brane_action_contract:

Kept GNLS medium:
  Bulk coordinates X^i=(x,y,z,w), i=1..4.
  psi := sqrt(rho) exp(i theta), rho=|psi|^2.
  rho is the T0 number density; bulk mass density is m rho.
  Quantum pressure density:
    QP := (hbar^2/(8 m rho)) (partial_i rho)(partial_i rho).
  U(rho) := (K/4) rho^5.
  P_eos(rho) := K rho^5.
  c_s^2(rho) := 5 K rho^4 / m.
  v_i := (hbar/m) partial_i theta - (q_star/m) A_i.
  Bulk circulation/vortices and the shear-free GNLS bulk are kept.

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

Postulated codim-1 shear surface:
  Surface:
    Sigma_0 := {w=0}.
    w_hat is an imposed normal axis and conceded inherited wall.
  Confinement profile:
    g_ell(w) := exp(-(w/ell_g)^2)/(sqrt(pi) ell_g), ell_g > 0.
    int_{-infty}^{infty} dw g_ell(w) = 1.
    No free-form g(w), kernel, extra shape parameter, or target-tuned profile is admitted.
  Brane displacement:
    u^a(x,y,z,t), a in {x,y,z}, [u]=L.
    u^a is a surface collective/material DOF of the same medium.
    u^a is carried for material rotation/traction bookkeeping but tangentially free-slip from bulk mass transport:
      partial_t u^a != v^a.
    u_w(x,y,z,t) is the out-of-plane bending scalar, [u_w]=L.
    rho_br is an independent brane inertia density, [rho_br]=M L^-3.
  Interface traction retained for Gate L:
    T_na := T_wa + (T_ww delta_ab - T_ab) partial_b u_w.
    T_wa := m rho v_w v_a.

Light-sector package on the flat brane:
  In-plane curl:
    Omega_u^a := (nabla_parallel x u)^a.
  MacCullagh term:
    L_Mac :=
        (1/2) rho_br (partial_t u^a)(partial_t u^a)
      - (1/2) mu_R Omega_u^a Omega_u^a.
    mu_R > 0 is an independent postulated modulus.

  P mass/gap baseline and alternates:
    Active baseline: massless T0 spin-wave modes around |P|=1, with the T0 radial amplitude retained.
    Named inactive alternate: gapped P spin waves.
    Named inactive alternate: slaved-rigid P_parallel = w_hat x Omega_u, carrying a k^4 correction to be examined only at Gate L.

  Couple-stress reservoir:
    No independent micro-rotation field is added.
    The Cosserat micro-rotation reservoir reuses the T0 polar field P^i.
    The T0 coefficients m rho a^2, m rho c_s^2(rho) a^2, and m rho c_s^2(rho)
    are the inherited inertia, Frank stiffness, and soft-spin radial scale.

  Parity-repaired P-u coupling:
    P_parallel^a is the projection of P^i to the tangent brane.
    varpi_a := (w_hat x P_parallel)_a.
    L_Pu := - lambda_Pu varpi_a Omega_u^a
          = - lambda_Pu w_hat . (P_parallel x (nabla_parallel x u)).
    lambda_Pu is an independent coupling normalization.
    The direct P_parallel . Omega_u operator is excluded as parity-odd.
    The retained operator depends on w_hat and re-admits the epsilon-contracted/chiral class excluded by T0.

  Longitudinal/C5 decision:
    No scalar-potential analog phi is present in the active baseline.
    No phi = u_w identification is present.
    No raw nabla_parallel . u = 0 projector or longitudinal constraint is present.
    Inactive consistency record: phi = u_w would require a massless scalar potential and tension with the u_w gap;
    a new phi would add a new field and variational term.

  Out-of-plane gap term:
    L_uw :=
        (1/2) rho_br (partial_t u_w)(partial_t u_w)
      - (1/2) rho_br Omega_w^2 u_w^2.
    Omega_w > 0 is an independent gap scale.

  Brane-localized action representation:
    S_brane = int dt d^4X g_ell(w) [L_Mac + L_Pu + L_uw].

Flat-brane linear bookkeeping used for G0.5 only:
  u-operator:
    O_u_ab := rho_br omega^2 delta_ab - mu_R (k^2 delta_ab - k_a k_b).
    rank(k^2 delta_ab - k_a k_b) = 2 for k_a != 0, with one longitudinal component.
    c_gamma^2 := mu_R / rho_br.
  T0 P spin-wave expression:
    omega_P^2 := c_s^2 k^2.
  T0 P radial expression:
    omega_radial^2 := c_s^2 k^2 + 2 c_s^2/a^2.
  u_w bare expression:
    omega_uw,bare^2 := Omega_w^2.

Full G0 action grammar:
  S_G0 =
    S_GNLS_existing
    + int dt d^4X L_pol
    + int dt d^4X g_ell(w) [L_Mac + L_Pu + L_uw].
```

## Freeze Hash

Extraction command from project root:

```sh
awk '/^```freeze-action$/ {on=1; next} /^```$/ && on {exit} on {print}' software/stage1_solver/reports/pathA_35_G0_freeze.md | sha256sum
```

Frozen-action SHA-256:

```text
d9520d3819c3f718290f9d0be57138c07d5bf02d2237106478e17b6a1e389ac3  -
```

Short frozen-action hash: `d9520d3819c3`

Byte range: content bytes excluding fence lines are `0-based [8110,13310)`, equivalently `1-based inclusive 8111-13310`, length `5200` bytes in this report revision.

## Dimensional Firewall and Script Results

Scripts:

- `software/stage1_solver/tools/pathA_35_G0_sympy.py`
- `software/stage1_solver/tools/pathA_35_G0.wl`

Run commands from `/var/projects/toy_physics`:

```text
timeout 600 python3 software/stage1_solver/tools/pathA_35_G0_sympy.py
timeout 600 math -script software/stage1_solver/tools/pathA_35_G0.wl
timeout 600 python3 software/stage1_solver/tools/pathA_35_G0_sympy.py --compare
```

Result: all three commands exited `0`.

The scripts check, inline as each expression is constructed: kept GNLS dimensions, T0 `L_pol`, confinement profile and normalization measure, brane inertia terms, MacCullagh term, reused T0 couple-stress coefficients, parity-repaired `P-u` term, `u_w` gap term, full projected traction including `T_wa=m rho v_w v_a`, action-measure dimensions, every G0.5 linearization quantity (`O_u`, `c_gamma^2`, `omega_P^2`, `omega_radial^2`, `omega_uw,bare^2`), and the computed flat-brane DOF count from quadratic-form ranks.

Able-to-fail ablations:

- `drop_m_from_T_wa`: corrupted `rho v_w v_a` reduced to `L^-2 T^-2` instead of the required `M L^-2 T^-2`; FIRED.
- `MacCullagh_without_curl`: corrupted `mu_R u^2` reduced to `M L T^-2` instead of the required `M L^-1 T^-2`; FIRED.
- `drop_u_w_gap_term`: computed flat-brane total changed from `8` to `7`; FIRED.
- `drop_P_soft_spin_radial_term`: computed flat-brane total changed from `8` to `7`; FIRED.
- `zero_u_longitudinal_component`: computed flat-brane total changed from `8` to `7`; FIRED.

Engine agreement: `ENGINE_AGREE`; the SymPy `agreement_payload` and Mathematica `agreement_payload` matched exactly.

## Machine Ledger

Machine-readable ledger: `software/stage1_solver/reports/pathA_35_G0_results.yaml`.

---

## POST-FREEZE COMMENTARY

No post-freeze commentary yet.

---

## ⚠ ERRATUM (2026-07-04, appended by NG5 / `pathA_41` — append-only, no frozen content altered)

**The `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)` count (line 1) is NOT inflated by a `ρ_br` overcount — this freeze's count stands.** (An earlier
draft of this note, based on a GLM tertiary catch, claimed the count was inflated because `ρ_br` here matched a *derived* functional
`varrho_br[ρ]=∫_layer dn·m·ρ` in `pathA_25`. The `pathA_41` framing adjudication (`_scratch/pathA_41_framing_codex.md`) SUPERSEDES that:
`pathA_25`'s `varrho_br` belongs to the **CLOSED density-smectic candidate** (`FAIL_NOT_CODIM1`), a *different* structure than this
active **shear-surface** brane. They are not the same active object; `varrho_br` is `OUT_OF_ACTIVE_NG5`.)

**Corrected finding:** this freeze's `ρ_br` (and `μ_R`) is a **genuine postulated shear-surface brane inertia/modulus with a
registered-but-unsolved reduction** — the `pathA_40` cone-lock **Route-A** (R3: "derive BOTH `ρ_br` and `μ_R` from one bulk/throat
integral", currently `ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT`). It is neither B4-derived nor a B4 same-object overcount. So the
count of 11 is not reduced. **The honest cross-sector drift** (per `pathA_41`) is instead `{ρ_B0, χ_c, C_hu}` — the on-brane compression +
embedding-mixing parameters whose 4D→3D reduction is not yet registered. Erratum NOTE only; no frozen content / earned physics / dimcheck
altered. Authoritative adjudication lives in `pathA_41`.
