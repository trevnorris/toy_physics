# 5PN progress notes — explicit overlap model and weak-axisymmetric transport

## What was added

### Stage 3 — explicit isotropic finite-throat overlap model

Files:
- `5pn_stage3_isotropic_overlap_model.py`
- `5pn_stage3_isotropic_overlap_model_output.txt`

Model choice:
- finite interval `s in [0,L]`
- wall profile
  `chi_eta(s) = sqrt(2/L) sin(pi s / L)`
- D/N-like half-wave
  `phi_DN(s) = sqrt(2/L) sin(pi s / (2L))`
- simplest brane-like gauge profile `u = chi_eta`
- simplest mixed profile `w = phi_DN`

Exact overlaps:
- `I_eta_u = 1`
- `I_eta_phi = I_eta_w = I_uw = 8/(3 pi)`

Using one BdG mode and one conservative Maxwell/mixed pair, the script constructs
- `B0,B2,B4`
- `Z0,Z2,Z4`
- `N0,N2,N4`
- `D0,D2,D4`
- `u2,u4`
- `P0,P2,P4`

and verifies exact grouped isotropy when the three grouped lanes share the same radial/axial data.

### Most useful Stage-3 result

In this explicit prototype, two important requirements both fix the same static wall stiffness `K`:

1. the 5PN conservative one-pole identity
   `u4 = 4 u2^2`
2. the universal 2.5PN / 4PN normalization target for `P0`

So simultaneous success reduces to one explicit compatibility equation:

`N0 / P0_target = 3 (M + B2 + Z2)^2 / (B4 + Z4)`

This is the first concrete radial/axial compatibility surface I have for the 5PN program.
It is the explicit prototype version of “the PDE must make the conservative and outgoing calibrations agree on one branch.”

## Stage 4 — weak-axisymmetric transport and microscopic grouped obstructions

Files:
- `5pn_stage4_axisymmetric_transport.py`
- `5pn_stage4_axisymmetric_transport_output.txt`

This stage starts from the exact weak-axisymmetric grouped signature
- `lambda_(20) = 1`
- `lambda_(21) = 1/2`
- `lambda_(22) = -1`

and carries the grouped operator / transfer slopes through the full Stage-154/155/156 logic.

### Exact microscopic slope decomposition

Using
- wall slopes `(K1_wall, M1_wall)`
- BdG slopes `(B01, B21, B41)`
- conservative Maxwell/mixed slopes `(Z01, Z21, Z41)`
- outgoing transfer slope `N01`

it verifies

`K_1 = D21 + D01/9 = W1 - Bcal1 - Zcal1`

and

`G_1 = N01 - P0 D01 = -P0 K1_wall + P0 B01 + Nbundle1`

with
- `D01 = K1_wall - B01 - Z01`
- `D21 = -(M1_wall + B21 + Z21)`
- `D41 = -(B41 + Z41)`

### Physical weak-axisymmetric slopes

The script derives

`u2^(1) = -(D21 + u2 D01)/D0`

`P1/P0 = N01/N0 - D01/D0`

and on the canonical compensated branch verifies the hidden-even relation

`u4^(1) = (8/9) u2^(1)`

iff

`D41 = (2/3) D21 + D01/27`

### Most useful Stage-4 result

On the even-preserving branch `u2^(1)=0`, the conservative grouped response collapses to

- `D21 = -D01/9`
- `D41 = -D01/27`

and the remaining linear grouped `2.5`PN defect becomes one scalar loading mismatch

`Xi_load = N01/N0 - D01/D0 = P1/P0`

with fixed grouped-lane signature

- `(20,21,22) ~ (1, 1/2, -1)`
- equivalently `b = 3 a`

So after Stage 4, the next theorem gate is no longer “solve all grouped anisotropies.”
It is much narrower:

> compute `D01` and `N01` — or directly `Xi_load` — from a primitive weak-axisymmetric deformation of the explicit finite-throat overlap model.

## Run status

Both new scripts completed here without running out of time or memory.

## Stage 5 — primitive deformation and exact compensation surfaces

Files:
- `5pn_stage5_primitive_deformation_compensation.py`
- `5pn_stage5_primitive_deformation_compensation_output.txt`

This stage takes the explicit isotropic overlap model from Stage 3 and adds a primitive weak-axisymmetric microscopic deformation through the slopes

- `dK`, `dM`
- `d(lambda_B)`, `d(varpi)`
- `d(lambda_U)`, `d(lambda_W)`, `d(lambda_R)`
- `d(Omega_U)`, `d(Omega_W)`

Then it computes the induced grouped-lane slope data

- `D01`, `D21`, `D41`
- `N01`

and from them the three key combinations

- `K1 = D21 + D01/9`
- `G1 = N01 - P0 D01`
- `Xi_load = N01/N0 - D01/D0 = G1/N0`

### Most useful Stage-5 results

1. **Even-preserving compensation** `K1 = 0` is algebraic and fixes the inertia-side slope `dM` exactly.
2. **Odd/normalization-preserving compensation** `Xi_load = 0` is algebraic and fixes the static loading slope `dK` exactly.
3. The remaining explicit **5PN hidden-even residual** is
   `D41 - (2/3) D21 - D01/27`.

So after Stage 5, the next theorem gate is now extremely sharp:

> choose a concrete primitive anisotropy mechanism — wall-only, BdG-only, Maxwell/mixed-only, or a mixed combination — substitute its slopes into Stage 5, and test
> `K1 = 0`, `Xi_load = 0`, and `D41 - (2/3) D21 - D01/27 = 0`.

## Run status

All three new scripts completed here without running out of time or memory.
