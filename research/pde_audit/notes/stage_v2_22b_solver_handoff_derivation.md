# Stage V2-22B — Numerical Branch Profile Schema and Solver Handoff

## 0. Purpose

Stage V2-22A created a profile-to-coefficient adapter:

\[
\text{profiles}\rightarrow \text{overlaps}\rightarrow B_n,Z_n,N_n\rightarrow D_n,u_n,P_n.
\]

Stage V2-22B freezes the **upstream solver export contract**.  Its job is to define exactly what a future moving-throat PDE solver must output before the V2-22A adapter and V2-21 extraction fixture are allowed to evaluate target residuals.

This is an anti-overfitting and branch-realization gate.  The solver must provide a **target-blind, pre-target-frozen, open-throat branch packet** with enough sampled data to compute the wall/BdG/Maxwell/mixed overlaps.  If the packet is a hard-cap branch, omits the mixed sector, has unstable primitive signs, or fails the grid/profile checks, it is rejected before any quadrupole target residual is computed.

## 1. Solver export schema

The frozen solver output schema is:

```text
stage_v2_22b_solver_handoff/v1
```

The required top-level fields are:

```text
schema
branch_id
freeze
geometry
constants
grid
wall
profiles
bdg_modes
mixed_ports
```

The required freeze block is:

```json
{
  "pre_target_freeze": true,
  "target_blind": true,
  "parent_action": "declared GNLS + localized Maxwell + wall/action convention",
  "gauge_convention": "declared gauge convention",
  "boundary_protocol": "open_impedance_AC_reflecting_DC_leaking"
}
```

The required geometry block is:

```json
{
  "L": "positive finite throat length",
  "R_mouth": "positive mouth radius",
  "R_exit": "strictly positive open exit radius",
  "boundary_class": "open_impedance",
  "Y_L_limit": "finite load-admittance parameter",
  "exit_model": "organ-pipe / impedance mismatch open exit model"
}
```

The validator explicitly forbids:

\[
R_{\rm exit}\le 0,
\qquad
\texttt{boundary\_class}=\texttt{hard\_cap}.
\]

So the old capped-tube interpretation cannot enter the branch-realization packet.

## 2. Mathematical checks inside the validator

The script verifies the D/N support profile convention:

\[
q(s)=\sin(ks),
\qquad
q(0)=0,
\qquad
q_s(L)=0
\quad\text{for}\quad
k=\frac{\pi}{2L}.
\]

It also checks the open-impedance endpoint law:

\[
T_w q_w(L,\omega)+Y_L(\omega)q(L,\omega)=0.
\]

In the AC support-mode reflection limit,

\[
Y_L\to0,
\]

this becomes

\[
q_w(L)=0,
\]

which is the Neumann end of the D/N ladder.  This is compatible with the V2-04 open-organ-pipe patch: AC support modes can reflect through impedance mismatch while the DC/background flow exits through the open finite-radius conduit.

The script also verifies the exact D/N overlap prototype:

\[
\chi_\eta(s)=\sqrt{\frac2L}\sin\frac{\pi s}{L},
\qquad
\phi_{DN}(s)=\sqrt{\frac2L}\sin\frac{\pi s}{2L},
\]

with

\[
\int_0^L \chi_\eta^2\,ds=1,
\qquad
\int_0^L \phi_{DN}^2\,ds=1,
\qquad
\int_0^L \chi_\eta\phi_{DN}\,ds=\frac{8}{3\pi}.
\]

For weak-axisymmetric grouped data, the script confirms the fixed signature:

\[
(20,21,22)\sim\left(1,\frac12,-1\right),
\qquad
b=3a.
\]

For the Maxwell/mixed block, it verifies that the leading outgoing transfer factor is a perfect square:

\[
N_0
=
\frac{(\Omega_U^2 g_W+Rg_U)^2}{(\Omega_U^2\Omega_W^2-R^2)^2}.
\]

So, once

\[
\Delta=\Omega_U^2\Omega_W^2-R^2>0,
\]

is enforced, the outgoing transfer has nonnegative static weight.

## 3. Runtime validation gates

For a solver output packet, the validator checks:

1. required schema keys and declared freeze metadata;
2. target-blind pre-target freeze status;
3. open-throat geometry;
4. monotone grid from \(0\) to \(L\);
5. finite sampled profiles with lengths matching the grid;
6. wall inertia/stiffness positivity;
7. BdG mode positivity \(\varpi>0\);
8. at least one mixed Maxwell/\(A_w\) port;
9. effective mixed-block positivity;
10. availability of the constants needed by V2-21/V2-22A.

The profile diagnostics compute, for example,

\[
\int_0^L \mu_s(s)\chi_\eta(s)^2\,ds
\]

and each BdG-mode norm.  Norm deviations generate warnings, not automatic failures, because a real PDE solver might output nonunit profiles plus explicit normalization conventions.  Structural issues such as hard-cap geometry or nonpositive \(\Delta\) are errors.

## 4. Conversion into the V2-22A profile manifest

If validation passes, the script converts the solver packet into the V2-22A profile manifest schema:

```text
stage_v2_22a_profile_adapter/v1
```

It maps:

```text
solver profiles.weight       -> profile "weight"
solver profiles.wall_chi_eta -> profile "wall"
solver bdg_modes[*]          -> sampled BdG profiles
solver mixed_ports[*].u/w    -> sampled U/W mixed profiles
```

and preserves the reduction coefficients:

\[
K,
\quad
M,
\quad
\lambda_B,
\quad
\varpi,
\quad
\lambda_U,
\quad
\Omega_U,
\quad
\lambda_W,
\quad
\Omega_W,
\quad
\lambda_R.
\]

The generated profile manifest is then immediately consumable by the V2-22A adapter.

## 5. Built-in valid sample result

The built-in sample solver export is an open D/N finite-throat prototype with:

```text
R_exit = 0.35
boundary_class = open_impedance
pre_target_freeze = true
target_blind = true
```

The V2-22B validation result was:

```text
symbolic_checks: 9/9 passed
validation_pass: True
error_count: 0
warning_count: 0
packet_hash: 959ba8b19d5c8afc683c006e0214c0fa2794f873ae7b611bf6d8a527df66b7d1
generated_profile_manifest_hash: 46da7fc420e753dd749e3e590f24614b82705399497fd79a3e227505dcba0a9c
```

A smoke test through V2-22A then produced:

```text
I_eta_phi::bdg_halfwave = 0.8488255450333207
I_eta_u::one_mixed_port = 1
I_eta_w::one_mixed_port = 0.8488255450333207
I_u_w::one_mixed_port = 0.8488255450333207
```

and the downstream observable packet remained open/stable but target-failing, as expected for an untuned prototype:

```text
open_gate_pass: True
stability_gate_pass: True
target_packet_pass: False
D0_bar: 1.560684600911219
N0_bar: 0.0008248594058984478
P0_bar: 0.0005285240883499759
P0_target: 10.8
R_norm: -10.79947147591165
```

This is the intended behavior.  V2-22B and V2-22A are not supposed to rescue a branch; they only validate, convert, and extract.

## 6. Built-in invalid hard-cap rejection

The invalid sample changes the exit to:

```text
R_exit = 0
boundary_class = hard_cap
boundary_protocol = hard_cap_DC_blocked
```

The validator rejected it:

```text
validation_pass: False
error_count: 4
```

The specific errors were:

```text
freeze.boundary_protocol: must declare open impedance AC/DC split
geometry.boundary_class: boundary_class must be open_impedance
geometry.R_exit: R_exit must be strictly positive; hard cap is forbidden
geometry.boundary_class: hard_cap geometry is forbidden by V2-04 patch
```

This confirms that the capped-tube mass-continuity bug cannot silently re-enter the PDE branch pipeline.

## 7. Status

V2-22B passes as a solver-handoff layer:

\[
\boxed{\text{PDE solver output}\rightarrow\text{validated frozen branch packet}\rightarrow\text{V2-22A profile manifest}.}
\]

It does not solve the moving-throat PDE.  It fixes the exact data contract required before an actual numerical branch can be judged.

## 8. Immediate continuation

The next natural stage is V2-22C:

```text
profile-to-fixture end-to-end smoke pipeline and tolerance budget
```

That stage should run the generated V2-22A manifest automatically through the V2-22A adapter and V2-21 extraction fixture, then emit a single end-to-end branch-realization report with:

```text
solver packet hash
profile manifest hash
V2-21 manifest hash
observable packet hash
open/stability gates
isotropy gates
target residuals
profile normalization diagnostics
```

That would turn the current validation stack into a one-command branch-realization pipeline.
