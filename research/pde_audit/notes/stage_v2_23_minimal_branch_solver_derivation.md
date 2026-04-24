# Stage V2-23 — Minimal open-throat branch solver and first real residual extraction

## 0. Purpose

This stage is the first reduced branch-realization prototype after the V2-22C
handoff pipeline.

It is **not** a full nonlinear GNLS/Maxwell/moving-wall PDE solution. It is a
target-blind, frozen, one-dimensional open-throat branch solve that moves beyond
mock profiles by solving:

1. a stationary finite-radius open throat profile,
2. a wall `l=2` Sturm-Liouville support profile,
3. a stable BdG support profile,
4. a brane-like localized gauge profile,
5. a mixed `A_w/F_{\mu w}/J^w` profile,

then computing overlap integrals and passing the resulting branch through the
same isotropic full-bundle residual packet used in V2-19 through V2-22C.

The branch is intentionally simple. Its job is to give Codex/local solvers a
concrete executable template for the next full PDE branch export.

## 1. Frozen branch definition

```text
branch_id: v2_23_minimal_open_throat_frozen_demo
branch_freeze_hash: 0feee821ffec9b23f79d80bb4e176139f01fbffb3f1358607dad6597032592ae
pre_target_freeze: true
target_blind: true
no_post_residual_refit: true
boundary_class: open_impedance
boundary_protocol: open_impedance_AC_reflecting_DC_leaking
```

The geometry is an open finite throat, not a capped tube:

\[
R(0)=a,\qquad R(L)>0.
\]

The stationary profile is obtained by minimizing

\[
E[R]
=
\frac12\int_0^L T_R (R')^2\,ds
+
\frac12\int_0^L K_R(R-R_{\rm pref})^2\,ds
+
\frac12Y_{\rm exit}(R(L)-R_{\rm exit,pref})^2,
\]

with fixed mouth radius \(R(0)=a\). The exit condition is a finite-radius open
penalty rather than a cap.

The solved branch gives

\[
R_{\rm mouth}=1,\qquad
R_{\rm exit}=0.452938042901,\qquad
R_{\rm min}=0.452938042901.
\]

The open gate passes because \(R_{\rm exit}>0\).

## 2. Axial measure

The solver uses the geometry-derived effective axial measure

\[
\mu_s(s)=R_0(s)^2\sqrt{1+R_0'(s)^2},
\]

renormalized so that

\[
\int_0^L\mu_s(s)\,ds=L.
\]

This remains a reduced one-dimensional model, but the overlap integrals now
depend on the solved open-throat geometry rather than on a flat hand-inserted
measure.

## 3. Solved Sturm-Liouville problems

Each profile is obtained from a finite-element Sturm-Liouville problem

\[
-\partial_s(T(s)\partial_s q)+V(s)q=\lambda\mu_s(s)q,
\]

with mouth Dirichlet condition

\[
q(0)=0,
\]

and open-end impedance condition

\[
T(L)q'(L)+Y_Lq(L)=0.
\]

For this first branch prototype the AC-reflecting organ-pipe limit is used:

\[
Y_L=0,
\]

so the exit is Neumann-like for AC support modes while remaining geometrically
open for DC/background flow.

The solved eigenvalues are:

\[
K=\lambda_{{\rm wall},l=2}=2.2393180779,
\]

\[
\varpi^2=3.16919632623,
\]

\[
\Omega_U^2=3.78509395487,\qquad
\Omega_W^2=4.05724378721.
\]

The FEM residuals are reported in `stage_v2_23_tolerance_report.json`.

## 4. Overlap integrals

The reduced overlap integrals are

\[
I_{\eta\phi}=\int_0^L\mu_s\chi_\eta\phi\,ds,
\]

\[
I_{\eta U}=\int_0^L\mu_s\chi_\eta u\,ds,
\]

\[
I_{\eta W}=\int_0^L\mu_s\chi_\eta w\,ds,
\]

\[
I_{UW}=\int_0^L\mu_s u w\,ds.
\]

For this frozen branch,

\[
I_{\eta\phi}=0.999946613757,
\qquad
I_{\eta U}=0.999180452256,
\]

\[
I_{\eta W}=0.994783990782,
\qquad
I_{UW}=0.998096444041.
\]

The reduced couplings are

\[
c_B=\lambda_B I_{\eta\phi},
\qquad
g_U=\lambda_U I_{\eta U},
\qquad
g_W=\lambda_W I_{\eta W},
\qquad
R=\lambda_R I_{UW}.
\]

## 5. Reduced coefficients

The stable BdG contribution is

\[
B_0=\frac{c_B^2}{\varpi^2},
\qquad
B_2=\frac{c_B^2}{\varpi^4},
\qquad
B_4=\frac{c_B^2}{\varpi^6}.
\]

The conservative Maxwell/mixed block uses

\[
\Delta=\Omega_U^2\Omega_W^2-R^2,
\]

\[
Q=g_U^2\Omega_W^2+2g_Ug_WR+g_W^2\Omega_U^2,
\]

\[
H=g_U^2+g_W^2,\qquad
S=\Omega_U^2+\Omega_W^2.
\]

Then

\[
Z_0=\frac Q\Delta,
\]

\[
Z_2=\frac{QS-H\Delta}{\Delta^2},
\]

\[
Z_4=\frac{Q(S^2-\Delta)-SH\Delta}{\Delta^3}.
\]

The outgoing-transfer moments are

\[
N_0=\frac{(\Omega_U^2g_W+Rg_U)^2}{\Delta^2},
\]

\[
N_2=
\frac{2P(P S-\Delta g_W)}{\Delta^3},
\qquad
P=\Omega_U^2g_W+Rg_U,
\]

\[
N_4=
\frac{\Delta^2g_W^2-2\Delta P^2-4\Delta PSg_W+3P^2S^2}{\Delta^4}.
\]

The total conservative operator is

\[
D_0=K-B_0-Z_0,
\]

\[
D_2=-(M+B_2+Z_2),
\]

\[
D_4=-(B_4+Z_4).
\]

The solved branch gives

\[
D_0=1.89448908712,\qquad
D_2=-1.09827027208,\qquad
D_4=-0.0282005090292.
\]

## 6. Observable packet

The normalized response and outgoing prefactor are

\[
u_2=-\frac{D_2}{D_0},
\]

\[
u_4=\frac{D_2^2-D_0D_4}{D_0^2},
\]

\[
P_0=\frac{N_0}{D_0},
\]

\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]

\[
P_4=
\frac{D_0^2N_4-2D_0(D_2N_2+D_4N_0)+3D_2^2N_0}{D_0^3}.
\]

For this branch,

\[
u_2=0.579718447339,\qquad
u_4=0.350959026601,
\]

\[
P_0=0.0197851167073,\qquad
P_2=0.0332589572512,\qquad
P_4=0.0365510452457.
\]

## 7. Target residuals

The residual packet is

\[
R_{\rm pole}=D_0(B_4+Z_4)-3(M+B_2+Z_2)^2,
\]

\[
R_{\rm norm}
=
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
-
\frac{54Gc_s^5}{5a^5c^5},
\]

\[
R_{P2}=P_2,\qquad R_{P4}=P_4,
\]

\[
R_{\rm tail}=\Theta_{\rm tail}\left(\frac c{c_s}\right)^3-1.
\]

The branch output is:

```text
R_pole = -3.56516721502
R_norm = -10.7802148833
R_P2   = 0.0332589572512
R_P4   = 0.0365510452457
R_tail = 0
```

The target packet therefore fails:

```text
target_packet_pass = False
```

The failure is expected. This first reduced branch is open and stable, but it was
not tuned to the isotropic full-bundle target surface.

The normalization ratio is

\[
\frac{P_0}{P_0^{\rm target}}=0.00183195525067.
\]

The one-pole ratio is

\[
\frac{D_0(B_4+Z_4)}{3(M+B_2+Z_2)^2}=0.0147641804366.
\]

So this branch undershoots the outgoing normalization and misses the conservative
one-pole surface.

## 8. Verdict

```text
open_gate_pass: True
stability_gate_pass: True
outgoing_transfer_gate_pass: True
target_packet_pass: False
```

This is the desired first-real-branch behavior:

- the open-throat geometry is valid,
- the reduced conservative branch is stable,
- the mixed/outgoing transfer is non-dark,
- but the branch fails the target residuals honestly.

That means the V2-22C pipeline is now usable for real reduced branch exports:
it does not merely process mock data, and it does not rescue a failing branch.

## 9. Next handoff to Codex

The next Codex/local continuation should replace the toy 1D operators with a
true moving-throat export while preserving this schema:

1. solve the stationary open-throat branch,
2. export `R0(s)`, wall coefficients, BdG profiles, gauge/mixed profiles,
3. compute `B_n,Z_n,N_n`,
4. run the same residual packet,
5. do not change the branch after reading residuals.

The most important data to improve are:

- `P0/P0_target`, currently far too small,
- `D0(B4+Z4)/(3A^2)`, currently far from the one-pole surface,
- the constant-prefactor conditions `P2=0`, `P4=0`.
