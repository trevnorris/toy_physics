# pathA_23 Stage 3b over-count and curvature localization

## VERDICT

`OVER_COUNT_CONFIRMED_CURVATURE_LOCALIZED`

CONDITIONAL flag: `CONDITIONAL`. The constitutive law remains the user-postulated
`ROTATIONAL_POSTULATED` MacCullagh/rotational brane law.

One-line summary: varying the full action with respect to the bulk fields gives
\(\delta S/\delta q_{\rm bulk}=\delta S_\psi^{\rm keep}/\delta q_{\rm bulk}
+\delta S_{\rm cpl}/\delta q_{\rm bulk}\); the brane stress
\(\sigma^R_{ab}\) is an internal brane restoring stress, not a bulk source.
Flat, density-preserving in-plane light free-slips, and the remaining leak
channels are defect/curvature localized.

Subtokens:

- `SIGMA_R_NOT_A_BULK_SOURCE`
- `LIGHT_FREE_SLIPS_NO_LEAK`

The adversarial defect-independent reading of Stage 3 is retired. Stage 3's
verdict token remains `LEAK_BOUNDED_CONDITIONAL`, but the bounded leak is now
classified as defect/curvature/throat data, not an intrinsic leak for every
photon.

## 1. Bulk EOM variation

Let \(q_\alpha\in\{\rho,\theta,v_x,v_y,v_z,v_w\}\). In the regulated
finite-thickness representation,
\[
\mathcal L_{\rm tot}^{(4)}
=\mathcal L_\psi(\rho,\theta,v_i)
+B_\ell\left[
\mathcal L_{\rm brane}(u_a,u_w)
+u_w\Pi_n(\rho,\theta,v_i;\Sigma)
+\alpha_uJ^a u_a-J^0\Phi_{\rm br}
\right].
\]
The rotational energy \(\mathcal U_R=\frac12\mu_R r_ar^a\) depends on
\(D_au_b\), not on \(\rho,\theta,v_i\). Therefore
\[
\frac{\delta S}{\delta q_\alpha}
=\frac{\delta S_\psi^{\rm keep}}{\delta q_\alpha}
+\int dw\,B_\ell u_w
\frac{\delta\Pi_n}{\delta q_\alpha}
\]
up to any later Stage-6 dependence of \(J^\mu_{\rm br}\) or \(\Phi_{\rm br}\)
on bulk/throat fields. There is no term proportional to
\(D_b\sigma^R_{ba}\), \(\mu_R\), or \(\mathcal U_R\) in the bulk-field
variation.

Equivalently, the hydrodynamic bulk equations have the schematic form
\[
\partial_t\rho+\partial_i(\rho v_i)=C_\theta[B_\ell u_w\Pi_n],
\]
\[
\partial_t(m\rho v_i)+\partial_jT^\psi_{ji}
+\rho\,\partial_iV_{\rm conf}
=f_i^{\rm cpl}[B_\ell u_w\Pi_n],
\]
\[
\partial_t\omega_{\rm bulk}
=\cdots+\nabla\times\left(f^{\rm cpl}/(m\rho)\right),
\]
where \(T^\psi_{ij}=m\rho v_iv_j+\delta_{ij}P+\sigma^Q_{ij}\). The brane
Euler equation does contain \(D_b\sigma^R_{ba}\); the bulk Euler/vorticity
equations do not.

Able-to-fail controls were included. If a forbidden term
\(B_\ell\beta v_aD_b\sigma^R_{ba}\) or \(B_\ell\eta\rho\mathcal U_R\) is
added, the bulk velocity or density variation immediately contains \(\mu_R\)
and \(D_b\sigma^R_{ba}\). The declared Stage-0 action does not contain those
couplings.

## 2. Tangential traction

For a flat surface with \(n=\hat w\),
\[
t_a=n_iT^\psi_{ia}=T^\psi_{wa}
=m\rho v_wv_a+\underbrace{P\delta_{wa}}_{0}+\sigma^Q_{wa}.
\]
Thus the ideal pressure is purely normal. The ideal-fluid tangential traction
is convective only, \(m\rho(v\cdot n)v_a=m\rho v_wv_a\).

The quantum stress was kept explicitly:
\[
\sigma^Q_{wa}
=\frac{\hbar^2}{4m}
\left(\frac{\partial_w\rho\,\partial_a\rho}{\rho}
-\partial_w\partial_a\rho\right).
\]
For \(\rho=\rho_0(w)+\delta\rho\,e^{ik\cdot x}\),
\[
\delta\sigma^Q_{wa}
=\frac{\hbar^2}{4m}\,ik_a
\left(\frac{\rho_0'}{\rho_0}\delta\rho-\partial_w\delta\rho\right).
\]
This is nonzero for a generic density perturbation. It vanishes for the
density-preserving transverse in-plane shear used for light:
\[
\delta\rho_{\rm br}=-i\rho_3 k_au^a=0,\qquad k\cdot u=0.
\]

## 3. Flat-brane light

For in-plane transverse light on a flat brane,
\[
u_w=0,\qquad s_a=D_au_w=0,\qquad k\cdot u=0.
\]
The kinematic normal relation is
\[
v_n=v_w-v^as_a=\dot u_w,
\]
so flat in-plane light gives
\[
v_w=0.
\]
Then
\[
T_{wa}^{\rm light}
=m\rho v_wv_a+\sigma^Q_{wa}=0.
\]
The flat-light token is therefore `LIGHT_FREE_SLIPS_NO_LEAK`. The scripts also
checked the failure controls: an independent \(v_w\ne0\) gives a convective
flat leak, and a compressive density perturbation gives a quantum-stress leak.
Those are not the transverse density-preserving light branch.

## 4. Curvature localization and scaling

Stage 1's projected stress relation remains the correct defect/geometry
channel:
\[
T_{na}=T_{wa}+(T_{ww}\delta_{ab}-T_{ab})s_b+O(s^2),
\qquad s_a=D_au_w.
\]
The Stage 3b refinement is that this is a bulk stress/projection channel, not
the brane's own \(\sigma^R\) force.

A constant slope is not an invariant leak: in a local tangent frame, pressure
has \(n\cdot e_a=0\), and tangent ideal flow has \(n\cdot v=0\), so its
convective tangential traction also vanishes. Physical localization begins
with slope variation across the interaction region,
\[
\Delta s_a\simeq K_{ab}\Delta\xi^b,\qquad K_{ab}=D_aD_bu_w.
\]
For a curvature/support scale \(L_{\rm mix}\),
\[
\Delta s\sim |K|L_{\rm mix}.
\]
The surviving transverse source has the leading force-amplitude estimate
\[
S_T^{\rm defect/curv}
\sim P_T\left[
\delta T_w^{\rm defect}
+A_{\rm mix}\,\Delta s
+\delta\sigma^Q_w{}^{\rm defect}
\right],
\qquad
A_{\rm mix}\sim T_{ww}I-T_\parallel.
\]
Therefore, in the Stage-3 force-ratio convention,
\[
\epsilon_{\rm leak}^{\rm force}
\sim
\frac{\|P_TA_{\rm mix}\|}{F_{\rm ref}}\,|K|L_{\rm mix}
+\frac{\|P_T\delta T_w^{\rm defect}\|}{F_{\rm ref}}
+\frac{\|P_T\delta\sigma^Q_w{}^{\rm defect}\|}{F_{\rm ref}}.
\]
If one instead reports leaked power/probability, the geometric mixing factor
is squared, \((|K|L_{\rm mix})^2\). With Stage 3's force-density definition,
the leading curvature amplitude is linear in \(|K|L_{\rm mix}\). In either
case it vanishes in the flat far field and turns on at bends/throats.

## Channel split and D1 meaning

| channel | classification after Stage 3b | survives as bulk transverse leak? |
| --- | --- | --- |
| \(D_b\sigma^R_{ba}\) | intrinsic brane restoring force | No. Brane EOM only; bulk-source inclusion was the over-count. |
| ideal pressure \(P\delta_{ia}\) | ideal bulk normal pressure | No tangential traction for \(n=\hat w\) or a local tangent surface. |
| convective \(m\rho v_wv_a\) | bulk normal-flow channel | No for flat light; yes for defect/throat normal flow. |
| \(\sigma^Q_{wa}\) | density-gradient quantum stress | No for density-preserving flat shear; yes for generic mixed density gradients. |
| \((T_{ww}\delta_{ab}-T_{ab})D_bu_w\) | geometry/projection channel | Vanishes flat; survives where slope/curvature and anisotropic/defect stress exist. |
| \(t_A^a\) mouth traction | defect boundary data | Defect-localized; throat solve deferred. |
| \(\alpha_uJ_a\) | brane/defect current drive | Not a direct bulk-velocity source under Stage-0 independence; defect/current channel. |
| \(v_n\to T_{wa}\) feedback | normal-flow feedback | Not a flat intrinsic light leak; survives with defect/throat normal flow. |

D1/Magnus consequence: the free-space half passes this audit. There is no
defect-independent transverse bulk forcing attached to every photon. D1 is not
failed by Stage 3's \(D_b\sigma^R_{ba}\) term. The remaining risk is
localized at curvature/defect/throat data and remains bounded/deferred under
Stage 3's existing `LEAK_BOUNDED_CONDITIONAL` verdict.

## Script pointers

Primary Mathematica:

- `software/stage1_solver/tools/pathA_23_stage3b_overcount_and_curvature.wl`
- Output: `software/stage1_solver/_scratch/pathA_23_stage3b_overcount_and_curvature_mathematica.json`
- Run: `timeout 600 math -script software/stage1_solver/tools/pathA_23_stage3b_overcount_and_curvature.wl`
- Result: 30/30 checks, exit 0, token `OVER_COUNT_CONFIRMED_CURVATURE_LOCALIZED`.

Independent SymPy cross-check:

- `software/stage1_solver/tools/pathA_23_stage3b_overcount_and_curvature_sympy.py`
- Output: `software/stage1_solver/_scratch/pathA_23_stage3b_overcount_and_curvature_sympy.json`
- Run: `timeout 600 python3 software/stage1_solver/tools/pathA_23_stage3b_overcount_and_curvature_sympy.py`
- Result: 29/29 checks, exit 0, same token.

Load-bearing checks: full-action bulk variation equality with the brane
sector removed; negative controls that make \(\sigma^R\) source the bulk when
forbidden couplings are injected; flat tangential traction decomposition;
quantum-stress density-preserving versus compressive controls; flat-light
\(v_w=0,T_{wa}=0\); local tangent-frame free-slip; curvature \(K\)-localized
mixing.

Bookkeeping checks: dimensional homogeneity, JSON output, check counting, and
the explicit string export of symbolic expressions.

Stage 3b stops here. No Stage 4 spectrum, \(u_w\) confinement, Maxwell,
charge, cone, paper, or `decisions/*` edit is made.
