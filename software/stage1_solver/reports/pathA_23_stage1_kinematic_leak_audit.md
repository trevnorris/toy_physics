# pathA_23 Stage 1 kinematic leak audit

## VERDICT

`LEAK_CONDITIONS_DEFERRED` - deriving the scalar-stress projection removes the old no-leak claim. The direct scalar/normal variations of `S_cpl` are zero or longitudinal at linear order, but the actual interface traction
\[
S_a^{\rm stress}=B_\ell T_{na}
\]
has a generically transverse part unless the selected scalar stress profile satisfies named no-leak conditions. The mouth tangential traction \(t_A^a\) is also a Stage-1-reachable transverse channel, deferred because the throat/mouth boundary solve is not available in Stage 1.

Stage 1 stops here. This is not a Stage-2 constitutive result and not a Stage-3 closure result.

## Derived Coupling Variation

The audited coupling is
\[
S_{\rm cpl}
=\int dt\,d^3\xi\,\sqrt\gamma\,
\left[
u_w\Pi_n+J_{\rm br}^aA_a^{\rm br}-J_{\rm br}^0\Phi_{\rm br}
\right]+S_{\rm mouth},
\qquad A_a^{\rm br}=\alpha_u u_a .
\]
The finite-thickness bulk representation keeps \(B_\ell(w)\):
\[
\mathcal L_{\rm cpl}^{(4)}=
B_\ell\left[u_wT_{nn}+J_{\rm br}^a\alpha_u u_a-J_{\rm br}^0\Phi_{\rm br}\right].
\]

For the kept scalar sector, the spatial momentum-flux/stress representative from the Noether balance is
\[
T_{ij}^{\psi}
=m_{\rm GNLS}\rho v_i v_j+\delta_{ij}P(\rho)+\sigma^Q_{ij},
\qquad P(\rho)=K\rho^5,
\]
\[
\sigma^Q_{ij}
=\frac{\hbar^2}{4m_{\rm GNLS}}
\left(\frac{\partial_i\rho\,\partial_j\rho}{\rho}-\partial_i\partial_j\rho\right).
\]
The explicit confinement term contributes a body force \(-\rho\partial_iV_{\rm conf}\), not a hidden tangential stress.

With small brane slope \(s_a=D_a u_w\),
\[
n_i=(-s_a,1)+O(s^2),
\qquad e_a^i=\delta_a^i+s_a\delta_w^i+O(s^2).
\]
The scripts derive, rather than assume,
\[
T_{nn}=T_{ww}-2s_bT_{wb}+O(s^2),
\]
\[
T_{na}=T_{wa}+(T_{ww}\delta_{ab}-T_{ab})s_b+O(s^2).
\]
For a Fourier bending mode \(s_a=i k_a u_w\),
\[
S_a^{\rm stress}
=B_\ell\left[T_{wa}+iu_w(T_{ww}\delta_{ab}-T_{ab})k_b\right]+O(s^2).
\]

The scalar variations of \(B_\ell u_wT_{nn}\) give scalar sources. If represented in the Euler equation, their in-plane force is a gradient:
\[
S_a^{(q)}=iB_\ell u_w\frac{\partial T_{nn}}{\partial q}k_a ,
\qquad q\in\{\rho,\theta,\ldots\}.
\]
The normal scalar work similarly gives
\[
S_a^{(nn)}=iB_\ell T_{nn}u_w k_a .
\]
For a possible scalar normal-flow dependence \(v_n=n_iv_i=v_w-v^aD_au_w+O(s^2)\),
\[
\delta v_n=\delta v_w-s_a\delta v_a-v_aD_a\delta u_w+O(s^2),
\]
so, with \(\Lambda_n=\partial\mathcal L_{\rm cpl}/\partial v_n\),
\[
S_a^{(v_n)}=-iB_\ell\Lambda_n u_w k_a
\]
for constant \(\Lambda_n\) at this linearized level.

The source-current term has no direct Stage-1 variation with respect to the bulk scalar variables:
\[
\frac{\delta}{\delta(\rho,\theta,v_i)}
B_\ell(J_{\rm br}^a\alpha_u u_a-J_{\rm br}^0\Phi_{\rm br})=0,
\]
under the Stage-1 premise that \(J_{\rm br}^\mu\) and \(\Phi_{\rm br}\) are brane/throat variables not yet derived as functionals of \(\rho,\theta,v_i\). Its brane variation is nonzero:
\[
\frac{\delta\mathcal L_{\rm cpl}^{(4)}}{\delta u_a}=B_\ell\alpha_uJ_{\rm br}^a .
\]

## Longitudinal/Transverse Split

Use the numerator projector
\[
P_T^{ab}=k^2\delta^{ab}-k^ak^b .
\]
The scalar/normal pieces above are longitudinal:
\[
P_T^{ab}S_b^{(q)}=P_T^{ab}S_b^{(nn)}=P_T^{ab}S_b^{(v_n)}=0,
\qquad
\epsilon^{abc}k_bS_c=0 .
\]

The derived stress-traction channel is different:
\[
P_T^{ab}S_b^{\rm stress}
=B_\ell\left[P_T^{ab}T_{wb}
-iu_wP_T^{ab}T_{bc}k_c\right],
\]
\[
\epsilon^{abc}k_bS_c^{\rm stress}
=B_\ell\left[(k\times T_w)^a-iu_w(k\times(T_\parallel k))^a\right].
\]
Therefore no-leak is not automatic. It holds at this stage only if
\[
P_T^{ab}T_{wb}=0,\qquad P_T^{ab}T_{bc}k_c=0,
\]
and the equivalent curl conditions hold on the selected scalar/throat branch.

One sufficient special case is \(T_{wa}=0\) and \(T_{ab}=T_\parallel\delta_{ab}\), where
\[
T_{na}=iu_w(T_{ww}-T_\parallel)k_a
\]
is longitudinal. But the scalar stress actually contains
\[
T_{wa}=m_{\rm GNLS}\rho v_wv_a+\sigma^Q_{wa},
\qquad
T_{ab}=m_{\rm GNLS}\rho v_av_b+\delta_{ab}P+\sigma^Q_{ab},
\]
so normal-tangent flow, anisotropic convective stress, or quantum-stress anisotropy can carry transverse/curl content. Stage 1 does not know the selected throat profile well enough to set these terms to zero or to compute a definite nonzero magnitude.

## Pi_n/T_na Result

The old scalar-\(\Pi_n\) premise is removed. \(\Pi_n\) is
\[
\Pi_n=T_{nn},
\]
and the associated in-plane interface reaction is
\[
T_{na}=T_{wa}+(T_{ww}\delta_{ab}-T_{ab})D_bu_w+O(s^2).
\]
The projector/curl of this derived source is generically nonzero. The honest Stage-1 status is:

No-leak iff the selected scalar branch satisfies \(P_TT_w=0\), \(P_T(T_\parallel k)=0\), and the corresponding curl conditions. Otherwise the reachable token is `FAIL_COUPLING_SOURCES_TRANSVERSE`. The actual branch data needed to decide between those cases is deferred to the throat/profile and Stage-3 closure solve.

## S_mouth

Stage 0 includes mouth work
\[
\delta W_A=t_A^a\delta u_a+f_A\delta u_w+m_A^{ab}\delta K_{ab}.
\]
The tangential mouth term contributes an in-plane source localized at each mouth:
\[
S_{a,A}^{\rm mouth}=t_{A,a}.
\]
Its transverse and curl tests are
\[
P_T^{ab}t_{A,b},\qquad \epsilon^{abc}k_bt_{A,c}.
\]
These are not computable from the Stage-0 action alone because \(t_A^a\) depends on the throat boundary condition, defect-current matching, anchoring-vs-traction choice, and later constitutive closure. It is deferred explicitly as a Stage-1-reachable leak channel, not dropped. No-leak additionally requires \(P_Tt_A=0\) and \(k\times t_A=0\) at each mouth.

## C2 Decomposition

The static source coupling is nonzero:
\[
S_J=-\int dt\,d^3\xi\,\sqrt\gamma\,J_{\rm br}^0\Phi_{\rm br}.
\]
So the C2 check is not the trivial zero-coupling limit. Under the Stage-1 independence premise,
\[
\delta_{\rho,\theta,v_i}S_J=0,
\]
and the direct bulk transverse/curl source from the Coulomb/static term is zero. This is a premise-limited bulk-zero result: if Stage 6a later derives \(J_{\rm br}^\mu\) or \(\Phi_{\rm br}\) as bulk/throat functionals, the zero-variation half must be revisited.

## u_w Channel And Scope

At linear order in brane slope, the normal scalar pieces from \(u_wT_{nn}\) and scalar \(v_n\) dependence are longitudinal for constant coefficients. The derived traction \(T_{na}\), however, carries the bending-mediated anisotropic-stress channel above.

Deferred scope:

- \(O(s^2)\) bending and curvature corrections from the normalized normal, tangent vectors, and \(K_{ab}=D_aD_bu_w\).
- \(u_w\,\partial T_{nn}/\partial v_a\), which is \(O(u_wD_au_w)\) for the convective stress and absent only in the strictly linear brane audit.
- \(v_n\to\) bulk-vorticity feedback through the Euler/vorticity equations. The direct \(v_n\) scalar source is longitudinal at linear order, but the evolution feedback into the Magnus reservoir is not settled here.
- Actual scalar stress components \(T_{wa}\), anisotropic \(T_{ab}\), and mouth traction \(t_A^a\).

## Script Pointers

Primary Mathematica:

- `software/stage1_solver/tools/pathA_23_stage1_kinematic_leak_audit.wl`
- Output: `software/stage1_solver/_scratch/pathA_23_stage1_kinematic_leak_audit_mathematica.json`
- Run: `timeout 600 math -script software/stage1_solver/tools/pathA_23_stage1_kinematic_leak_audit.wl`
- Result: 42/42 checks, exit 0, token `LEAK_CONDITIONS_DEFERRED`.

Independent SymPy cross-check:

- `software/stage1_solver/tools/pathA_23_stage1_kinematic_leak_audit_sympy.py`
- Output: `software/stage1_solver/_scratch/pathA_23_stage1_kinematic_leak_audit_sympy.json`
- Run: `timeout 600 python3 software/stage1_solver/tools/pathA_23_stage1_kinematic_leak_audit_sympy.py`
- Result: 42/42 checks, exit 0, token `LEAK_CONDITIONS_DEFERRED`.

What the scripts verify:

- \(T_{nn}\) and \(T_{na}\) are derived from the projected scalar stress, not hardcoded.
- The Fourier factor \(i\), the finite-thickness factor \(B_\ell\), and units for \(B_\ell u_wT_{nn}\), \(B_\ell T_{na}\), \(v_n\), and the source-current term are kept.
- Scalar/normal variations are longitudinal or direct-bulk-zero at linear order.
- The derived \(T_{na}\) projector/curl is generically nonzero, while the isotropic/no-normal-tangent special case is longitudinal.
- Negative controls for arbitrary tangential source and forbidden direct \(J^a\to\) bulk-velocity coupling are detected.
- `FAIL_COUPLING_SOURCES_TRANSVERSE` is reachable when the derived stress or mouth channel is known transverse.

What the scripts do not verify:

- The actual throat/profile values of \(T_{wa}\), anisotropic \(T_{ab}\), or \(t_A^a\).
- Constitutive brane traction/couple-stress closure.
- \(O(s^2)\) bending/curvature corrections.
- \(v_n\to\) bulk-vorticity feedback and comparison to the Magnus/gravity-flow terms.
- Maxwell/gauge structure, physical current derivation, charge normalization, spectrum, cone lock, or paper integration.

Honesty note: 42/42 certifies the algebra, dimensions, projector/curl machinery, and controls. It does not certify no-leak. The no-leak verdict was removed because the derived interface traction contains unresolved transverse channels.
