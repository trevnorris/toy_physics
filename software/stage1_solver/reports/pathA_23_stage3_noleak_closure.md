# pathA_23 Stage 3 constitutive no-leak closure

## VERDICT

`LEAK_BOUNDED_CONDITIONAL(epsilon_leak << 1 + transverse-cancellation/impedance price; otherwise FAIL_CONSTITUTIVE_TRACTION_LEAK)`

CONDITIONAL flag: `CONDITIONAL`. The constitutive law used here is the user-elected postulate
`ROTATIONAL_POSTULATED`, not a medium-derived shear law.

One-line summary: the postulated rotational package can be made locally angular-momentum consistent with a minimal gyrostatic
spin reservoir, but the resulting interface exchange has a generically nonzero transverse projector/curl once the MacCullagh
force, Stage-1 `T_na`, mouth traction, current drive, and `v_n -> T_wa -> bulk-vorticity` feedback are kept. A clean no-leak
closure is not earned.

Stage-3b clarification: `software/stage1_solver/reports/pathA_23_stage3b_overcount_and_curvature.md` re-varied the full action
with respect to the bulk fields and found `SIGMA_R_NOT_A_BULK_SOURCE`. The term \(D_b\sigma^R_{ba}\) below is a brane-internal
restoring force in the brane Euler equation, not a direct bulk Euler/vorticity source. The Stage-3 verdict token is unchanged,
but the fatal reading "every free photon drives bulk vorticity through \(D_b\sigma^R_{ba}\)" is an over-count. The surviving
bulk leak channels are the Stage-1 stress/projection, density-gradient, normal-flow, mouth, and current/throat channels.

If the named bound below is not imposed, the same algebra is the concept-fatal branch
`FAIL_CONSTITUTIVE_TRACTION_LEAK` = `FAIL_LEAK_BREAKS_MAGNUS`.

## Postulated Package

Adopted provenance label:

| ingredient | status |
| --- | --- |
| Rotational in-plane energy | `ROTATIONAL_POSTULATED` |
| Brane-bulk coupling | carried from Stage 0, `SYMMETRY_ALLOWED_POSTULATED` |
| Minimal spin closure | `GYROSTATIC_SPIN_POSTULATED_GAP`; not derived from substructure |
| Nonzero couple-stress sector | not added; would be an extra Cosserat/gyrostat ingredient |

The fixed in-plane energy is
\[
\mathcal U_R=\frac12\mu_R r_ar^a,\qquad
r_a=\epsilon_{abc}D^bu^c .
\]

The force stress derived from this energy is
\[
\sigma^R_{ab}=\frac{\partial\mathcal U_R}{\partial(D_a u_b)}
=\mu_R(D_a u_b-D_bu_a)=\mu_R\epsilon_{abc}r_c .
\]
It is antisymmetric. Without a spin/couple completion the local angular-momentum residual is nonzero:
\[
\sigma^R_{ab}-\sigma^R_{ba}\ne0 .
\]

The minimal closure used here is spin-only:
\[
M_{cab}=0,\qquad
\dot s_{ab}=-(\sigma^R_{ab}-\sigma^R_{ba}),
\]
so
\[
\sigma^R_{ab}-\sigma^R_{ba}+D_cM_{cab}+\dot s_{ab}=0 .
\]
Equivalently, in axial form, \(\dot s_c=-2\mu_R r_c\). This closes total angular momentum algebraically, but it is an
acknowledged gyrostatic reservoir postulate. If that reservoir is rejected, Stage 3 has the obstruction
`FAIL_ANTISYMMETRIC_STRESS_NO_SPIN_CLOSURE` before the leak audit.

## Interface Conditions

For an in-brane region \(M\),
\[
\delta U_R
=-\int_M d^3\xi\,D_b\sigma^R_{ba}\,\delta u_a
+\int_{\partial M}d^2S\,\nu_b\sigma^R_{ba}\,\delta u_a .
\]
Thus the rotational boundary traction is
\[
t_a^R=\nu_b\sigma^R_{ba},
\]
and the minimal couple traction is
\[
m_{\nu ab}=\nu_cM_{cab}=0 .
\]

With the Stage-0 coupling retained, Stage 3 originally assembled the following tangential residual in the sign convention of
the Stage-1 source audit:
\[
R_a^{(3)}
=D_b\sigma^R_{ba}
+T_{na}
+\alpha_u J_a
+\mathcal t_{A,a}
+\delta T_{wa}^{(v_n)} .
\]
Stage 3b splits this residual by variational ownership. The brane-internal piece is
\[
R_{a,\rm brane}^{(3)}=D_b\sigma^R_{ba}+\alpha_uJ_a+\cdots ,
\]
where \(D_b\sigma^R_{ba}\) appears in the \(u_a\) equation. It is not part of the bulk-field source obtained by varying
\((\rho,\theta,v_i)\). The bulk/projection side is instead
\[
S_{a,\rm bulk}^{(4)}=B_\ell(w)S_{a,\rm bulk}^{(3b)},
\]
\[
S_{a,\rm bulk}^{(3b)}
=T_{na}
+\mathcal t_{A,a}
+\delta T_{wa}^{(v_n)}
+\text{bulk-dependent pieces of }S_{\rm cpl},
\]
with no \(D_b\sigma^R_{ba}\) term unless an additional forbidden coupling such as \(v_aD_b\sigma^R_{ba}\) is added.
\(\mathcal t_{A,a}=t_{A,a}\delta_{\partial M_A}\) is the mouth traction converted to an effective brane-volume force density.
The Stage-1 stress channel is kept:
\[
T_{na}=T_{wa}+(T_{ww}\delta_{ab}-T_{ab})D_bu_w+O(s^2),
\]
\[
T_{wa}=m_{\rm GNLS}\rho v_wv_a+\sigma^Q_{wa}.
\]

The normal channel remains
\[
S_w^{(4)}=B_\ell\left(T_{nn}-\frac{\delta U_w}{\delta u_w}+f_A+\cdots\right),
\]
with the \(u_w\) stability/gap question deferred to Stage 4. Stage 3 uses it only through the already-flagged
\(v_n=v^w-v^aD_au_w\) feedback channel.

Units are restored:
\[
[B_\ell]=L^{-1},\quad [\mu_R]=[\sigma^R]=E/L^3,
\]
\[
[D_b\sigma^R_{ba}]=[T_{na}]=[\alpha_uJ_a]=[\mathcal t_A]=E/L^4,
\]
\[
[B_\ell R_a^{(3)}]=[B_\ell S_{a,\rm bulk}^{(3b)}]=E/L^5=M L^{-3}T^{-2}.
\]
The raw boundary/mouth traction \(t_A^a\) has units \(E/L^3\); the localized \(\mathcal t_A^a\) has one additional inverse
length from the mouth support delta. No equality among \(B_\ell\), \(Z(w)\), and \(W(w)\) is used.

## Transverse Source

The rotational algebra in this subsection remains correct for the brane Euler residual. Stage 3b changes its interpretation:
the nonzero transverse projector of \(D_b\sigma^R_{ba}\) is not, by itself, a bulk source.

For a Fourier mode,
\[
D_b\sigma^R_{ba}
=\mu_R\left[k_a(k\cdot u)-k^2u_a\right].
\]
This vanishes for a longitudinal displacement \(u_a=k_a\phi\), but for a transverse brane mode \(k\cdot u=0\) it is
\[
D_b\sigma^R_{ba}=-\mu_R k^2u_a ,
\]
so its transverse projector is nonzero.

Stage 3b result:
\[
P_TD_b\sigma^R_{ba}\ne0\quad\text{in the brane equation, but}\quad
\frac{\delta S}{\delta(\rho,\theta,v_i)}\not\supset D_b\sigma^R_{ba}.
\]
Thus this term is removed from the intrinsic-to-light bulk-vorticity source.

Using \(P_T^{ab}=k^2\delta^{ab}-k^ak^b\), the aggregate residual tested by the original Stage-3 scripts was
\[
P_T^{ab}R_b^{(3)}
=P_T^{ab}\Big[
\mu_R(k_b(k\cdot u)-k^2u_b)
+T_{wb}
+iu_w(T_{ww}\delta_{bc}-T_{bc})k_c
+\alpha_uJ_b
+\mathcal t_{A,b}
+m_{\rm GNLS}\rho_3\,\delta v_w\,V_b
\Big].
\]
The scripts verify two sides of the test:

- The special branch \(u_a=k_a\phi\), \(P_TT_w=0\), isotropic \(T_{ab}\), no transverse mouth/current force, and no transverse
  \(v_n\) feedback gives zero transverse source.
- A transverse rotational brane mode, generic `T_wa`, or generic \(v_n\) feedback gives a nonzero transverse aggregate residual and
  nonzero curl.

After Stage 3b, the intrinsic-to-light bulk source excludes the first term. The flat-light bulk source has
\[
v_w=0,\qquad T_{wa}^{\rm light}=m_{\rm GNLS}\rho v_wv_a+\sigma^Q_{wa}=0
\]
for density-preserving transverse shear \(k\cdot u=0\). The surviving bulk channels are defect/curvature channels:
\[
P_T^{ab}S_{b,\rm bulk}^{(3b)}
=P_T^{ab}\Big[
T_{wb}
+iu_w(T_{ww}\delta_{bc}-T_{bc})k_c
+\mathcal t_{A,b}
+\delta T_{wb}^{(v_n)}
+\text{bulk-dependent coupling pieces}
\Big],
\]
with the \(T_{na}\) slope/projection term scaling as \(s_a=D_au_w\) in a fixed chart and as
\(\Delta s_a\simeq K_{ab}\Delta\xi^b\) in a local tangent frame.

The \(v_n\) priority channel is not a direct constant-coefficient leak:
\[
S_a^{(v_n,{\rm direct})}=-i\Lambda_n u_w k_a,\qquad
P_TS^{(v_n,{\rm direct})}=0,\qquad k\times S^{(v_n,{\rm direct})}=0 .
\]
The feedback leak appears one step later through the scalar stress:
\[
\delta T_{wa}^{(v_n)}=m_{\rm GNLS}\rho_3\,\delta v_w\,V_a+\delta\sigma^Q_{wa}.
\]
For the convective part,
\[
P_T\delta T_w=m_{\rm GNLS}\rho_3\delta v_w\,P_TV,\qquad
k\times\delta T_w=m_{\rm GNLS}\rho_3\delta v_w\,(k\times V),
\]
which is nonzero unless the in-plane background flow is parallel to the perturbation wavevector or \(\delta v_w=0\). This is
the channel that feeds the bulk vorticity equation:
\[
\partial_t\omega_{\rm bulk}\supset \nabla\times\left(S^{(4)}/\rho_4\right).
\]

Therefore no-leak is not structurally closed for defect/throat data. The free flat-light intrinsic channel is closed by Stage 3b,
but the curvature/projection, normal-flow, quantum-density-gradient, mouth, and current/throat channels still require the named
bounded/deferred treatment.

## Order Comparison

Original Stage-3 named small parameter for the aggregate residual:
\[
\epsilon_{\rm leak}^{\rm agg}(k)
=\frac{\|P_TR^{(3)}\|}
{\max(F_{\rm Mag},F_g)} .
\]
After the Stage-3b split, the bulk-leak version replaces \(R^{(3)}\) by \(S_{\rm bulk}^{(3b)}\).

Use the brane-level reference force densities
\[
F_{\rm Mag}\sim\rho_3\Gamma n_vV_{\rm rel},\qquad
F_g\sim\rho_3\frac{v_r^2}{L_g}.
\]
Both have units \(E/L^4=M L^{-2}T^{-2}\). The finite-thickness comparison uses
\(B_\ell F\sim F/\ell\) on both numerator and denominator, so the common layer factor cancels only if the same support is
being compared.

The leading pieces in the numerator are
\[
F_R\sim\mu_R k^2U_T\qquad\text{(brane-internal after Stage 3b; not a direct bulk source)},
\]
\[
F_{\rm Stage1}\sim
\left\|P_T\left[T_w+iu_w(T_{ww}I-T_\parallel)k\right]\right\|,
\]
\[
F_{v_n}\sim m_{\rm GNLS}\rho_3|\delta v_w|\,|V_\parallel|\sin\theta_{kV},
\qquad
F_A\sim\|P_T\mathcal t_A\|,
\qquad
F_J\sim\|\alpha_uP_TJ\|.
\]

Stage 3b removes \(F_R\) from the intrinsic-to-light bulk-leak numerator. For flat, density-preserving in-plane light,
\[
F_R\ne0\ \text{in the brane equation},\qquad
F_{\rm bulk}^{\rm flat\,light}=0 .
\]
The remaining bulk numerator is
\[
F_{\rm bulk}^{(3b)}
\sim F_{\rm Stage1}+F_{v_n}+F_A+F_J+\|\delta\sigma^Q_w{}^{\rm defect}\| ,
\]
where the geometric projection part obeys
\[
F_{\rm Stage1}^{\rm geom}
\sim \|P_TA_{\rm mix}\|\,|K|L_{\rm mix},
\qquad A_{\rm mix}\sim T_{ww}I-T_\parallel ,
\]
in a local tangent frame. If the leakage measure is an energy/power fraction rather than the Stage-3 force-density ratio, the
geometric mixing factor is squared, \((|K|L_{\rm mix})^2\).

Bounded survival condition:
\[
\epsilon_{\rm leak}\ll1
\]
for the realized throat/mouth/curvature branch and for the brane modes being allowed. After Stage 3b this is not an arbitrary
free-space impedance imposed on every photon; it is a defect/curvature bound to be supplied by the throat/profile solve. If
\(\epsilon_{\rm leak}=O(1)\) for that localized branch, the verdict becomes
`FAIL_CONSTITUTIVE_TRACTION_LEAK`.

## Stage 3b Channel Split

Cross-reference: `software/stage1_solver/reports/pathA_23_stage3b_overcount_and_curvature.md`.

| term/channel | Stage-3b classification | bulk transverse status |
| --- | --- | --- |
| \(D_b\sigma^R_{ba}\) | intrinsic brane restoring force | Not a bulk source. This was the Stage-3 over-count. |
| flat ideal pressure \(P\delta_{ia}\) | ideal bulk normal pressure | No tangential traction. |
| \(m_{\rm GNLS}\rho v_wv_a\) | convective normal-flow channel | Zero for flat in-plane light with \(v_w=0\); survives at defect/throat normal flow. |
| \(\sigma^Q_{wa}\) | quantum density-gradient channel | Zero for density-preserving transverse shear; survives for generic mixed density gradients. |
| \((T_{ww}\delta_{ab}-T_{ab})D_bu_w\) | geometric projection channel | Vanishes flat; turns on with slope/curvature and anisotropic defect stress. |
| \(\mathcal t_{A,a}\) | mouth boundary channel | Defect-localized and deferred to throat data. |
| \(\alpha_uJ_a\) | brane/defect current drive | Not a direct bulk-velocity source under Stage-0 independence; current/throat channel. |
| \(\delta T_{wa}^{(v_n)}\) | normal-flow feedback | Not a flat intrinsic light leak; survives with defect/throat normal flow. |

Stage 3b tokens: `SIGMA_R_NOT_A_BULK_SOURCE`, `LIGHT_FREE_SLIPS_NO_LEAK`, and
`OVER_COUNT_CONFIRMED_CURVATURE_LOCALIZED`. Thus D1 is not failed by an intrinsic free-photon leak. D1 remains exposed only to
the localized defect/curvature channels above.

## C2 Decomposition

The static source coupling remains nonzero:
\[
S_J^{\rm static}=-\int dt\,d^3\xi\,J^0_{\rm br}\Phi_{\rm br}.
\]
Under the Stage-0 independence premise, the direct variation of this static term with respect to
\((\rho,\theta,v_i)\) is zero, and it supplies no direct transverse bulk force. Thus the longitudinal/static Coulomb sector can
still be sourced without a direct static transverse-shear source.

That C2 result is limited. After Stage 3b, \(D_b\sigma^R_{ba}\) is not a direct bulk source, but the full back-reaction still has
defect/curvature transverse channels through `T_wa`, anisotropic `T_ab`, mixed density-gradient `sigma^Q_wa`, mouth traction,
\(J^a\), and \(v_n\) feedback. If Stage 6 later derives \(J^\mu_{\rm br}\) or \(\Phi_{\rm br}\) as bulk/throat functionals, the
direct-zero half of this C2 decomposition must be revisited.

## Decided vs Deferred

Decided:

- The postulated rotational stress is antisymmetric and needs the postulated gyrostatic spin reservoir to close angular
  momentum.
- With that closure, the aggregate brane-side residual is generically transverse/curl-bearing.
- Stage 3b split: \(D_b\sigma^R_{ba}\) is not a bulk source, and flat density-preserving light has \(v_w=0\), \(T_{wa}=0\).
- `NO_LEAK_CLOSED` is not earned for defect/curvature/throat data.
- The bounded path is the explicit \(\epsilon_{\rm leak}\ll1\) condition for the localized curvature/defect branch.
- The direct static \(J^0\Phi_{\rm br}\) coupling can coexist with no direct static transverse bulk source under the Stage-0
  independence premise.

Deferred:

- The actual throat/profile values of `T_wa`, anisotropic `T_ab`, \(t_A^a\), and \(\delta v_w\).
- Any derivation of the gyrostatic spin reservoir or a nonzero couple-stress sector.
- Stage 4 spectrum/`u_w`, Stage 5 Maxwell/gauge, and Stage 6 charge/cone tests.

## Script Pointers

Primary Mathematica:

- `software/stage1_solver/tools/pathA_23_stage3_noleak_closure.wl`
- Output: `software/stage1_solver/_scratch/pathA_23_stage3_noleak_closure_mathematica.json`
- Run: `timeout 600 math -script software/stage1_solver/tools/pathA_23_stage3_noleak_closure.wl`
- Result: 36/36 checks, exit 0, token
  `LEAK_BOUNDED_CONDITIONAL(epsilon_leak<<1 + transverse-cancellation/impedance price; otherwise FAIL_CONSTITUTIVE_TRACTION_LEAK)`.

Independent SymPy cross-check:

- `software/stage1_solver/tools/pathA_23_stage3_noleak_closure_sympy.py`
- Output: `software/stage1_solver/_scratch/pathA_23_stage3_noleak_closure_sympy.json`
- Run: `timeout 600 python3 software/stage1_solver/tools/pathA_23_stage3_noleak_closure_sympy.py`
- Result: 36/36 checks, exit 0, same token.

Load-bearing checks: derived \(\sigma^R\) from \(\mathcal U_R\); showed angular closure fails without spin and closes with the
minimal spin-rate postulate; computed \(D_b\sigma^R_{ba}\), the Stage-1 \(T_{na}\) projector/curl, mouth/current slots, and the
\(v_n\) feedback projector/curl; verified a special no-leak control and multiple nonzero leak controls. Dimensional checks are
bookkeeping except for the equal-dimension \(\epsilon_{\rm leak}\) comparison.
Stage 3b supersedes only the direct-bulk-source interpretation of the \(D_b\sigma^R_{ba}\) part of that aggregate residual.

Stage 3 stops here.
