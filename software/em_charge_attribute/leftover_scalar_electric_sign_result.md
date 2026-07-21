# Leftover-scalar electric-sign result

**Scope.** This is the target-blind analytic far-monopole result for the supplied G0 scalar closure, not the assembled R1/R3 throat solve. Both engines transcribe

\[
S_{Lh}=\int dt\,d^3x\left[\tfrac12A_{\rm eff}\dot u_L^2+\tfrac12M_h\dot h^2-\tfrac12B_{\rm eff}|\nabla u_L|^2-\tfrac12K_h|\nabla h|^2-C_{hu}\nabla u_L\!\cdot\!\nabla h\right]
\]

and \(\sum_i\int\eta_i(\tfrac12k_mh^2-g_{\chi h}s_i h)\). Here \(T_L=B_{\rm eff}=\rho_{B0}^2/\chi_c>0\), \(K_h=M_hc_E^2\), \(D=B_{\rm eff}K_h-C_{hu}^2>0\), and G0 has \((B_{\rm eff},K_h,C_{hu},D)=(2,1,1/2,7/4)\). Native conditions are continuity of \(u_L\) and \(B_{\rm eff}\partial_nu_L+C_{hu}\partial_nh\), continuity of \(h\) and its mixed conormal away from the mouth source, and \(u_L\to0\) in the IR. All 13 G0 zero-ledger rows used here retain their **[POSTULATE]** tags.

## Q1 — realizability

**Q1a: `NO_NATIVE_CLAMP`.** A common barrier evaluator tests whether each candidate contains the signed datum, has it as a stationary point, has positive holding curvature, and has protection. The actual reduced mouth energies have curvatures \(B_{\rm eff}>0\) or \(B_{\rm eff}-C_{hu}^2/K_h=7/4>0\): their unique minima relax, and a nonzero signed target is not stationary. The \(h\) source plus mixing shifts that one minimum but does not hold it. IR decay fixes zero, not a signed mouth value. Conditional \(q_L\propto s\) likewise shifts one convex minimum and is absent from G0. The parsed zero ledger removes wall–\(u_L\), drain–\(u_L\), and geon–\(u_L\) holding terms; `S_hold` fixes only \(r_B\). PathA_24 gives connected \(S^3\) and zero wall-orientation unwinding barrier. The mutation control supplies a double well plus the explicit odd selector \(-\mu s[3u_L/(2u_0)-u_L^3/(2u_0^3)]\); for \(0<\mu<2\lambda u_0^4/3\), both wells remain protected and the selected minimum is \(u_L=su_0\). The same evaluator then returns `NATIVE_CLAMP_EXISTS(injected_s_uL_hold)`.

**Q1b: `BOLT_ON_DEFERRED_TO_R1`.** No added action was supplied, so none is invented. It would need a protected/disconnected \(\pm w\) sector and a direct \(s\leftrightarrow u_L\) datum coupling, with an exact coefficient domain.

## Q2 — coupled far-monopole reaction

The profiles are explicit: for V, \(f_V=1\) on the spherical mouth \(A\) (unit area-average); for M/J, \(f_b=1/|A|\) on \(A\), zero off it, so \(\int_Af_b\,dS=1\). Thus M holds conormal \(sqf_b\), while J adds \(-sJ_0\int_Af_bu_L\,dS\).

Completing the gradient square with \(w=u_L+(C_{hu}/B_{\rm eff})h\) gives \(\kappa=D/B_{\rm eff}\) and the live Robin operator \(L_h=-\kappa\nabla^2+k_m\eta\). Define its exact vertices
\(z_g=1-k_m\langle\eta,L_h^{-1}\eta\rangle>0\) (the positive maximum-principle/Schur escape factor) and
\(z_b=1-k_m\langle\eta,L_h^{-1}f_b\rangle\). Then

\[
m=Z^T\!\begin{pmatrix}B_{\rm eff}^{-1}&0\\0&\kappa^{-1}\end{pmatrix}\!Z,
\quad Z=\begin{pmatrix}1&0\\-(C_{hu}/B_{\rm eff})z_b&z_g\end{pmatrix},
\quad \det m=\frac{z_g^2}{D}>0.
\]

Write the entries as \(m_{uu},m_{ug},m_{gg}\). The isolated one-body response defines \(u_A=S_{uu}Q+S_{ug}g\), with \(S_{uu}>0\). For V, \(Q_0=(u_0-S_{ug}g)/S_{uu}\) is the physical reaction and
\(v_{hh}=m_{gg}-2(S_{ug}/S_{uu})m_{ug}+(S_{ug}/S_{uu})^2m_{uu}>0\). Each row uses one normalization and \(U=s_1s_2A/(4\pi R)\), \(F_{\rm out}=s_1s_2A/(4\pi R^2)\) when \(A\ne0\).

| ensemble; conjugate functional | \(h\)-only \(A\) | \(u_L\)-only \(A\) | interference \(A\) | computed total sign |
|---|---:|---:|---:|---|
| **V**, \(u_L|_M=su_0\); fixed-potential \(E_0-gh-Qu\) | \(-v_{hh}g^2\) | \(+m_{uu}u_0^2/S_{uu}^2\) | \(0\) | negative, zero, or positive |
| **M**, \(\oint(B_{\rm eff}\partial_nu_L+C_{hu}\partial_nh)=sq\); stored energy with only the held-\(h\) source work removed | \(-m_{gg}g^2\) | \(+m_{uu}q^2\) | \(0\) | negative, zero, or positive |
| **J**, \(J=sJ_0\); \(E_0-Ju-gh\) | \(-m_{gg}g^2\) | \(-m_{uu}J_0^2\) | \(-2m_{ug}J_0g\) | strictly negative for a nonzero source vector |

The four readouts reconstruct every total exactly. Wrong-functional controls recompute the bare/nonconjugate member and change each V/M/J expression. The \(1/R\) Green behavior follows from the decaying harmonic branch in three dimensions, not from a target exponent.

## Q3 and §5 landing

Tier A passes: \(D>0\); the transverse block stays decoupled; the zero ledger is preserved; the charge datum excludes separate mass source \(q_M\); and the supplied ensemble makes no selection claim over U2’s 144/144 unresolved cells. Tier B remains deferred: assembled gravity/drain response, curved-sleeve stability, momentum/return closure, and the dressed nonperturbative two-body response.

Only here is the real-EM comparison made. V and M can yield the target sign, its opposite, or a null leading term; J yields the opposite sign. With no native holder, total precedence gives exactly:

`NO_NATIVE_CLAMP; ALGEBRAIC_SIGN={V:range(attract_1/R^2|null_leading|repel_1/R^2),M:range(attract_1/R^2|null_leading|repel_1/R^2),J:attract_1/R^2}; BOLT_ON_DEFERRED_TO_R1`

**DECIDED:** G0 has no native nonrelaxable signed \(u_L\) datum; the three conjugate ensemble structures, algebraic sign ranges, \(1/R^2\) force falloff, and Tier-A result are fixed above. **R1/R3 must settle:** whether an explicit added holder exists and survives consistency, which mouth ensemble and normalization a throat realizes, the full dressed monopole/nonlinear pair response, and every Tier-B item. Analysis does not earn a bolt-on or a physical electric sign.
