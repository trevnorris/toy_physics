# 5PN GR program — final wrap-up

## Bottom line

The symbolic narrowing phase is complete.

Within the current closure hierarchy, the full 5PN question has been reduced to two exact tests:

1. **Branch packet**
   
   \[
   \Delta_{\rm branch}
   =
   \bigl(
   a_2,\ b_2,\ a_4,\ b_4,\ a_{P_0},\ b_{P_0},\ \Delta_{\rm pole},\ \Delta_{\rm norm}
   \bigr),
   \]
   with
   \[
   \Delta_{\rm pole}=\bar u_4-4\bar u_2^{\,2},
   \qquad
   \Delta_{\rm norm}=\hat m_0^{\,2}\bar P_0-\frac{54Gc_s^5}{5a^5c^5}.
   \]

2. **Orbit packet**
   
   Any one of the equivalent forms
   \[
   (m_T,m_K,m_\mu),
   \qquad
   (R_{\rm tr},R_{\rm nt},R_\eta),
   \qquad
   (q_{\rm tr},q_{\rm nt},q_\eta),
   \]
   with exact orbit lock iff
   \[
   m_T=m_K=m_\mu=1
   \iff
   q_{\rm tr}=q_{\rm nt}=q_\eta=0.
   \]

Everything else in the long symbolic chain now compiles into those two packets.

## What we established

### 1. Conservative grouped-`P2` front end

The full conservative grouped real `P2` problem is reduced to the isotropy defects and the single-pole test:

\[
 a_2=b_2=0,
 \qquad
 a_4=b_4=0,
 \qquad
 u_4=4u_2^2.
\]

The exact isotropic full-bundle target is
\[
D_0(B_4+Z_4)=3(M+B_2+Z_2)^2,
\]
with
\[
P_0=\frac{N_0}{D_0},
\qquad
\hat m_0^{\,2}P_0=\frac{54Gc_s^5}{5a^5c^5}.
\]

### 2. Support/source side

The support/source branch was compressed all the way down through the explicit Family-1 construction. On the minimal isotropic contact-plus-pole conservative module,
\[
Y_Q^{\rm cons}(\omega)=\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2},
\]
the loading data are fixed to
\[
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43 C_{\rm mix},
\]
and the explicit Family-1 branch already clears that threshold at `Pe_req = 0`.
So the support/source side is no longer the active bottleneck.

### 3. Isotropic outgoing quadrupole branch

The exact outgoing `l=2` DtN fingerprint on the canonical compact branch is
\[
\widehat Y_2^{\rm out}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\frac{a^5\omega^5}{27c_s^5}+O(\omega^6),
\]
which fixes the canonical isotropic outgoing normalization to
\[
\chi_Q=1.
\]
The remaining reduced retarded gap is therefore not a generic outgoing ambiguity; it is whether the actual moving-throat branch realizes this canonical branch, or an explicitly deformed branch parameterized by the Stage-90 variables.

### 4. Mouth/core and compensation chain

Stages 91–150 showed that the isotropic outlet/core side has a nontrivial compensated surface, but once the lower compensated branch and its parent transport laws are imposed, arbitrary first-order isotropic bundle drift in
\[
(\Theta_w,K_s,K_q,P_0)
\]
is tangent to the exact compensated Family-1 parent family. So the first-order isotropic branch transport problem is closed.

### 5. Weak-axisymmetric / orbit-lock side

Stages 151–170 compressed the whole weak-axisymmetric grouped problem to preservation of three quotient coordinates
\[
(\mathfrak C_{{\rm tr},*},\ \mathfrak C_{{\rm nt},*},\ \epsilon_\eta),
\]
or equivalently tangency to a single exact similarity orbit `G_*`.

Stages 171–198 then sharpened that further:
- under the **adiabatic wall** condition `δln Θ_w = 0`, one failure channel is removed,
- under the stronger **adiabatic-elastic** condition
  \[
  \delta\ln\Theta_w=0,
  \qquad
  \varepsilon_L=\varepsilon_v=\varepsilon_T=0,
  \]
  the remaining first-order branch question is exactly orbit lock,
- and the exact finite orbit-lock verdict is reduced to the dependent residual triple
  \[
  (m_T,m_K,m_\mu).
  \]

The pairwise/reference-free version is now exact:
- two branch states are on the same orbit iff
  \[
  m_T=m_K=m_\mu=1,
  \]
- the equivalent invariant packet is
  \[
  (R_{\rm tr},R_{\rm nt},R_\eta),
  \]
- and the equivalent additive packet is
  \[
  (q_{\rm tr},q_{\rm nt},q_\eta),
  \]
  with orbit-distance
  \[
  D^2=q_{\rm tr}^2+q_{\rm nt}^2+q_\eta^2.
  \]

## Final theorem shape

Within the current hierarchy, the 5PN endgame is now:

1. Compute **Packet A** from the actual moving-throat branch:
   \[
   (D_{A0},D_{A2},D_{A4},N_{A0},N_{A2},N_{A4}),
   \qquad A\in\{20,21,22\},
   \]
   together with `mhat_0`.

2. Compute **Packet B** from the actual branch:
   either
   \[
   (m_T,m_K,m_\mu),
   \quad
   (R_{\rm tr},R_{\rm nt},R_\eta),
   \quad\text{or}\quad
   (q_{\rm tr},q_{\rm nt},q_\eta).
   \]

3. Evaluate
   \[
   \Delta_{\rm branch},
   \qquad
   \Delta_{\rm orbit}=(q_{\rm tr},q_{\rm nt},q_\eta).
   \]

4. Verdict:
   - if `Δ_branch ≠ 0`, the reduced GR-like 5PN branch fails;
   - if `Δ_branch = 0` but `Δ_orbit ≠ 0`, the branch is isotropic but off the exact similarity orbit;
   - only if both vanish is the reduced 5PN closure complete inside the current hierarchy.

## Most likely next real calculation

The symbolic stage is finished. The next nontrivial calculation is to extract actual PDE branch estimates for:

- the grouped-lane operator and transfer moments,
- and the dependent orbit-lock triple,

and feed them directly into the Stage-201 compilers.
