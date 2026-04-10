# Stage 151–160 continuation notes

This batch fills the missing bridge between the Stage 150 tangent-compensation theorem and the later outgoing-load / similarity-orbit chain.

The key logical contraction is:

1. **Stage 151**: first-order off-bundle motion is not a large microscopic vector; it collapses to three scalar slippages \((\varepsilon_L,\varepsilon_v,\varepsilon_T)\), and then to one weighted scalar \(\varepsilon_\perp\) with \(\delta_\perp=-\varepsilon_\perp\).

2. **Stage 152**: pure grouped real `P2` anisotropy cannot linearly feed those scalar slippages. Its scalar feed-down begins only at quadratic order through the grouped invariant
   \[
   \mathcal I[x,y]=\frac15\,\delta x^T G_{\rm grp}\,\delta y
   =4a_x a_y+\frac45 b_x b_y.
   \]
   On the weak-axisymmetric branch this becomes
   \[
   \mathcal I[x,y]=\frac{7}{10}\epsilon^2 x^{(1)}y^{(1)}.
   \]

3. **Stage 153**: the remaining **linear** grouped problem therefore lives only in the direct outlet coefficients. It collapses to
   \[
   \mathcal K_A=\delta D_{A,2}+\frac19\,\delta D_{A,0},
   \qquad
   \mathcal G_A=\delta N_{A,0}-P_0\,\delta D_{A,0},
   \]
   together with the hidden-even compatibility relation
   \[
   \delta D_{A,4}=\frac23\,\delta D_{A,2}+\frac1{27}\,\delta D_{A,0}.
   \]

4. **Stage 154**: those two grouped obstructions are not sourced by the whole microscopic bundle independently. They decompose into:
   \[
   \mathcal K_A=\mathcal W_A-\mathcal B_A-\mathcal Z_A,
   \qquad
   \mathcal G_A=-P_0\delta K_A+P_0\delta B_{A,0}+\mathcal N_A.
   \]
   So the linear grouped-even/odd problem is driven only by wall, BdG, conservative Maxwell/mixed, and outgoing-transfer anisotropies.

5. **Stage 155**: on the canonical compensated branch, the microscopic obstruction pair is just the physical slope pair
   \[
   \mathfrak K_1=-D_0 u_2^{(1)},
   \qquad
   \mathfrak G_1=D_0 P_1.
   \]
   The direct outlet amplitudes become
   \[
   \kappa_1=-\frac{3(1-\sigma_*)}{\sigma_*}u_2^{(1)},
   \qquad
   \gamma_1=-\frac{1-\sigma_*}{9\sigma_*}\frac{P_1}{P_0}.
   \]

6. **Stage 156**: expanding the actual grouped operator moments on the weak-axisymmetric branch yields
   \[
   u_2^{(1)}=-\frac{D_{21}+D_{01}/9}{D_0},
   \qquad
   \frac{P_1}{P_0}=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
   \]
   On the even-preserving branch,
   \[
   D_{21}=-\frac{D_{01}}{9},
   \qquad
   D_{41}=-\frac{D_{01}}{27},
   \]
   so the whole remaining linear grouped normalization defect is one static loading mismatch
   \[
   \Xi_{\rm load}:=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
   \]

7. **Stage 157**: that mismatch is the weighted failure of static self-similarity relative to the wall baseline:
   \[
   \Xi_{\rm load}
   =(\delta_N-\delta_K)+\omega_B(\delta_B-\delta_K)+\omega_Z(\delta_Z-\delta_K).
   \]

8. **Stage 158**: factoring by the wall baseline yields wall-normalized shape/load variables
   \[
   B_{0,\alpha}=K\chi_\alpha^2,
   \qquad
   Z_0^{(r)}=K\Upsilon_r,
   \qquad
   N_0^{(r)}=\Lambda_r^2,
   \]
   and the outgoing-load theorem
   \[
   2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r=\delta_K
   \]
   on conservative-shape-preserving branches.

9. **Stage 159**: the outgoing load factor splits exactly into
   \[
   \frac{\Lambda_r^2}{K}
   =
   \mathcal M_r^2\frac{(1+\mathcal I_r)^2}{(1-\mathcal H_r)^2},
   \]
   with
   \[
   \mathcal M_r=\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K},
   \quad
   \mathcal I_r=\frac{R_rG_{U,r}}{\Omega_{U,r}^2G_{W,r}},
   \quad
   \mathcal H_r=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2}.
   \]
   Under interference/hybridization rigidity, zero defect is the square-root mixed-leg law
   \[
   \frac{G_{W,r}}{\Omega_{W,r}^2}\propto \sqrt K.
   \]

10. **Stage 160**: on the weak-axisymmetric branch the whole outgoing-slippage bundle collapses to one scalar amplitude
    \[
    \Xi_1
    =
    \sum_r \rho_r^{(N)}
    \left[
      2\mathfrak m_r
      +\frac{2\mathcal I_r}{1+\mathcal I_r}\mathfrak i_r
      +\frac{2\mathcal H_r}{1-\mathcal H_r}\mathfrak h_r
    \right],
    \]
    with the grouped pattern
    \[
    \Xi_{\rm load}^{(20)}=\epsilon\,\Xi_1,
    \qquad
    \Xi_{\rm load}^{(21)}=\frac{\epsilon}{2}\Xi_1,
    \qquad
    \Xi_{\rm load}^{(22)}=-\epsilon\,\Xi_1,
    \]
    and the exact physical identification
    \[
    \Xi_1=\frac{P_1}{P_0}.
    \]

So after Stage 160, the remaining weak-axisymmetric grouped `2.5`PN bottleneck is no longer a diffuse outlet-bundle problem. It is just the single microscopic amplitude \(\Xi_1 = P_1/P_0\) on the actual moving-throat branch.
