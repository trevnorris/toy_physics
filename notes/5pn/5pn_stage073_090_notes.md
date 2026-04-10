
# 5PN / Moving-Throat Continuation — Stages 73–90

This batch continues directly from the Stage 72 result that the explicit Family-1 support/source side is no longer the active bottleneck. The live question is whether the actual grouped-\(P_2\)/geometry branch realizes the minimal isotropic contact-plus-pole conservative quadrupole module, and then whether the passive/outgoing \(l=2\) branch carries the canonical outgoing normalization.

## Stage 73
Recasts the post-72 status in theorem language: the explicit Family-1 support/source branch already succeeds on the minimal isotropic branch with
\[
\rho_\alpha=\frac43,\qquad \zeta_{\rm req}=\frac13,\qquad Pe_{\rm req}=0,
\]
so the remaining reduced theorem gap is no longer on the support/source side.

## Stage 74
Derives the `3/4 + 1/4` conservative module directly from the 3PN conservative split.
If the isotropic grouped-\(P_2\) branch is carried by one effective pole and the geometry lane is static through \(O(\omega^4)\), then the minimal isotropic branch identity forces
\[
K_{\rm pole}=\frac{K_0}{4},\qquad K_{\rm geom}=\frac{3K_0}{4},
\]
hence
\[
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
\]

## Stage 75
Allows the geometry lane to carry dynamic even moments and derives the exact obstruction formula
\[
c_{\rm pole}
=
\frac{1+\epsilon_4}{4(1+\epsilon_2)^2},
\qquad
\epsilon_2=\frac{\Omega_Q^2 K_{(g,2)}}{K_{\rm pole}},
\qquad
\epsilon_4=\frac{\Omega_Q^4 K_{(g,4)}}{K_{\rm pole}}.
\]
So the `3/4 + 1/4` split is exact iff the geometry lane is static at the relevant orders.

## Stage 76
Freezes the updated reduced status: the only remaining reduced ambiguity is whether the actual moving-throat geometry lane is dynamically inert through \(O(\omega^4)\).

## Stage 77
Proves the isotropic geometry-decoupling theorem. In the exact isotropic quadratic wall theory, the \(l=0\) scalar/geometry lane is orthogonal to the grouped real \(l=2\) bundle, so no dynamic \(O(\omega^2)\) or \(O(\omega^4)\) geometry contamination can enter the grouped quadrupole carrier at linear order:
\[
K_{(g,2)}=K_{(g,4)}=0,\qquad \epsilon_2=\epsilon_4=0.
\]

## Stage 78
Shows that if weak anisotropy later induces an \(l=0\leftrightarrow l=2\) mixing source, the first nonzero geometry contamination appears only at second order in that mixing:
\[
\epsilon_2,\epsilon_4 = O(\chi^2).
\]
So the Stage-74 split is perturbatively stable.

## Stage 79
Records the actual reduced-branch verdict: on the natural isotropic branch,
\[
\epsilon_2=\epsilon_4=0,
\]
hence
\[
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2},
\qquad
\rho_\alpha=\frac43,\qquad
\zeta_{\rm req}=\frac13.
\]

## Stage 80
Shows that once the actual isotropic grouped-\(P_2\) one-pole branch is accepted, the full reduced 2.5PN normalization problem collapses to one scalar defect
\[
N_Q:=\frac{\overline K_0}{\overline K_0^{\rm target}}.
\]
All low-frequency coefficients scale by the same factor on that branch:
\[
\overline K_2,\ \overline K_4,\ \overline\Gamma_5 \propto N_Q.
\]

## Stage 81
Proves the explicit Family-1 support/source theorem is automatic on the actual isotropic branch. Since
\[
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}^{(\rm act)}(\epsilon_{\rm blk})=\frac{1}{3-2\epsilon_{\rm blk}},
\]
any explicit family with \(\zeta_{\max}>1\) already passes throughout the admissible blocked regime. Family-1 has
\[
\zeta_{\max}^{(F1)}\approx 2.46752922945601 > 1,
\]
so it is automatically safe.

## Stage 82
Freezes the reduced finish line: the grouped-\(P_2\) conservative branch is geometry-clean, the support/source side is automatic, and the only remaining reduced theorem gap is the single normalization defect \(N_Q-1\).

## Stage 83
Factors the last 2.5PN obstruction into a conservative and an outgoing piece. Introducing one outgoing-normalization factor \(\chi_Q\),
\[
\widehat Y_Q^{\rm ret}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2-i\chi_Q \sigma_Q^{\rm can}\omega^5}+O(\omega^6),
\]
one gets
\[
\frac{\overline\Gamma_5}{\overline\Gamma_5^{\rm target}} = \chi_Q N_Q,
\]
and the full odd normalization condition is
\[
\hat m_0^{\,2}\chi_Q N_Q = 1.
\]

## Stage 84
Uses the natural point-particle source map \(\hat m_0\to 1\) to show the remaining reduced obstruction is purely outgoing:
\[
N_Q=\frac{1}{\chi_Q}.
\]
So all remaining reduced uncertainty is now in the outgoing normalization factor \(\chi_Q\).

## Stage 85
Shows higher odd retarded data beginning at \(O(\omega^7)\) are irrelevant to the 2.5PN theorem. The only live retarded obstruction at 2.5PN is the leading \(\omega^5\) normalization factor \(\chi_Q\).

## Stage 86
States the conditional reduced closure:
if \(\chi_Q=1\), the reduced GR-like point-particle 2.5PN theorem is closed; if not, the whole remaining failure is measured by \(\Delta_Q=\chi_Q-1\).

## Stage 87
Computes the exact outgoing spherical \(l=2\) DtN fingerprint:
\[
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}-\frac{2z^6}{27}-i\frac{z^7}{27}+O(z^8),
\]
and therefore
\[
\widehat Y_2^{\rm out}(z)
=
1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}-\frac{11z^6}{729}-i\frac{z^7}{243}+O(z^8).
\]

## Stage 88
Matches the Stage-87 DtN fingerprint to the retarded grouped-\(P_2\) one-pole-plus-contact module and fixes
\[
\chi_Q=1
\]
on the canonical compact passive/outgoing branch. A deformed DtN branch with
\[
\Lambda_2^{\rm def}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\xi_Q\frac{z^5}{9}+O(z^6)
\]
simply gives \(\chi_Q=\xi_Q\).

## Stage 89
Inserts \(\chi_Q=1\) into the reduced normalization stack and closes the reduced 2.5PN theorem on the canonical outgoing branch. In the strict point-particle limit,
\[
N_Q=1,
\]
and the canonical invariant coefficients are exactly
\[
\overline K_0=\frac{54Gc_s^5}{5a^5c^5},\qquad
\overline K_2=\frac{6Gc_s^3}{5a^3c^5},\qquad
\overline K_4=\frac{8Gc_s}{15ac^5},\qquad
\overline\Gamma_5=\frac{2G}{5c^5}.
\]

## Stage 90
Builds the first explicit isotropic DtN deformation algebra. For
\[
\Lambda_2^{\rm def}(z)
=
S\,\Lambda_2^{\rm out}(\beta z)
+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5+O(z^6),
\]
the normalized outgoing factor is
\[
\chi_Q=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0},
\]
with \(\Sigma_2,\Sigma_4\) fixed by the requirement that the canonical even moments remain unchanged. So the only isotropic branch data that can shift \(\chi_Q\) are \(\beta,\Sigma_0,\Sigma_5\).

## Net result of Stages 73–90

Inside the present reduced hierarchy, the conservative grouped-\(P_2\)/geometry branch now reaches the same minimal isotropic `3/4 + 1/4` module that the 2.5PN audit wanted, the explicit Family-1 support/source side is automatic, and the outgoing \(l=2\) DtN branch fixes the last reduced scalar to
\[
\chi_Q=1
\]
on the canonical compact branch. So the reduced nonspinning point-particle 2.5PN theorem is closed on that canonical branch; what remains genuinely PDE-facing is branch realization and explicit isotropic DtN deformation data.
