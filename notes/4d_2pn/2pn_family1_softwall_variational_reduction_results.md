# 2PN Family-1 soft-wall variational reduction: current result

## 1. What this step adds

This step turns the previously fitted static axisymmetric wall operator into the first **reduced variational wall-Hessian** that reproduces the passed 2PN support data exactly, while also tying its angular structure back to the frozen **Family-1** flared soft wall.

The correction to the earlier isotropic 4D-ball unit-test series does **not** affect this step. This derivation depends only on the already-passed **static channel data** and on exact axisymmetric angular moments, not on the corrected \(z^8\) coefficients.

---

## 2. Exact reduced modal wall-Hessian

For the static support channels
\[
K_{00}=\frac4{45},\qquad
K_{1\perp}=\frac27,\qquad
K_{10}=\frac14,\qquad
K_{20}=\frac49,\qquad
K_{21}=\frac23,\qquad
K_{22}=\frac83,
\]
a compact reduced wall law closes **exactly**:

\[
K_{00}=K_{\rm mono},
\]
\[
K_{\ell m}
=
K_{\rm mono}
+
T_0\,\ell(\ell+1)
+
\big(A_0+A_1\,\ell(\ell+1)\big)\,\langle P_2\rangle_{\ell m}
+
\big(B_0+B_1\,\ell(\ell+1)\big)\,\langle P_2^2\rangle_{\ell m},
\qquad \ell=1,2.
\]

The exact angular moments are
\[
\langle P_2\rangle_{1\perp}=-\frac15,
\quad
\langle P_2\rangle_{10}=\frac25,
\quad
\langle P_2\rangle_{20}=\frac27,
\quad
\langle P_2\rangle_{21}=\frac17,
\quad
\langle P_2\rangle_{22}=-\frac27,
\]
\[
\langle P_2^2\rangle_{1\perp}=\frac17,
\quad
\langle P_2^2\rangle_{10}=\frac{11}{35},
\quad
\langle P_2^2\rangle_{20}=\frac37,
\quad
\langle P_2^2\rangle_{21}=\frac17,
\quad
\langle P_2^2\rangle_{22}=\frac17.
\]

Solving against the passed channel data gives
\[
K_{\rm mono}=\frac4{45},
\qquad
T_0=\frac{23}{135},
\]
\[
A_0=\frac{9095}{3528},
\qquad
A_1=-\frac{25559}{21168},
\]
\[
B_0=-\frac{109}{56},
\qquad
B_1=\frac{1765}{3024}.
\]

Substituting these values back reproduces all six passed channel stiffnesses with zero residual.

So the static axisymmetric support operator is no longer just a fitted Robin table.
It is the Hessian of a concrete **reduced variational wall law** on the truncated \(\ell=0,1,2\) mouth basis.

---

## 3. Family-1 flare automatically generates the needed \(\{1,P_2,P_2^2\}\) basis

The frozen Family-1 throat profile uses a flared radius
\[
a(z)=a_0\Big(1+\beta e^{-z^2/z_m^2}\Big).
\]

On the mouth sphere, define
\[
\xi \equiv \frac{a_0^2}{z_m^2}.
\]

Expanding the flare to second order gives
\[
\frac{\delta a(\theta)}{a_0}
=
\beta\Big[
1-\frac{\xi}{3}+\frac{\xi^2}{18}
\Big]
+
\beta\Big[
-\frac{2\xi}{3}+\frac{2\xi^2}{9}
\Big]P_2
+
\beta\Big[
\frac{2\xi^2}{9}
\Big]P_2^2
+O(\xi^3).
\]

This is the key structural result:

> the actual frozen Family-1 flare automatically produces the exact axisymmetric basis \(\{1,P_2,P_2^2\}\) that the reduced wall-Hessian needs.

So the earlier \(P_2^2\) term is no longer ad hoc.
At the reduced wall-law level, it is the generic second-order imprint of the flared Family-1 soft wall.

---

## 4. Minimal linear-gradient interpretation

There is also a useful sharp sub-result.

If the **anisotropic gradient block** of the reduced wall law is assumed to be **linear** in the Family-1 flare profile, then
\[
\frac{B_1}{A_1}=rac{D_2}{D_1}=rac{\xi}{\xi-3},
\]
where
\[
D_1=-\frac{2\xi}{3}+\frac{2\xi^2}{9},
\qquad
D_2=\frac{2\xi^2}{9}.
\]

Using the exact solved coefficients gives
\[
\frac{B_1}{A_1}=-\frac{12355}{25559},
\]
so
\[
\xi=\frac{12355}{12638}\approx 0.9776072163,
\qquad
\frac{z_m}{a_0}=\frac{1}{\sqrt{\xi}}\approx 1.0113880097.
\]

So the minimal linear-gradient reading points to a flare width of order the throat radius,
\[
z_m\sim a_0.
\]

---

## 5. What this means

The constructive state of play is now:

1. the conservative 2PN comparable-mass cross sector is solved algebraically,
2. that solution was rebuilt as a small mouth-port / support-minus-closure EFT,
3. the static axisymmetric support data now come from an explicit **reduced variational wall-Hessian**,
4. and the exact angular basis of that Hessian is now traced back to the actual **Family-1** flared soft wall rather than inserted by hand.

So the remaining “physics gap” is narrower than before.

The next derivation target is no longer “invent the 2PN wall law.”
It is:

- derive the coefficients \((K_{\rm mono},T_0,A_0,A_1,B_0,B_1)\) from a more explicit soft-wall or traction-balance reduction,
- then fold that wall-Hessian into the inner-throat DtN / port solver so the support coefficients emerge from the PDE side directly.

That is a real step toward closing the full 2PN derivation.
