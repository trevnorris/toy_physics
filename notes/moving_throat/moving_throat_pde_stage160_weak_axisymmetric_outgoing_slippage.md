# Moving-Throat PDE — Stage 160: Weak-Axisymmetric Collapse of the Outgoing-Slippage Bundle to One Scalar Amplitude

## Purpose

Stage 159 factored the remaining linear grouped weak-axisymmetric defect into three microscopic port slippages,
\[
\mathcal M_r=\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K},
\qquad
\mathcal I_r=\frac{R_rG_{U,r}}{\Omega_{U,r}^2G_{W,r}},
\qquad
\mathcal H_r=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2},
\]
and showed that, on conservative-shape-preserving branches,
\[
\Xi_{\rm load}=
\sum_r \rho_r^{(N)}
\Big[
2\,\delta\ln\mathcal M_r
+2\,\delta\ln(1+\mathcal I_r)
-2\,\delta\ln(1-\mathcal H_r)
\Big].
\]

That is already sharp, but it is still phrased port by port.
The next honest step is therefore not to reopen the full grouped outlet problem again, but to ask:

> what does this factorization look like on the actual weak-axisymmetric grouped branch?

This stage answers that question exactly.

The main result is that, on the weak-axisymmetric branch,
all three microscopic outgoing slippages inherit the same grouped signature
\[
\lambda_{20}=1,
\qquad
\lambda_{21}=\frac12,
\qquad
\lambda_{22}=-1,
\]
so the full remaining grouped defect collapses to **one scalar amplitude**
\[
\boxed{
\Xi_{\rm load}^{(A)}=\epsilon\,\lambda_A\,\Xi_1.
}
\]
Equivalently,
\[
\boxed{
\Xi_1
=
\sum_r \rho_r^{(N)}
\left[
2\,\mathfrak m_r
+\frac{2\mathcal I_r}{1+\mathcal I_r}\,\mathfrak i_r
+\frac{2\mathcal H_r}{1-\mathcal H_r}\,\mathfrak h_r
\right].
}
\]
So the whole Stage-159 port bundle has now collapsed one step further:

- not to every microscopic grouped drift separately,
- but to one weighted scalar combination \(\Xi_1\).

On the same conservative-shape-preserving branch already used in Stages 157–159, this scalar is exactly the logarithmic outgoing-prefactor slope:
\[
\boxed{
\Xi_1=\frac{P_1}{P_0}.
}
\]

So the remaining linear grouped `2.5`PN bottleneck is now smaller still.
It is no longer “compute all outgoing slippages lane by lane.”
It is just:

> compute the single weak-axisymmetric amplitude \(\Xi_1\), or equivalently \(P_1/P_0\), from the microscopic moving-throat branch.

---

## 1. Weak-axisymmetric microscopic slope bookkeeping

Write every grouped outgoing microscopic quantity in logarithmic weak-axisymmetric form:
\[
\delta\ln K_A=\epsilon\lambda_A\,\kappa_1,
\]
\[
\delta\ln G_{W,A,r}=\epsilon\lambda_A\,\mathfrak g_{W,r},
\qquad
\delta\ln G_{U,A,r}=\epsilon\lambda_A\,\mathfrak g_{U,r},
\]
\[
\delta\ln R_{A,r}=\epsilon\lambda_A\,\mathfrak r_r,
\qquad
\delta\ln\Omega_{U,A,r}^2=\epsilon\lambda_A\,\mathfrak o_{U,r},
\qquad
\delta\ln\Omega_{W,A,r}^2=\epsilon\lambda_A\,\mathfrak o_{W,r}.
\]

So the grouped weak-axisymmetric branch is parameterized by the scalar wall-baseline slope
\(\kappa_1\)
and the portwise logarithmic microscopic slopes
\(
\mathfrak g_{W,r},\mathfrak g_{U,r},\mathfrak r_r,\mathfrak o_{U,r},\mathfrak o_{W,r}
\).

This is the natural Stage-159 continuation because the outgoing-load theorem is already logarithmic.

---

## 2. Exact weak-axisymmetric transport of the three port slippages

By definition,
\[
\mathcal M_r=\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K},
\qquad
\mathcal I_r=\frac{R_rG_{U,r}}{\Omega_{U,r}^2G_{W,r}},
\qquad
\mathcal H_r=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2}.
\]
Therefore the weak-axisymmetric grouped slopes are
\[
\boxed{
\delta\ln\mathcal M_{A,r}=\epsilon\lambda_A\,\mathfrak m_r,
\qquad
\mathfrak m_r:=\mathfrak g_{W,r}-\mathfrak o_{W,r}-\frac12\kappa_1,
}
\]
\[
\boxed{
\delta\ln\mathcal I_{A,r}=\epsilon\lambda_A\,\mathfrak i_r,
\qquad
\mathfrak i_r:=\mathfrak r_r+\mathfrak g_{U,r}-\mathfrak o_{U,r}-\mathfrak g_{W,r},
}
\]
\[
\boxed{
\delta\ln\mathcal H_{A,r}=\epsilon\lambda_A\,\mathfrak h_r,
\qquad
\mathfrak h_r:=2\mathfrak r_r-\mathfrak o_{U,r}-\mathfrak o_{W,r}.
}
\]

So the three Stage-159 microscopic outgoing slippages automatically inherit the same grouped signature.
That is the first exact collapse of the stage.

---

## 3. Exact portwise outgoing-defect amplitude

Stage 159 already gave the first-order defect field
\[
\Sigma_r^{(N)}
=
2\,\delta\ln\mathcal M_r
+\frac{2\mathcal I_r}{1+\mathcal I_r}\,\delta\ln\mathcal I_r
+\frac{2\mathcal H_r}{1-\mathcal H_r}\,\delta\ln\mathcal H_r.
\]
Substituting the weak-axisymmetric slopes above gives
\[
\boxed{
\Sigma_{A,r}^{(N)}=\epsilon\lambda_A\,\sigma_r,
}
\]
with the exact port amplitude
\[
\boxed{
\sigma_r
:=
2\,\mathfrak m_r
+\frac{2\mathcal I_r}{1+\mathcal I_r}\,\mathfrak i_r
+\frac{2\mathcal H_r}{1-\mathcal H_r}\,\mathfrak h_r.
}
\]
So at the level of each outgoing port, the whole Stage-159 defect bundle collapses to one weak-axisymmetric amplitude \(\sigma_r\).

---

## 4. Exact grouped collapse of the full remaining defect

On the conservative-shape-preserving branch,
\[
\Xi_{\rm load}=
\sum_r \rho_r^{(N)}\Sigma_r^{(N)}.
\]
Therefore, on the weak-axisymmetric grouped branch,
\[
\boxed{
\Xi_{\rm load}^{(A)}
=
\epsilon\lambda_A\,\Xi_1,
\qquad
\Xi_1:=\sum_r \rho_r^{(N)}\sigma_r.
}
\]
Equivalently,
\[
\boxed{
\Xi_1
=
\sum_r \rho_r^{(N)}
\left[
2\,\mathfrak m_r
+\frac{2\mathcal I_r}{1+\mathcal I_r}\,\mathfrak i_r
+\frac{2\mathcal H_r}{1-\mathcal H_r}\,\mathfrak h_r
\right].
}
\]
This is the stage’s main structural theorem.

The grouped-lane pattern is now completely fixed:
\[
\boxed{
\Xi_{\rm load}^{(20)}=\epsilon\,\Xi_1,
\qquad
\Xi_{\rm load}^{(21)}=\frac{\epsilon}{2}\,\Xi_1,
\qquad
\Xi_{\rm load}^{(22)}=-\epsilon\,\Xi_1.
}
\]
So the remaining linear grouped defect is once again one-dimensional.

---

## 5. Grouped anisotropy defects of the outgoing-load theorem

Because the grouped weak-axisymmetric signature is fixed,
\[
\lambda_{20}=1,
\qquad
\lambda_{21}=\frac12,
\qquad
\lambda_{22}=-1,
\]
the grouped trace/anomaly variables of the remaining defect are
\[
\bar\Xi=0,
\qquad
 a_\Xi=\frac{\epsilon}{4}\,\Xi_1,
\qquad
 b_\Xi=\frac{3\epsilon}{4}\,\Xi_1.
\]
So the exact weak-axisymmetric relation is
\[
\boxed{b_\Xi=3a_\Xi.}
\]
This is the outgoing-load analog of the earlier grouped weak-axisymmetric transport laws.
It says the first anisotropy of the full remaining grouped defect is never arbitrary; it always lies on the one-dimensional axisymmetric fingerprint.

---

## 6. Identification with the physical outgoing-prefactor slope

Stage 156 already showed that on the weak-axisymmetric branch
\[
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}
=\Xi_{\rm load}.
\]
Combining that with the present Stage-160 collapse gives the exact identification
\[
\boxed{
\frac{P_1}{P_0}=\Xi_1.
}
\]
So the Stage-159 outgoing-slippage factorization and the Stage-156 physical prefactor slope are now the **same object**.

This is the cleanest interpretation reached so far:

> the remaining weak-axisymmetric grouped `2.5`PN defect is exactly the logarithmic slope of the outgoing prefactor, and its microscopic content is the weighted port combination \(\Xi_1\) above.

That is the second main theorem of the stage.

---

## 7. Quadrupole-normalization defect on the even-preserving branch

Stage 156 also showed that on the even-preserving branch the remaining linear grouped normalization defect is
\[
\Delta_Q^{(20)}=\epsilon\,\Xi_{\rm load},
\qquad
\Delta_Q^{(21)}=\frac{\epsilon}{2}\,\Xi_{\rm load},
\qquad
\Delta_Q^{(22)}=-\epsilon\,\Xi_{\rm load}.
\]
So the present Stage-160 result immediately gives
\[
\boxed{
\Delta_Q^{(20)}=\epsilon\,\Xi_1,
\qquad
\Delta_Q^{(21)}=\frac{\epsilon}{2}\,\Xi_1,
\qquad
\Delta_Q^{(22)}=-\epsilon\,\Xi_1,
}
\]
with
\[
\boxed{
\Xi_1=\frac{P_1}{P_0}
=
\sum_r \rho_r^{(N)}
\left[
2\,\mathfrak m_r
+\frac{2\mathcal I_r}{1+\mathcal I_r}\,\mathfrak i_r
+\frac{2\mathcal H_r}{1-\mathcal H_r}\,\mathfrak h_r
\right].
}
\]
So on the even-preserving branch the full remaining linear grouped `2.5`PN theorem problem has now collapsed to one scalar amplitude.

---

## 8. Rigidity and dominant-port corollaries

### 8.1 Interference/hybridization rigidity

If the interference and hybridization ratios are rigid,
\[
\mathfrak i_r=0,
\qquad
\mathfrak h_r=0,
\]
then
\[
\boxed{
\Xi_1=2\sum_r \rho_r^{(N)}\,\mathfrak m_r.
}
\]
So the weak-axisymmetric outgoing-prefactor slope becomes the weighted square-root mixed-leg slope.

### 8.2 Dominant-port limit

If one outgoing port dominates,
\[
\rho_{r_*}^{(N)}\approx 1,
\qquad
\rho_{r\neq r_*}^{(N)}\approx 0,
\]
then
\[
\boxed{
\Xi_1\approx
2\,\mathfrak m_{r_*}
+\frac{2\mathcal I_{r_*}}{1+\mathcal I_{r_*}}\,\mathfrak i_{r_*}
+\frac{2\mathcal H_{r_*}}{1-\mathcal H_{r_*}}\,\mathfrak h_{r_*}.
}
\]
Under interference/hybridization rigidity this reduces further to
\[
\boxed{
\Xi_1\approx 2\,\mathfrak m_{r_*}.
}
\]
So in the dominant-port rigid branch, zero defect is exactly the square-root mixed-leg slope condition
\[
\mathfrak m_{r_*}=0
\qquad\Longleftrightarrow\qquad
\delta\ln\!\left(\frac{G_{W,r_*}}{\Omega_{W,r_*}^2}\right)=\frac12\,\delta\ln K.
\]

---

## 9. What Stage 160 changes

Stage 159 reduced the remaining grouped weak-axisymmetric defect to three microscopic slippages per outgoing port:
\(
\mathcal M_r,\mathcal I_r,\mathcal H_r
\).
That was already a big collapse, but it still left the grouped theorem problem looking like a port-by-port vector problem.

Stage 160 sharpens the problem one more step.

It shows that on the weak-axisymmetric grouped branch:

1. each port slippage inherits the same grouped signature,
2. each port therefore contributes only one scalar amplitude \(\sigma_r\),
3. the full remaining grouped defect collapses to one weighted scalar \(\Xi_1\),
4. that scalar is exactly the physical outgoing-prefactor slope \(P_1/P_0\),
5. and on the even-preserving branch it is exactly the remaining linear grouped `2.5`PN normalization defect.

So the next honest theorem gate is now narrower again.
It is no longer

> “compute the grouped outgoing slippage bundle somehow.”

It is simply:

> compute the single weak-axisymmetric microscopic amplitude
> \[
> \Xi_1=\frac{P_1}{P_0}
> \]
> from the real moving-throat outgoing ports.

That is the direct continuation point.
