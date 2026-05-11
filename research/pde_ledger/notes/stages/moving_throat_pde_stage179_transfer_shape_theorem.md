# Moving-Throat PDE — Stage 179: Wall-Normalized Transfer-Shape Theorem and the Final Collapse of the Port Co-Loading Gate

## Purpose

Stage 246 reduced the remaining weak-axisymmetric grouped `2.5`PN defect to

a single scalar,
\[
\Xi_1=\bar\nu_N-\kappa_1,
\]
where
\[
\nu_r=\frac{\delta\ln N_{A,0}^{(r)}}{\epsilon\lambda_A}
\]
is the weak-axisymmetric logarithmic slope of the actual static outgoing-transfer coefficient
\[
N_{A,0}^{(r)}=\frac{P_{A,r}^2}{\Delta_{A,r}^2}.
\]

That was already sharp, but it still left the theorem gate phrased in terms of the raw pair
\((P_r,\Delta_r)\).
The next honest step is therefore to normalize the actual outgoing ports directly by the wall baseline and ask whether the remaining defect is really nothing more than the drift of one dimensionless transfer shape.

That is what this stage does.

The main result is that each outgoing port admits an exact wall-normalized factorization
\[
\boxed{
N_{A,0}^{(r)} = K_A\,\mathcal T_{A,r}^2,
}
\]
with one dimensionless transfer shape
\[
\boxed{
\mathcal T_{A,r}
=
\frac{\widehat G_{W,A,r}+\widehat R_{A,r}\,\widehat G_{U,A,r}}{1-\widehat R_{A,r}^2}.
}
\]
So the remaining grouped defect collapses one step further:
\[
\boxed{
\Xi_1 = 2\sum_r \rho_r^{(N)}\,\tau_r,
\qquad
\tau_r:=\frac{\delta\ln\mathcal T_{A,r}}{\epsilon\lambda_A}.
}
\]

Equivalently, the Stage 246 co-loading gate is now exactly this:

> does the wall-normalized transfer shape of the actual moving-throat outgoing ports stay weak-axisymmetrically rigid?

That is the direct continuation point.

---

## 1. Exact wall-normalized factorization of the actual static outgoing-transfer coefficient

Recall the actual static outgoing-port data from Stage 246:
\[
P_r=\Omega_{U,r}^2 G_{W,r}+R_r G_{U,r},
\qquad
\Delta_r=\Omega_{U,r}^2\Omega_{W,r}^2-R_r^2,
\qquad
N_0^{(r)}=\frac{P_r^2}{\Delta_r^2}.
\]

Introduce the wall-normalized dimensionless port variables
\[
\boxed{
\widehat G_{W,r}:=\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K},
\qquad
\widehat G_{U,r}:=\frac{G_{U,r}}{\Omega_{U,r}\Omega_{W,r}\sqrt K},
\qquad
\widehat R_r:=\frac{R_r}{\Omega_{U,r}\Omega_{W,r}}.
}
\]
Here
\(\Omega_{U,r}:=\sqrt{\Omega_{U,r}^2}\) and \(\Omega_{W,r}:=\sqrt{\Omega_{W,r}^2}\).

Then the actual static numerator and detuning factor exactly as
\[
P_r
=
\sqrt K\,\Omega_{U,r}^2\Omega_{W,r}^2
\bigl(\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}\bigr),
\]
\[
\Delta_r
=
\Omega_{U,r}^2\Omega_{W,r}^2(1-\widehat R_r^2).
\]
Therefore
\[
\boxed{
\frac{N_0^{(r)}}{K}
=
\left[
\frac{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}}{1-\widehat R_r^2}
\right]^2.
}
\]
This is the exact wall-normalized factorization of the actual static outgoing-transfer coefficient.
Define the dimensionless transfer shape
\[
\boxed{
\mathcal T_r
:=
\frac{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}}{1-\widehat R_r^2},
}
\]
so that
\[
\boxed{
N_0^{(r)}=K\,\mathcal T_r^2.
}
\]
Thus the actual static outgoing-transfer coefficient is the wall baseline times one dimensionless transfer shape squared.

---

## 2. Weak-axisymmetric transport of the wall-normalized port variables

On the weak-axisymmetric grouped branch,
\[
\lambda_{20}=1,
\qquad
\lambda_{21}=\frac12,
\qquad
\lambda_{22}=-1,
\]
and the wall baseline obeys
\[
\delta\ln K_A=\epsilon\lambda_A\,\kappa_1.
\]

Now define the weak-axisymmetric slopes of the three wall-normalized port variables by
\[
\boxed{
\delta\ln\widehat G_{W,A,r}=\epsilon\lambda_A\,\mathfrak w_r,
}
\]
\[
\boxed{
\delta\ln\widehat G_{U,A,r}=\epsilon\lambda_A\,\mathfrak u_r,
}
\]
\[
\boxed{
\delta\ln\widehat R_{A,r}=\epsilon\lambda_A\,\mathfrak c_r.
}
\]
In terms of the primitive port slopes already used in Stages 228–229,
\[
\boxed{
\mathfrak w_r
=
\mathfrak g_{W,r}-\mathfrak o_{W,r}-\frac12\kappa_1,
}
\]
\[
\boxed{
\mathfrak u_r
=
\mathfrak g_{U,r}-\frac12\mathfrak o_{U,r}-\frac12\mathfrak o_{W,r}-\frac12\kappa_1,
}
\]
\[
\boxed{
\mathfrak c_r
=
\mathfrak r_r-\frac12\mathfrak o_{U,r}-\frac12\mathfrak o_{W,r}.
}
\]
So the three wall-normalized port variables already isolate the wall-baseline scaling from the genuine dimensionless port-shape drift.

---

## 3. Exact transfer-shape slope and the identity `\nu_r = \kappa_1 + 2\tau_r`

Because
\[
\mathcal T_r
=
\frac{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}}{1-\widehat R_r^2},
\]
its weak-axisymmetric logarithmic slope is
\[
\boxed{
\delta\ln\mathcal T_{A,r}=\epsilon\lambda_A\,\tau_r,
}
\]
with
\[
\boxed{
\tau_r
=
\widehat\alpha_r\,\mathfrak w_r
+
\widehat\beta_r\,(\mathfrak u_r+\mathfrak c_r)
+
\frac{2\widehat R_r^2}{1-\widehat R_r^2}\,\mathfrak c_r,
}
\]
where
\[
\boxed{
\widehat\alpha_r:=\frac{\widehat G_{W,r}}{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}},
\qquad
\widehat\beta_r:=\frac{\widehat R_r\widehat G_{U,r}}{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}},
\qquad
\widehat\alpha_r+\widehat\beta_r=1.
}
\]

Now use the exact factorization
\[
N_{A,0}^{(r)}=K_A\mathcal T_{A,r}^2.
\]
Taking the weak-axisymmetric logarithmic slope gives
\[
\boxed{
\nu_r
=
\frac{\delta\ln N_{A,0}^{(r)}}{\epsilon\lambda_A}
=
\kappa_1+2\tau_r.
}
\]
So the Stage 246 port slope is just the wall-baseline slope plus twice the transfer-shape slope.

This is the central identity of the stage.

---

## 4. Exact collapse of the remaining grouped defect

Stage 246 already gave
\[
\Xi_1
=
\sum_r \rho_r^{(N)}(\nu_r-\kappa_1),
\qquad
\sum_r \rho_r^{(N)}=1.
\]
Substituting the transfer-shape identity above yields
\[
\boxed{
\Xi_1
=
2\sum_r \rho_r^{(N)}\tau_r.
}
\]
Equivalently, the remaining weak-axisymmetric grouped defect is twice the outgoing-weighted average transfer-shape slope.

So the exact zero-defect condition becomes
\[
\boxed{
\sum_r \rho_r^{(N)}\tau_r=0.
}
\]
This is the sharpest reduced theorem gate reached so far.

---

## 5. Exact equivalence to the Stage 244/228/229 slippage languages

Stage 244 used
\[
\mathcal M_r=\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K},
\qquad
\mathcal I_r=\frac{R_r G_{U,r}}{\Omega_{U,r}^2 G_{W,r}},
\qquad
\mathcal H_r=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2}.
\]
These are related to the present variables by
\[
\boxed{
\mathcal M_r=\widehat G_{W,r},
\qquad
\mathcal I_r=\frac{\widehat R_r\widehat G_{U,r}}{\widehat G_{W,r}},
\qquad
\mathcal H_r=\widehat R_r^2.
}
\]
Likewise the Stage 245/229 microscopic slopes satisfy
\[
\boxed{
\mathfrak m_r=\mathfrak w_r,
\qquad
\mathfrak i_r=(\mathfrak u_r+\mathfrak c_r)-\mathfrak w_r,
\qquad
\mathfrak h_r=2\mathfrak c_r.
}
\]
With these substitutions, the present transfer-shape slope becomes exactly
\[
\boxed{
\tau_r
=
\mathfrak m_r
+
\frac{\mathcal I_r}{1+\mathcal I_r}\,\mathfrak i_r
+
\frac{\mathcal H_r}{1-\mathcal H_r}\,\mathfrak h_r.
}
\]
So Stage 245’s port amplitude
\(
\sigma_r
\)
is simply
\[
\boxed{
\sigma_r=2\tau_r.
}
\]
This proves that the present theorem is not a different branch of the algebra.
It is the exact compressed form of the Stage 244/228/229 slippage theorems.

---

## 6. Port-transfer-shape rigidity theorem

The exact zero-defect condition is
\[
\sum_r \rho_r^{(N)}\tau_r=0.
\]
A stronger per-port sufficient condition is
\[
\boxed{
\tau_r=0
\qquad\text{for every active outgoing port }r.
}
\]
Because
\(
\tau_r = \delta\ln\mathcal T_{A,r}/(\epsilon\lambda_A)
\),
this is equivalent to
\[
\boxed{
\delta\ln\mathcal T_{A,r}=0
\qquad\text{for every active outgoing port }r.
}
\]
So the exact reduced meaning of outgoing-port co-loading is now:

> each wall-normalized transfer shape must be weak-axisymmetrically rigid.

When this holds,
\[
N_{A,0}^{(r)}=K_A\mathcal T_r^2+O(\epsilon^2),
\]
so every static outgoing-transfer coefficient co-loads with the wall baseline individually, and therefore
\[
\boxed{\Xi_1=0.}
\]

---

## 7. Corollaries

### 7.1 Dominant-port limit

If one outgoing port dominates,
\[
\rho_{r_*}^{(N)}\approx1,
\qquad
\rho_{r\neq r_*}^{(N)}\approx0,
\]
then
\[
\boxed{
\Xi_1\approx 2\tau_{r_*}.
}
\]
So in the dominant-port regime the last linear grouped defect is just twice the dominant port’s transfer-shape slope.

### 7.2 Recovery of the Stage 244 square-root mixed-leg law

If the wall-normalized upper-leg shape and the normalized coupling shape are rigid,
\[
\mathfrak u_r=0,
\qquad
\mathfrak c_r=0,
\]
then
\[
\boxed{
\tau_r=\mathfrak w_r=
\mathfrak g_{W,r}-\mathfrak o_{W,r}-\frac12\kappa_1.
}
\]
So the present theorem collapses exactly to the Stage 244 square-root mixed-leg law,
\[
\frac{G_{W,r}}{\Omega_{W,r}^2}\propto \sqrt K.
\]
Thus Stage 247 is a strict generalization of Stage 244 rather than a competing condition.

### 7.3 Common normalized-leg co-scaling

A slightly more symmetric sufficient condition is
\[
\mathfrak c_r=0,
\qquad
\mathfrak u_r=\mathfrak w_r.
\]
Then the numerator shape
\(
\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}
\)
co-scales as one rigid object, so again
\[
\boxed{
\tau_r=\mathfrak w_r.
}
\]
Therefore the square-root mixed-leg law is recovered whenever the normalized upper and mixed numerator legs drift together while the normalized coupling shape stays rigid.

---

## 8. What Stage 247 changes

Stage 246 already reduced the remaining grouped weak-axisymmetric defect to the mismatch
between the outgoing-weighted static port slope and the wall-baseline slope,
\[
\Xi_1=\bar\nu_N-\kappa_1.
\]
Stage 247 sharpens that result one more step.

It shows that each actual static outgoing-transfer coefficient factors exactly as
\[
N_{0}^{(r)}=K\,\mathcal T_r^2,
\]
with one dimensionless wall-normalized transfer shape \(\mathcal T_r\).
Therefore the whole remaining theorem problem collapses to one question:

> are the wall-normalized transfer shapes of the actual moving-throat outgoing ports weak-axisymmetrically rigid?

Equivalently, the next honest theorem gate is no longer

> “compute the raw static outgoing-transfer slopes \(\nu_r\).”

It is simply

> compute the weak-axisymmetric transfer-shape slopes \(\tau_r\), and test whether their outgoing-weighted average vanishes.

That is the direct continuation point.
