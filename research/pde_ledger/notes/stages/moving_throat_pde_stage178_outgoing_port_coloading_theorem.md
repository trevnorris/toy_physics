# Moving-Throat PDE — Stage 178: Outgoing-Port Co-Loading Theorem and Direct Portwise Formula for `\Xi_1`

## Purpose

Stage 245 reduced the remaining weak-axisymmetric grouped `2.5`PN defect to one scalar,
\[
\Xi_1=\frac{P_1}{P_0},
\]
and expressed it in the port-slippage variables
\[
\mathfrak m_r,\qquad \mathfrak i_r,\qquad \mathfrak h_r.
\]

That was already a major collapse, but it still left the final theorem gate looking more “internal” than “physical.”
The next honest step is therefore:

> rewrite the remaining defect directly in terms of the weak-axisymmetric slopes of the **actual outgoing port data**
> \[
> P_r,\qquad \Delta_r,\qquad N_{0}^{(r)}=\frac{P_r^2}{\Delta_r^2},
> \]
> and identify the exact co-loading law those ports must satisfy relative to the wall baseline.

This stage does exactly that.

The main result is that on the same conservative-shape-preserving branch used in Stages 226–228, the full remaining grouped defect is
\[
\boxed{
\Xi_1
=
\sum_r \rho_r^{(N)}\bigl(\nu_r-\kappa_1\bigr)
=
\bar\nu_N-\kappa_1,
}
\]
where

- \(\kappa_1\) is the weak-axisymmetric wall-baseline slope,
- \(\nu_r\) is the weak-axisymmetric logarithmic slope of the **actual static outgoing-transfer coefficient**
  \[
  N_{A,0}^{(r)}=\frac{P_{A,r}^2}{\Delta_{A,r}^2},
  \]
- and
  \[
  \bar\nu_N:=\sum_r \rho_r^{(N)}\,\nu_r=\frac{N_{01}}{N_0}
  \]
  is the outgoing-weighted average of those port slopes.

So the remaining linear grouped `2.5`PN theorem problem is no longer “compute every microscopic outgoing slippage.”
It is just:

> does the weighted static outgoing-port load slope co-load with the wall baseline?

That is the direct port-level version of the Stage 245 result.

---

## 1. Weak-axisymmetric slopes of the actual outgoing-port data

Write the weak-axisymmetric grouped branch as usual,
\[
\lambda_{20}=1,
\qquad
\lambda_{21}=\frac12,
\qquad
\lambda_{22}=-1,
\]
and define the wall-baseline slope by
\[
\boxed{
\delta\ln K_A=\epsilon\lambda_A\,\kappa_1.
}
\]

For each outgoing port \(r\), define the weak-axisymmetric logarithmic slopes of the **actual static numerator** \(P_r\) and **actual static detuning** \(\Delta_r\):
\[
\boxed{
\delta\ln P_{A,r}=\epsilon\lambda_A\,\mathfrak p_r,
\qquad
\delta\ln \Delta_{A,r}=\epsilon\lambda_A\,\mathfrak d_r.
}
\]

Now the actual static outgoing-transfer coefficient of that port is
\[
N_{A,0}^{(r)}=\frac{P_{A,r}^2}{\Delta_{A,r}^2}.
\]
So its weak-axisymmetric logarithmic slope is
\[
\boxed{
\delta\ln N_{A,0}^{(r)}=\epsilon\lambda_A\,\nu_r,
\qquad
\nu_r:=2(\mathfrak p_r-\mathfrak d_r).
}
\]

This is the first exact portwise identity of the stage.
It says the static outgoing-transfer slope is simply the doubled numerator slope minus the doubled detuning slope.

---

## 2. Exact portwise rewrite of the remaining grouped defect

Stage 243 already showed that on conservative-shape-preserving branches the remaining grouped defect is carried only by the outgoing load factor:
\[
\Xi_{\rm load}^{(A)}
=
\sum_r \rho_r^{(N)}\,
\delta\ln\!\left(\frac{\Lambda_{A,r}^2}{K_A}\right),
\qquad
\Lambda_{A,r}=\frac{P_{A,r}}{\Delta_{A,r}}.
\]
But \(\Lambda_{A,r}^2=N_{A,0}^{(r)}\), so
\[
\delta\ln\!\left(\frac{\Lambda_{A,r}^2}{K_A}\right)
=
\delta\ln N_{A,0}^{(r)}-\delta\ln K_A.
\]
Using the weak-axisymmetric definitions above gives
\[
\boxed{
\delta\ln\!\left(\frac{\Lambda_{A,r}^2}{K_A}\right)
=
\epsilon\lambda_A\,(\nu_r-\kappa_1).
}
\]

Therefore the full remaining grouped defect collapses exactly to
\[
\boxed{
\Xi_{\rm load}^{(A)}
=
\epsilon\lambda_A\,\Xi_1,
\qquad
\Xi_1
=
\sum_r \rho_r^{(N)}(\nu_r-\kappa_1).
}
\]

Because the outgoing weights satisfy
\[
\sum_r \rho_r^{(N)}=1,
\]
this can be rewritten as
\[
\boxed{
\Xi_1
=
\bar\nu_N-\kappa_1,
\qquad
\bar\nu_N:=\sum_r \rho_r^{(N)}\,\nu_r.
}
\]

And since Stage 241/228 already identified the same scalar with the physical outgoing-prefactor slope,
\[
\boxed{
\frac{P_1}{P_0}=\Xi_1=\bar\nu_N-\kappa_1.
}
\]

So the whole remaining linear grouped normalization defect is exactly the mismatch between

- the outgoing-weighted static transfer slope \(\bar\nu_N\),
- and the wall-baseline slope \(\kappa_1\).

This is the stage’s main theorem.

---

## 3. Exact formulas for the actual port slopes `\mathfrak p_r` and `\mathfrak d_r`

The static outgoing numerator and detuning are, lane by lane,
\[
P_r=\Omega_{U,r}^2 G_{W,r}+R_r G_{U,r},
\]
\[
\Delta_r=\Omega_{U,r}^2\Omega_{W,r}^2-R_r^2.
\]

### 3.1 Numerator slope

Define the positive numerator weights
\[
\boxed{
\alpha_r:=\frac{\Omega_{U,r}^2 G_{W,r}}{P_r},
\qquad
\beta_r:=\frac{R_r G_{U,r}}{P_r},
\qquad
\alpha_r+\beta_r=1.
}
\]
Then the weak-axisymmetric numerator slope is
\[
\boxed{
\mathfrak p_r
=
\alpha_r(\mathfrak o_{U,r}+\mathfrak g_{W,r})
+\beta_r(\mathfrak r_r+\mathfrak g_{U,r}).
}
\]

So \(\mathfrak p_r\) is the convex average of the two static numerator channels:

- the brane-like contribution \(\Omega_{U,r}^2 G_{W,r}\),
- the mixed contribution \(R_r G_{U,r}\).

### 3.2 Detuning slope

Define the detuning weights
\[
\boxed{
\chi_r:=\frac{\Omega_{U,r}^2\Omega_{W,r}^2}{\Delta_r}
=\frac{1}{1-\mathcal H_r},
\qquad
\zeta_r:=\frac{R_r^2}{\Delta_r}
=\frac{\mathcal H_r}{1-\mathcal H_r}.
}
\]
Then the weak-axisymmetric detuning slope is
\[
\boxed{
\mathfrak d_r
=
\chi_r(\mathfrak o_{U,r}+\mathfrak o_{W,r})
-
2\zeta_r\,\mathfrak r_r.
}
\]

So \(\mathfrak d_r\) measures how the static port detuning is changed by:

- the combined frequency drift of the internal brane-like and mixed legs,
- and the coupling drift \(R_r\).

### 3.3 Static outgoing-transfer slope in actual port variables

Combining the last two formulas with \(\nu_r=2(\mathfrak p_r-\mathfrak d_r)\) gives
\[
\boxed{
\nu_r
=
2\alpha_r(\mathfrak o_{U,r}+\mathfrak g_{W,r})
+
2\beta_r(\mathfrak r_r+\mathfrak g_{U,r})
-
2\chi_r(\mathfrak o_{U,r}+\mathfrak o_{W,r})
+
4\zeta_r\,\mathfrak r_r.
}
\]

This is the direct expression of the static outgoing-transfer slope in the real moving-throat port variables.

---

## 4. Exact link back to the Stage 245 slippage formula

Stage 245 used the port slippages
\[
\mathfrak m_r:=\mathfrak g_{W,r}-\mathfrak o_{W,r}-\frac12\kappa_1,
\]
\[
\mathfrak i_r:=\mathfrak r_r+\mathfrak g_{U,r}-\mathfrak o_{U,r}-\mathfrak g_{W,r},
\]
\[
\mathfrak h_r:=2\mathfrak r_r-\mathfrak o_{U,r}-\mathfrak o_{W,r},
\]
together with
\[
\mathcal I_r=\frac{R_r G_{U,r}}{\Omega_{U,r}^2 G_{W,r}},
\qquad
\mathcal H_r=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2}.
\]

Since
\[
\alpha_r=\frac{1}{1+\mathcal I_r},
\qquad
\beta_r=\frac{\mathcal I_r}{1+\mathcal I_r},
\qquad
\chi_r=\frac{1}{1-\mathcal H_r},
\qquad
\zeta_r=\frac{\mathcal H_r}{1-\mathcal H_r},
\]
the port-slope formula above is exactly equivalent to
\[
\boxed{
\nu_r
=
\kappa_1
+
2\,\mathfrak m_r
+\frac{2\mathcal I_r}{1+\mathcal I_r}\,\mathfrak i_r
+\frac{2\mathcal H_r}{1-\mathcal H_r}\,\mathfrak h_r.
}
\]
So the Stage 245 port amplitude \(\sigma_r\) is simply
\[
\boxed{
\sigma_r=\nu_r-\kappa_1.
}
\]

Therefore
\[
\boxed{
\Xi_1
=
\sum_r \rho_r^{(N)}\,\sigma_r
=
\sum_r \rho_r^{(N)}(\nu_r-\kappa_1)
=
\bar\nu_N-\kappa_1.
}
\]

This is the exact equivalence between the Stage 245 slippage language and the present direct outgoing-port language.

---

## 5. Outgoing-port co-loading theorem

The previous formulas make the real theorem gate extremely sharp.

### 5.1 Exact zero-defect condition

The remaining linear grouped defect vanishes if and only if
\[
\boxed{
\bar\nu_N=\kappa_1.
}
\]
That is, the outgoing-weighted static transfer slope must match the wall-baseline slope.

Equivalently,
\[
\boxed{
\frac{P_1}{P_0}=0
\qquad\Longleftrightarrow\qquad
\bar\nu_N=\kappa_1.
}
\]

### 5.2 Strong per-port sufficient condition

A stronger sufficient condition is
\[
\boxed{
\nu_r=\kappa_1
\quad\text{for every active outgoing port }r.
}
\]
Because \(\nu_r=\delta\ln N_{A,0}^{(r)}/(\epsilon\lambda_A)\), this is equivalent to
\[
\boxed{
N_{A,0}^{(r)}\propto K_A
\quad\text{lane by lane for every active }r.
}
\]
So if each static outgoing-transfer coefficient co-loads with the wall baseline individually, then
\[
\Xi_1=0.
\]

This is the direct outgoing-port version of the earlier wall-referenced self-similarity theorem.

### 5.3 Dominant-port limit

If one outgoing port dominates,
\[
\rho_{r_*}^{(N)}\approx1,
\qquad
\rho_{r\neq r_*}^{(N)}\approx0,
\]
then
\[
\boxed{
\Xi_1\approx \nu_{r_*}-\kappa_1
=2(\mathfrak p_{r_*}-\mathfrak d_{r_*})-\kappa_1.
}
\]
So in the dominant-port regime the last linear grouped defect is just the mismatch between

- the dominant port’s static outgoing-transfer slope,
- and the wall-baseline slope.

### 5.4 Naive rigidity is not enough

If the actual outgoing ports are rigid,
\[
\mathfrak p_r=\mathfrak d_r=0
\qquad\Longrightarrow\qquad
\nu_r=0,
\]
then
\[
\boxed{
\Xi_1=-\kappa_1.
}
\]
So the zero-defect condition is **not** achieved by freezing the outgoing ports.
They must actively co-load with the wall baseline.
This is the direct outgoing-port restatement of the earlier “naive common self-similarity fails” theorem.

---

## 6. What Stage 246 changes

Stage 245 already showed that the remaining weak-axisymmetric grouped defect is one scalar,
\[
\Xi_1=\frac{P_1}{P_0},
\]
and that its microscopic content is the weighted slippage combination
\[
\sum_r \rho_r^{(N)}
\left[
2\,\mathfrak m_r
+\frac{2\mathcal I_r}{1+\mathcal I_r}\,\mathfrak i_r
+\frac{2\mathcal H_r}{1-\mathcal H_r}\,\mathfrak h_r
\right].
\]

Stage 246 sharpens that result in a more physical language.

It shows that the same scalar is exactly the mismatch between

- the outgoing-weighted static transfer slope of the **actual** moving-throat outgoing ports,
- and the wall-baseline slope.

So the next honest theorem gate is now smaller again.

It is no longer

> “compute the microscopic slippage bundle.”

It is simply:

> compute the weak-axisymmetric static outgoing-transfer slopes
> \[
> \nu_r=\frac{\delta\ln N_{A,0}^{(r)}}{\epsilon\lambda_A}
> \]
> of the real moving-throat outgoing ports, and check whether their outgoing-weighted average equals \(\kappa_1\).

That is the direct continuation point.
