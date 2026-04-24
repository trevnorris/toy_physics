# Stage V2-08 — BdG–Wall Schur Complement and Stability/Softening Audit

## Purpose

This stage audits the first matter-support closure in the Volume 2 roadmap.

The target is the reduced **stable BdG normal-mode** block coupled to a moving wall/worldtube mode. The question is not whether the full nonlinear moving-throat PDE is solved. The narrower question is:

> If the linearized GNLS/BdG sector is reduced to positive-energy stable normal coordinates, does integrating those modes out produce a consistent conservative wall kernel, and what stability constraints does that impose?

The answer is a conditional **yes**:

\[
\boxed{
\text{Stable positive-norm BdG support modes give a valid Schur-complement softening of the wall.}
}
\]

But the audit also finds three important restrictions:

1. the closure requires positive-Krein/positive-energy support modes;
2. near-softening must preserve \(D_0>0\);
3. the low-frequency BdG moments are constrained Stieltjes moments, not freely tunable coefficients.

---

## 1. One wall mode coupled to one stable support mode

Start with the reduced quadratic Lagrangian

\[
L
=
\frac12 M\dot q^2
-\frac12 Kq^2
+\frac12 \dot X^2
-\frac12 \varpi^2X^2
+gqX.
\]

The Euler–Lagrange equations are

\[
M\ddot q+Kq-gX=0,
\]

\[
\ddot X+\varpi^2X-gq=0.
\]

Using the frequency convention \(e^{-i\omega t}\),

\[
(K-M\omega^2)q-gX=0,
\]

\[
(\varpi^2-\omega^2)X-gq=0.
\]

Solving the second equation gives

\[
X=\frac{g}{\varpi^2-\omega^2}q.
\]

Substitution gives the exact effective wall operator

\[
\boxed{
D_{\rm eff}(\omega)
=
K-M\omega^2
-
\frac{g^2}{\varpi^2-\omega^2}.
}
\]

This is the one-mode Schur complement.

---

## 2. Low-frequency expansion

Expanding for \(|\omega|<\varpi\),

\[
\frac{1}{\varpi^2-\omega^2}
=
\frac1{\varpi^2}
+
\frac{\omega^2}{\varpi^4}
+
\frac{\omega^4}{\varpi^6}
+
O(\omega^6).
\]

Therefore

\[
D_{\rm eff}(\omega)
=
\left(K-\frac{g^2}{\varpi^2}\right)
-
\omega^2
\left(M+\frac{g^2}{\varpi^4}\right)
-
\omega^4
\left(\frac{g^2}{\varpi^6}\right)
+
O(\omega^6).
\]

Define

\[
K_{\rm eff}=K-\frac{g^2}{\varpi^2},
\]

\[
M_{\rm eff}=M+\frac{g^2}{\varpi^4},
\]

\[
N_{\rm eff}=\frac{g^2}{\varpi^6}.
\]

Then

\[
\boxed{
D_{\rm eff}(\omega)
=
K_{\rm eff}
-
M_{\rm eff}\omega^2
-
N_{\rm eff}\omega^4
+
O(\omega^6).
}
\]

The support mode therefore:

\[
\boxed{
\text{lowers static stiffness, raises inertia, and adds a positive higher even moment.}
}
\]

This is exactly the desired conservative support behavior.

---

## 3. Static positivity and softening bound

The static potential is

\[
V(q,X)
=
\frac12Kq^2
+
\frac12\varpi^2X^2
-
gqX.
\]

Completing the square,

\[
V(q,X)
=
\frac12\varpi^2
\left(
X-\frac{g}{\varpi^2}q
\right)^2
+
\frac12
\left(
K-\frac{g^2}{\varpi^2}
\right)q^2.
\]

So the static stability gate is

\[
\boxed{
K_{\rm eff}=K-\frac{g^2}{\varpi^2}>0.
}
\]

Equivalently,

\[
\boxed{
K\varpi^2-g^2>0.
}
\]

This is the exact Schur-complement positivity condition for the static Hessian

\[
H_{\rm stat}
=
\begin{pmatrix}
K & -g\\
-g & \varpi^2
\end{pmatrix}.
\]

The determinant is

\[
\det H_{\rm stat}=K\varpi^2-g^2.
\]

So conservative support softening is allowed only up to the stability boundary. A branch that uses softening to raise downstream normalization must keep

\[
D_0=K_{\rm eff}>0.
\]

For many support modes in a one-wall lane,

\[
D_0
=
K-\sum_\alpha\frac{g_\alpha^2}{\varpi_\alpha^2}.
\]

It is useful to define

\[
\epsilon_B
=
\frac{1}{K}
\sum_\alpha
\frac{g_\alpha^2}{\varpi_\alpha^2}.
\]

Then

\[
\boxed{
D_0=K(1-\epsilon_B),
\qquad
\epsilon_B<1.
}
\]

---

## 4. Dynamic pole stability

The exact dispersion relation is

\[
(K-M\omega^2)(\varpi^2-\omega^2)-g^2=0.
\]

Let

\[
\Omega_\eta^2=\frac{K}{M}.
\]

Then the two squared frequencies are

\[
\boxed{
\omega_\pm^2
=
\frac{
\Omega_\eta^2+\varpi^2
\pm
\sqrt{(\Omega_\eta^2-\varpi^2)^2+4g^2/M}
}{2}.
}
\]

Both poles are real. They are both positive exactly when

\[
\Omega_\eta^2\varpi^2-\frac{g^2}{M}>0,
\]

which is again

\[
\boxed{
K\varpi^2-g^2>0.
}
\]

So the same Schur-complement condition protects both static and dynamic stability.

If the support mode lies above the wall mode,

\[
\varpi^2=\Omega_\eta^2+\Delta,
\qquad
\Delta>0,
\]

and

\[
h=\frac{g^2}{M},
\]

then the wall-like pole shifts as

\[
\boxed{
\omega_-^2
=
\Omega_\eta^2
-
\frac{h}{\Delta}
+
\frac{h^2}{\Delta^3}
+
O(h^3).
}
\]

So an above-wall support mode softens the wall-like pole at weak coupling.

---

## 5. Matrix Schur complement

For wall coordinates \(Q^A\) and stable support coordinates \(X_\alpha\),

\[
L
=
\frac12\dot Q^TM\dot Q
-
\frac12 Q^TKQ
+
\frac12\dot X^T\dot X
-
\frac12 X^T\Omega^2X
+
Q^TCX.
\]

Eliminating \(X\) gives

\[
\boxed{
D_{\rm eff}(\omega)
=
K-\omega^2M
-
C(\Omega^2-\omega^2I)^{-1}C^T.
}
\]

The low-frequency expansion is

\[
\boxed{
D_{\rm eff}(\omega)
=
K_{\rm eff}
-\omega^2M_{\rm eff}
-\omega^4N_{\rm eff}
+O(\omega^6),
}
\]

with

\[
\boxed{
K_{\rm eff}=K-C\Omega^{-2}C^T,
}
\]

\[
\boxed{
M_{\rm eff}=M+C\Omega^{-4}C^T,
}
\]

\[
\boxed{
N_{\rm eff}=C\Omega^{-6}C^T.
}
\]

The positivity gate is the block-Hessian condition

\[
\begin{pmatrix}
K & -C\\
-C^T & \Omega^2
\end{pmatrix}>0.
\]

Since \(\Omega^2>0\), this is equivalent to

\[
\boxed{
K-C\Omega^{-2}C^T>0.
}
\]

That is the general matrix stability gate.

---

## 6. Negative-Krein / ghost warning

The reduced closure only works for positive-energy stable support modes.

If a support coordinate has negative kinetic/potential sign,

\[
L_X
=
-\frac12\dot X^2
+
\frac12\varpi^2X^2
+
gqX,
\]

then the static Hessian is effectively

\[
H_{\rm ghost}
=
\begin{pmatrix}
K & -g\\
-g & -\varpi^2
\end{pmatrix}.
\]

Its determinant is

\[
\det H_{\rm ghost}
=
-K\varpi^2-g^2<0.
\]

Thus the static energy is indefinite for \(K>0,\varpi^2>0\). The SymPy audit marks this as a failure mode.

So the Volume 2 statement must be:

\[
\boxed{
\text{V2-08 passes only for the positive-Krein stable normal-mode sector.}
}
\]

A later full BdG/Krein audit should verify that the modes retained in the reduced branch really have positive energy.

---

## 7. Moment constraints: \(B_0,B_2,B_4\) are not free

For scalar lane data, define

\[
B_0=\sum_\alpha\frac{g_\alpha^2}{\lambda_\alpha},
\]

\[
B_2=\sum_\alpha\frac{g_\alpha^2}{\lambda_\alpha^2},
\]

\[
B_4=\sum_\alpha\frac{g_\alpha^2}{\lambda_\alpha^3},
\]

where

\[
\lambda_\alpha=\varpi_\alpha^2>0.
\]

For two modes with positive weights \(w_i=g_i^2\),

\[
B_0B_4-B_2^2
=
\boxed{
\frac{
w_1w_2(\lambda_1-\lambda_2)^2
}{
\lambda_1^3\lambda_2^3
}
\ge0.
}
\]

This is a Stieltjes/Hankel moment constraint. It means the stable BdG front end cannot arbitrarily fit the conservative even coefficients. Any proposed branch data must satisfy these positivity/moment relations.

In particular, equality occurs only when the support effectively has one spectral scale.

---

## 8. Grouped real \(P_2\) isotropy

For each grouped real lane \(A\in\{20,21,22\}\),

\[
D_A(\omega)
=
D_{A0}+D_{A2}\omega^2+D_{A4}\omega^4+\cdots.
\]

The normalized response coefficient is

\[
u_2^{(A)}
=
-\frac{D_{A2}}{D_{A0}}.
\]

Grouped trace and anisotropy variables are

\[
\bar u_2
=
\frac{
u_2^{(20)}
+2u_2^{(21)}
+2u_2^{(22)}
}{5},
\]

\[
a_2
=
\frac{
2u_2^{(20)}
-u_2^{(21)}
-u_2^{(22)}
}{10},
\]

\[
b_2
=
\frac{
u_2^{(21)}
-u_2^{(22)}
}{2}.
\]

If all grouped lanes share identical stable support data,

\[
D_{20,n}=D_{21,n}=D_{22,n},
\]

then

\[
\boxed{
a_2=b_2=0.
}
\]

So the stable BdG Schur complement preserves grouped \(P_2\) isotropy if the support spectrum and overlap data are isotropic.

---

## 9. Verdict

The V2-08 audit passes as a **controlled reduced-sector gate**.

\[
\boxed{
\text{The Schur complement is algebraically exact and physically stable under positive-mode conditions.}
}
\]

The gate does **not** prove the full moving-throat branch. It establishes the conditions that any branch realization must satisfy:

\[
\boxed{
\Omega^2>0,\qquad
M>0,\qquad
K-C\Omega^{-2}C^T>0.
}
\]

It also establishes that downstream coefficient matching cannot treat the BdG moments as free:

\[
\boxed{
B_0B_4-B_2^2\ge0
}
\]

in the scalar two-mode prototype, with the corresponding positive-moment constraints in the multimode case.

## Carry-forward patch

Use this wording in Volume 2:

> The BdG-wall reduced closure is valid after projecting the linearized GNLS/BdG system onto positive-energy stable normal coordinates. Integrating those modes out gives an exact Schur-complement wall kernel. Stable support modes lower static stiffness, raise effective inertia, and generate positive higher even moments. The branch must satisfy \(K-C\Omega^{-2}C^T>0\); near-softened normalization enhancement is allowed only below that stability boundary. Negative-Krein modes are outside this closure and require a separate full BdG signature audit.
