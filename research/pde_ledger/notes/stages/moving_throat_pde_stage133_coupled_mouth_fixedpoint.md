
# Moving-Throat PDE — Stage 133: Full Coupled Mouth-Layer Fixed-Point Law

## Goal

Replace the effective parent datum
\[
\Pi_m=\frac{L}{\Theta_\sigma}
\left(
\left.\partial_z\delta V_{{\rm conf}}\right|_{\rm m}
-
q_*\left.\partial_zA_0\right|_{\rm m}
\right)
\]
by the first explicit **coupled mouth-layer solve** in which the mechanical and
localized-Maxwell throat channels are solved together and the source law
\(\sigma_\Pi\) closes self-consistently.

The point is not yet to solve the full moving-throat PDE. The point is to reduce
the last mouth-bias ambiguity to a small set of explicit dimensionless gains.

---

## 1. Dimensionless coupled mouth-layer system

Let
\[
x:=\frac{z}{L}\in[0,1],
\qquad
\Sigma_\Pi(x)=\frac{\Pi e^{-\Pi x}}{1-e^{-\Pi}}
\]
be the normalized mouth-source law from Stage 129.

Introduce the two-component boundary-layer field
\[
\mathbf U(x)=
\begin{pmatrix}
U_s(x)\\
U_q(x)
\end{pmatrix},
\]
where:

- \(U_s\) is the shell / confinement-response field,
- \(U_q\) is the localized-Maxwell / mixed-channel field.

The first coupled linear mouth-layer model is
\[
\boxed{
\left[-\partial_x^2\,\mathbf I+\mathbf K\right]\mathbf U(x)
=
\Sigma_\Pi(x)\,\mathbf G,
}
\]
with D/N-type mouth/bottom conditions
\[
\boxed{
\mathbf U(0)=0,
\qquad
\mathbf U'(1)=0.
}
\]

Here \(\mathbf K\) is a constant positive symmetric \(2\times2\) stiffness matrix
and \(\mathbf G\) is the source-coupling vector.

The parent mouth bias is read from a linear projection of the mouth derivative,
\[
\boxed{
V_1 = \mathbf H^{\!T}\mathbf U'(0),
\qquad
\Pi=\frac{L}{\Theta_\sigma}V_1,
}
\]
for some channel-weight vector \(\mathbf H\).

So the source parameter \(\Pi\) is no longer an input. It is the fixed point of a
coupled mechanical/electromagnetic boundary-layer problem.

---

## 2. Exact diagonalization

Diagonalize
\[
\mathbf K = \mathbf R
\begin{pmatrix}
\kappa_+^2 & 0\\
0 & \kappa_-^2
\end{pmatrix}
\mathbf R^{T},
\qquad
\mathbf R\in O(2).
\]
Define eigenbasis couplings
\[
\mathbf G_\pm := \mathbf R^T\mathbf G,
\qquad
\mathbf H_\pm := \mathbf R^T\mathbf H.
\]

Then each eigenchannel \(u_\alpha(x)\) satisfies
\[
\left(-\partial_x^2+\kappa_\alpha^2\right)u_\alpha(x)
=
G_\alpha\,\Sigma_\Pi(x),
\qquad
u_\alpha(0)=0,\quad
u_\alpha'(1)=0,
\qquad
\alpha\in\{+,-\}.
\]

So the fully coupled two-channel mouth-layer problem reduces exactly to two scalar
D/N response problems.

---

## 3. Exact scalar D/N response to the exponential source

For one channel with stiffness \(\kappa\), solve
\[
\left(-\partial_x^2+\kappa^2\right)u(x)=G\,\Sigma_\Pi(x),
\qquad
u(0)=0,\quad u'(1)=0.
\]

A direct elementary solve gives
\[
u(x)=A\sinh(\kappa x)-C\cosh(\kappa x)+C e^{-\Pi x},
\]
with
\[
C=\frac{G\Pi}{(1-e^{-\Pi})(\kappa^2-\Pi^2)},
\qquad
A=
\frac{C\left(\kappa\sinh\kappa+\Pi e^{-\Pi}\right)}{\kappa\cosh\kappa}.
\]

Therefore the exact mouth derivative is
\[
u'(0)=G\,\mathcal S(\Pi,\kappa),
\]
where
\[
\boxed{
\mathcal S(\Pi,\kappa)
=
\frac{\Pi\left[\kappa\tanh\kappa+\Pi\left(e^{-\Pi}\operatorname{sech}\kappa-1\right)\right]}
{(1-e^{-\Pi})(\kappa^2-\Pi^2)}.
}
\]

This is the exact D/N mouth-response kernel for the exponential source branch.

A useful exact limit is
\[
\boxed{
\mathcal S(\Pi,0)=1,
}
\]
so a purely static shell-compliance lane contributes a constant unit slope factor.

---

## 4. Full fixed-point law

Returning to the coupled system,
\[
V_1
=
\sum_{\alpha=\pm} H_\alpha u_\alpha'(0)
=
\sum_{\alpha=\pm} H_\alpha G_\alpha\,\mathcal S(\Pi,\kappa_\alpha).
\]

So the explicit self-consistency law for the mouth-source bias is
\[
\boxed{
\Pi
=
\sum_{\alpha=\pm}
M_\alpha\,\mathcal S(\Pi,\kappa_\alpha),
\qquad
M_\alpha:=\frac{L\,H_\alpha G_\alpha}{\Theta_\sigma}.
}
\]

This is the first honest coupled replacement for the Stage 129 effective parent
datum.

The open mouth-source problem has now shrunk to:

1. the two D/N eigen-stiffnesses \(\kappa_\pm\),
2. the two dimensionless channel gains \(M_\pm\).

Everything else is fixed by the exact source law and the exact D/N response.

---

## Result

The actual mouth-layer selection problem is no longer “what profile \(\sigma(z)\)?”
and no longer “what effective slope \(V_1\)?”

It is now the explicit coupled fixed-point equation
\[
\boxed{
\Pi
=
M_+\mathcal S(\Pi,\kappa_+)+M_-\mathcal S(\Pi,\kappa_-).
}
\]

So the remaining ambiguity is now just a **small dimensionless gain/stiffness
quadruple** \((M_+,M_-,\kappa_+,\kappa_-)\).
