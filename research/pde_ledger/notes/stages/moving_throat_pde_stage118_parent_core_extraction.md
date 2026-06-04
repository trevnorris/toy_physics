
# Moving-Throat PDE — Stage 118: Parent-Action Extraction of Core Parameters

## Goal

Replace the reduced core parameters
\[
(K_s,\;K_q,\;\lambda,\;g_s,\;g_q)
\]
by explicit overlap formulas from one concrete **GNLS + localized-Maxwell throat-core ansatz**.

The point is not yet to prove the full moving-throat PDE. The point is to reduce the surviving outlet-core ambiguity to parent-action quantities that can be computed on an explicit throat branch.

## Parent sectors used

The parent 4D stack already gives:

- a gauged 4D GNLS matter sector with confinement \(V_{\rm conf}(\mathbf X;a,L)\) and frozen \(n=5\) EOS, so the bulk matter Hamiltonian has both gradient and compressional terms, fileciteturn33file4
- a localized \(4+1\) Maxwell sector with real mixed channels \((A_w,F_{\mu w},J^w)\) outside the strict far-field zero-mode reduction, fileciteturn33file1turn33file3
- and explicit mixed-sector EM energy channels such as \(E_w\), \(C_a=F_{aw}\), and \(S_{\rm EM}^w\), so a concrete throat-core mixed outlet is on-model rather than ad hoc. fileciteturn34file7turn34file12

## Core ansatz

Take one grouped-real \(P_2\) lane and introduce two core amplitudes:

- \(s(t)\): shell / compliance mode,
- \(q(t)\): mixed \(A_w/F_{\mu w}\)-type side-channel mode.

Use the factorized parent ansatz
\[
\delta\rho(\mathbf X,t)=s(t)\,\varrho_s(\mathbf X),
\qquad
\delta A_w(\mathbf X,t)=q(t)\,\mathcal A_q(\mathbf X),
\]
with all other fluctuating fields omitted at this first closure level.

For the shell mode, use the carried thin-wall GNLS profile
\[
\varrho_s(\mathbf X)=\rho_w\,\chi_s(y)\,Y_{2m}(\Omega),
\qquad
y:=\frac{r-a}{\ell},
\qquad
\chi_s(y)=f'(y),
\]
and on the canonical tanh wall,
\[
f(y)=\frac{1+\tanh y}{2},
\qquad
\chi_s(y)=\frac12\operatorname{sech}^2 y.
\]

For the mixed side-channel, use the first finite D/N half-wave on an auxiliary tube \(z\in[0,L_W]\):
\[
\mathcal A_q(\mathbf X)=\Xi_q(\Omega,w)\,\chi_{1/2}(z),
\qquad
\chi_{1/2}(z)=\sqrt{\frac{2}{L_W}}\sin\!\left(\frac{\pi z}{2L_W}\right).
\]

## 1. Shell stiffness \(K_s\)

Expanding the GNLS matter Hamiltonian to quadratic order in the shell mode gives
\[
\delta^2 H_s=\frac12 K_s s^2,
\]
with
\[
\boxed{
K_s
=
\int d^4X\left[
H_w\,\chi_s^2
+\frac{\hbar^2}{4m_\psi\rho_w}\,|\nabla\chi_s|^2
\right],
}
\]
where
\[
H_w=h'(\rho_w)=\frac{m_\psi c_{s,w}^2}{\rho_w}
\]
is the local compressional stiffness.

On the canonical thin wall, using the already carried moments
\[
I_f=\int_{-\infty}^{\infty}\chi_s^2\,dy=\frac13,
\qquad
I_g=\int_{-\infty}^{\infty}(\chi_s')^2\,dy=\frac4{15},
\]
this becomes
\[
\boxed{
K_s
=
4\pi a^2\left(
\frac{H_w\ell}{3}
+
\frac{\hbar^2}{15m_\psi\rho_w\,\ell}
\right).
}
\]

If the wall thickness is locked by the static GNLS healing relation
\[
\ell=\frac{\hbar}{2m_\psi c_{s,w}},
\]
then
\[
\boxed{
K_s
=
\frac{3\pi a^2\hbar^2}{5m_\psi\rho_w\,\ell}.
}
\]

## 2. Mixed side-channel stiffness \(K_q\)

Reducing the localized Maxwell action on the D/N half-wave gives
\[
\delta^2 H_q=\frac12 K_q q^2
\]
with
\[
\boxed{
K_q
=
\frac{\mathcal Z_q}{\mu_0}\int_0^{L_W} (\chi'_{1/2})^2\,dz,
}
\]
where
\[
\mathcal Z_q:=\int d\Omega\,dw\,Z(w)\,\Xi_q(\Omega,w)^2
\]
is the transverse localization norm of the mixed channel.

Since
\[
\int_0^{L_W}\chi_{1/2}^2\,dz=1,
\qquad
\int_0^{L_W}(\chi'_{1/2})^2\,dz=\frac{\pi^2}{4L_W^2},
\]
one gets
\[
\boxed{
K_q
=
\frac{\mathcal Z_q}{\mu_0}\,\frac{\pi^2 c_s^2}{4L_W^2}.
}
\]

So \(K_q\) is not arbitrary: it is fixed by one localization norm and the D/N tube length.

## 3. Static/mixed hybridization \(\lambda\)

The GNLS kinetic energy in Madelung form contains
\[
\mathcal H_{\rm kin}=\frac{m_\psi}{2}\rho\,v_Av_A.
\]
Take a static background mixed flow \(v_{w0}\neq0\), hold the phase fixed in the \(q\)-channel so that
\[
\delta v_w = -\frac{q_*}{m_\psi}\,\delta A_w,
\]
and expand
\[
\rho=\rho_0+s\,\varrho_s,
\qquad
A_w=A_{w0}+q\,\mathcal A_q.
\]
The bilinear \(sq\) term is exactly
\[
\delta^2 H_{sq}
=
- q_* v_{w0}\int d^4X\,\varrho_s\,\mathcal A_q\;s q.
\]
So the concrete hybridization is
\[
\boxed{
\lambda
=
- q_* v_{w0}\,\mathcal I_{sq},
\qquad
\mathcal I_{sq}:=\int d^4X\,\varrho_s\,\mathcal A_q.
}
\]

In the simplest uniform-core overlap closure,
\[
\mathcal I_{sq}=J_s I_q,
\qquad
J_s:=4\pi a^2\ell I_f=\frac{4\pi a^2\ell}{3},
\qquad
I_q:=\int_0^{L_W}\chi_{1/2}(z)\,dz=\frac{2\sqrt{2L_W}}{\pi},
\]
so
\[
\boxed{
\lambda
=
-\frac{8\sqrt2}{3}\,q_* v_{w0}\,a^2\ell\,\sqrt{L_W}.
}
\]

This is the first explicit parent-level source of the Stage-114 bilinear shell/mixed coupling.

## 4. Mouth couplings \(g_s\) and \(g_q\)

Use one external mouth amplitude \(u\) and let it couple to:

- shell traction at the mouth,
- mixed flux at the mouth.

### Shell coupling

If the mouth forcing is a lip traction \(\mathcal T_m u\), then the shell overlap gives
\[
\delta H_{\rm mouth}^{(s)}=-u g_s s,
\qquad
\boxed{
g_s=\mathcal T_m J_s=\mathcal T_m \frac{4\pi a^2\ell}{3}.
}
\]

### Mixed coupling

For the D/N mixed tube, the natural mouth observable is the derivative at the Dirichlet mouth:
\[
\delta H_{\rm mouth}^{(q)}=-u g_q q,
\qquad
g_q=\frac{\mathcal Z_q}{\mu_0}\chi'_{1/2}(0).
\]
Since
\[
\chi'_{1/2}(0)=\frac{\pi}{\sqrt2\,L_W^{3/2}},
\]
one gets
\[
\boxed{
g_q=\frac{\mathcal Z_q}{\mu_0}\,\frac{\pi}{\sqrt2\,L_W^{3/2}}.
}
\]

So every parameter in the Stage-114 core matrix now has a concrete parent-action interpretation.

## Result

The effective core matrix
\[
\begin{pmatrix}
K_s & \lambda\\
\lambda & -K_q D_W^{\rm bare}(z)
\end{pmatrix}
\binom{s}{q}
=
u\binom{g_s}{g_q}
\]
is no longer just a reduced model.

In the first explicit throat-core closure its entries are:

\[
\boxed{
K_s
=
4\pi a^2\left(
\frac{H_w\ell}{3}
+
\frac{\hbar^2}{15m_\psi\rho_w\,\ell}
\right),
}
\]
\[
\boxed{
K_q
=
\frac{\mathcal Z_q}{\mu_0}\,\frac{\pi^2 c_s^2}{4L_W^2},
}
\]
\[
\boxed{
\lambda
=
- q_* v_{w0}\,\mathcal I_{sq},
}
\]
\[
\boxed{
g_s=\mathcal T_m \frac{4\pi a^2\ell}{3},
\qquad
g_q=\frac{\mathcal Z_q}{\mu_0}\,\frac{\pi}{\sqrt2\,L_W^{3/2}}.
}
\]

The next step is to rewrite the compensation surface entirely in terms of the parent overlap ratios these formulas define.
