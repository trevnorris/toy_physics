# Moving-Throat PDE — Stage 140: Self-Matched Mouth Susceptibility Closure

## Goal

Push the explicit Family-1 gain formulas one step further by replacing the remaining
source susceptibility `\Theta_\sigma` with the first self-matched mouth-layer closure.

This does **not** solve the full PDE. It gives the first explicit parent formula for the
overall shell gain scale `\Sigma_0=M_s`.

---

## 1. Self-matched mouth susceptibility

Take the mouth source to live on the same active shell layer as the shell/compliance mode.
Then the entropic/source susceptibility is the Stage-43 type quantity
\[
\Theta_\sigma = h'(\rho_w) N_{\sigma\sigma}^{(m)}.
\]
On the self-matched mouth layer,
\[
N_{\sigma\sigma}^{(m)} = J_s = \frac{4\pi a^2\ell}{3},
\qquad
h'(\rho_w)=H_w=\frac{m_\psi c_{s,w}^2}{\rho_w}.
\]
So
\[
\boxed{
\Theta_\sigma = H_w J_s.
}
\]
This is a same-layer closure: it removes an otherwise free susceptibility scale by
identifying the source channel with the already active shell layer, rather than
fitting a new branch-dependent coefficient by hand. If the actual mouth source
lives on a materially different layer, the numerical prefactor derived below would
change with that microscopic choice.

Using the Stage 118 shell coupling
\[
g_s=\mathcal T_m J_s,
\]
and the healing-locked shell stiffness
\[
K_s=\frac{3\pi a^2\hbar^2}{5m_\psi\rho_w\ell},
\]
Stage 138 gives
\[
\Sigma_0=M_s=\frac{L g_s^2}{K_s\Theta_\sigma}.
\]
A direct simplification yields
\[
\boxed{
\Sigma_0
=
\frac{20L\ell^2\rho_w^2\mathcal T_m^2}{9\hbar^2 c_{s,w}^2}.
}
\]

So the overall mouth shell gain is no longer abstract: it is fixed by one explicit
traction scale within this self-matched mouth-layer closure.

---

## 2. Convenient normalized traction parameter

Define the self-matched normalized mouth-traction amplitude
\[
\boxed{
\widehat T_m:=\frac{\rho_w\ell\sqrt{L}\,\mathcal T_m}{\hbar c_{s,w}}.
}
\]
Then
\[
\boxed{
\Sigma_0 = \frac{20}{9}\widehat T_m^2.
}
\]

So the full actual Family-1 gains are
\[
M_s = \frac{20}{9}\widehat T_m^2,
\qquad
M_q = -R_q\frac{20}{9}\widehat T_m^2.
\]

---

## 3. Required traction on the Family-1 branch

From Stage 139:

- natural equal-normalized branch requires
  \[
  M_s^{\rm nat,*}\approx 1.66854252965624,
  \qquad
  M_q^{\rm nat,*}\approx -0.242696939724365;
  \]
- exact compensated branch requires
  \[
  M_s^{\rm comp,*}\approx 1.80594111095636,
  \qquad
  M_q^{\rm comp,*}\approx -0.451485277739090.
  \]

Therefore the required normalized traction amplitudes are
\[
\boxed{
\widehat T_m^{\rm nat,*}
=\sqrt{\frac{9M_s^{\rm nat,*}}{20}}
\approx 0.866512630228382,
}
\]
\[
\boxed{
\widehat T_m^{\rm comp,*}
=\sqrt{\frac{9M_s^{\rm comp,*}}{20}}
\approx 0.901484054174206.
}
\]

So the exact outlet-consistent compensated branch requires only
\[
\boxed{
\frac{\widehat T_m^{\rm comp,*}}{\widehat T_m^{\rm nat,*}}-1
\approx 0.0403588161624,
}
\]
that is, only about a `4.04%` enhancement in the normalized mouth traction relative
to the natural equal-normalized branch.

---

## Result

Under the first self-matched mouth susceptibility closure, the overall shell gain is
fully explicit:
\[
\boxed{
M_s=\Sigma_0=\frac{20}{9}\widehat T_m^2.
}
\]

On the explicit Family-1 branch, the natural and exact-compensated mouth closures differ
by only about `4%` in the normalized traction amplitude.

That is the cleanest explicit parent-level narrowing of the mouth-gain problem so far.
