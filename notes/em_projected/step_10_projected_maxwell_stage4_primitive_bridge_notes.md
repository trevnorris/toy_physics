# Projected Maxwell Primitive One-Port Bridge Notes

## What this adds beyond the grouped-bundle bridge

The first bridge note worked at the **bundle-moment** level:
\[
Z_n\to Z_n+\varepsilon z_n,\qquad
N_n\to N_n+\varepsilon n_n.
\]

That was enough to show where projected-Maxwell near-throat physics enters the
moving-throat PDE targets, but it still treated the moment shifts
\[
z_0,\ z_2,\ z_4,\ n_0,\ n_2,\ n_4
\]
as abstract.

This follow-up goes one step deeper and starts from the **Stage-4/Stage-5 one-port
primitive data**:
\[
Z_0=\frac{Q}{\Delta},\qquad
Z_2=\frac{QS_2-H\Delta}{\Delta^2},\qquad
Z_4=\frac{Q(S_2^2-\Delta)-S_2H\Delta}{\Delta^3},
\]
\[
N_0=\frac{P^2}{\Delta^2},
\]
\[
N_2=\frac{2P(PS_2-\Delta G_W)}{\Delta^3},
\]
\[
N_4=\frac{\Delta^2G_W^2-2\Delta P^2-4\Delta PS_2G_W+3P^2S_2^2}{\Delta^4}.
\]

Then it introduces a **mouth-local projected-Maxwell slippage packet**
\[
(Q,S_2,H,\Delta,P,G_W)\to
(Q+\ell q_1,\ S_2+\ell s_1,\ H+\ell h_1,\ \Delta+\ell d_1,\ P+\ell p_1,\ G_W+\ell g_1),
\]
and computes the induced
\[
z_0,\ z_2,\ z_4,\ n_0,\ n_2,\ n_4
\]
exactly at first order.

So this note is the first real answer to the question:

> If the projected-Maxwell near-throat correction is missing from the PDE, which
> microscopic one-port quantities would have to move to repair the bundle data?

---

## 1. The induced primitive-to-bundle map

The SymPy derivation gives

\[
z_0=\frac{\Delta q_1-Qd_1}{\Delta^2},
\]

\[
z_2=
\frac{-\Delta^2 h_1+\Delta(Hd_1+Qs_1+S_2q_1)-2QS_2d_1}{\Delta^3},
\]

\[
z_4=
\frac{-\Delta^2Hs_1-\Delta^2S_2h_1-\Delta^2q_1
+2\Delta HS_2d_1+2\Delta QS_2s_1+2\Delta Qd_1
+\Delta S_2^2q_1-3QS_2^2d_1}{\Delta^4},
\]

\[
n_0=\frac{2P(\Delta p_1-Pd_1)}{\Delta^3}.
\]

Two immediate structural points follow.

First,
\[
z_0
\]
depends only on the primitive conservative source slippage \(q_1\) and
primitive denominator slippage \(d_1\).

Second,
\[
n_0
\]
depends only on the primitive transfer-leg slippage \(p_1\) and the same
denominator slippage \(d_1\).

So the first isotropic normalization mover is already very constrained.

---

## 2. Static prefactor transport and \(\Xi_1\)

The grouped-bundle bridge showed that a weak-axisymmetric projected-Maxwell static
slippage contributes to the actual prefactor bottleneck by
\[
\Xi_1^{(\mathrm{proj})}=\frac{n_0}{N_0}+\frac{z_0}{D_0}.
\]

After substituting the primitive one-port formulas, this becomes
\[
\Xi_1^{(\mathrm{primitive})}
=
2\frac{p_1}{P}
-2\frac{d_1}{\Delta}
+\frac{q_1}{D_0\Delta}
-\frac{Qd_1}{D_0\Delta^2}.
\]

That is useful because it splits the static projected-Maxwell contribution into
three clearly different microscopic pieces:

1. **transfer-leg slippage** \(p_1/P\),
2. **denominator slippage** \(d_1/\Delta\),
3. **direct conservative source slippage** \(q_1\) through \(z_0/D_0\).

So if the actual weak-axisymmetric bottleneck is
\[
\Xi_1=\frac{P_1}{P_0},
\]
this one-port bridge says the first projected-Maxwell place to look is **not**
all of Maxwell space. It is the small primitive packet
\[
(q_1,\ d_1,\ p_1).
\]

That is much sharper than “include more mixed-sector electromagnetism somehow.”

---

## 3. Isotropic compatibility transport

The grouped-bundle bridge also showed that the isotropic compatibility shift
\[
\delta\mathcal C
=
\frac{n_0}{P_{0,\mathrm{target}}}
-\frac{6S z_2}{T}
+\frac{3S^2 z_4}{T^2}
\]
depends only on
\[
n_0,\ z_2,\ z_4,
\]
so it is blind to \(z_0\).

The primitive substitution confirms that the first one-port contributors are then

- \(q_1\),
- \(s_1\),
- \(h_1\),
- \(d_1\),
- \(p_1\),

while the explicit mixed-channel slope
\[
g_1
\]
does **not** appear until \(n_2\) and \(n_4\), meaning it first matters for the
constant-prefactor branch transport rather than for the isotropic compatibility
surface itself.
The updated script now reconstructs \(\delta\mathcal C\) in two independent
ways after primitive substitution:

1. from the competing one-pole and normalization \(K\)-surfaces,
2. from the directly eliminated compatibility transport
   \[
   \mathcal C(\ell)=
   \frac{N_0+\ell n_0}{P_{0,\mathrm{target}}}
   -
   \frac{3(S+\ell z_2)^2}{T+\ell z_4}.
   \]

The script verifies that these two routes agree exactly at first order. Since
the eliminated surface itself carries no \(z_0\) slot, the loss of the \(z_0\)
channel is now a genuine algebraic cancellation rather than an omitted symbol.
It also introduces an independent \(z_0\) probe and checks that both
compatibility shifts are insensitive to that probe, while the normalization
\(K\)-surface remains explicitly sensitive before compatibility elimination.

It now also checks the stronger ratio-form target transport
\[
P_{0,\mathrm{target}}(\ell)
=
\frac{N_0+\ell n_0}{D_{0,\mathrm{target}}},
\]
with an explicit transported denominator slot \(D_{0,\mathrm{target}}\). In
that variant the normalization solve returns
\[
K=B_0+Z_0+\ell z_0 + D_{0,\mathrm{target}},
\]
and the resulting compatibility shift loses both the \(z_0\) and \(n_0\)
channels, leaving only the one-pole geometry transport from \(z_2\) and
\(z_4\).

That gives a useful hierarchy:

- **static normalization movers:** \(q_1,d_1,p_1\),
- **one-pole / compatibility movers:** \(q_1,s_1,h_1,d_1,p_1\),
- **constant-prefactor movers:** \(q_1,s_1,h_1,d_1,p_1,g_1\).

---

## 4. Why this is helpful

This one-port primitive bridge does not solve the moving-throat PDE. But it tells us
what kind of **concrete throat-local projected ansatz** is worth trying next.

If the goal is the weak-axisymmetric prefactor slope \(\Xi_1\), the first good
ansatz should expose:

- the primitive transfer slippage \(p_1\),
- the denominator slippage \(d_1\),
- and the conservative source slippage \(q_1\).

If the goal is the isotropic compatibility surface, the first good ansatz must also
resolve:

- \(s_1\), which changes the \(S_2\) combination,
- \(h_1\), which changes the conservative \(H\)-slot.

If the goal is the constant-prefactor outgoing branch, then the first good ansatz
must also resolve:

- \(g_1\), the primitive slope of the explicit mixed outgoing leg.

So the next move after this note is not “invent more algebra.” It is to choose a
projected near-mouth throat ansatz rich enough to produce these specific primitive
slippages.
