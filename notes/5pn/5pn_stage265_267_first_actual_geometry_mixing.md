# 5PN continuation — Stages 265–267

## Goal

Take the Stage-78 weak-anisotropy statement

\[
\epsilon_2,\epsilon_4 = O(\chi^2)
\]

and turn it into the first concrete induced-geometry mixing law for the actual grouped-`P2` branch. The point is to stop talking about a generic `O(chi^2)` correction and instead extract the exact coefficients

\[
\epsilon_2 = a_2 \chi^2,
\qquad
\epsilon_4 = a_4 \chi^2,
\]

for the first explicit `l=0 \leftrightarrow l=2` mechanism, then test those coefficients against the exact twin-safe and Family-1 corridor conditions from the later grouped-`P2` support analysis.

---

## Stage 265 — exact coefficient extraction from the Stage-78 Schur complement

Start from the exact Stage-78 reduced mixed scalar-geometry / grouped-`P2` model,

\[
D_q(\omega)=K_{\rm stat}+\frac{K_{\rm pole}}{1-\omega^2/\Omega_Q^2},
\qquad
D_g(\omega)=G_0+G_2\omega^2+G_4\omega^4+O(\omega^6),
\]

with weak anisotropy generating a bilinear mixing term

\[
\chi M_0\,qg.
\]

Integrating out the scalar/geometry lane gives

\[
D_{\rm eff}(\omega)
=
D_q(\omega)-\frac{\chi^2 M_0^2}{D_g(\omega)}.
\]

Expanding through `O(omega^4)` yields the exact contamination moments

\[
K_{(g,2)}^{\rm eff} = \chi^2\frac{M_0^2 G_2}{G_0^2},
\qquad
K_{(g,4)}^{\rm eff} = \chi^2\frac{M_0^2(G_0G_4-G_2^2)}{G_0^3}.
\]

Therefore the grouped-`P2` obstruction coefficients are exactly

\[
a_2 = \frac{\Omega_Q^2 M_0^2 G_2}{K_{\rm pole} G_0^2},
\qquad
a_4 = \frac{\Omega_Q^4 M_0^2 (G_0G_4-G_2^2)}{K_{\rm pole} G_0^3}.
\]

So the Stage-78 `O(chi^2)` statement already has a fully explicit coefficient-level form.

### First concrete scalar-lane specialization

Take the scalar/geometry lane itself to be one effective contact plus one scalar pole,

\[
D_g(\omega)=G_c+\frac{G_p}{1-\omega^2/\Omega_g^2}.
\]

Define

\[
G_0 := G_c+G_p,
\qquad
c_g := \frac{G_p}{G_0},
\qquad
r := \frac{\Omega_Q^2}{\Omega_g^2},
\qquad
\mu_{\rm mix}:=\frac{M_0^2}{K_{\rm pole}G_0}.
\]

Then the exact induced-mixing coefficients collapse to

\[
a_2 = \mu_{\rm mix}\,r\,c_g,
\qquad
a_4 = \mu_{\rm mix}\,r^2 c_g(1-c_g).
\]

The initial grouped-`P2` support-demand drift is therefore controlled by the single combination

\[
a_4-2a_2
=
\mu_{\rm mix}\,r\,c_g\,[r(1-c_g)-2].
\]

So the first actual `l=0 \leftrightarrow l=2` mechanism is already much more rigid than a generic perturbative correction: it is governed by one amplitude `mu_mix`, one pole-ratio `r`, and one scalar-lane pole fraction `c_g`.

---

## Stage 266 — exact corridor theorem for the one-pole scalar lane

Introduce the compact control variable

\[
u := r(1-c_g).
\]

Then the exact weak-anisotropy corridor tests become

\[
M_{\rm twin}(y)=1+\alpha(4-u)y+2\alpha^2 y^2,
\qquad
\alpha:=\mu_{\rm mix}rc_g,
\qquad y:=\chi^2,
\]

and, for a general Family-1 ceiling `c_*`,

\[
M_{F1}(y)=(4c_*-1)+\alpha(8c_*-u)y+4c_*\alpha^2y^2.
\]

This gives a complete exact threshold ladder.

### Initial-drift thresholds

\[
u > 2
\quad\Longleftrightarrow\quad
\text{the actual grouped-`P2` support demand grows initially},
\]

\[
u > 4
\quad\Longleftrightarrow\quad
\text{the universal twin-safe margin shrinks initially},
\]

\[
u > 8c_*
\quad\Longleftrightarrow\quad
\text{the exact Family-1 margin shrinks initially}.
\]

### Actual failure thresholds

The exact twin discriminant is

\[
\Delta_{\rm twin}=\alpha^2[(u-4)^2-8],
\]

so a positive twin-boundary root exists only if

\[
u \ge 4+2\sqrt2.
\]

Likewise the exact Family-1 discriminant is

\[
\Delta_{F1}=\alpha^2[(u-8c_*)^2-16c_*(4c_*-1)],
\]

so an actual Family-1 failure requires

\[
u \ge 8c_* + 4\sqrt{c_*(4c_*-1)}.
\]

That is a strong theorem. The first induced geometry-mixing mechanism threatens the isotropic branch only when the scalar lane is both sufficiently pole-light and sufficiently faster than the grouped quadrupole pole.

### Pure-pole scalar lane

If `c_g=1`, then

\[
a_4=0,
\qquad
u=0,
\]

so the branch sits far below every drift/failure threshold. A pure scalar pole can soften the demand, but by itself it can never drive the actual grouped-`P2` branch out of either the universal twin-safe strip or the exact Family-1 corridor.

---

## Stage 267 — exact numerical thresholds on the Lambda_EM-refreshed Family-1 branch

Using the exact refreshed unblocked Family-1 ceiling,

\[
c_* = c_{\rm pole,max}^{(F1)}(0) \approx 0.7116102605,
\]

the hard threshold values are

\[
u_{\rm twin}^{\rm fail} = 4+2\sqrt2 \approx 6.8284271247,
\]

\[
u_{F1}^{\rm fail} = 8c_* + 4\sqrt{c_*(4c_*-1)} \approx 10.2779821110.
\]

Since `u = r(1-c_g)`, the required pole-ratio thresholds scale as

\[
r_{\rm crit} = \frac{u_{\rm crit}}{1-c_g}.
\]

Representative values on the exact Lambda_EM-refreshed hard ceiling are:

- `c_g = 0.25`:
  \[
  r_{\rm twin}^{\rm fail} \approx 9.10457,
  \qquad
  r_{F1}^{\rm fail} \approx 13.70398;
  \]
- `c_g = 0.50`:
  \[
  r_{\rm twin}^{\rm fail} \approx 13.65685,
  \qquad
  r_{F1}^{\rm fail} \approx 20.55596;
  \]
- `c_g = 0.75`:
  \[
  r_{\rm twin}^{\rm fail} \approx 27.31371,
  \qquad
  r_{F1}^{\rm fail} \approx 41.11193.
  \]

So unless the scalar geometry pole is *much* faster than the grouped quadrupole pole, the first actual `l=0 \leftrightarrow l=2` mixing mechanism does not by itself drive the real isotropic branch out of the universal twin-safe or exact Family-1 strips.

---

## Best current reading after Stages 265–267

The next theorem gate is now much narrower.

The first concrete induced geometry-mixing mechanism is no longer an unspecified `O(chi^2)` nuisance. It has one exact coefficient pair `(a_2,a_4)` and one exact danger variable

\[
u = r(1-c_g).
\]

The resulting message is surprisingly strong:

1. pure-pole scalar-lane mixing is automatically safe;
2. mixed contact+pole scalar lanes become dangerous only at rather large quadrupole/scalar pole ratios;
3. on the exact Lambda_EM-refreshed Family-1 branch, actual twin/F1 failure requires `u` above about `6.83` / `10.28`, respectively;
4. therefore this first induced `l=0 \leftrightarrow l=2` mechanism looks much more like a controlled correction than like the actual source of 5PN failure.

So the next honest continuation is no longer “some mixing exists.” It is:

> derive from the actual moving-throat reduced operator whether the first nontrivial scalar/geometry lane is closer to the safe pure-pole end, the safe contact-dominated end, or the genuinely dangerous fast mixed contact+pole regime.
