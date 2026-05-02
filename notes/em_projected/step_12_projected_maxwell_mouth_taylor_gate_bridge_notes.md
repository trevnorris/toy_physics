# Projected Maxwell — Mouth-Taylor Gate Bridge Notes

This continuation takes the abstract projected-Maxwell primitive slippages
\[
(q_1,\ s_1,\ h_1,\ d_1,\ p_1,\ g_1)
\]
and ties them to a **concrete near-throat ansatz**.

The goal is to see, in the cleanest possible way, how a projection-first
near-mouth Maxwell correction could feed the actual moving-throat PDE
bottlenecks rather than just the bundle moments in isolation.

---
## 1. Concrete near-throat ansatz

Let \(s\ge 0\) be a local throat coordinate measured inward from the mouth, and
let the projected one-port primitive packet be
\[
X(s)\in\{Q(s),\,S_2(s),\,H_{\rm port}(s),\,\Delta(s),\,P(s),\,G_W(s)\}.
\]

Use a normalized one-sided mouth kernel
\[
W_\ell(s)=\frac1\ell\,w(s/\ell),\qquad
\int_0^\infty w(u)\,du=1,
\]
with first moment
\[
\mu_1=\int_0^\infty u\,w(u)\,du.
\]

Then every primitive quantity has the near-mouth expansion
\[
X_{\rm proj}(\ell)=X(0)+\ell\,\mu_1\,X'(0)+O(\ell^2).
\]

So the earlier abstract slippages become
\[
q_1=\mu_1 Q'(0),\qquad
s_1=\mu_1 S_2'(0),\qquad
h_1=\mu_1 H'_{\rm port}(0),
\]
\[
d_1=\mu_1 \Delta'(0),\qquad
p_1=\mu_1 P'(0),\qquad
g_1=\mu_1 G_W'(0).
\]

This makes the earlier “mouth-local” language completely concrete:

- a **mouth-anchored** kernel has \(\mu_1\neq 0\), so the first correction is \(O(\ell)\),
- a **symmetric interior** kernel has \(\mu_1=0\), so this whole first layer vanishes.

---
## 2. Primitive-to-bundle derivative map

Dividing out the common factor \(\mu_1\), the induced projected-Maxwell bundle
shifts are

\[
\frac{z_0}{\mu_1}
=
\frac{\Delta Q' - Q\Delta'}{\Delta^2},
\]

\[
\frac{z_2}{\mu_1}
=
-\frac{\Delta^2 H'_{\rm port}-\Delta H_{\rm port}\Delta'-\Delta Q S_2'
-\Delta S_2 Q'+2Q S_2 \Delta'}{\Delta^3},
\]

\[
\frac{z_4}{\mu_1}
=
-\frac{
\Delta^2 H_{\rm port}S_2'
+\Delta^2 S_2 H'_{\rm port}
+\Delta^2 Q'
-2\Delta H_{\rm port}S_2\Delta'
-2\Delta Q\Delta'
-2\Delta Q S_2 S_2'
-\Delta S_2^2 Q'
+3Q S_2^2\Delta'
}{\Delta^4},
\]

\[
\frac{n_0}{\mu_1}
=
\frac{2P(\Delta P'-P\Delta')}{\Delta^3},
\]

\[
\frac{n_2}{\mu_1}
=
-\frac{2\left(
\Delta^2 G_W P'
+\Delta^2 P G_W'
-2\Delta G_W P \Delta'
-\Delta P^2 S_2'
-2\Delta P P' S_2
+3P^2 S_2 \Delta'
\right)}{\Delta^4},
\]

\[
\frac{n_4}{\mu_1}
=
\frac{2\left(
\Delta^3 G_W G_W'
-\Delta^2 G_W^2\Delta'
-2\Delta^2 G_W P S_2'
-2\Delta^2 G_W P' S_2
-2\Delta^2 P G_W' S_2
-2\Delta^2 P P'
+6\Delta G_W P S_2 \Delta'
+3\Delta P^2 \Delta'
+3\Delta P^2 S_2 S_2'
+3\Delta P P' S_2^2
-6P^2 S_2^2 \Delta'
\right)}{\Delta^5}.
\]

So the concrete near-mouth ansatz does exactly what we wanted: it turns the
projected-Maxwell one-port packet into explicit local derivative data.

---
## 3. Direct bridge to the actual 5PN / moving-throat bottlenecks

Use the existing weak-axisymmetric bottlenecks
\[
\Xi_{\rm load}=\frac{P_1}{P_0},\qquad
K_1=D_{21}+\frac{D_{01}}{9},\qquad
H_{\rm even}=D_{41}-\frac23 D_{21}-\frac{D_{01}}{27}.
\]

The projected-Maxwell mouth layer contributes

\[
\frac{\Xi_{\rm load}^{(\rm proj)}}{\mu_1}
=
2\frac{P'}{P}
-2\frac{\Delta'}{\Delta}
+\frac{Q'}{D_0\Delta}
-\frac{Q\Delta'}{D_0\Delta^2},
\]

\[
\frac{K_1^{(\rm proj)}}{\mu_1}
=
\frac{
9\Delta^2 H'_{\rm port}
-\Delta^2 Q'
-9\Delta H_{\rm port}\Delta'
+\Delta Q\Delta'
-9\Delta Q S_2'
-9\Delta S_2 Q'
+18Q S_2 \Delta'
}{9\Delta^3},
\]

\[
\frac{H_{\rm even}^{(\rm proj)}}{\mu_1}
=
-\frac{
18\Delta^3 H'_{\rm port}
+\Delta^3 Q'
-18\Delta^2 H_{\rm port}\Delta'
-\Delta^2 Q\Delta'
-27\Delta^2 H_{\rm port} S_2'
-27\Delta^2 S_2 H'_{\rm port}
-18\Delta^2 Q S_2'
-18\Delta^2 S_2 Q'
-27\Delta^2 Q'
+54\Delta H_{\rm port}S_2\Delta'
+36\Delta Q S_2 \Delta'
+54\Delta Q\Delta'
+54\Delta Q S_2 S_2'
+27\Delta S_2^2 Q'
-81Q S_2^2\Delta'
}{27\Delta^4}.
\]

The immediate factorization is the most useful result:

- \(\Xi_{\rm load}^{(\rm proj)}\) depends only on \((Q',\Delta',P')\),
- \(K_1^{(\rm proj)}\) and \(H_{\rm even}^{(\rm proj)}\) depend only on \((Q',\Delta',S_2',H'_{\rm port})\),
- \(G_W'\) does **not** enter \(\Xi_{\rm load},K_1,H_{\rm even}\),
- \(G_W'\) first enters only through the constant-prefactor transport slots.

So the projected-Maxwell throat packet naturally splits into separate jobs.

---
## 4. Exact mechanism sieve

Two “simple” possibilities fail immediately.

### 4.1 Source/denominator only

If
\[
S_2'=H'_{\rm port}=0,
\]
then solving
\[
K_1^{(\rm proj)}=H_{\rm even}^{(\rm proj)}=0
\]
gives only
\[
Q'=\Delta'=0.
\]

So a pure \((Q',\Delta')\) projected-Maxwell mouth correction cannot close the
even gates nontrivially.

### 4.2 Spectral only

If
\[
Q'=\Delta'=0,
\]
then solving the same equations gives only
\[
S_2'=H'_{\rm port}=0.
\]

So a pure \((S_2',H'_{\rm port})\) correction also cannot close the even gates
nontrivially.

---
## 5. Exact mixed compensation surface

A nontrivial projected-Maxwell repair of the conservative even gates therefore
requires a **mixed** correction involving both source/denominator and spectral
slots.

Solving
\[
K_1^{(\rm proj)}=H_{\rm even}^{(\rm proj)}=0
\]
for \((H'_{\rm port},S_2')\) gives an exact compensation surface:

\[
H'_{\rm port}
=
\frac{
\Delta^3 H_{\rm port}Q'
+\Delta^3 Q Q'
+9\Delta^2 H_{\rm port}^2\Delta'
-\Delta^2 H_{\rm port}Q\Delta'
-\Delta^2 Q^2\Delta'
+9\Delta^2 H_{\rm port}S_2 Q'
-2\Delta^2 Q S_2 Q'
-9\Delta^2 Q Q'
-18\Delta H_{\rm port}Q S_2\Delta'
+2\Delta Q^2 S_2\Delta'
+18\Delta Q^2\Delta'
-9\Delta Q S_2^2 Q'
+9Q^2 S_2^2\Delta'
}{9\Delta^2(\Delta H_{\rm port}-Q S_2)},
\]

\[
S_2'
=
\frac{
\Delta^3 Q'
-\Delta^2 Q\Delta'
-\Delta^2 S_2 Q'
-9\Delta^2 Q'
+9\Delta H_{\rm port}S_2\Delta'
+\Delta Q S_2\Delta'
+18\Delta Q\Delta'
-9Q S_2^2\Delta'
}{9\Delta(\Delta H_{\rm port}-Q S_2)}.
\]

The denominator
\[
\Delta H_{\rm port}-Q S_2
\]
is exactly the primitive combination behind the conservative
\[
Z_2=\frac{Q S_2-H_{\rm port}\Delta}{\Delta^2}
\]
slot. So the compensation surface becomes singular on
\[
\Delta H_{\rm port}=Q S_2
\qquad\Leftrightarrow\qquad
Z_2=0.
\]

That is an unexpectedly clean structural result.

---
## 6. Constant-prefactor transport

On the constant-prefactor branch, the projected-Maxwell mouth layer gives
derivative-level shifts \(\delta P_2,\delta P_4\) that are linear in \(G_W'\).

The simplest coefficients are

\[
\frac{\partial \Xi_{\rm load}^{(\rm proj)}}{\partial P'}=\frac{2}{P},
\]

\[
\frac{\partial (\delta P_2^{(\rm proj)})}{\partial G_W'}
=
-\frac{2P}{D_0\Delta^2},
\]

\[
\frac{\partial (\delta P_4^{(\rm proj)})}{\partial G_W'}
=
\frac{2\left(D_0\Delta G_W - 2D_0 P S_2 + 2D_2\Delta P\right)}{D_0^2\Delta^3}.
\]

So at derivative level:

- \(P'\) is the clean direct prefactor mover,
- \(G_W'\) is the clean outgoing-leg mover for the constant-prefactor branch.

---
## 7. Practical interpretation

The concrete mouth-Taylor ansatz says that a real near-throat projected-Maxwell
missing piece would naturally have **three** separate jobs:

1. \((Q',\Delta',S_2',H'_{\rm port})\) repair the conservative even gates,
2. \(P'\) tunes the weak-axisymmetric prefactor slope \(\Xi_1=P_1/P_0\),
3. \(G_W'\) tunes the constant-prefactor transport.

That is much more structured than “one extra EM correction near the mouth.”
It looks exactly like the kind of multi-channel missing ingredient a real
moving-throat PDE would generate.

---
## 8. Why this is useful for the PDE program

This continuation does **not** solve the full moving-throat branch. But it does
sharpen the next theorem target:

- a single-channel projected-Maxwell correction is not enough,
- the even-gate repair must live on a mixed compensation surface,
- the static prefactor bottleneck can be tuned without disturbing that surface,
- and the outgoing-leg slope only appears at the next transport layer.

That is a strong sign that projection-first Maxwell near the throat is not just
another side correction. It has the right internal structure to be a genuine
missing PDE-side ingredient.
