# 2PN mouth-operator / DtN reduction: current preliminary result

## Core structural result

Take the low-frequency monopole mouth operator in the form
\[
Z_{00}(\omega;u)=Z_2(u)\,\omega^2+Z_4(u)\,\omega^4+O(\omega^6),
\qquad u\equiv U/c^2.
\]

Define the normalized DtN invariants
\[
C_s(u)\equiv \frac{Z_4/Z_2^3}{(Z_4/Z_2^3)_0},
\qquad
G(u)\equiv \frac{Z_4/Z_2^2}{(Z_4/Z_2^2)_0}.
\]

For the cylinder / Neumann-bottom unit-test branch,
\[
Z_2=-\frac{L}{c_s^2},
\qquad
Z_4=-\frac{L^3}{3c_s^4},
\]
so these become exactly
\[
C_s(u)=\frac{c_s^2(u)}{c_{s0}^2},
\qquad
G(u)=\frac{L(u)}{L_0}.
\]

Under the current 2PN freeze:
- exact Bernoulli gives \(c_s^2/c^2 = 1-4u\),
- the reduced throat closure gives \(L/L_0 = a(u)\),
- and the known 1PN breathing slope implies
  \[
  a(u)=1+\frac{57}{64}u+O(u^2).
  \]

## Unique minimal one-body closure that preserves the frozen 1PN slot

If the corrected one-body denominator is built from the DtN geometry invariant as
\[
D_{\rm eff}(u)=C_s(u)\left[1+\alpha\,(G-1)+\beta\,(G-1)^2\right],
\]
then preserving the already-fixed 1PN coefficient forces
\[
\alpha=0.
\]
So the first nontrivial allowed correction is quadratic:
\[
D_{\rm eff}(u)=C_s(u)\left[1+\mu\,(G-1)^2\right].
\]

At 2PN only the linear geometry slope matters. Writing
\[
G(u)=1+g_1 u+g_2 u^2+\cdots,
\]
one gets
\[
D_{\rm eff}(u)=1-4u+\mu g_1^2 u^2+O(u^3).
\]
Matching the exact isotropic Schwarzschild test-mass target through 2PN therefore fixes
\[
\mu g_1^2 = 8.
\]
Using the already-fixed throat slope
\[
g_1=\frac{57}{64},
\]
this gives
\[
\mu = \frac{8}{g_1^2} = \frac{32768}{3249} \approx 10.0855647892.
\]

So the final one-body denominator is
\[
D_{\rm eff}(u)=C_s(u)\left[1+\frac{32768}{3249}(G(u)-1)^2\right].
\]

With the current closure this expands to
\[
D_{\rm eff}(u)=1-4u+8u^2+O(u^3),
\]
which exactly closes the isotropic one-body 2PN target.

## Explicit current-series values

From the current throat closure:
\[
a(u)=1+\frac{57}{64}u+\frac{298821}{131072}u^2+O(u^3).
\]

Therefore the cylinder DtN coefficients become
\[
\frac{Z_2(u)}{Z_2(0)} = 1+\frac{313}{64}u+\frac{2862917}{131072}u^2+O(u^3),
\]
\[
\frac{Z_4(u)}{Z_4(0)} = 1+\frac{683}{64}u+\frac{10301487}{131072}u^2+O(u^3).
\]

The extracted invariants are then
\[
C_s(u)=1-4u,
\qquad
G(u)=1+\frac{57}{64}u+\frac{298821}{131072}u^2+O(u^3).
\]

## Relation to the earlier raw resonance-proxy fit

The raw resonance proxy derived from the same DtN data is
\[
D_{\rm raw}(u)=\frac{[Z_2/Z_4](u)}{[Z_2/Z_4](0)}
=\frac{1-4u}{G(u)^2}
=1-\frac{185}{32}u+\frac{324075}{65536}u^2+O(u^3).
\]

The multiplicative factor needed to convert it into the exact target denominator is
\[
P_{\rm port}(u)=G(u)^2\left[1+\mu (G(u)-1)^2\right]
=1+\frac{57}{32}u+\frac{875093}{65536}u^2+O(u^3).
\]

So the earlier port fit was not arbitrary; it factorizes cleanly as
\[
D_{\rm eff}(u)=D_{\rm raw}(u)\,P_{\rm port}(u).
\]

## Current status

This closes the **one-body missing 2PN response slot** under a minimal DtN-invariant conservative ansatz.

It does **not** yet derive that ansatz from the full inner PDE, and it does **not** yet close the comparable-mass 2PN wake/cross sector.
