# Stage V2-14 — Compact outgoing `l=2` fingerprint derivation

## Purpose

This stage rederives the compact outgoing quadrupole fingerprint

\[
\widehat Y_2^{\rm out}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6)
\]

from the outgoing spherical Hankel solution.  The main convention fixed here is:

\[
\boxed{
\widehat Y_2^{\rm out}
\text{ is the normalized inverse DtN compliance, not the raw DtN admittance.}
}
\]

If one uses the raw normalized Dirichlet-to-Neumann map, the signs of the
low-frequency even terms and the leading odd term are reversed before inversion.

---

## 1. Setup

Use the time convention

\[
e^{-i\omega t},
\]

so the outgoing radial solution in three spatial dimensions is the spherical
Hankel function

\[
h_l^{(1)}(kr).
\]

Let

\[
z=ka=\frac{a\omega}{c_s}.
\]

For \(l=2\), the exact outgoing function can be written as

\[
h_2^{(1)}(z)
=
i\,e^{iz}\frac{z^2+3iz-3}{z^3}.
\]

Its static singular behavior is

\[
h_2^{(1)}(z)\sim -\frac{3i}{z^3},
\]

which is the outgoing continuation of the exterior static quadrupole field
\(\sim r^{-3}\).

---

## 2. Raw normalized DtN map

For a boundary datum at \(r=a\), the outgoing exterior field is proportional to

\[
f(r)=h_2^{(1)}(kr).
\]

The raw DtN ratio is

\[
a\frac{f'(a)}{f(a)}
=
z\frac{h_2'(z)}{h_2(z)}.
\]

Since the static \(l=2\) exterior field scales as \(r^{-3}\), the normalized raw
DtN is

\[
\widehat D_2^{\rm raw}(z)
=
-\frac{z}{3}\frac{h_2'(z)}{h_2(z)}.
\]

SymPy gives

\[
\widehat D_2^{\rm raw}(z)
=
1-\frac{z^2}{9}
-\frac{z^4}{27}
-i\frac{z^5}{27}
+O(z^6).
\]

So the quoted bridge fingerprint is **not** this raw DtN object.

---

## 3. Normalized inverse DtN compliance

The compact outgoing bridge uses the inverse normalized compliance

\[
\widehat Y_2^{\rm out}(z)
=
\left(\widehat D_2^{\rm raw}(z)\right)^{-1}
=
-\,\frac{3h_2^{(1)}(z)}{z\,h_2^{(1)\prime}(z)}.
\]

The exact rational form is

\[
\widehat Y_2^{\rm out}(z)
=
\frac{
3(3-3iz-z^2)
}{
9-4z^2+i(z^3-9z)
}.
\]

Expanding at low frequency,

\[
\widehat Y_2^{\rm out}(z)
=
1+\frac{z^2}{9}
+\frac{4z^4}{81}
+i\frac{z^5}{27}
-\frac{11z^6}{729}
-i\frac{z^7}{243}
+O(z^8).
\]

Therefore, through the order needed by the 2.5PN and 4PN bridge,

\[
\boxed{
\widehat Y_2^{\rm out}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
}
\]

The leading odd coefficient is

\[
\boxed{
\Gamma_5^{\rm port}
=
\frac{a^5}{27c_s^5}.
}
\]

---

## 4. Channel interpretation

This result is a theorem for a **three-spatial-dimensional exterior outgoing
worldtube/STF quadrupole port**.

It is not automatically the result for an unrestricted \(4\)-spatial-dimensional
bulk radiation channel.  In \(d\) spatial dimensions, the small-\(z\) absorptive
power of a clean outgoing \(l\)-partial wave scales as

\[
z^{2l+d-2}.
\]

Thus

\[
d=3,\ l=2
\quad\Rightarrow\quad
z^5,
\]

while

\[
d=4,\ l=2
\quad\Rightarrow\quad
z^6
\]

(up to the usual integer-Bessel logarithmic subtleties).  Therefore the
\(i\omega^5\) branch is specifically the brane/worldtube \(3\)-space STF
quadrupole port.  A future full-PDE realization stage must show that the actual
moving-throat branch projects onto this port rather than an unrestricted
\(4\)-space bulk outgoing law.

---

## 5. Audit result

The SymPy script verifies:

1. the exact \(h_2^{(1)}\) closed form,
2. the static normalization,
3. the raw normalized DtN expansion,
4. the inverse-compliance expansion,
5. the coefficients \(1/9\), \(4/81\), and \(i/27\),
6. the leading coefficient \(\Gamma_5=a^5/(27c_s^5)\),
7. and the channel discriminator \(2l+d-2\).

The stage passes with one wording patch:

\[
\boxed{
\text{call the object a normalized outgoing compliance / inverse DtN response.}
}
\]

Calling the same expansion the raw DtN admittance would be incorrect.
