# Moving-Throat PDE — Stage 104: Exact Outgoing `l=2` DtN Fingerprint

## Goal

Compute the canonical compact passive/outgoing quadrupole normalization directly from an explicit outgoing `l=2` Dirichlet-to-Neumann model instead of carrying the coefficient `chi_Q` symbolically.

## Exact outgoing DtN model

Let
\[
z := \frac{a\omega}{c_s},
\]
and let the outgoing partial wave be the spherical Hankel mode
\[
h_2^{(1)}(z)=j_2(z)+i\,y_2(z).
\]

The exact `l=2` outgoing DtN operator is
\[
\Lambda_2^{\rm out}(z)
=
z\,\frac{d}{dz}\ln h_2^{(1)}(z)
=
z\,\frac{h_2^{(1)\prime}(z)}{h_2^{(1)}(z)}.
\]

Its small-\(z\) expansion is
\[
\boxed{
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}
+i\,\frac{z^5}{9}
-\frac{2z^6}{27}
-i\,\frac{z^7}{27}
+O(z^8).
}
\]

So the exact static slot is \(\Lambda_2^{\rm out}(0)=-3\).

## Normalized outgoing quadrupole admittance

Define the normalized outgoing branch by
\[
\widehat Y_2^{\rm out}(z)
:=
\frac{\Lambda_2^{\rm out}(0)}{\Lambda_2^{\rm out}(z)}
=
-\frac{3}{\Lambda_2^{\rm out}(z)}.
\]

Then
\[
\boxed{
\widehat Y_2^{\rm out}(z)
=
1+\frac{z^2}{9}
+\frac{4z^4}{81}
+i\,\frac{z^5}{27}
-\frac{11z^6}{729}
-i\,\frac{z^7}{243}
+O(z^8).
}
\]

Restoring \(\omega\),
\[
z=\frac{a\omega}{c_s},
\]
this becomes
\[
\boxed{
\widehat Y_2^{\rm out}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\,\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
}
\]

This is exactly the outgoing compact `l=2` fingerprint used earlier in the reduced 2.5PN program, but now derived directly from the explicit DtN model.

## Consequence

The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch, the leading odd quadrupole coefficient is fixed to
\[
\Gamma_{5,\rm can}^{\rm DtN}=\frac{a^5}{27c_s^5}.
\]

So any later normalization mismatch must come from branch selection or source normalization, not from ambiguity in the canonical outgoing `l=2` DtN model itself.
