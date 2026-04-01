# Moving-Throat PDE — Stage 88: Exact Fixing of `chi_Q`

## Goal

Use the explicit DtN fingerprint to determine the last reduced 2.5PN scalar
\[
\chi_Q
\]
in the canonical retarded grouped-`P2` one-pole-plus-contact module.

## Retarded minimal isotropic grouped-`P2` module

Write the retarded normalized quadrupole module as
\[
\widehat Y_Q^{\rm ret}(\omega)
=
\frac34
+
\frac14\,
\frac{1}{
1-\omega^2/\Omega_Q^2
-i\,\chi_Q\,\sigma_Q^{\rm can}\,\omega^5
}
+O(\omega^6).
\]

The conservative moment match already fixed
\[
\Omega_Q=\frac{3c_s}{2a}.
\]

Define
\[
\sigma_Q^{\rm can}:=\frac{9}{8\Omega_Q^5}.
\]
Using \(\Omega_Q=3c_s/(2a)\),
\[
\boxed{
\sigma_Q^{\rm can}
=
\frac{4a^5}{27c_s^5}.
}
\]

## Low-frequency expansion

Expanding through \(O(\omega^5)\) gives
\[
\widehat Y_Q^{\rm ret}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\,\chi_Q\,\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
\]

From Stage 87, the explicit outgoing DtN branch gives
\[
\widehat Y_2^{\rm out}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\,\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
\]

Matching the \(O(\omega^5)\) coefficient yields
\[
\boxed{\chi_Q=1.}
\]

So on the canonical compact passive/outgoing grouped-`P2` DtN branch, the last reduced normalization scalar is fixed exactly.

## General deformed DtN obstruction

A useful way to parametrize the only remaining PDE-facing freedom is to deform the outgoing DtN operator by
\[
\Lambda_2^{\rm def}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}
+i\,\xi_Q\,\frac{z^5}{9}+O(z^6).
\]

Then
\[
\widehat Y_2^{\rm def}(z)
=
1+\frac{z^2}{9}
+\frac{4z^4}{81}
+i\,\xi_Q\,\frac{z^5}{27}
+O(z^6),
\]
so
\[
\boxed{\chi_Q=\xi_Q.}
\]

That means the only reduced 2.5PN obstruction left after the present calculation is a deviation of the actual moving-throat DtN branch from the canonical outgoing `l=2` coefficient \(\xi_Q=1\).
