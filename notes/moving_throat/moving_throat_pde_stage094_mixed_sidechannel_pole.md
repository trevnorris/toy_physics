# Moving-Throat PDE — Stage 94: Explicit Mixed `A_w/F_{\mu w}`-Type Side-Channel Pole

## Goal

Build the first explicit isotropic hidden mixed-sector side-channel model and check whether it can preserve the already-fixed conservative even `l=2` fingerprint.

## Hidden pole model

Take the canonical outgoing branch and add the first isotropic Schur-complement-style mixed side-channel
\[
\boxed{
\Lambda_2^{\rm mix}(z)
=
\Lambda_2^{\rm out}(z)
-
\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}
+O(z^6),
}
\]
with
- `sigma_W > 0` the static mixed loading,
- `kappa_W > 0` the hidden even dispersion scale,
- `gamma_W > 0` the hidden odd outgoing coefficient.

Expanding,
\[
\Lambda_2^{\rm mix}(z)
=
-(3+\sigma_W)
+
\left(\frac13-\sigma_W\kappa_W\right)z^2
+
\left(\frac19-\sigma_W\kappa_W^2\right)z^4
+
 i\left(\frac19-\sigma_W\gamma_W\right)z^5
+O(z^6).
\]

## Exact even-branch no-go

Demand that this branch preserve the canonical conservative even fingerprint. The `z^2` condition gives
\[
-\frac{L_2}{L_0}=\frac19
\quad\Longrightarrow\quad
\boxed{\kappa_W=-\frac19.}
\]
This is already incompatible with a standard positive passive pole parameter `\kappa_W>0`.

Even if one formally inserts that value, the `z^4` condition becomes
\[
\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}=\frac{4}{81}
\quad\Longrightarrow\quad
\boxed{\sigma_W=0.}
\]
So a standalone isotropic hidden pole of this type cannot sit on the already-fixed canonical even branch unless it is absent.

## Outgoing-normalization shift

The raw branch normalization factor is
\[
\boxed{
\chi_Q^{\rm mix}
=
\frac{3(1-9\sigma_W\gamma_W)}{3+\sigma_W}.
}
\]
For small loading,
\[
\chi_Q^{\rm mix}
=
1-
\sigma_W\left(\frac13+9\gamma_W\right)
+O(\sigma_W^2).
\]
So the linearized branch-selection triple is
\[
\boxed{(b,a_0,a_5)=(0,-\sigma_W,-\sigma_W\gamma_W).}
\]

## Consequence

A naive passive mixed `A_w/F_{\mu w}` side-channel pole is **too rigid**. It generically shifts the outgoing normalization and, more importantly, it cannot preserve the already-fixed canonical even `l=2` branch. If a mixed sector survives on the actual branch, it must appear in a more structured, compensated outlet law.
