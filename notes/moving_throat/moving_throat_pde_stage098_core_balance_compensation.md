
# Moving-Throat PDE — Stage 98: Exact Core-Balance Compensation Theorem

## Goal

Determine when the concrete two-channel core model of Stage 97 lands exactly on the nontrivial compensated canonical branch found algebraically in Stage 95.

## Canonical branch conditions

The Stage-95 nontrivial compensated outlet is
\[
\delta\Lambda_{\rm can}(z)
=
4\sigma_*-\frac{\sigma_*}{1-z^2/3-i z^5/9},
\]
which is exactly the branch that preserves the canonical outgoing quadrupole fingerprint.

For the concrete core model, this requires
\[
\boxed{
\rho_c=4\sigma_c,
\qquad
\kappa_c=\frac13,
\qquad
\gamma_c=\frac19.
}
\]

## Exact coupling-balance law

Using the Stage-97 identifications,
\[
\rho_c=\frac{g_s^2}{K_s},
\qquad
\sigma_c=
\frac{(K_s g_q-\lambda g_s)^2}
{K_s^2K_q(1+r_c)},
\qquad
r_c=\frac{\lambda^2}{K_sK_q},
\]
the nontrivial compensation condition becomes
\[
\boxed{
g_s^2\bigl(K_sK_q+\lambda^2\bigr)
=
4\bigl(K_s g_q-\lambda g_s\bigr)^2.
}
\]
Solving for the mixed coupling gives the exact two-branch law
\[
\boxed{
g_q=
\frac{g_s}{2K_s}
\left(
2\lambda \pm \sqrt{K_sK_q+\lambda^2}
\right).
}
\]

So the required Robin–mixed balance is not a mystery coefficient fit. It is an explicit codimension-one surface in the concrete core-coupling space.

## Bare mixed-channel normalization conditions

The even/odd preservation conditions become
\[
\boxed{
\kappa_0=\frac{1+r_c}{3},
\qquad
\gamma_0=\frac{1+r_c}{9}.
}
\]
So the bare mixed side-channel must itself be a scale-deformed copy of the canonical compact outgoing branch.

## Exact collapse to the Stage-95 branch

On the coupling-balance surface together with the bare-channel conditions above, the concrete core outlet reduces identically to
\[
\boxed{
\delta\Lambda_{\rm core}(z)
=
4\sigma_*-\frac{\sigma_*}{1-z^2/3-i z^5/9},
\qquad
\sigma_*=\frac{g_s^2}{4K_s}.
}
\]
Adding this to the canonical exterior branch
\[
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}+O(z^6)
\]
gives a normalized response with exactly the same outgoing fingerprint:
\[
\widehat Y_2(z)=1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}+O(z^6).
\]

So the concrete core model does not merely mimic the reduced algebra. It reproduces the exact nontrivial compensated branch.

## Interpretation

This is the first explicit throat-core theorem in the outlet program:

- the Stage-95 compensation law is not just an algebraic accident,
- it is realized by a concrete two-channel core model,
- and the surviving free data are sharply reduced to one coupling-balance surface plus one scale-deformed bare mixed outlet.
