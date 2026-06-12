# Moving-Throat PDE — Stage 137: Explicit Core-to-Mouth Gain Map

## Goal

Derive the actual coupled mouth-layer gains `(M_s,M_q)` from an explicit
GNLS + localized-Maxwell throat-core ansatz instead of leaving them as abstract
fixed-point coefficients.

The key idea is to use the already derived concrete two-channel core response from
Stages 114–116 and embed it directly into the mouth electrochemical free energy.

---

## 1. Explicit electrochemical mouth free energy

Take the positive normalized mouth-source density `\sigma(z)` and let the brane-facing
mouth electrochemical potential be sourced by two scalar channels:

- a shell/compliance channel `U_s(z)`,
- a mixed localized-Maxwell channel `U_q(z)`.

Use the explicit free-energy density
\[
F_{\rm mouth}[\sigma,U_s,U_q]
=
\int_0^L dz\,
\Big[
\Theta_\sigma\,\sigma\big(\ln(\sigma/\sigma_*)-1\big)
+\sigma\big(\rho_c U_s(z)-\sigma_c U_q(z)\big)
\Big].
\]

Here the shell contribution enters with positive sign, while the mixed
localized-Maxwell contribution enters with the opposite sign, exactly as expected
for the parent electrochemical combination
\[
\partial_z\delta V_{\rm conf}-q_*\partial_z A_0.
\]

Varying with respect to `\sigma` gives the same exponential source law as before,
but with a self-consistent slope parameter
\[
\Pi
=
\frac{L}{\Theta_\sigma}
\Big[
\rho_c U_s'(0)-\sigma_c U_q'(0)
\Big].
\]

So the gains are now set by the *same* static and mixed coefficients that appear in
the explicit throat-core outlet.

---

## 2. Core coefficients from the exact Schur complement

Stage 114 already gave the concrete two-channel core outlet in the form
\[
\delta\Lambda_{\rm core}(z)
=
\rho_c-
\frac{\sigma_c}{1-\kappa_c z^2-i\gamma_c z^5},
\]
with exact coefficients
\[
\rho_c=\frac{g_s^2}{K_s},
\qquad
\sigma_c=
\frac{(K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)}.
\]

Provenance note. The SymPy audit names the unnormalized denominator
\(D_W^{\rm bare}(z)=1-\kappa_c z^2-i\gamma_c z^5\) before the static normalization \(D_W^{\rm bare}(0)=1\). This is a
script-side name for the same Schur-complement denominator, not an extra parameter.

So the actual mouth-layer gains are
\[
\boxed{
M_s=
\frac{L}{\Theta_\sigma}\,\rho_c
=
\frac{L g_s^2}{K_s\Theta_\sigma},
}
\]
\[
\boxed{
M_q=
-\frac{L}{\Theta_\sigma}\,\sigma_c
=
-\frac{L (K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)\Theta_\sigma}.
}
\]

This is the first direct derivation of `(M_s,M_q)` from the explicit parent throat-core
ansatz.

The Stage 134 Family-1 fixed-point law therefore becomes
\[
\boxed{
\Pi=M_s+M_q\,\mathcal S_q(\Pi),
\qquad
\mathcal S_q(\Pi)=\mathcal S\!\left(\Pi,\frac{\pi}{2}\right).
}
\]

---

## 3. Immediate consequence

The abstract outlet-consistent closure of Stage 135 is no longer an extra guess.
It is the direct consequence of using the same Schur-complement core weights in the
mouth electrochemical source law.

The only remaining mouth-core data are now:

1. the shell stiffness `K_s`,
2. the mixed stiffness `K_q`,
3. the hybridization `\lambda`,
4. the mouth couplings `(g_s,g_q)`,
5. the source susceptibility `\Theta_\sigma`.

Once these are known on a branch, `(M_s,M_q)` are fixed.

---

## Result

The actual coupled mouth-layer gains are
\[
\boxed{
M_s=
\frac{L g_s^2}{K_s\Theta_\sigma},
\qquad
M_q=
-\frac{L (K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)\Theta_\sigma}.
}
\]

So the mouth fixed-point ambiguity has now collapsed from an abstract gain pair to
one explicit set of parent core quantities.
