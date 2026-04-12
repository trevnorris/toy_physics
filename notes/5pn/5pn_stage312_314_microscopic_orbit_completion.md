# 5PN continuation notes — Stages 312–314

These stages continue directly from the Stage-309–311 coherent branch defect normal form.

At Stage 311 the first coherent weak-axisymmetric theorem had already collapsed to
three exact branch-adapted scalars,
\[
(\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta),
\]
and equivalently to invariance of three direct branch composites.

The next honest step was to remove the remaining dependence on those intermediate
composites and rewrite the whole zero-defect test directly in the microscopic
coherent-kernel variables. That is what Stages 312–314 do.

## Stage 312 — direct microscopic monomial coordinates

Using the coherent-kernel definitions
\[
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\delta_U=\frac{\pi^2 T_U}{L^2K_U},
\qquad
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}},
\]
\[
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}},
\qquad
\frac{Z_W}{\Omega_W^2}=\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2},
\]
with frozen reference-branch exponents
\[
E_*=
\frac{2\epsilon_{W,*}(11+9\delta_{U,*})}{11(1-\epsilon_*)(1+\delta_{U,*})},
\qquad
F_*=
\frac{2\chi_{0,*}}{1+\delta_{U,*}}
+
\frac{4\epsilon_{W,*}\delta_{U,*}}{11(1-\epsilon_*)(1+\delta_{U,*})^2},
\]
the three direct microscopic monomials are
\[
\boxed{
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}}
}
\]
\[
\boxed{
\mathfrak C_{{\rm nt},*}
=
\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}
\left(
\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}}
\right)^{E_*}
\left(
\frac{\pi^2T_U}{L^2K_U}
\right)^{-F_*}
}
\]
\[
\boxed{
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}}.
}
\]

With microscopic grouped weak-axisymmetric log-drifts
\[
(\lambda_1,c_1,\gamma_1,\kappa_U,\kappa_\eta,\kappa_W,\mu_1,\tau_1),
\]
the exact logarithmic derivatives are
\[
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt},
\qquad
\delta\ln\epsilon_\eta=\Sigma_\eta.
\]
So the coherent first-order zero-defect theorem becomes
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln \mathfrak C_{{\rm tr},*}=0,
\quad
\delta\ln \mathfrak C_{{\rm nt},*}=0,
\quad
\delta\ln\epsilon_\eta=0.
}
\]

That is the smallest exact microscopic formulation reached so far.

## Stage 313 — exact microscopic compatibility ledger

Once the three direct monomials are used, the zero-defect theorem is equivalent
to the explicit microscopic compatibility system
\[
\boxed{
(1+\chi_{0,*})(\tau_1-\kappa_U)
+
(1+\delta_{U,*})(\gamma_1+c_1-\kappa_U)=0,
}
\]
\[
\boxed{
2c_1-\kappa_U-\kappa_\eta=0,
}
\]
\[
\boxed{
2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W
+
E_*(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W)
-
F_*(\tau_1-\kappa_U)=0.
}
\]

These solve exactly for the three dependent drifts:
\[
\boxed{
\tau_1
=
\kappa_U
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}
(\gamma_1+c_1-\kappa_U),
}
\]
\[
\boxed{
\kappa_\eta = 2c_1-\kappa_U,
}
\]
\[
\boxed{
\mu_1
=
2c_1-\kappa_U+2\kappa_W-2\lambda_1
-
E_*(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W)
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\gamma_1+c_1-\kappa_U).
}
\]

So the zero-defect branch is no longer an abstract slippage statement. It is a
codimension-3 microscopic rigidity ledger inside the eight-drift kernel space.

## Stage 314 — exact microscopic similarity orbit

The compatibility ledger is not an isolated fine-tuning condition. It is the
logarithmic tangent form of a finite five-parameter multiplicative similarity orbit.

Take free logarithmic orbit parameters
\[
(u_\lambda,u_c,u_\gamma,u_{K_U},u_{K_W}).
\]
The dependent exponents are fixed by exact monomial preservation:
\[
\boxed{
u_{K_\eta}=2u_c-u_{K_U},}
\]
\[
\boxed{
u_{T_U}=u_{K_U}-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(u_\gamma+u_c-u_{K_U}),}
\]
\[
\boxed{
u_{\mu_W}=2u_c-u_{K_U}+2u_{K_W}-2u_\lambda
-
E_*(2u_\gamma+2u_\lambda-u_{K_U}-u_{K_W})
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(u_\gamma+u_c-u_{K_U}).}
\]

The finite orbit acts multiplicatively by
\[
\lambda_W\mapsto e^{u_\lambda}\lambda_W,
\qquad
c_{\eta U}\mapsto e^{u_c}c_{\eta U},
\qquad
\gamma\mapsto e^{u_\gamma}\gamma,
\]
\[
K_U\mapsto e^{u_{K_U}}K_U,
\qquad
K_\eta^{(\mathrm{eff})}\mapsto e^{u_{K_\eta}}K_\eta^{(\mathrm{eff})},
\qquad
K_W^{(\mathrm{eff})}\mapsto e^{u_{K_W}}K_W^{(\mathrm{eff})},
\]
\[
\mu_W\mapsto e^{u_{\mu_W}}\mu_W,
\qquad
T_U\mapsto e^{u_{T_U}}T_U.
\]

On this finite family,
\[
\mathfrak C_{{\rm tr},*},
\qquad
\mathfrak C_{{\rm nt},*},
\qquad
\epsilon_\eta
\]
are preserved **exactly**.

The tangent matrix of this orbit has rank `5`, while the monomial-drift matrix has
rank `3`, so
\[
\dim\ker M_*=5,
\]
and the orbit tangent space is precisely the zero-defect tangent space.

Therefore the coherent weak-axisymmetric theorem can now be read geometrically:
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\mathbf x\in T_{\mathrm{id}}\mathcal G_*,
}
\]
where `\mathcal G_*` is the exact five-parameter monomial-preserving similarity orbit.

## Net result after Stage 314

The algebraic side of the coherent weak-axisymmetric branch is now effectively finished.

1. The remaining defect is measured by the drift of three exact microscopic monomials.
2. The zero-defect conditions are an explicit three-equation compatibility ledger.
3. That ledger is exactly the tangent-space condition for a finite five-parameter similarity orbit.

So the next honest theorem gate is no longer algebraic compression.
It is the actual branch-realization question:

> does the true moving-throat grouped weak-axisymmetric branch move tangent to the exact monomial-preserving similarity orbit `\mathcal G_*`?

That is the cleanest continuation point we have right now.
