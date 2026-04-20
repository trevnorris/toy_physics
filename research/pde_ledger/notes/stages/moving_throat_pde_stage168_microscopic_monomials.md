
# Moving-Throat PDE — Stage 168: Direct Microscopic Monomial Coordinates and the First-Order Zero-Defect Compatibility Ledger

## Purpose

Stage 167 reduced the coherent weak-axisymmetric continuation point to the first grouped drifts of three exact branch composites,
\[
R_{\rm tr},
\qquad
\mathcal T^2 R_{\rm tr}^{B_*},
\qquad
\epsilon_\eta.
\]
That was already a sharp compression, but the branch variables still mixed direct observables with derived ratios.

The next honest step is even smaller:

> push the Stage-167 branch coordinates all the way down to direct **microscopic monomials** of the coherent-kernel couplings, and rewrite the zero-defect condition as an explicit linear compatibility ledger for the grouped weak-axisymmetric microscopic drifts.

That is what this stage does.

The main result is that, at first grouped weak-axisymmetric/reference-branch order, the three branch-adapted defect coordinates are exactly the first logarithmic drifts of three direct microscopic monomials:
\[
\boxed{
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt},
\qquad
\delta\ln \epsilon_\eta=\Sigma_\eta.
}
\]
Here
\[
\boxed{
\mathfrak C_{{\rm tr},*}:=\chi_0^{\,1+\delta_{U,*}}\,
\delta_U^{\,1+\chi_{0,*}},
}
\]
and
\[
\boxed{
\mathfrak C_{{\rm nt},*}
:=
\frac{Z_W}{\Omega_W^2}\,
\epsilon_W^{E_*}\,
\delta_U^{-F_*},
}
\]
with
\[
\boxed{
E_*:=\frac{2\epsilon_{W,*}}{1-\epsilon_*}\,
\frac{11+9\delta_{U,*}}{11(1+\delta_{U,*})},
}
\]
\[
\boxed{
F_*:=\frac{2\chi_{0,*}}{1+\delta_{U,*}}
+\frac{4\epsilon_{W,*}\delta_{U,*}}
{11(1-\epsilon_*)(1+\delta_{U,*})^2}.
}
\]

So the coherent weak-axisymmetric zero-defect theorem at this same linearized reference-branch order is now equivalent to the first-order invariance of three direct microscopic objects:
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln \mathfrak C_{{\rm tr},*}=0,
\quad
\delta\ln \mathfrak C_{{\rm nt},*}=0,
\quad
\delta\ln \epsilon_\eta=0.
}
\]

This is a genuine sharpening of Stage 167:

- the **tracking** coordinate becomes a monomial in \((\chi_0,\delta_U)\),
- the **nontracking** coordinate becomes a monomial in \((Z_W/\Omega_W^2,\epsilon_W,\delta_U)\),
- and the **dressing** coordinate is already the microscopic ratio \(\epsilon_\eta\).

The continuation point is therefore no longer “compute drifts of three branch composites.” It is now the smaller theorem gate:

> determine whether the true grouped weak-axisymmetric branch preserves, to first order about the coherent reference branch, the three direct microscopic monomials
> \(\mathfrak C_{{\rm tr},*}\), \(\mathfrak C_{{\rm nt},*}\), and \(\epsilon_\eta\).

---

## 1. Reference-branch microscopic log coordinates

Stage 165 already gave the microscopic drift variables
\[
\Sigma_\chi=\delta\ln\chi_0,
\qquad
\Sigma_\delta=\delta\ln\delta_U,
\qquad
\Sigma_Z=\delta\ln\!\left(\frac{Z_W}{\Omega_W^2}\right),
\]
\[
\Sigma_\epsilon=\delta\ln\epsilon_W,
\qquad
\Sigma_\eta=\delta\ln\epsilon_\eta.
\]

Using the microscopic coherent-kernel definitions,
\[
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\delta_U=\frac{\pi^2T_U}{L^2K_U},
\qquad
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}},
\]
\[
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}},
\qquad
\frac{Z_W}{\Omega_W^2}
=
\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2},
\]
these are explicitly
\[
\Sigma_\chi=\gamma_1+c_1-\kappa_U,
\qquad
\Sigma_\delta=\tau_1-\kappa_U,
\]
\[
\Sigma_Z=2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W,
\qquad
\Sigma_\epsilon=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W,
\]
\[
\Sigma_\eta=2c_1-\kappa_U-\kappa_\eta.
\]

So the Stage-166 branch-adapted coordinates already sit inside a five-dimensional logarithmic microscopic ledger.

---

## 2. Direct microscopic tracking monomial

Stage 166 defined the tracking coordinate by
\[
\Sigma_{\rm tr}
=
(1+\chi_0)\Sigma_\delta
+
(1+\delta_U)\Sigma_\chi.
\]
At first grouped weak-axisymmetric order on the coherent reference branch, the coefficients are frozen at
\[
1+\chi_{0,*},
\qquad
1+\delta_{U,*}.
\]
Therefore define the direct microscopic tracking monomial
\[
\boxed{
\mathfrak C_{{\rm tr},*}
:=
\chi_0^{\,1+\delta_{U,*}}\,
\delta_U^{\,1+\chi_{0,*}}.
}
\]
Its first logarithmic drift is
\[
\delta\ln \mathfrak C_{{\rm tr},*}
=
(1+\delta_{U,*})\,\delta\ln\chi_0
+
(1+\chi_{0,*})\,\delta\ln\delta_U
=
\Sigma_{\rm tr}.
\]
So the tracking branch-adapted coordinate is no longer indirect:
\[
\boxed{
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr}.
}
\]

### Explicit microscopic form

Substituting the coherent-kernel definitions gives
\[
\boxed{
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}}.
}
\]

So the first zero-defect condition can be read directly as invariance of one monomial in the microscopic interference and split-\(U\) ratios.

---

## 3. Direct microscopic nontracking monomial

Stage 166 defined the genuine nontracking transfer-shape slippage by
\[
\Sigma_{\rm nt}
=
\Sigma_Z
+
\frac{2\epsilon_W}{1-\epsilon}
\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-
\left[
\frac{2\chi_0}{1+\delta_U}
+
\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
\]
Freeze the reference-branch coefficients:
\[
\boxed{
E_*:=
\frac{2\epsilon_{W,*}}{1-\epsilon_*}\,
\frac{11+9\delta_{U,*}}{11(1+\delta_{U,*})},
}
\]
\[
\boxed{
F_*:=
\frac{2\chi_{0,*}}{1+\delta_{U,*}}
+
\frac{4\epsilon_{W,*}\delta_{U,*}}
{11(1-\epsilon_*)(1+\delta_{U,*})^2}.
}
\]
Then the exact first-order nontracking coordinate is the logarithmic drift of the direct microscopic monomial
\[
\boxed{
\mathfrak C_{{\rm nt},*}
:=
\frac{Z_W}{\Omega_W^2}\,
\epsilon_W^{E_*}\,
\delta_U^{-F_*}.
}
\]
Indeed,
\[
\delta\ln \mathfrak C_{{\rm nt},*}
=
\delta\ln\!\left(\frac{Z_W}{\Omega_W^2}\right)
+
E_*\,\delta\ln\epsilon_W
-
F_*\,\delta\ln\delta_U
=
\Sigma_{\rm nt}.
\]
So
\[
\boxed{
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt}.
}
\]

### Explicit microscopic form

Using the coherent-kernel definitions,
\[
\boxed{
\mathfrak C_{{\rm nt},*}
=
\frac{\lambda_W^2\mu_W}
{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}
\left(
\frac{\gamma^2\lambda_W^2\sigma}
{K_U K_W^{(\mathrm{eff})}}
\right)^{E_*}
\left(
\frac{\pi^2T_U}{L^2K_U}
\right)^{-F_*}.
}
\]

This is a sharper result than Stage 167 gave:
the nontracking coordinate has now lost all explicit dependence on the support variables
\(\lambda_\phi, K_\phi^{(\mathrm{eff})}\),
and it also carries no direct \(\Sigma_\chi\) term.
After exact tracking-feed-through subtraction, the nontracking coordinate is controlled only by the three microscopic ratios
\[
\frac{Z_W}{\Omega_W^2},
\qquad
\epsilon_W,
\qquad
\delta_U.
\]

---

## 4. Direct microscopic dressing coordinate

The dressing coordinate was already direct:
\[
\boxed{
\mathfrak D:=\epsilon_\eta
=
\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}}.
}
\]
So
\[
\boxed{
\delta\ln \mathfrak D
=
\delta\ln\epsilon_\eta
=
\Sigma_\eta.
}
\]

This is the simplest of the three branch-adapted coordinates.

---

## 5. Observable triangular law in microscopic monomials

Stage 166 already proved the exact triangular observable structure
\[
\Theta_1=-C_{{\rm tr},*}\Sigma_{\rm tr},
\qquad
\Xi_1=A_{{\rm tr},*}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,\Sigma_\eta,
\]
with
\[
C_{{\rm tr},*}
=
\frac{\chi_{0,*}\delta_{U,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})},
\]
\[
A_{{\rm tr},*}
=
\frac{2\chi_{0,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})}.
\]

Using the direct microscopic monomials above, this becomes
\[
\boxed{
\Theta_1
=
-C_{{\rm tr},*}\,\delta\ln \mathfrak C_{{\rm tr},*},
}
\]
\[
\boxed{
\Xi_1
=
A_{{\rm tr},*}\,\delta\ln \mathfrak C_{{\rm tr},*}
+
\delta\ln \mathfrak C_{{\rm nt},*},
}
\]
\[
\boxed{
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,
\delta\ln \epsilon_\eta.
}
\]

So the direct weak-axisymmetric branch equations are now microscopic monomial equations.

---

## 6. First-order zero-defect microscopic compatibility ledger

Because
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\Sigma_{\rm tr}=\Sigma_{\rm nt}=\Sigma_\eta=0,
\]
Stage 168 turns the first-order zero-defect theorem on the coherent reference branch into the explicit microscopic compatibility system
\[
\boxed{
(1+\chi_{0,*})(\tau_1-\kappa_U)
+
(1+\delta_{U,*})(\gamma_1+c_1-\kappa_U)
=0,
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
E_*\,(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W)
-
F_*\;(\tau_1-\kappa_U)
=0.
}
\]

These can be solved directly for the three dependent drifts \((\tau_1,\kappa_\eta,\mu_1)\):

### 6.1 Tracking compatibility
\[
\boxed{
\tau_1
=
\kappa_U
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}
\bigl(\gamma_1+c_1-\kappa_U\bigr).
}
\]

### 6.2 Dressing compatibility
\[
\boxed{
\kappa_\eta = 2c_1-\kappa_U.
}
\]

### 6.3 Nontracking compatibility
\[
\boxed{
\mu_1
=
\kappa_\eta+2\kappa_W-2\lambda_1
-
E_*\,(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W)
+
F_*(\tau_1-\kappa_U).
}
\]
Equivalently, after eliminating \(\tau_1\) and \(\kappa_\eta\),
\[
\boxed{
\mu_1
=
2c_1-\kappa_U+2\kappa_W-2\lambda_1
-
E_*\,(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W)
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}
\bigl(\gamma_1+c_1-\kappa_U\bigr).
}
\]

So the zero-defect branch is now a three-equation microscopic rigidity ledger rather than an abstract slippage statement.

---

## 7. What Stage 168 changes

Stage 167 reduced the continuation point to the drifts of three exact branch composites. That was already a major simplification.

Stage 168 sharpens the target one more time:

1. the tracking coordinate becomes a single monomial in \((\chi_0,\delta_U)\);
2. the nontracking coordinate becomes a single monomial in \((Z_W/\Omega_W^2,\epsilon_W,\delta_U)\);
3. the dressing coordinate remains \(\epsilon_\eta\);
4. and the full zero-defect theorem becomes an explicit linear compatibility ledger for the microscopic grouped drifts
   \[
   (\lambda_1,c_1,\gamma_1,\kappa_U,\kappa_\eta,\kappa_W,\mu_1,\tau_1).
   \]

So the direct next theorem gate is now smaller again.

It is no longer

> compute the first drifts of \(R_{\rm tr}\), \(\mathcal T^2 R_{\rm tr}^{B_*}\), and \(\epsilon_\eta\).

It is now

> determine whether the true grouped weak-axisymmetric branch preserves, to first order about the coherent reference branch, the three direct microscopic monomials
> \[
> \mathfrak C_{{\rm tr},*},
> \qquad
> \mathfrak C_{{\rm nt},*},
> \qquad
> \epsilon_\eta,
> \]
> or equivalently satisfies the three microscopic compatibility equations above.

That is the clean continuation point.

---

## 8. Best current theorem statement after Stage 168

On the coherent local D/N tracking branch, the first grouped weak-axisymmetric defect is no longer most naturally described by the observable composites \((R_{\rm tr},\mathcal T^2R_{\rm tr}^{B_*},\epsilon_\eta)\). At the same linearized reference-branch order it is exactly equivalent to the logarithmic drift of three direct microscopic monomials:
\[
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt},
\qquad
\delta\ln\epsilon_\eta=\Sigma_\eta.
\]
Here
\[
\mathfrak C_{{\rm tr},*}
=
\chi_0^{\,1+\delta_{U,*}}\delta_U^{\,1+\chi_{0,*}},
\]
\[
\mathfrak C_{{\rm nt},*}
=
\frac{Z_W}{\Omega_W^2}\,
\epsilon_W^{E_*}\,
\delta_U^{-F_*},
\]
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}},
\]
with \(E_*,F_*\) as above.

Therefore the full coherent weak-axisymmetric zero-defect theorem is now equivalent to the linearized invariance of three direct microscopic kernel monomials:
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

Equivalently, the grouped weak-axisymmetric branch at this same linearized reference-branch order must satisfy a three-equation microscopic compatibility ledger for the drifts of
\[
(\lambda_W,c_{\eta U},\gamma,K_U,K_\eta^{(\mathrm{eff})},K_W^{(\mathrm{eff})},\mu_W,T_U).
\]

So the continuation point is now smaller than Stage 167 left it:

> determine whether the true grouped weak-axisymmetric branch is microscopically monomial-rigid, at this same first-order reference-branch level, in the three direct coordinates
> \(\mathfrak C_{{\rm tr},*}\), \(\mathfrak C_{{\rm nt},*}\), and \(\epsilon_\eta\).
