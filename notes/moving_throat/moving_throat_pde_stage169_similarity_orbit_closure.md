
# Moving-Throat PDE — Stage 169: Exact Microscopic Similarity Orbit and the Final Weak-Axisymmetric Closure Theorem

## Purpose

Stage 168 reduced the coherent weak-axisymmetric problem to the logarithmic drift of three direct microscopic monomials,
\[
\mathfrak C_{{\rm tr},*},
\qquad
\mathfrak C_{{\rm nt},*},
\qquad
\epsilon_\eta,
\]
or equivalently to the three exact compatibility equations
\[
\Sigma_{\rm tr}=0,
\qquad
\Sigma_{\rm nt}=0,
\qquad
\Sigma_\eta=0.
\]

That was already very sharp, but it still left the closure in the form of three linear relations among eight microscopic grouped drifts. There is one more natural step:

> identify the exact finite microscopic transformation family whose tangent space is cut out by those three equations.

That is what this stage does.

The main result is that the full coherent weak-axisymmetric zero-defect problem is **exactly** the tangent problem of a finite five-parameter multiplicative similarity orbit acting on the microscopic kernel couplings
\[
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U).
\]

Equivalently, the Stage-168 compatibility equations are not an isolated fine-tuning ledger. They are the tangent-space equations of an exact codimension-3 similarity family that preserves the three direct monomials
\[
\mathfrak C_{{\rm tr},*},
\qquad
\mathfrak C_{{\rm nt},*},
\qquad
\epsilon_\eta.
\]

So within the coherent reduced hierarchy, the weak-axisymmetric defect is now completely characterized by one geometric statement:

\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\mathbf x \in T_{\mathrm{id}}\mathcal G_*,
}
\]
where \(\mathcal G_*\) is the exact five-parameter monomial-preserving similarity orbit defined below and \(\delta\mathbf x\) is the microscopic grouped weak-axisymmetric log-drift vector.

This is the cleanest “bring-it-home” closure reached so far.

---

## 1. Direct microscopic monomials from Stage 168

Stage 168 already reduced the three branch-adapted coordinates to the logarithmic drifts
\[
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt},
\qquad
\delta\ln \epsilon_\eta=\Sigma_\eta,
\]
with
\[
\mathfrak C_{{\rm tr},*}
=
\chi_0^{\,1+\delta_{U,*}}\,
\delta_U^{\,1+\chi_{0,*}},
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
and
\[
E_*=\frac{2\epsilon_{W,*}}{1-\epsilon_*}\frac{11+9\delta_{U,*}}{11(1+\delta_{U,*})},
\qquad
F_*=\frac{2\chi_{0,*}}{1+\delta_{U,*}}
+\frac{4\epsilon_{W,*}\delta_{U,*}}{11(1-\epsilon_*)(1+\delta_{U,*})^2}.
\]

Using the coherent branch definitions
\[
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\delta_U=\frac{\pi^2 T_U}{L^2 K_U},
\qquad
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}},
\]
\[
Z_W=\frac{\lambda_W^2}{K_\eta^{(\mathrm{eff})}K_W^{(\mathrm{eff})}},
\qquad
\Omega_W^2=\frac{K_W^{(\mathrm{eff})}}{\mu_W},
\]
the direct monomials can be written explicitly as
\[
\boxed{
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}},
}
\]
\[
\boxed{
\mathfrak C_{{\rm nt},*}
=
\left(\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}\right)
\left(\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}}\right)^{E_*}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{-F_*}.
}
\]
Since the grouped weak-axisymmetric branch keeps the scalar axial scale \(L\) fixed at first order in the present reduction, \(L\) is inert here.

---

## 2. The exact monomial-drift map as a rank-3 linear operator

Introduce the microscopic grouped weak-axisymmetric log-drift vector
\[
\delta\mathbf x
=
\begin{pmatrix}
\lambda_1\\
c_1\\
\gamma_1\\
\kappa_U\\
\kappa_\eta\\
\kappa_W\\
\mu_1\\
\tau_1
\end{pmatrix}
=
\begin{pmatrix}
\delta\ln\lambda_W\\
\delta\ln c_{\eta U}\\
\delta\ln\gamma\\
\delta\ln K_U\\
\delta\ln K_\eta^{(\mathrm{eff})}\\
\delta\ln K_W^{(\mathrm{eff})}\\
\delta\ln\mu_W\\
\delta\ln T_U
\end{pmatrix}_{\!\rm grp}.
\]

Then the three direct monomial drifts are exactly
\[
\begin{pmatrix}
\delta\ln\mathfrak C_{{\rm tr},*}\\[2pt]
\delta\ln\mathfrak C_{{\rm nt},*}\\[2pt]
\delta\ln\epsilon_\eta
\end{pmatrix}
=
M_*
\,
\delta\mathbf x,
\]
with
\[
\boxed{
M_*=
\begin{pmatrix}
0 & 1+\delta_{U,*} & 1+\delta_{U,*} & -(2+\chi_{0,*}+\delta_{U,*}) & 0 & 0 & 0 & 1+\chi_{0,*}\\[4pt]
2(1+E_*) & 0 & 2E_* & F_*-E_* & -1 & -(2+E_*) & 1 & -F_*\\[4pt]
0 & 2 & 0 & -1 & -1 & 0 & 0 & 0
\end{pmatrix}.
}
\]

The three rows are simply the logarithmic exponent vectors of
\(\mathfrak C_{{\rm tr},*}\), \(\mathfrak C_{{\rm nt},*}\), and \(\epsilon_\eta\).

A convenient \(3\times3\) minor is obtained by selecting the dependent variables
\[
(\tau_1,\kappa_\eta,\mu_1).
\]
Its determinant is
\[
\boxed{
\det M_*^{(\tau,\kappa_\eta,\mu_1)} = 1+\chi_{0,*} > 0.
}
\]
So the monomial-drift map has exact rank three on the constructive coherent branch \(\chi_{0,*}>0\). Therefore
\[
\boxed{
\dim\ker M_* = 8-3 = 5.
}
\]

This is the exact linear reason the zero-defect branch is a codimension-3 manifold rather than an isolated point.

---

## 3. Exact finite microscopic similarity orbit

Choose five free logarithmic co-scaling parameters
\[
(\Lambda,\ C,\ \Gamma,\ U,\ W)
\]
for
\[
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})}).
\]
Define the finite multiplicative transformation
\[
\lambda_W \mapsto e^{\Lambda}\lambda_W,
\qquad
c_{\eta U} \mapsto e^{C}c_{\eta U},
\qquad
\gamma \mapsto e^{\Gamma}\gamma,
\qquad
K_U \mapsto e^{U}K_U,
\qquad
K_W^{(\mathrm{eff})} \mapsto e^{W}K_W^{(\mathrm{eff})},
\]
and determine the remaining three scalings by exact monomial preservation:
\[
\boxed{
K_\eta^{(\mathrm{eff})}
\mapsto
e^{\,2C-U}\,K_\eta^{(\mathrm{eff})},
}
\]
\[
\boxed{
T_U
\mapsto
e^{\,U-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U)}\,T_U,
}
\]
\[
\boxed{
\mu_W
\mapsto
e^{\,M(\Lambda,C,\Gamma,U,W)}\,\mu_W,
}
\]
with
\[
\boxed{
M(\Lambda,C,\Gamma,U,W)
=
2C-U+2W-2\Lambda
-
E_*(2\Gamma+2\Lambda-U-W)
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U).
}
\]

Call the resulting five-parameter multiplicative family \(\mathcal G_*\).

---

## 4. Exact preservation of the three direct monomials

By construction, the orbit \(\mathcal G_*\) preserves the three direct branch monomials exactly:

### 4.1 Tracking monomial
\[
\delta_{\mathcal G_*}\ln\mathfrak C_{{\rm tr},*}
=
(1+\delta_{U,*})(\Gamma+C-U)
+
(1+\chi_{0,*})\!\left[
U-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U)-U
\right]
=0.
\]

### 4.2 Dressing monomial
\[
\delta_{\mathcal G_*}\ln\epsilon_\eta
=
2C-U-(2C-U)=0.
\]

### 4.3 Nontracking monomial
Using
\[
\mathfrak C_{{\rm nt},*}
=
\left(\frac{Z_W}{\Omega_W^2}\right)\epsilon_W^{E_*}\delta_U^{-F_*},
\]
the chosen exponent \(M(\Lambda,C,\Gamma,U,W)\) gives
\[
\delta_{\mathcal G_*}\ln\mathfrak C_{{\rm nt},*}=0.
\]

So \(\mathcal G_*\) is an exact monomial-preserving similarity orbit, not merely a first-order ansatz.

---

## 5. Tangent-space equivalence to the Stage-168 compatibility ledger

Linearizing the finite orbit at the identity,
\[
(\Lambda,C,\Gamma,U,W) \mapsto 0,
\]
gives the tangent-space relations
\[
\boxed{
\kappa_\eta = 2c_1-\kappa_U,
}
\]
\[
\boxed{
\tau_1
=
\kappa_U
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\gamma_1+c_1-\kappa_U),
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

These are exactly the Stage-168 microscopic compatibility formulas. Therefore
\[
\boxed{
\ker M_* = T_{\mathrm{id}}\mathcal G_*.
}
\]

So the Stage-168 linear zero-defect ledger is not an isolated set of three equations. It is the tangent-space equation of the exact five-parameter similarity orbit \(\mathcal G_*\).

---

## 6. Final weak-axisymmetric closure theorem

Stage 166 already proved
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\Sigma_{\rm tr}=\Sigma_{\rm nt}=\Sigma_\eta=0,
\]
and Stage 168 proved
\[
\Sigma_{\rm tr}
=
\delta\ln\mathfrak C_{{\rm tr},*},
\qquad
\Sigma_{\rm nt}
=
\delta\ln\mathfrak C_{{\rm nt},*},
\qquad
\Sigma_\eta
=
\delta\ln\epsilon_\eta.
\]
Combining these with the rank-3 monomial-drift map above gives the exact closure theorem:

\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
M_*\,\delta\mathbf x=0
\iff
\delta\mathbf x \in T_{\mathrm{id}}\mathcal G_*.
}
\]

Equivalently, the coherent weak-axisymmetric defect vanishes exactly when the microscopic grouped drift is tangent to the finite monomial-preserving similarity orbit \(\mathcal G_*\).

This is the sharpest final formulation of the reduced theorem reached so far.

---

## 7. Interpretation

The importance of this stage is conceptual as much as algebraic.

### 7.1 The zero-defect branch is not fine-tuning

Because \(\ker M_*\) has dimension five, the reduced coherent weak-axisymmetric zero-defect set is a codimension-3 similarity manifold, not an isolated point. So within the coherent reduced hierarchy, exact first-order defect cancellation is a substantial orbit of allowed co-scalings rather than a one-shot numerical accident.

### 7.2 The last real question is geometric, not algebraic

All purely algebraic compression is now finished. What remains is no longer “solve another drift identity.” It is:

> does the actual moving-throat branch follow a tangent direction inside the exact similarity orbit \(\mathcal G_*\)?

That is a genuine microscopic branch-selection question.

### 7.3 The direct microscopic meaning of failure is now explicit

If the actual grouped weak-axisymmetric moving-throat branch leaves \(T_{\mathrm{id}}\mathcal G_*\), then the failure is immediately localized:

- motion normal to the first row of \(M_*\) drives tracking failure,
- motion normal to the second row drives nontracking transfer-shape failure,
- motion normal to the third row drives dressing failure.

So the reduced defect geometry is completely transparent.

---

## 8. Best current theorem statement after Stage 169

On the coherent local D/N tracking branch, the entire first grouped weak-axisymmetric normalization problem is equivalent to the action of an exact five-parameter multiplicative similarity orbit \(\mathcal G_*\) on the microscopic kernel couplings
\[
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U).
\]
This orbit preserves the three direct branch monomials
\[
\mathfrak C_{{\rm tr},*},
\qquad
\mathfrak C_{{\rm nt},*},
\qquad
\epsilon_\eta
\]
exactly, and its tangent space is precisely the Stage-168 compatibility ledger. Therefore the full coherent weak-axisymmetric zero-defect theorem can be written in the compact geometric form
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\mathbf x \in T_{\mathrm{id}}\mathcal G_*.
}
\]

So the continuation point is now as small and as clean as it can be without the completed moving-throat PDE itself:

> determine whether the true moving-throat grouped weak-axisymmetric branch is tangent to the exact monomial-preserving similarity orbit \(\mathcal G_*\).

If it is, the coherent first-order grouped weak-axisymmetric defect vanishes automatically.

