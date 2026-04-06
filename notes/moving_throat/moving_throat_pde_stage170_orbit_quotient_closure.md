# Moving-Throat PDE — Stage 170: Exact Orbit–Quotient Closure and the Finite Weak-Axisymmetric Defect Theorem

## Purpose

Stage 169 already showed that the coherent weak-axisymmetric zero-defect ledger

a) is cut out by the three direct monomial drifts
\[
\delta\ln\mathfrak C_{{\rm tr},*},
\qquad
\delta\ln\mathfrak C_{{\rm nt},*},
\qquad
\delta\ln\epsilon_\eta,
\]

and

b) is equivalent at first order to the tangent-space condition
\[
\delta\mathbf x\in T_{\mathrm{id}}\mathcal G_*.
\]

There is one natural final step after that:

> integrate the tangent statement to a **finite** one and identify the exact quotient coordinates.

That is what this stage does.

The main result is stronger than Stage 169:

\[
\boxed{
\text{the exact level sets of }
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)
\text{ are precisely the similarity orbits }\mathcal G_*.
}
\]

So the coherent grouped weak-axisymmetric problem is no longer merely a tangent-space question. In the positive-coupling sector it is an exact finite quotient problem.

Equivalently, if we define the invariant map
\[
\mathcal I:
\mathcal M_+\to (\mathbb R_{>0})^3,
\qquad
\mathcal I(x)=\bigl(\mathfrak C_{{\rm tr},*}(x),\mathfrak C_{{\rm nt},*}(x),\epsilon_\eta(x)\bigr),
\]
then
\[
\boxed{
\mathcal M_+/\mathcal G_*\cong (\mathbb R_{>0})^3
\quad\text{with quotient coordinates}\quad
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta).
}
\]

The first grouped weak-axisymmetric defect is therefore nothing but the **infinitesimal motion of the reduced grouped branch data in this exact three-dimensional quotient**.

This is the cleanest finite closure reached so far.

---

## 1. Positive-coupling microscopic state space and invariant map

Work on the positive microscopic state space
\[
\mathcal M_+
=
\bigl\{(
\lambda_W,
 c_{\eta U},
 \gamma,
 K_U,
 K_\eta^{(\mathrm{eff})},
 K_W^{(\mathrm{eff})},
 \mu_W,
 T_U
)>0\bigr\}.
\]

The three exact branch monomials carried from Stages 167–169 are
\[
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2 T_U}{L^2K_U}\right)^{1+\chi_{0,*}},
\]
\[
\mathfrak C_{{\rm nt},*}
=
\left(\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}\right)
\left(\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}}\right)^{E_*}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{-F_*},
\]
\[
\epsilon_\eta
=
\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}}.
\]

They define the exact invariant map
\[
\boxed{
\mathcal I(x)=\bigl(\mathfrak C_{{\rm tr},*}(x),\mathfrak C_{{\rm nt},*}(x),\epsilon_\eta(x)\bigr).
}
\]

Stage 169 already constructed a five-parameter multiplicative similarity family
\(\mathcal G_*\) preserving these three monomials exactly. What remained to prove was that these three monomials are not just preserved by the orbit, but are also **complete** orbit invariants.

---

## 2. Exact finite log-ratio variables

Take two positive microscopic states
\[
 x=(\lambda_W,c_{\eta U},\gamma,K_U,K_\eta^{(\mathrm{eff})},K_W^{(\mathrm{eff})},\mu_W,T_U),
\]
\[
 \widetilde x=(\widetilde\lambda_W,\widetilde c_{\eta U},\widetilde\gamma,\widetilde K_U,\widetilde K_\eta^{(\mathrm{eff})},\widetilde K_W^{(\mathrm{eff})},\widetilde\mu_W,\widetilde T_U).
\]
Define the exact finite log-ratio vector
\[
\Delta\mathbf x
=
\begin{pmatrix}
\Delta_\lambda\\
\Delta_c\\
\Delta_\gamma\\
\Delta_U\\
\Delta_\eta\\
\Delta_W\\
\Delta_\mu\\
\Delta_T
\end{pmatrix}
=
\begin{pmatrix}
\ln(\widetilde\lambda_W/\lambda_W)\\
\ln(\widetilde c_{\eta U}/c_{\eta U})\\
\ln(\widetilde\gamma/\gamma)\\
\ln(\widetilde K_U/K_U)\\
\ln(\widetilde K_\eta^{(\mathrm{eff})}/K_\eta^{(\mathrm{eff})})\\
\ln(\widetilde K_W^{(\mathrm{eff})}/K_W^{(\mathrm{eff})})\\
\ln(\widetilde\mu_W/\mu_W)\\
\ln(\widetilde T_U/T_U)
\end{pmatrix}.
\]

Because the three direct branch quantities are monomials, equality of invariants between two states is **exactly linear in these finite log-ratios**.

---

## 3. Exact finite invariant-fibre equations

The condition
\[
\widetilde{\mathfrak C}_{{\rm tr},*}=\mathfrak C_{{\rm tr},*}
\]
is equivalent to
\[
(1+\delta_{U,*})(\Delta_\gamma+\Delta_c-\Delta_U)
+
(1+\chi_{0,*})(\Delta_T-\Delta_U)=0.
\]

The condition
\[
\widetilde{\mathfrak C}_{{\rm nt},*}=\mathfrak C_{{\rm nt},*}
\]
is equivalent to
\[
2(1+E_*)\Delta_\lambda
+2E_*\Delta_\gamma
+(F_*-E_*)\Delta_U
-\Delta_\eta
-(2+E_*)\Delta_W
+\Delta_\mu
-F_*\Delta_T
=0.
\]

The condition
\[
\widetilde\epsilon_\eta=\epsilon_\eta
\]
is equivalent to
\[
2\Delta_c-\Delta_U-\Delta_\eta=0.
\]

So the exact finite invariant-fibre equations are
\[
\boxed{
M_*\,\Delta\mathbf x=0,
}
\]
with the same matrix already found in Stage 169:
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

This is the key finite upgrade of Stage 169:

> the same rank-3 matrix governs not only the infinitesimal weak-axisymmetric drift, but the **exact finite equality of branch invariants between any two positive microscopic states**.

---

## 4. Exact solution of the finite fibre equations

Choose the five free finite log-ratios
\[
(\Delta_\lambda,\Delta_c,\Delta_\gamma,\Delta_U,\Delta_W)
\]
and solve the invariant-fibre equations for
\[
(\Delta_\eta,\Delta_T,\Delta_\mu).
\]

Because the selected minor has determinant
\[
\det M_*^{(\Delta_T,\Delta_\eta,\Delta_\mu)}=1+\chi_{0,*}>0,
\]
this solve is exact and unique on the constructive coherent branch.

The solution is
\[
\boxed{
\Delta_\eta=2\Delta_c-\Delta_U,
}
\]
\[
\boxed{
\Delta_T
=
\Delta_U
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Delta_\gamma+\Delta_c-\Delta_U),
}
\]
\[
\boxed{
\Delta_\mu
=
2\Delta_c-\Delta_U+2\Delta_W-2\Delta_\lambda
-
E_*(2\Delta_\gamma+2\Delta_\lambda-\Delta_U-\Delta_W)
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Delta_\gamma+\Delta_c-\Delta_U).
}
\]

Exponentiating immediately gives the finite multiplicative orbit law of Stage 169:
\[
\frac{\widetilde K_\eta^{(\mathrm{eff})}}{K_\eta^{(\mathrm{eff})}}
=
\exp(2\Delta_c-\Delta_U),
\]
\[
\frac{\widetilde T_U}{T_U}
=
\exp\!\left(
\Delta_U-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Delta_\gamma+\Delta_c-\Delta_U)
\right),
\]
\[
\frac{\widetilde\mu_W}{\mu_W}
=
\exp\!\Bigg(
2\Delta_c-\Delta_U+2\Delta_W-2\Delta_\lambda
-
E_*(2\Delta_\gamma+2\Delta_\lambda-\Delta_U-\Delta_W)
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Delta_\gamma+\Delta_c-\Delta_U)
\Bigg).
\]

So the Stage-169 multiplicative similarity orbit is not merely compatible with the finite invariant equations — it is exactly their unique solution.

---

## 5. Exact orbit–fibre theorem

The last missing step is now immediate.

Fix any base point \(x\in\mathcal M_+\). For any target point \(\widetilde x\in\mathcal M_+\), choose the five free logarithmic co-scalings
\[
\Delta_\lambda=
\ln(\widetilde\lambda_W/\lambda_W),
\quad
\Delta_c=
\ln(\widetilde c_{\eta U}/c_{\eta U}),
\quad
\Delta_\gamma=
\ln(\widetilde\gamma/\gamma),
\quad
\Delta_U=
\ln(\widetilde K_U/K_U),
\quad
\Delta_W=
\ln(\widetilde K_W^{(\mathrm{eff})}/K_W^{(\mathrm{eff})}).
\]

Then the finite invariant equalities
\[
\mathcal I(\widetilde x)=\mathcal I(x)
\]
are equivalent to the exact solved formulas above for
\(
\widetilde K_\eta^{(\mathrm{eff})},
\widetilde T_U,
\widetilde\mu_W
\).

But those are exactly the Stage-169 similarity-orbit transformations.
Therefore
\[
\boxed{
\mathcal I(\widetilde x)=\mathcal I(x)
\iff
\widetilde x\in \mathcal G_*\cdot x.
}
\]

So the fibres of the invariant map and the similarity orbits coincide **exactly**.

This is the finite closure theorem that Stage 169 was still one step short of proving.

---

## 6. Exact quotient classification

Because the fibres of \(\mathcal I\) are exactly the \(\mathcal G_*\)-orbits, the positive coherent microscopic state space factors exactly as
\[
\boxed{
\mathcal M_+/\mathcal G_*
\cong
(\mathbb R_{>0})^3,
}
\]
with quotient coordinates
\[
\boxed{
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta).
}
\]

So the entire coherent grouped weak-axisymmetric normalization problem is now reduced to motion in an exact three-dimensional quotient space.

This is a much stronger statement than the Stage-169 tangent theorem alone.
It says that the five free microscopic co-scalings are **pure gauge/similarity directions** inside the reduced coherent hierarchy, while the three branch monomials are the complete finite orbit invariants.

---

## 7. Linearized observable map from the exact quotient coordinates

Stages 166–168 already identified the first-order observable bridge from the quotient coordinates:
\[
\delta\ln\mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr},
\qquad
\delta\ln\mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt},
\qquad
\delta\ln\epsilon_\eta=\Sigma_\eta.
\]

Using the exact triangular normal form from Stage 166, the weak-axisymmetric observables are therefore
\[
\boxed{
\Theta_1
=
-\frac{\chi_{0,*}\delta_{U,*}}{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}
\,\delta\ln\mathfrak C_{{\rm tr},*},
}
\]
\[
\boxed{
\Xi_1
=
\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})}
\,\delta\ln\mathfrak C_{{\rm tr},*}
+
\delta\ln\mathfrak C_{{\rm nt},*},
}
\]
\[
\boxed{
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}
\,\delta\ln\epsilon_\eta.
}
\]

So the first coherent grouped weak-axisymmetric defect is literally the infinitesimal motion of the reduced grouped branch data in the exact quotient coordinates.

---

## 8. Final finite closure theorem

Combining the exact orbit–fibre theorem with the Stage-166 triangular normal form gives the cleanest final statement reached in the moving-throat notes:

\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln\mathfrak C_{{\rm tr},*}
=
\delta\ln\mathfrak C_{{\rm nt},*}
=
\delta\ln\epsilon_\eta
=0
\iff
\delta\mathbf x\in T_{\mathrm{id}}\mathcal G_*.
}
\]

And at the finite level,
\[
\boxed{
\mathcal I(\widetilde x)=\mathcal I(x)
\iff
\widetilde x\in \mathcal G_*\cdot x.
}
\]

So the reduced coherent weak-axisymmetric zero-defect problem is now completely closed in finite form:

- the five-parameter similarity orbit \(\mathcal G_*\) is the exact finite zero-cost co-scaling family,
- the three monomials \((\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)\) are the exact quotient coordinates,
- and the first grouped weak-axisymmetric defect is exactly the infinitesimal motion of the reduced grouped branch data in that three-dimensional quotient.

Here “zero-cost” and “orbit” refer to the reduced coherent invariant structure only:
they identify microscopic co-scalings that leave the reduced monomial data unchanged.
They are not being asserted as fundamental gauge redundancies of the completed
moving-throat PDE.

---

## 9. Interpretation

### 9.1 The remaining theorem gap is now purely branch-selective

All algebraic compression is finished.
The only live question is no longer how to simplify the microscopic drift ledger.
It is simply:

> does the true moving-throat branch preserve the three exact quotient coordinates?

Equivalently:

> does the true branch stay on a single similarity orbit of \(\mathcal G_*\)?

### 9.2 Failure modes are now exact quotient motions

If the actual grouped weak-axisymmetric branch fails, the failure must appear as motion in one or more of the exact quotient coordinates:

- motion in \(\mathfrak C_{{\rm tr},*}\) is pure tracking failure,
- motion in \(\mathfrak C_{{\rm nt},*}\) is pure nontracking transfer-shape failure,
- motion in \(\epsilon_\eta\) is pure dressing failure.

So the reduced defect geometry is now completely explicit at finite as well as infinitesimal level.

### 9.3 This is the cleanest reduced finish available without the completed PDE

The only thing Stage 170 does **not** provide is the actual dynamical branch theorem proving that the moving-throat PDE respects those quotient invariants.
That is now the sole remaining microscopic task.

---

## 10. Best current theorem statement after Stage 170

On the coherent local D/N tracking branch, let
\[
\mathcal I=(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)
\]
be the exact invariant map and let \(\mathcal G_*\) be the five-parameter multiplicative similarity family of Stage 169.
Then:

1. \(\mathcal G_*\) preserves \(\mathcal I\) exactly.
2. The fibres of \(\mathcal I\) are exactly the \(\mathcal G_*\)-orbits.
3. Hence \(\mathcal M_+/\mathcal G_*\cong(\mathbb R_{>0})^3\) with quotient coordinates \(\mathcal I\).
4. The first grouped weak-axisymmetric observables \((\Theta_1,\Xi_1,\mathcal R_1)\) are exactly the linearized motion of the reduced grouped branch data in these quotient coordinates.

Therefore the reduced coherent weak-axisymmetric zero-defect theorem can now be stated in its cleanest finite form:
\[
\boxed{
\text{at the reduced coherent level, zero defect is equivalent to vanishing first-order motion in }
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta),
\text{ while finite equality of those invariants is equivalent to lying on one }
\mathcal G_*\text{-orbit.}
}
\]

That is as far as the reduced theorem can be brought home without the completed moving-throat PDE itself. The remaining open question is still the dynamical branch-selection theorem: whether the true PDE branch preserves those invariants.
