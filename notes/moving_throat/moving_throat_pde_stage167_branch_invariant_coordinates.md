# Moving-Throat PDE — Stage 167: Exact Branch-Invariant Coordinates for the Coherent Weak-Axisymmetric Defect

## Purpose

Stage 166 compressed the coherent weak-axisymmetric defect to the three branch-adapted scalars
\[
(\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta),
\]
but it still presented them through the observable drift ledger
\[
(\Theta_1,\Xi_1,\mathcal R_1).
\]
That is already useful, but it is not yet the cleanest possible continuation point.

The next honest step is now extremely small:

> derive the three branch-adapted coordinates **directly from the exact coherent branch equations themselves**, so that the remaining theorem gate is phrased in terms of exact branch composites rather than reconstructed drift variables.

That is what this stage does.

The main result is that the full coherent weak-axisymmetric problem is equivalent to the linearized invariance of three exact branch composites:
\[
R_{\rm tr},
\qquad
\mathfrak N_*,
\qquad
\epsilon_\eta,
\]
where
\[
R_{\rm tr}=\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
\]
\[
\mathfrak N_*:=\mathcal T^2\,R_{\rm tr}^{B_*},
\qquad
B_*:=\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}},
\]
and the starred quantities are the coherent reference-branch values about which the weak-axisymmetric grouped perturbation is taken.

Equivalently, after one harmless positive normalization of the tracking factor,
\[
\mathfrak T_*:=R_{\rm tr}^{-C_*},
\qquad
C_*:=\frac{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}{\chi_{0,*}\delta_{U,*}},
\]
we obtain the exact logarithmic defect coordinates
\[
\boxed{\delta\ln \mathfrak T_* = \Sigma_{\rm tr},}
\qquad
\boxed{\delta\ln \mathfrak N_* = \Sigma_{\rm nt},}
\qquad
\boxed{\delta\ln \epsilon_\eta = \Sigma_\eta.}
\]

So the whole continuation point is now smaller than Stage 166 left it:

the actual moving-throat branch no longer needs to supply five microscopic slippages, or even three reconstructed slippages. It only needs the first grouped drifts of three exact branch composites.

---

## 1. Exact coherent branch equations to be used

The coherent local D/N branch already supplies three exact structural relations.

### 1.1 Tracking factor

The exact coherent tracking factor is
\[
\boxed{
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}
=
\frac{1+\chi_0+\delta_U}{(1+\chi_0)(1+\delta_U)}.
}
\]
By definition of the weak-axisymmetric tracking drift,
\[
\boxed{
\delta\ln R_{{\rm tr},A}=\epsilon\lambda_A\,\Theta_1.
}
\]
So the tracking factor itself is already one exact branch observable.

### 1.2 Effective transfer shape

The exact coherent transfer shape is
\[
\boxed{
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
\qquad
\epsilon=\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right].
}
\]
By definition of the grouped weak-axisymmetric defect,
\[
\boxed{
\delta\ln \mathcal T_A^2=\epsilon\lambda_A\,\Xi_1.
}
\]
So the transfer shape is the exact branch observable whose logarithmic drift is the direct grouped defect.

### 1.3 Selected-branch demand identity

The selected-branch and direct-transfer descriptions meet through the exact identity
\[
\boxed{
R_{\rm target}\,\mathcal T^2
=
\Lambda_0(1-\epsilon_\eta),
\qquad
\Lambda_0:=\frac{27\pi^2Gc_s^5}{20a^5c^5},
}
\]
where the front factor \(\Lambda_0\) is inert at linear grouped order on the isotropic branch.

Therefore
\[
\boxed{
\delta\ln(R_{{\rm target},A}\mathcal T_A^2)
=
\delta\ln(1-\epsilon_{\eta,A})
=
\epsilon\lambda_A\,(\mathcal R_1+\Xi_1).
}
\]
Since
\[
\delta\ln \epsilon_{\eta,A}=\epsilon\lambda_A\,\Sigma_\eta,
\]
this immediately reproduces the Stage-166 dressing relation
\[
\boxed{
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\Sigma_\eta.
}
\]

So the coherent branch already contains one exact direct dressing coordinate:
\[
\boxed{\delta\ln\epsilon_\eta=\Sigma_\eta.}
\]

---

## 2. Direct branch-invariant tracking coordinate

Stage 166 showed that
\[
\Sigma_{\rm tr}
=
-\frac{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}{\chi_0\delta_U}\,\Theta_1.
\]
On the coherent reference branch, the prefactor is a positive constant. Define
\[
\boxed{
C_*:=\frac{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}{\chi_{0,*}\delta_{U,*}}>0,
}
\]
and the normalized tracking invariant
\[
\boxed{
\mathfrak T_*:=R_{\rm tr}^{-C_*}.
}
\]
Then
\[
\delta\ln \mathfrak T_*
=
-C_*\,\delta\ln R_{\rm tr}
=
-C_*\Theta_1
=
\Sigma_{\rm tr}.
\]
So the tracking coordinate is now a direct branch-invariant logarithmic drift:
\[
\boxed{
\delta\ln \mathfrak T_* = \Sigma_{\rm tr}.
}
\]

### Interpretation

Nothing essentially new had to be added. Stage 166 already showed that \(\Sigma_{\rm tr}\) is just a positively normalized version of the tracking-factor drift. This stage only packages that fact into an exact branch-invariant coordinate.

Because \(C_*>0\), the equivalent zero-tracking statements are
\[
\boxed{
\Sigma_{\rm tr}=0
\iff
\delta\ln \mathfrak T_*=0
\iff
\delta\ln R_{\rm tr}=0
\iff
\Theta_1=0.
}
\]

So the first branch equation the actual moving-throat solution must preserve is simply the tracking factor itself.

---

## 3. Direct branch-invariant nontracking coordinate

Stage 166 also showed that
\[
\Sigma_{\rm nt}
=
\Xi_1+\frac{2(1+\chi_0+\delta_U)}{\delta_U}\,\Theta_1.
\]
Freeze the background coefficient at the coherent reference branch,
\[
\boxed{
B_*:=\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}}.
}
\]
Now define the exact branch composite
\[
\boxed{
\mathfrak N_*:=\mathcal T^2\,R_{\rm tr}^{B_*}.
}
\]
Its linearized logarithmic drift is
\[
\delta\ln\mathfrak N_*
=
\delta\ln\mathcal T^2 + B_*\,\delta\ln R_{\rm tr}
=
\Xi_1+B_*\Theta_1
=
\Sigma_{\rm nt}.
\]
So the second branch-adapted coordinate is also a direct exact branch-invariant drift:
\[
\boxed{
\delta\ln\mathfrak N_* = \Sigma_{\rm nt}.
}
\]

### Explicit microscopic form

Using the coherent branch formulas,
\[
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
\qquad
R_{\rm tr}=
\frac{1+\chi_0+\delta_U}{(1+\chi_0)(1+\delta_U)},
\]
this composite can be written directly as
\[
\boxed{
\mathfrak N_*
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2}
\left[
\frac{1+\chi_0+\delta_U}{(1+\chi_0)(1+\delta_U)}
\right]^{B_*}.
}
\]
This is support-blind for the same reason \(\mathcal T^2\) itself is support-blind: no coherent-support factor \(\zeta\) appears.

### Interpretation

The direct grouped defect \(\Xi_1\) is not yet the genuine nontracking branch coordinate, because it still carries a universal tracking feed-through. The composite \(\mathfrak N_*\) subtracts that feed-through exactly. So the second branch equation the completed moving-throat branch must preserve is not just the transfer shape \(\mathcal T^2\), but the corrected nontracking composite \(\mathfrak N_*\).

---

## 4. Direct branch-invariant dressing coordinate

The third branch-adapted coordinate is already direct at the microscopic level:
\[
\boxed{
\Sigma_\eta = \delta\ln\epsilon_\eta.
}
\]
So we may simply take the third branch invariant to be
\[
\boxed{
\mathfrak D:=\epsilon_\eta.
}
\]
Then
\[
\boxed{
\delta\ln\mathfrak D = \Sigma_\eta.
}
\]

It is also useful to define the exact selected-branch composite
\[
\boxed{
\mathfrak E:=\frac{R_{\rm target}\,\mathcal T^2}{\Lambda_0}=1-\epsilon_\eta.
}
\]
Its logarithmic drift is
\[
\boxed{
\delta\ln\mathfrak E
=
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\Sigma_\eta.
}
\]
So the direct dressing coordinate can be read in two equivalent ways:

- as the logarithmic drift of the microscopic ratio \(\epsilon_\eta\), or
- as the logarithmic drift of the selected-branch composite \(\mathfrak E=1-\epsilon_\eta\).

### Interpretation

The important gain is conceptual: the selected-branch residual is no longer some independent extra datum. It is just the complementary drift of the same exact microscopic dressing variable.

---

## 5. Composite rigidity theorem in direct branch variables

We now have three exact branch-invariant coordinates:
\[
\boxed{
\delta\ln \mathfrak T_* = \Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak N_* = \Sigma_{\rm nt},
\qquad
\delta\ln \mathfrak D = \Sigma_\eta.
}
\]
Since Stage 166 already proved
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\Sigma_{\rm tr}=\Sigma_{\rm nt}=\Sigma_\eta=0,
\]
we can rewrite the full triple-rigidity condition as
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln \mathfrak T_*=
\delta\ln \mathfrak N_*=
\delta\ln \mathfrak D=0.
}
\]
Equivalently, because \(C_*>0\), this is the same as
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln R_{\rm tr}=0,
\qquad
\delta\ln \mathfrak N_*=0,
\qquad
\delta\ln \epsilon_\eta=0.
}
\]

This is the sharpest form of the continuation point so far.

The actual moving-throat branch no longer needs to be tested against a broad set of microscopic slippages. It only needs to answer three exact first-order questions:

1. **tracking:** is the coherent tracking factor \(R_{\rm tr}\) invariant?
2. **nontracking transfer shape:** is the corrected composite \(\mathfrak N_*\) invariant?
3. **dressing:** is the microscopic dressing ratio \(\epsilon_\eta\) invariant?

If the answer to all three is yes, then the coherent weak-axisymmetric defect vanishes completely at first order.

---

## 6. Best current theorem statement after Stage 167

On the coherent local D/N tracking branch, the exact branch equations already contain three direct composite quantities whose first grouped weak-axisymmetric drifts are exactly the three Stage-166 branch-adapted defect coordinates:
\[
\delta\ln \mathfrak T_* = \Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak N_* = \Sigma_{\rm nt},
\qquad
\delta\ln \epsilon_\eta = \Sigma_\eta.
\]
Here
\[
\mathfrak T_*:=R_{\rm tr}^{-C_*},
\qquad
\mathfrak N_*:=\mathcal T^2 R_{\rm tr}^{B_*},
\qquad
C_*>0,
\qquad
B_*=rac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}}.
\]
Moreover,
\[
\frac{R_{\rm target}\mathcal T^2}{\Lambda_0}=1-\epsilon_\eta,
\qquad
\Lambda_0=\frac{27\pi^2Gc_s^5}{20a^5c^5},
\]
so the selected-branch residual is just the complementary drift of the same dressing variable.

Therefore the full coherent weak-axisymmetric zero-defect theorem is now equivalent to the linearized invariance of three exact branch composites:
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln R_{\rm tr}=0,
\quad
\delta\ln(\mathcal T^2 R_{\rm tr}^{B_*})=0,
\quad
\delta\ln\epsilon_\eta=0.
}
\]

So the direct next theorem gate is now even smaller than Stage 166 suggested.

It is no longer

> compute the three branch-adapted coordinates \((\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta)\).

It is now

> compute the first grouped weak-axisymmetric drifts of the exact branch composites
> \[
> R_{\rm tr},
> \qquad
> \mathcal T^2 R_{\rm tr}^{B_*},
> \qquad
> \epsilon_\eta,
> \]
> on the actual moving-throat branch.

That is the direct continuation point.
