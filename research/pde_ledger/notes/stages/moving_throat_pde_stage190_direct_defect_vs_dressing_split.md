# Moving-Throat PDE — Stage 190: Direct Transfer-Shape Defect vs. Dressing Residual, Support-Blindness, and the Scalar No-Go Filters

## Status

**Exact within the carried coherent local D/N weak-axisymmetric reference-branch closure, together with the exact reduced grouped-lane symmetry theorem for pure real `P2` anisotropy.**

This stage does **not** enlarge the endgame packet from Stage 189.
It sharpens the internal logic of that packet by separating:

1. the **direct transfer-shape defect**,
2. the **tracking feed-through** inside that defect,
3. the **selected-branch dressing residual**,
4. and the two main obstruction theorems that prevent the remaining gap from being blamed on the wrong microscopic channels.

---

## Purpose

Stage 188 completed the first-order branch-observable compiler in the direct observables
\[
R_{\rm tr},
\qquad
\mathfrak N_*:=\mathcal T^2 R_{\rm tr}^{B_*},
\qquad
\epsilon_\eta,
\]
and Stage 189 carried that packet all the way to the isotropic grouped outgoing-prefactor compiler.

At that point the remaining algebra was no longer the issue.
The actual missing clarification was structural:

> which microscopic slippages belong to the **direct transfer-shape defect**, which belong only to the **selected-branch dressing residual**, and which channels are ruled out as first-order rescue mechanisms by exact symmetry?

That is what this stage fixes.

The main outputs are:

1. an exact **direct-defect / dressing** split of the coherent weak-axisymmetric problem,
2. an exact **support-blindness theorem** for the direct transfer shape,
3. the exact branch-adapted direct-defect coordinates
   \[
   \Sigma_{\rm tr},\qquad \Sigma_{\rm nt},\qquad \Sigma_\eta,
   \]
   with
   \[
   \Theta_1=-C_{\rm tr}\Sigma_{\rm tr},
   \qquad
   \Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
   \qquad
   \mathcal R_1+\Xi_1=-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta,
   \]
4. exact inverse reconstruction formulas,
5. and the exact no-go filter saying that a **pure grouped real `P2` anisotropy cannot linearly feed the scalar off-bundle slippages**.

So this stage is the clean separation layer between the Stage 189 transfer compiler and the later finite packet theorem gate.

---

## 1. Carry-forward exact transfer-shape and selected-branch identities

### 1.1 Direct coherent transfer shape

Carry forward the exact coherent transfer shape
\[
\boxed{
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2}.
}
\]
The corrected nontracking branch observable from Stage 188 is
\[
\boxed{
\mathfrak N_*:=\mathcal T^2 R_{\rm tr}^{B_*},
\qquad
R_{\rm tr}:=\frac{1+\chi_0+\delta_U}{(1+\chi_0)(1+\delta_U)},
\qquad
B_*:=\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}}.
}
\]
By Stage 189,
\[
\boxed{
\delta\ln\mathcal T^2=\Xi_1,
\qquad
\delta\ln\mathfrak N_*=\Sigma_{\rm nt}.
}
\]
So the direct grouped defect is already the logarithmic drift of one exact transfer shape, while the corrected nontracking coordinate is the logarithmic drift of a feed-through-subtracted branch composite.

### 1.2 Selected-branch identity

Carry forward the exact selected-branch relation
\[
\boxed{
R_{\rm target}\,\mathcal T^2
=
\Lambda_0(1-\epsilon_\eta),
\qquad
\Lambda_0:=\frac{27\pi^2Gc_s^5}{20a^5c^5}.
}
\]
Equivalently,
\[
\boxed{
\mathfrak E:=\frac{R_{\rm target}\,\mathcal T^2}{\Lambda_0}=1-\epsilon_\eta.
}
\]
So the complementary selected-branch observable is not independent. It is just the complementary dressing variable built from the same microscopic ratio \(\epsilon_\eta\).

At first grouped weak-axisymmetric order,
\[
\boxed{
\delta\ln(1-\epsilon_\eta)=\mathcal R_1+\Xi_1,
\qquad
\delta\ln R_{\rm target}=\mathcal R_1.
}
\]

### 1.3 Immediate structural lesson

The problem already separates at the level of exact observables:

- `\(\Xi_1\)` measures the direct transfer-shape drift,
- `\(\mathcal R_1+\Xi_1\)` measures only the complementary dressing drift,
- and `\(\mathcal R_1\)` is what remains after subtracting the direct defect from the selected-branch demand drift.

That separation will now be pushed down to the microscopic slippage level.

---

## 2. Exact support-blindness theorem for the direct transfer shape

### 2.1 Coherent support enhancement changes the baseline, not the transfer shape

On the coherent branch the support enhancement may be written as
\[
M_{\rm tr}=M_{\rm mix}\,S(\zeta;\epsilon),
\qquad
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon}.
\]
But the transfer shape itself is still
\[
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
\]
with no explicit support factor `\(\zeta\)`.
Therefore
\[
\boxed{
\frac{\partial\mathcal T^2}{\partial\zeta}=0,
\qquad
\frac{\partial\ln\mathcal T^2}{\partial\zeta}=0.
}
\]
The same is true for the selected-branch demand ratio,
\[
R_{\rm target}
=
\Lambda_0\,\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\]
so
\[
\boxed{
\frac{\partial R_{\rm target}}{\partial\zeta}=0,
\qquad
\frac{\partial\ln R_{\rm target}}{\partial\zeta}=0.
}
\]
Hence also
\[
\boxed{
\frac{\partial\mathfrak N_*}{\partial\zeta}=0,
\qquad
\frac{\partial\ln\mathfrak N_*}{\partial\zeta}=0,
\qquad
\partial_\zeta\Xi_1=0.
}
\]

### 2.2 Physical meaning

This is a sharp role separation:

- the coherent support lane can alter the **baseline loading** of the branch,
- but it cannot alter the **effective transfer shape**,
- so it can help satisfy steady normalization,
- yet it cannot cancel the weak-axisymmetric direct defect.

That is the first obstruction theorem promoted explicitly by this stage.

---

## 3. Exact microscopic slippage ledger

The coherent weak-axisymmetric branch depends on the five microscopic logarithmic slippages
\[
\boxed{
\Sigma_Z,
\qquad
\Sigma_\chi,
\qquad
\Sigma_\epsilon,
\qquad
\Sigma_\delta,
\qquad
\Sigma_\eta.
}
\]
They are defined by
\[
\boxed{
\Sigma_\chi=\delta\ln\chi_0,
\qquad
\Sigma_\delta=\delta\ln\delta_U,
\qquad
\Sigma_Z=\delta\ln\!\left(\frac{Z_W}{\Omega_W^2}\right),
}
\]
\[
\boxed{
\Sigma_\epsilon=\delta\ln\epsilon_W,
\qquad
\Sigma_\eta=\delta\ln\epsilon_\eta.
}
\]
Using the coherent-kernel definitions,
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
the microscopic drift forms are
\[
\boxed{
\Sigma_\chi=\gamma_1+c_1-\kappa_U,
\qquad
\Sigma_\delta=\tau_1-\kappa_U,
}
\]
\[
\boxed{
\Sigma_Z=2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W,
\qquad
\Sigma_\epsilon=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W,
}
\]
\[
\boxed{
\Sigma_\eta=2c_1-\kappa_U-\kappa_\eta.
}
\]

The important structural split is now visible already:

- the **direct transfer-shape defect** depends only on
  \[
  (\Sigma_Z,\Sigma_\chi,\Sigma_\epsilon,\Sigma_\delta),
  \]
- the **selected-branch dressing residual** introduces only the additional slippage
  \[
  \Sigma_\eta.
  \]

---

## 4. Exact direct transfer-shape defect law

### 4.1 Direct microscopic law before branch adaptation

The exact coherent direct-defect law is
\[
\boxed{
\Xi_1
=
\Sigma_Z
+\frac{2\chi_0}{1+\chi_0}\Sigma_\chi
+\frac{2\epsilon_W}{1-\epsilon}
\left[
\frac{11+9\delta_U}{11(1+\delta_U)}\Sigma_\epsilon
-\frac{2\delta_U}{11(1+\delta_U)^2}\Sigma_\delta
\right].
}
\]
So the direct grouped defect is already a scalar built purely from the four direct slippages.

### 4.2 Exact tracking combination

Define the exact tracking combination
\[
\boxed{
\Sigma_{\rm tr}:=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi.
}
\]
Then the tracking-factor drift is
\[
\boxed{
\Theta_1=-C_{\rm tr}\Sigma_{\rm tr},
\qquad
C_{\rm tr}:=
\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}.
}
\]

### 4.3 Exact nontracking direct-defect scalar

Now define the branch-adapted nontracking direct-defect coordinate
\[
\boxed{
\Sigma_{\rm nt}
:=
\Sigma_Z
+\frac{2\epsilon_W}{1-\epsilon}\frac{11+9\delta_U}{11(1+\delta_U)}\Sigma_\epsilon
-\left[
\frac{2\chi_0}{1+\delta_U}
+\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
}
\]
Then the exact direct grouped-defect law becomes
\[
\boxed{
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
A_{\rm tr}:=\frac{2\chi_0}{(1+\chi_0)(1+\delta_U)}.
}
\]

So the direct defect packet has a canonical split:

- `\(\Sigma_{\rm tr}\)` = universal tracking feed-through,
- `\(\Sigma_{\rm nt}\)` = genuine nontracking transfer-shape defect.

Two remarks are important.

First, after the exact feed-through subtraction, `\(\Sigma_{\rm nt}\)` carries **no explicit** `\(\Sigma_\chi\)` term.
The interference ratio only survives inside the universal tracking channel.

Second, `\(\Sigma_{\rm nt}\)` is support-blind for the same reason `\(\mathcal T^2\)` and `\(\mathfrak N_*\)` are support-blind: no coherent-support factor enters.

---

## 5. Exact dressing residual and the direct-defect / dressing split

The third branch-adapted coordinate is simply the microscopic dressing slippage
\[
\boxed{
\Sigma_\eta=\delta\ln\epsilon_\eta.
}
\]
The selected-branch residual then satisfies the exact complementary law
\[
\boxed{
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta.
}
\]
So the direct-defect / dressing split is exact:
\[
\boxed{
(\Theta_1,\Xi_1)
\text{ depend only on }
(\Sigma_{\rm tr},\Sigma_{\rm nt}),
\qquad
\mathcal R_1+\Xi_1
\text{ depends only on }
\Sigma_\eta.
}
\]
Equivalently, the observable compiler is block-triangular:
\[
\boxed{
\begin{pmatrix}
\Theta_1\\[2pt]
\Xi_1\\[2pt]
\mathcal R_1+\Xi_1
\end{pmatrix}
=
\begin{pmatrix}
-C_{\rm tr} & 0 & 0\\[2pt]
A_{\rm tr} & 1 & 0\\[2pt]
0 & 0 & -\dfrac{\epsilon_\eta}{1-\epsilon_\eta}
\end{pmatrix}
\begin{pmatrix}
\Sigma_{\rm tr}\\[2pt]
\Sigma_{\rm nt}\\[2pt]
\Sigma_\eta
\end{pmatrix}.
}
\]
The determinant is
\[
\boxed{
\det = \frac{C_{\rm tr}\epsilon_\eta}{1-\epsilon_\eta}>0
}
\]
on the constructive physical branch
\[
\chi_0>0,
\qquad
\delta_U>0,
\qquad
0<\epsilon_\eta<1.
\]
So the split is not heuristic. It is exact and invertible.

### Physical meaning

The selected-branch residual is **not** an additional direct transfer-shape error.
It is just the complementary drift of the same microscopic dressing variable once the direct transfer-shape defect has already been accounted for.

That is the second major structural clarification of this stage.

---

## 6. Exact inverse reconstruction and rigidity statements

Because the direct/dressing compiler is triangular, it inverts exactly.

### 6.1 Tracking coordinate

\[
\boxed{
\Sigma_{\rm tr}
=
-\frac{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}{\chi_0\delta_U}\,\Theta_1.
}
\]
So
\[
\boxed{
\Theta_1=0
\iff
\Sigma_{\rm tr}=0.
}
\]

### 6.2 Nontracking direct-defect coordinate

Using `\(\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt}\)`,
\[
\boxed{
\Sigma_{\rm nt}
=
\Xi_1+\frac{A_{\rm tr}}{C_{\rm tr}}\Theta_1
=
\Xi_1+\frac{2(1+\chi_0+\delta_U)}{\delta_U}\Theta_1.
}
\]
Therefore
\[
\boxed{
\Xi_1=0
\iff
\Sigma_{\rm nt}=-A_{\rm tr}\Sigma_{\rm tr}.
}
\]
So a vanishing grouped defect can still hide a compensating tracking slippage unless tracking is checked separately.
On the tracking-rigid branch,
\[
\boxed{
\Theta_1=0\ \text{and}\ \Xi_1=0
\iff
\Sigma_{\rm tr}=0,
\quad
\Sigma_{\rm nt}=0.
}
\]

### 6.3 Dressing coordinate

\[
\boxed{
\Sigma_\eta
=
-\frac{1-\epsilon_\eta}{\epsilon_\eta}(\mathcal R_1+\Xi_1).
}
\]
Hence
\[
\boxed{
\mathcal R_1+\Xi_1=0
\iff
\Sigma_\eta=0.
}
\]
And if the direct transfer-shape defect already vanishes,
\[
\Xi_1=0,
\]
then
\[
\boxed{
\mathcal R_1=0
\iff
\Sigma_\eta=0.
}
\]

### 6.4 Full triple-rigidity theorem

Combining the three blocks,
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\Sigma_{\rm tr}=\Sigma_{\rm nt}=\Sigma_\eta=0.
}
\]
So the first-order coherent problem has decomposed completely into:

1. a tracking sector,
2. a direct nontracking transfer-shape sector,
3. and a dressing sector.

---

## 7. Exact no-go filter: pure grouped real `P2` anisotropy cannot linearly feed the scalar off-bundle slippages

The second obstruction theorem carried into this stage is the pure grouped-lane scalar no-go.

Let
\[
x=(x_{20},x_{21},x_{22})^T
\]
be any grouped real `P2` bundle, with grouped metric
\[
G_{\rm grp}=\operatorname{diag}(1,2,2).
\]
The weighted trace/anomaly variables are
\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}{5},
\qquad
 a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\qquad
 b_x=\frac{x_{21}-x_{22}}{2}.
\]
The exact grouped quadratic invariant is
\[
\boxed{
\mathcal I[x,y]
=
\frac15\,\delta x^T G_{\rm grp}\,\delta y
=
4a_x a_y+\frac45 b_x b_y,
}
\]
where `\(\delta x:=x-\bar x(1,1,1)^T\)`.

On the weak axisymmetric `\(Y_{20}\)` branch,
\[
 x_{20}=x^{(0)}+\varepsilon x^{(1)},
\qquad
 x_{21}=x^{(0)}+\frac\varepsilon2 x^{(1)},
\qquad
 x_{22}=x^{(0)}-\varepsilon x^{(1)}.
\]
Then
\[
\boxed{
\bar x=x^{(0)},
\qquad
b_x=3a_x,
\qquad
\mathcal I[x,x]=\frac{7}{10}\,\varepsilon^2\,(x^{(1)})^2.
}
\]
So the first weak-axisymmetric grouped anisotropy is already **quadratic** as a scalar invariant.

The representation-theoretic consequence is the exact no-go statement:
\[
\boxed{
\delta^{(1)}_{P_2}\mathcal S=0
}
\]
for every rotational scalar observable `\(\mathcal S\)` extracted from the isotropic branch.
In particular, for the scalar off-bundle slippages of Stage 168,
\[
\boxed{
\varepsilon_L^{(1,P_2)}=0,
\qquad
\varepsilon_v^{(1,P_2)}=0,
\qquad
\varepsilon_T^{(1,P_2)}=0,
}
\]
and therefore
\[
\boxed{
\varepsilon_\perp^{(1,P_2)}=0,
\qquad
\delta_\perp^{(1,P_2)}=0.
}
\]
So a pure grouped real `P2` anisotropy cannot be the first linear source of the scalar normal/off-family defect.
Its scalar feed-down begins only at quadratic order through the grouped invariants `\(\mathcal I[X,Y]\)`.

### Why this matters here

This means the remaining **linear** first-order theorem gate cannot be blamed on a hidden scalar backreaction from pure grouped-lane anisotropy.
At linear order, the grouped-anisotropy problem must live in the direct outlet / transfer-shape data themselves.

That is the second obstruction theorem packaged explicitly by this stage.

---

## 8. Best current theorem statement after Stage 190

After Stages 188–190, the coherent weak-axisymmetric first-order problem is now split as sharply as it can be before the final finite packet compiler.

### 8.1 Direct defect vs dressing

The exact microscopic split is
\[
\boxed{
\Theta_1=-C_{\rm tr}\Sigma_{\rm tr},
\qquad
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
\mathcal R_1+\Xi_1=-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta.
}
\]
So the first-order coherent problem is exactly the zero-set test of
\[
\boxed{
(\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta).
}
\]

### 8.2 Two obstruction theorems are now explicit

1. **Support-blindness:** the direct transfer shape `\(\mathcal T^2\)` and the corrected nontracking composite `\(\mathfrak N_*\)` are blind to the coherent support-enhancement lane. So that lane cannot cancel the direct grouped defect.

2. **No linear scalar feed-down from pure grouped real `P2`:** pure grouped-lane anisotropy cannot linearly generate the scalar off-bundle slippages. Its scalar feed-down begins only at quadratic order.

So two common escape routes are now removed from the first-order theorem gate.

### 8.3 Smallest remaining first-order microscopic gate

The actual moving-throat branch must now answer only three first-order questions:

1. does the tracking coordinate `\(\Sigma_{\rm tr}\)` vanish?
2. does the genuine nontracking direct-defect coordinate `\(\Sigma_{\rm nt}\)` vanish?
3. does the dressing coordinate `\(\Sigma_\eta\)` vanish?

If yes, then
\[
\Theta_1=\Xi_1=\mathcal R_1=0.
\]
If not, the remaining mismatch is already classified as either

- tracking,
- direct transfer-shape,
- or dressing.

That is the cleanest microscopic separation reached so far.

---

## 9. Immediate next stage

The next clean move is now exactly the Stage 191 compiler step.

The purpose of that next stage should be:

1. take the direct/dressing split fixed here,
2. take the grouped operator / transfer data from Stage 189,
3. and compile the **smallest exact finite PDE data packet** that still has to be computed on the actual moving-throat branch.

So Stage 190 is the last structural-cleanup stage before the minimal PDE packet theorem gate.
