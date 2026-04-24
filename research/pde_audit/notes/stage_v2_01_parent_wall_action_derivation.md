# Stage V2-01 — Parent-Action and Throat-Action Audit

## Executive result

This stage tests whether the moving-throat variable is already a **dynamical parent field** or whether it is currently an **effective wall/interface closure** attached to the exact GNLS/Maxwell parent.

The result is:

\[
\boxed{
\text{The current stated parent action gives a wall force, not an autonomous wall PDE.}
}
\]

More precisely:

- If the parent action is only
  \[
  S_{\rm parent}=\int dt\,d^4X\,(\mathcal L_\psi+\mathcal L_{\rm EM}),
  \]
  with geometry entering through
  \[
  V_{\rm conf}(\mathbf X;\Sigma),
  \]
  then \(\Sigma\) or \(R\) appears as a coupling argument, but not as an independent dynamical field with its own inertia/stiffness.

- If the quadratic distributed wall action
  \[
  S_\eta^{(2)}
  =
  \frac12\int dt\,dw\,d\Omega\,
  \Big[
  \mu_\eta \eta_t^2
  -T_w\eta_w^2
  -T_\Omega\eta(-\Delta_{S^2})\eta
  -K_\eta\eta^2
  \Big]
  \]
  is included, then the advertised moving-wall PDE, boundary terms, modal split, and positivity gates are mathematically consistent.

So V2-01 has a **split verdict**:

\[
\boxed{
\text{Strict parent-level moving-throat claim: not yet passed.}
}
\]

\[
\boxed{
\text{Linear effective wall-closure algebra: passed.}
}
\]

This is not a fatal flaw. It is a status correction. The program should either:

1. promote \(S_\eta\), or better a nonlinear \(S_\Sigma\), into the declared parent action, or  
2. explicitly label the moving-wall PDE as an effective closure until the parent throat action is derived.

---

## Source anchors

The compact program currently lists the core fields as \(\psi,A_M,\Sigma\) or \(R\), but its exact parent action is written as GNLS plus localized Maxwell, with geometry entering through \(V_{\rm conf}(\mathbf X;\Sigma)\). The same compact ledger also tags the moving-throat lift as a closure hierarchy rather than a completed branch theorem.

The Stage-1 geometry lift explicitly introduces the distributed wall action as a **new ansatz** and says the existing parent theory has exact GNLS, exact localized Maxwell, and a geometry sector previously carried only by \(a(t),L(t)\) through \(V_{\rm conf}\). The full derivation ledger later shows that the distributed wall lift is compatible with the old \(a,L\) closure, but still treats the later matter/gauge/outgoing realizations as reduced theorem work rather than a solved nonlinear moving-throat theorem.

---

## 1. Audit target

The ideal parent-level statement would be

\[
S_{\rm total}=S_\psi+S_{\rm EM}+S_\Sigma,
\]

with \(S_\Sigma\) supplying the dynamics of the moving throat.

The current compact parent action instead has the schematic form

\[
S_{\rm current}
=
\int dt\,d^4X\,
\left[
\mathcal L_\psi(\psi,A;\Sigma)
+
\mathcal L_{\rm EM}(A)
\right],
\]

where

\[
\mathcal L_\psi
=
\frac{i\hbar}{2}
(\psi^\ast D_t\psi-\psi D_t\psi^\ast)
-\frac{\hbar^2}{2m}(D_i\psi)^\ast(D_i\psi)
-V_{\rm conf}(\mathbf X;\Sigma)\rho
-U(\rho).
\]

Thus \(\Sigma\) appears through \(V_{\rm conf}\), but there is no term of the form

\[
\Sigma_t^2,\quad
|\nabla\Sigma|^2,\quad
\eta_t^2,\quad
\eta_w^2,\quad
\eta(-\Delta_{S^2})\eta,
\]

unless the wall action is added separately.

---

## 2. Variation of the confinement-only parent term

Use the moving-surface representation

\[
\Sigma=r-R(\Omega,w,t),
\qquad
R=R_0+\eta.
\]

Take

\[
V_{\rm conf}
=
V_{\rm wall}\!\left(\frac{\Sigma}{\ell_c}\right).
\]

Linearizing around the background gives

\[
\delta\Sigma=-\eta,
\]

so

\[
\delta V_{\rm conf}
=
-\frac{V_{\rm wall}'(\Sigma_0/\ell_c)}{\ell_c}\eta.
\]

Since the matter Lagrangian contains

\[
-\rho V_{\rm conf},
\]

the linearized wall-dependent contribution is

\[
\delta\mathcal L_\psi
=
+\rho\frac{V_{\rm wall}'(\Sigma_0/\ell_c)}{\ell_c}\eta.
\]

Therefore the Euler derivative with respect to \(\eta\) is only

\[
\frac{\delta \mathcal L_\psi}{\delta\eta}
=
\rho\frac{V_{\rm wall}'(\Sigma_0/\ell_c)}{\ell_c}.
\]

There are no terms proportional to

\[
\eta_{tt},\qquad
\partial_w(T_w\eta_w),\qquad
\Delta_{S^2}\eta.
\]

So the confinement-only parent action gives a force/source term. It does **not** give an autonomous moving-throat PDE.

This is the core V2-01 result.

---

## 3. Variation after adding the quadratic wall action

Now add the modal wall Lagrangian

\[
L_{\ell m}
=
\frac12\mu_\eta q_{\ell m,t}^2
-\frac12T_w q_{\ell m,w}^2
-\frac12
\left[
K_\eta+\ell(\ell+1)T_\Omega
\right]q_{\ell m}^2
+S_{\ell m}q_{\ell m}.
\]

The Euler-Lagrange expression is

\[
\frac{\partial L}{\partial q}
-\partial_t\frac{\partial L}{\partial q_t}
-\partial_w\frac{\partial L}{\partial q_w}
=0.
\]

This gives

\[
-\left[
\mu_\eta q_{\ell m,tt}
-\partial_w(T_w q_{\ell m,w})
+
(K_\eta+\ell(\ell+1)T_\Omega)q_{\ell m}
-S_{\ell m}
\right]=0.
\]

Equivalently,

\[
\boxed{
\mu_\eta q_{\ell m,tt}
-\partial_w(T_w q_{\ell m,w})
+
\left[
K_\eta+\ell(\ell+1)T_\Omega
\right]q_{\ell m}
=
S_{\ell m}.
}
\]

Special cases:

For the scalar lane,

\[
\boxed{
\mu_\eta q_{00,tt}
-\partial_w(T_w q_{00,w})
+
K_\eta q_{00}
=
S_{00}.
}
\]

For the grouped real \(P_2\) lane, \(\ell=2\), so \(\ell(\ell+1)=6\):

\[
\boxed{
\mu_\eta q_{2m,tt}
-\partial_w(T_w q_{2m,w})
+
(K_\eta+6T_\Omega)q_{2m}
=
S_{2m}.
}
\]

This verifies that the advertised \(l=0\) and grouped \(P_2\) split follows from \(S_\eta^{(2)}\).

---

## 4. Boundary terms

For the wall action,

\[
\pi_q=\frac{\partial L}{\partial q_t}
=
\mu_\eta q_t,
\]

and

\[
p_w=\frac{\partial L}{\partial q_w}
=
-T_w q_w.
\]

The \(w\)-boundary part of the variation is

\[
\left[p_w\,\delta q\right]_{\partial I_w}
=
\left[-T_w q_w\,\delta q\right]_{\partial I_w}.
\]

Thus:

- Dirichlet mouth data use \(\delta q=0\).
- A free natural endpoint uses
  \[
  T_w q_w=0.
  \]
- A driven mouth/worldtube port should specify the traction conjugate to \(q\), with sign fixed by the outward normal convention.

This is important for later D/N and cap-regularity work: the boundary conditions must be derived from the variational boundary terms, not just imposed as an interval mnemonic.

---

## 5. Canonical energy and positivity gate

For the source-free quadratic wall action,

\[
\mathcal H_{\ell m}
=
\pi_q q_t-L_{\ell m}
\]

gives

\[
\boxed{
\mathcal H_{\ell m}
=
\frac12\mu_\eta q_t^2
+
\frac12T_w q_w^2
+
\frac12
\left[
K_\eta+\ell(\ell+1)T_\Omega
\right]q^2.
}
\]

The local quadratic positivity gate is therefore

\[
\boxed{
\mu_\eta>0,\qquad
T_w>0,\qquad
K_\eta+\ell(\ell+1)T_\Omega\ge 0.
}
\]

For the grouped \(P_2\) lane this becomes

\[
\boxed{
\mu_\eta>0,\qquad
T_w>0,\qquad
K_\eta+6T_\Omega\ge 0.
}
\]

These are only the local quadratic gates. Coupling to BdG, Maxwell, mixed channels, and outgoing ports adds Schur-complement stability conditions in later stages.

---

## 6. Recovery of the old \(a,L\) closure

Use the axisymmetric two-profile ansatz

\[
\eta_0(w,t)
=
2\sqrt{\pi}
\left[
\alpha_a(w)\delta a(t)
+
\alpha_L(w)\delta L(t)
\right].
\]

After integrating over the normalized \(Y_{00}\) angular mode, the reduced Lagrangian has the form

\[
L_{\rm red}^{(0)}
=
\frac12M_{AB}\dot Q^A\dot Q^B
-\frac12K_{AB}Q^AQ^B,
\qquad
Q^A=(\delta a,\delta L).
\]

The matrices are

\[
\boxed{
M_{AB}
=
4\pi\int dw\,
\mu_\eta(w)\alpha_A(w)\alpha_B(w).
}
\]

\[
\boxed{
K_{AB}
=
4\pi\int dw\,
\left[
T_w(w)\alpha_A'(w)\alpha_B'(w)
+
K_\eta(w)\alpha_A(w)\alpha_B(w)
\right].
}
\]

So the old \(a,L\) closure is recovered as a lowest-mode truncation of the distributed wall action, provided \(S_\eta\) is included.

---

## 7. One-mode grouped \(P_2\) reduction

For one real \(P_2\) component,

\[
\eta_{2m}
=
\beta_2(w)q_{2m}(t)Y_{2m}^{\rm real}(\Omega),
\]

with

\[
-\Delta_{S^2}Y_{2m}=6Y_{2m},
\]

the reduced one-mode action gives

\[
L_{2m}
=
\frac12M_2\dot q_{2m}^2
-\frac12K_2 q_{2m}^2,
\]

where

\[
\boxed{
M_2
=
\int dw\,\mu_\eta(w)\beta_2(w)^2.
}
\]

\[
\boxed{
K_2
=
\int dw\,
\left[
T_w(w)\beta_2'(w)^2
+
(K_\eta(w)+6T_\Omega(w))\beta_2(w)^2
\right].
}
\]

This confirms that the grouped real \(P_2\) lane is not an artificial add-on; it is the next harmonic family of the same wall action.

---

## 8. SymPy audit summary

The accompanying script verifies:

1. \(V_{\rm conf}(\Sigma)\) alone produces only an algebraic wall source.
2. The quadratic wall action produces the modal PDE.
3. The \(l=0\) and \(l=2\) modal specializations are correct.
4. The canonical momentum and \(w\)-boundary momentum are correct.
5. The source-free Hamiltonian density has the expected positive quadratic form.
6. The \(a,L\) two-profile reduction gives the \(4\pi\) overlap matrices.
7. The one-mode \(P_2\) reduction gives \(K_\eta+6T_\Omega\).

The script output ends with:

```text
STRICT_PARENT_DYNAMIC_WALL: FAIL unless S_eta/S_Sigma is included in S_total.
EFFECTIVE_LINEAR_WALL_CLOSURE: PASS; S_eta^(2) supplies a consistent linear wall PDE.
PATCH_REQUIRED: promote S_eta to parent status or relabel the moving wall as an effective closure.
```

---

## 9. Recommended patch to the program ledger

I recommend adding this status line to Volume 2:

> The compact program’s \(\Sigma/R\) variable is parent-level as a **confinement-coupling argument**. It is not parent-level as a **dynamical throat field** unless a throat action \(S_\Sigma\) is added to \(S_{\rm total}\). The quadratic \(S_\eta^{(2)}\) is a consistent linear effective wall action and may be promoted as the first approximation to \(S_\Sigma\), but its coefficients must be branch data rather than post-hoc fit parameters.

A clean parent-complete statement would be:

\[
\boxed{
S_{\rm total}
=
S_\psi[\psi,A,\Sigma]
+
S_{\rm EM}[A]
+
S_\Sigma[\Sigma;\mathcal C_\Sigma].
}
\]

Here \(\mathcal C_\Sigma\) denotes the constitutive wall data. At quadratic order around a stationary throat branch,

\[
S_\Sigma
\longrightarrow
S_\eta^{(2)}
+
O(\eta^3).
\]

This distinction should be kept explicit in all later stages.
