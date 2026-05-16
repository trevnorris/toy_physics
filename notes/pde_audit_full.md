stage_v2_01_parent_wall_action_derivation.md
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

stage_v2_02_maxwell_gauge_localization_derivation.md
# Stage V2-02 — Maxwell Gauge-Localization Audit

## Purpose

This stage audits the localized Maxwell sector used by the 4D / moving-throat program, focusing on the fact that the Maxwell kinetic term is localized by a transverse profile \(Z(w)\), while the displayed Lorenz-type gauge-fixing term is not localized.

The current action form is

\[
S_{\rm EM}
=
\int dt\,d^3x\,dw\,
\left[
-\frac{Z(w)}{4\mu_0}F_{MN}F^{MN}
-\frac{1}{2\xi\mu_0}(\partial\!\cdot\!A)^2
-A_MJ^M
\right].
\]

The audit asks whether this is safe for the bulk equations, the zero-mode brane reduction, and the moving-throat mixed-sector response stack.

---

## 1. General weighted gauge-fixing functional

To separate the issue cleanly, introduce a general gauge-fixing weight \(H(w)\):

\[
\mathcal L_{\rm gf}
=
-\frac{H(w)}{2\xi\mu_0}\,B^2,
\qquad
B\equiv\partial_MA^M.
\]

The current papers correspond to

\[
H(w)=1.
\]

The localized patch corresponds to something like

\[
H(w)=Z(w)
\]

or, more generally, any finite localized profile \(H_{\rm loc}(w)\) with finite integral.

Varying the gauge-fixing term gives

\[
\delta S_{\rm gf}
=
-\frac{1}{\xi\mu_0}
\int H B\,\partial_M\delta A^M
=
\frac{1}{\xi\mu_0}
\int \partial_M(HB)\,\delta A^M,
\]

so the Maxwell equation becomes

\[
\boxed{
\partial_M\bigl(ZF^{MN}\bigr)
+
\frac{1}{\xi}\partial^N\bigl(HB\bigr)
=
\mu_0J^N.
}
\]

For \(H=1\), this reduces to the current displayed equation:

\[
\partial_M\bigl(ZF^{MN}\bigr)
+
\frac{1}{\xi}\partial^N B
=
\mu_0J^N.
\]

For \(H=Z\), it becomes

\[
\partial_M\bigl(ZF^{MN}\bigr)
+
\frac{1}{\xi}\partial^N\bigl(ZB\bigr)
=
\mu_0J^N.
\]

---

## 2. Divergence consistency

Because \(ZF^{MN}\) is antisymmetric in \(M,N\), while \(\partial_M\partial_N\) is symmetric,

\[
\partial_N\partial_M\bigl(ZF^{MN}\bigr)=0.
\]

The script verifies this explicitly in a symbolic antisymmetric subblock.

Taking the divergence of the general equation gives

\[
\boxed{
\frac{1}{\xi}\Box_5\bigl(HB\bigr)
=
\mu_0\partial_NJ^N.
}
\]

So for conserved sources,

\[
\Box_5(HB)=0.
\]

For \(H=1\), this is

\[
\Box_5B=0.
\]

For \(H=Z\), it is

\[
\Box_5(ZB)=0.
\]

The two are not identical as gauge-propagation equations. Since \(H=Z\) depends on \(w\),

\[
\Box_5(HB)
=
H\Box_5B+2H'\partial_wB+H''B.
\]

The script verifies the identity

\[
\Box(HB)-\left(H\Box B+2H'B_w+H''B\right)=0.
\]

This is not a fatal issue. Gauge fixing is not a physical observable. But it means that switching from \(H=1\) to \(H=Z\) is not a purely cosmetic rewrite; it changes the gauge-condition propagation while preserving gauge-invariant observables when handled consistently.

---

## 3. Zero-mode brane reduction

Take the usual zero-mode ansatz

\[
A_\mu(x,w)\simeq a_\mu(x),
\qquad
A_w=0,
\qquad
\partial_wA_\mu\simeq0.
\]

Then

\[
F_{\mu\nu}(x,w)=f_{\mu\nu}(x),
\qquad
F_{\mu w}=0,
\qquad
B=\partial_\mu a^\mu.
\]

The Maxwell kinetic term reduces to

\[
S_{\rm kin}^{(0)}
=
-\frac{Z_{\rm int}}{4\mu_0}
\int d^4x\,f_{\mu\nu}f^{\mu\nu},
\qquad
Z_{\rm int}=\int_{-\infty}^{\infty}Z(w)\,dw.
\]

For the Gaussian profile

\[
Z(w)=e^{-w^2/\lambda^2},
\]

SymPy verifies

\[
\boxed{Z_{\rm int}=\sqrt\pi\,\lambda.}
\]

The effective coupling is therefore

\[
\mu_0^{\rm eff}=\frac{\mu_0}{Z_{\rm int}}.
\]

### 3.1 Unweighted gauge fixing

For the current unweighted term \(H=1\), the zero-mode gauge-fixing contribution on a finite regulator interval \([-R,R]\) is

\[
S_{\rm gf,unweighted}^{(0)}(R)
=
-\frac{2R}{2\xi\mu_0}
\int d^4x\,\left(\partial_\mu a^\mu\right)^2.
\]

Thus

\[
\frac{H_{\rm int}}{Z_{\rm int}}
=
\frac{2R}{\sqrt\pi\lambda}
\longrightarrow \infty
\qquad(R\to\infty).
\]

So the unweighted gauge-fixing action diverges for a noncompact zero mode unless

\[
\partial_\mu a^\mu=0
\]

is imposed before reduction.

Equivalently, if one tries to interpret the finite-box unweighted term as a 4D gauge-fixing term with

\[
-\frac{Z_{\rm int}}{2\xi_4\mu_0}(\partial\!\cdot a)^2,
\]

then matching gives

\[
\xi_4=\xi\frac{Z_{\rm int}}{2R}
\longrightarrow 0.
\]

So the noncompact unweighted zero-mode limit is a singular Landau-type gauge limit, not an ordinary finite 4D gauge-fixed action.

### 3.2 Weighted localized gauge fixing

For \(H=Z\), the zero-mode gauge-fixing term is instead

\[
S_{\rm gf,weighted}^{(0)}
=
-\frac{Z_{\rm int}}{2\xi\mu_0}
\int d^4x\,\left(\partial_\mu a^\mu\right)^2.
\]

This is finite and has the same \(Z_{\rm int}\) normalization as the Maxwell kinetic term.

For any finite gauge weight \(H_{\rm int}=\int H(w)dw\), the 4D gauge parameter obeys

\[
\boxed{
\xi_4=\xi\frac{Z_{\rm int}}{H_{\rm int}}.
}
\]

The script verifies the matching identity

\[
\frac{H_{\rm int}}{\xi}=\frac{Z_{\rm int}}{\xi_4}.
\]

For \(H=Z\), one gets \(\xi_4=\xi\).

---

## 4. Mixed-sector gauge invariants

The moving-throat program must keep the mixed sector alive:

\[
A_w,
\qquad
F_{\mu w},
\qquad
J^w,
\qquad
E_w,
\qquad
C_a.
\]

Using the convention

\[
A_0\to A_0-\partial_t\chi,
\qquad
A_a\to A_a+\partial_a\chi,
\qquad
A_w\to A_w+\partial_w\chi,
\]

the mixed fields are

\[
E_w=-\partial_tA_w-\partial_wA_0,
\]

\[
C_a=\partial_aA_w-\partial_wA_a.
\]

The script verifies

\[
\delta E_w=0,
\qquad
\delta C_a=0.
\]

So the mixed fields are exact gauge-invariant observables. They should not be treated as artifacts of the gauge-fixing choice.

---

## 5. Stage verdict

### 5.1 What passes

The bulk equation with the current unweighted \(H=1\) Lorenz term is mathematically consistent as a formal gauge-fixed bulk equation:

\[
\partial_M(ZF^{MN})+rac1\xi\partial^N(\partial\!\cdot A)=\mu_0J^N.
\]

Its divergence identity is also consistent when the total current is conserved:

\[
\Box_5(\partial\!\cdot A)=0.
\]

The mixed-sector observables \(E_w\) and \(C_a\) are gauge invariant and survive independently of this choice.

### 5.2 What fails

The unweighted gauge-fixing term fails as a noncompact zero-mode gauge-fixed action. For \(A_\mu(x,w)=a_\mu(x)\), the Maxwell kinetic coefficient is finite:

\[
Z_{\rm int}=\sqrt\pi\lambda,
\]

but the unweighted gauge-fixing coefficient is

\[
\int_{-\infty}^{\infty}dw=\infty.
\]

Therefore the current unweighted term cannot be naively reduced to a standard finite 3+1 gauge-fixed Maxwell action by replacing the integral with \(Z_{\rm int}\).

### 5.3 Safe interpretations

There are two safe ways to state the program.

#### Option A — keep the current unweighted gauge fixing as a bulk gauge device

Use \(H=1\) only at the bulk gauge-fixed equation level, impose Lorenz gauge

\[
\partial\!\cdot A=0
\]

before the zero-mode brane reduction, and then choose any convenient 3+1 gauge fixing after the reduced Maxwell theory has been obtained.

In this reading, the clean zero-mode Maxwell equation is obtained from the gauge-invariant part:

\[
\partial_\mu f^{\mu\nu}=\frac{\mu_0}{Z_{\rm int}}J_b^\nu.
\]

The unweighted gauge-fixing term is not part of the reduced finite action.

#### Option B — localize the gauge-fixing functional

Replace the gauge-fixing term by

\[
\mathcal L_{\rm gf}^{\rm loc}
=
-\frac{H_{\rm loc}(w)}{2\xi\mu_0}(\partial\!\cdot A)^2,
\]

with finite

\[
H_{\rm int}=\int H_{\rm loc}(w)dw.
\]

The simplest choice is

\[
H_{\rm loc}=Z.
\]

Then the zero-mode gauge-fixed action is finite and has the same localization normalization as the kinetic term. The cost is that the gauge-condition propagation becomes

\[
\Box_5(Z\partial\!\cdot A)=0
\]

for conserved sources.

---

## 6. Recommended Volume 2 patch

The cleanest patch for the derivation ledger is to state the Maxwell sector in one of these two forms:

### Conservative patch

Keep the already-published equation, but add the warning:

> The unweighted Lorenz term is a bulk gauge-fixing device. In the zero-mode brane reduction, the Lorenz condition is imposed before integrating over the noncompact transverse coordinate. A finite 3+1 gauge-fixing term may then be chosen in the reduced brane Maxwell theory. The unweighted five-dimensional gauge-fixing term should not be integrated over the noncompact zero mode unless \(\partial_\mu a^\mu=0\).

### Structural patch

Use the localized gauge-fixing action

\[
\boxed{
\mathcal L_{\rm gf}^{\rm loc}
=
-\frac{Z(w)}{2\xi\mu_0}(\partial\!\cdot A)^2.
}
\]

Then the brane zero mode reduces cleanly to a finite 3+1 gauge-fixed Maxwell action with \(\xi_4=\xi\).

My recommendation is to use the conservative patch immediately for consistency with the existing papers, then mark the structural patch as the preferred parent-action cleanup for any future fully gauge-fixed PDE derivation.

---

## 7. Final status line

\[
\boxed{
\text{The unweighted Maxwell gauge-fixing term is not fatal, but it is unsafe for naive zero-mode action reduction.}
}
\]

A correct Volume 2 statement is:

\[
\boxed{
\text{Either impose Lorenz gauge before reducing, or localize the gauge-fixing functional.}
}
\]

stage_v2_03_dimension_ledger_derivation.md
# Stage V2-03 — Dimension Ledger and Normalization Audit

## Theorem target

This stage audits the unit consistency of the normalization packets used by the moving-throat grouped `P2` and outgoing-quadrupole bridge. The main targets are

\[
P_0^{\rm target}=\frac{54G c_s^5}{5a^5c^5},
\]

\[
\Lambda_0=\frac{27\pi^2G c_s^5}{20a^5c^5},
\]

and the invariant outgoing-quadrupole condition

\[
\widehat m_0^2 P_0\frac{a^5}{27c_s^5}=\frac{2G}{5c^5}.
\]

The audit also checks the grouped-response conversion formulas, constant-prefactor branch equations, EM zero-mode coupling dimensions, support-wave mass theorem, and the raw Maxwell/mixed one-lane transfer factor.

## Base dimensions

The script uses the exact base vector

\[
(M,L,T,Q,O),
\]

where `O` is an abstract reduced wall/operator unit. The `O` axis is kept separate because the absolute dimension of the reduced wall operator depends on the chosen normalization of the wall amplitude. Physical observables should not retain an uncancelled `O` unless they are still operator-valued.

The primitive assignments are

\[
[G]=M^{-1}L^3T^{-2},
\qquad
[c]=[c_s]=LT^{-1},
\qquad
[a]=L,
\qquad
[\omega]=T^{-1},
\]

\[
[\hbar]=ML^2T^{-1},
\qquad
[E]=ML^2T^{-2}.
\]

For the EM zero-mode ledger the script assumes a dimensionless localization profile `Z(w)`, so

\[
[Z_{\rm int}]=L.
\]

Then

\[
\mu_0^{\rm eff}=\frac{\mu_0^{(5)}}{Z_{\rm int}}
\]

is consistent if

\[
[\mu_0^{(5)}]=[\mu_0^{\rm eff}]L.
\]

Similarly,

\[
q_{\rm eff}=\frac{q_\star}{\sqrt{Z_{\rm int}}}
\]

is consistent if

\[
[q_\star]=[q_{\rm eff}]L^{1/2}.
\]

## Main target dimensions

The target prefactor has dimension

\[
\left[\frac{G c_s^5}{a^5c^5}\right]
= M^{-1}L^{-2}T^{-2}.
\]

Therefore

\[
[P_0^{\rm target}]=[\Lambda_0]=M^{-1}L^{-2}T^{-2}.
\]

The outgoing GR quadrupole coefficient has dimension

\[
\left[\frac{G}{c^5}\right]=M^{-1}L^{-2}T^3.
\]

The compact outgoing `l=2` factor has

\[
\left[\frac{a^5}{c_s^5}\right]=T^5.
\]

Therefore

\[
\left[P_0^{\rm target}\frac{a^5}{c_s^5}\right]
= M^{-1}L^{-2}T^3
=\left[\frac{G}{c^5}\right].
\]

So the central 2.5PN/4PN normalization product is dimensionally consistent.

## Outgoing fingerprint

The normalized outgoing branch is

\[
\widehat Y_2^{\rm out}(\omega)
=1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\frac{a^5\omega^5}{27c_s^5}+O(\omega^6).
\]

The script verifies

\[
\left[\frac{a^2\omega^2}{c_s^2}\right]
=
\left[\frac{a^4\omega^4}{c_s^4}\right]
=
\left[\frac{a^5\omega^5}{c_s^5}\right]
=1.
\]

So the outgoing fingerprint is dimensionless, as required.

## Operator and response moments

Let

\[
D(\omega)=D_0+D_2\omega^2+D_4\omega^4+O(\omega^6).
\]

The script assigns

\[
[D_0]=O,
\qquad
[D_2]=OT^2,
\qquad
[D_4]=OT^4,
\]

so that every term in `D(omega)` has unit `O`.

The normalized response

\[
Y(\omega)=\frac{D_0}{D(\omega)}
=1+u_2\omega^2+u_4\omega^4+O(\omega^6)
\]

has

\[
[u_2]=T^2,
\qquad
[u_4]=T^4.
\]

SymPy verifies the exact coefficients

\[
u_2=-\frac{D_2}{D_0},
\]

\[
 u_4=\frac{D_2^2-D_0D_4}{D_0^2}.
\]

The one-pole condition

\[
 u_4=4u_2^2
\]

is also dimensionally consistent because both sides have unit `T^4`.

## Outgoing prefactor moments

The Stage-5 prefactor is

\[
\mathrm{Pref}(\omega)
=\frac{D_0N(\omega)}{D(\omega)^2},
\qquad
N(\omega)=N_0+N_2\omega^2+N_4\omega^4+O(\omega^6).
\]

The gravitationally normalized assignment is

\[
[N_0]=[P_0]O,
\qquad
[N_2]=[P_0]OT^2,
\qquad
[N_4]=[P_0]OT^4.
\]

Then

\[
[P_0]=M^{-1}L^{-2}T^{-2},
\qquad
[P_2]=M^{-1}L^{-2},
\qquad
[P_4]=M^{-1}L^{-2}T^2.
\]

SymPy verifies

\[
P_0=\frac{N_0}{D_0},
\]

\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]

\[
P_4=\frac{D_0^2N_4-2D_0(D_2N_2+D_4N_0)+3D_2^2N_0}{D_0^3}.
\]

The constant-prefactor branch equations are dimensionally consistent:

\[
N_2=\frac{2D_2N_0}{D_0},
\]

\[
N_4=\frac{2D_0(D_2N_2+D_4N_0)-3D_2^2N_0}{D_0^2}.
\]

## Isotropic one-pole surface

The isotropic full-bundle one-pole surface is

\[
D_0(B_4+Z_4)=3(M+B_2+Z_2)^2.
\]

With

\[
[B_4+Z_4]=OT^4,
\qquad
[M+B_2+Z_2]=OT^2,
\]

both sides have unit

\[
O^2T^4.
\]

So the one-pole target surface is dimensionally consistent.

## Raw Maxwell/mixed transfer warning

The only nontrivial warning appears when the Stage-4 one-lane Maxwell/mixed transfer factor is interpreted as a raw canonical mechanical transfer.

For canonical internal coordinates `U,W`, the reduced kinetic terms imply

\[
[U]=[W]=E^{1/2}T.
\]

If the wall displacement coordinate has unit `L`, then the wall operator has

\[
[D_{\rm mech}]=E/L^2=M.
\]

The raw one-lane transfer factor

\[
N_0^{\rm raw}
=\frac{(\Omega_U^2g_W+Rg_U)^2}{\Delta^2}
\]

has dimension

\[
[N_0^{\rm raw}]=[D_{\rm mech}]T^2.
\]

Therefore

\[
P_0^{\rm raw}=\frac{N_0^{\rm raw}}{D_{\rm mech}}
\]

has unit

\[
[P_0^{\rm raw}]=T^2,
\]

not

\[
[P_0^{\rm target}]=M^{-1}L^{-2}T^{-2}.
\]

This is not a contradiction in the reduced algebra. It means the raw canonical transfer factor is not yet the gravitationally normalized `P0` used in the 2.5PN/4PN target. A port/source normalization scale is required:

\[
P_0^{\rm grav}=\mathcal S_{\rm port}P_0^{\rm raw},
\]

with

\[
[\mathcal S_{\rm port}]
=M^{-1}L^{-2}T^{-4}.
\]

Equivalently,

\[
N_0^{\rm grav}=\mathcal S_{\rm port}N_0^{\rm raw}.
\]

Volume 2 should therefore avoid identifying raw canonical `N0/D0` with the gravitational outgoing prefactor unless this port/source normalization is explicitly included.

## SymPy audit summary

The script produced:

- 29 exact dimension checks,
- 29 passes,
- 8 SymPy series-identity checks,
- 8 passes.

Key output dimensions:

```text
P0_target_dim = M^-1 L^-2 T^-2
Gamma_GR_dim = M^-1 L^-2 T^3
D0_dim = O
N0_normalized_dim = M^-1 L^-2 T^-2 O
raw_mixed_P0_dim = T^2
required_raw_to_grav_bridge_scale_dim = M^-1 L^-2 T^-4
```

## Verdict

```text
PASS for published target identities and response-moment algebra.
WARN that the raw canonical Maxwell/mixed transfer factor is a mechanical T^2 object and must be multiplied by an explicit port/source normalization scale to become the gravitational P0 used in the 2.5PN/4PN target.
```

## Carry-forward patch

Introduce an explicit normalization bridge in the Volume 2 notation:

\[
\boxed{
P_0^{\rm grav}
=\mathcal S_{\rm port}\,P_0^{\rm raw},
\qquad
[\mathcal S_{\rm port}]=M^{-1}L^{-2}T^{-4}.
}
\]

Then the actual theorem target becomes

\[
\boxed{
\widehat m_0^2\mathcal S_{\rm port}P_0^{\rm raw}
=\frac{54Gc_s^5}{5a^5c^5}.
}
\]

If later branch work defines `N0` directly in the gravitationally normalized port convention, then `S_port` is implicitly already absorbed. The derivation should say which convention is being used.

stage_v2_04_open_junction_impedance_derivation.md
# Stage V2-04 — Open-Junction / Organ-Pipe Impedance Audit

## Executive verdict

This audit **accepts the architectural patch** from a capped throat to an open finite-radius conduit:

\[
R(L)>0,
\]

with the bottom of the throat opening into unconfined four-spatial-dimensional bulk. A hard cap,

\[
R(L)=0,
\]

is incompatible with a nonzero baseline/DC superfluid throughput, because the exit area vanishes and a fixed flux would require a divergent flux density.

The audit gives a split result for the D/N support ladder:

\[
\textbf{D/N survives only if the support variable is the free-end flow/displacement variable.}
\]

A sudden expansion into a large bulk reservoir is a **low-impedance open end**. For a pressure-like or potential-like scalar, that open end reflects with phase \(\pi\),

\[
R_p\to -1,
\]

which is Dirichlet-like. For the dual flow/displacement variable, the reflected amplitude has the opposite sign,

\[
R_q=-R_p\to +1,
\]

which is Neumann-like:

\[
\frac{q_w(L)}{q(L)}\to0.
\]

So the organ-pipe patch is viable, but the Volume 2 wording must be precise:

> The open throat dynamically produces the D/N ladder for the **flow/displacement support coordinate**, not for a generic pressure-like scalar field.

---

## 1. Junction model

Use the variable-area scalar wave action

\[
S_\psi^{(2)}
=
\frac12\int dt\,dx\,
A(x)\left[
\frac{1}{c_s^2}\psi_t^2-\psi_x^2
\right].
\]

The Euler-Lagrange equation is

\[
\frac{A}{c_s^2}\psi_{tt}
-
\partial_x(A\psi_x)
=0.
\]

In divided form,

\[
\frac{1}{c_s^2}\psi_{tt}
-
\frac1A\partial_x(A\psi_x)
=0.
\]

For the 1D throat tube,

\[
A(x)=A_t=\text{constant},
\]

so

\[
\psi_{tt}=c_s^2\psi_{xx}.
\]

For the open 4D spatial bulk, the radial area scales as

\[
A_{\rm bulk}(r)=S_3r^3,
\qquad
S_3=2\pi^2,
\]

so the radial wave equation is

\[
\frac{1}{c_s^2}\Psi_{tt}
-
\frac{1}{S_3r^3}\partial_r(S_3r^3\Psi_r)
=0,
\]

or

\[
\frac{1}{c_s^2}\Psi_{tt}
-
\Psi_{rr}
-
\frac3r\Psi_r
=0.
\]

The conservative pressure/potential-like scalar matching conditions are continuity of the field and continuity of the weighted normal derivative:

\[
\psi_{\rm tube}(L)=\Psi_{\rm bulk}(r_e),
\]

\[
A_t\partial_w\psi_{\rm tube}(L)
=
A_{\rm bulk}(r_e)\partial_r\Psi_{\rm bulk}(r_e).
\]

For acoustic variables it is often clearer to use effort/flow impedance:

\[
p=Z_LQ
\]

at the load.

---

## 2. Reflection and transmission coefficients

Let

\[
Z_t
\]

be the tube characteristic impedance and

\[
Z_L
\]

the effective load impedance of the open bulk. Define

\[
\epsilon=\frac{Z_L}{Z_t}.
\]

For a pressure-like amplitude,

\[
R_p
=
\frac{Z_L-Z_t}{Z_L+Z_t}
=
\frac{\epsilon-1}{\epsilon+1},
\]

\[
T_p=1+R_p
=
\frac{2\epsilon}{1+\epsilon}.
\]

For the dual flow/displacement amplitude,

\[
R_q=-R_p
=
\frac{1-\epsilon}{1+\epsilon},
\]

\[
T_q=1+R_q
=
\frac{2}{1+\epsilon}.
\]

The energy reflection and transmission coefficients are

\[
\mathcal R_E
=
\left(\frac{\epsilon-1}{\epsilon+1}\right)^2,
\]

\[
\mathcal T_E
=
\frac{4\epsilon}{(1+\epsilon)^2}.
\]

A sudden expansion has

\[
\epsilon\ll1.
\]

Then

\[
R_p=-1+2\epsilon-2\epsilon^2+O(\epsilon^3),
\]

\[
R_q=1-2\epsilon+2\epsilon^2+O(\epsilon^3),
\]

\[
\mathcal T_E=4\epsilon+O(\epsilon^2).
\]

So the AC wave is strongly reflected, with only \(O(\epsilon)\) energy leakage per encounter.

If the effective open area is written as

\[
\chi=\frac{A_{\rm eff}}{A_t}\gg1,
\]

then

\[
\epsilon\simeq\chi^{-1},
\]

and

\[
R_p=\frac{1-\chi}{1+\chi}\to -1,
\]

\[
R_q=\frac{\chi-1}{\chi+1}\to +1.
\]

---

## 3. D/N validation

For a tube mode near the exit,

\[
f(w)
=
e^{ik(w-L)}
+
R e^{-ik(w-L)}.
\]

At \(w=L\),

\[
\frac{f_w}{f}
=
ik\frac{1-R}{1+R}.
\]

### Pressure/potential-like scalar

Using

\[
R_p=\frac{\epsilon-1}{\epsilon+1},
\]

the derivative ratio becomes

\[
\frac{p_w}{p}
=
\frac{ik}{\epsilon}.
\]

As \(\epsilon\to0\),

\[
p(L)=1+R_p\to0.
\]

So the open end is Dirichlet-like for \(p\):

\[
p(L)\simeq0.
\]

If the mouth is also Dirichlet, this gives a D/D ladder,

\[
k_j=\frac{(j+1)\pi}{L}.
\]

This **does not** reproduce the desired half-shifted D/N branch.

### Flow/displacement-like support variable

For the dual support variable,

\[
R_q=\frac{1-\epsilon}{1+\epsilon},
\]

and

\[
\frac{q_w}{q}
=
ik\epsilon.
\]

As \(\epsilon\to0\),

\[
\frac{q_w(L)}{q(L)}\to0,
\]

so

\[
q_w(L)\simeq0.
\]

With the mouth condition

\[
q(0)=0,
\]

the support ladder is

\[
q_j(w)=A\sin(k_j w),
\]

\[
q_j'(L)=0,
\]

which gives

\[
\cos(k_jL)=0,
\]

and therefore

\[
\boxed{
k_j=\frac{\pi}{L}\left(j+\frac12\right).
}
\]

So the D/N support ladder is validated by the open-junction mechanism only after the support coordinate is identified as the free-end flow/displacement variable, or an equivalent conjugate variable.

---

## 4. DC leakage and the zero-mode

For a steady throughput \(\Phi\), continuity in the 4D bulk gives

\[
A_{\rm bulk}(r)J_r(r)=\Phi.
\]

Since

\[
A_{\rm bulk}(r)=S_3r^3=2\pi^2r^3,
\]

the steady flux density is

\[
J_r(r)
=
\frac{\Phi}{2\pi^2r^3}.
\]

The divergence check is exact:

\[
\frac1{A_{\rm bulk}}\partial_r(A_{\rm bulk}J_r)=0.
\]

At a finite tube exit with effective 3D-ball cross-sectional area

\[
A_{\rm exit}=\frac43\pi a_{\rm exit}^3,
\]

the tube flux density is finite:

\[
J_{\rm tube}
=
\frac{\Phi}{A_{\rm exit}}
=
\frac{3\Phi}{4\pi a_{\rm exit}^3}.
\]

But for a hard cap,

\[
A_{\rm exit}\to0,
\]

and

\[
\frac{\Phi}{A_{\rm exit}}\to\infty.
\]

Therefore a capped geometry is incompatible with a finite nonzero DC superfluid throughput.

This is the main reason the old cap language must be retired.

---

## 5. Pass/fail table

| Gate | Status | Reason |
|---|---:|---|
| Hard cap with nonzero DC flux | **Fail** | \(\Phi/A_{\rm exit}\to\infty\) as \(A_{\rm exit}\to0\). |
| Open finite-radius exit for DC flow | **Pass** | \(J_{\rm tube}\) finite and \(J_r=\Phi/(2\pi^2r^3)\) satisfies radial continuity. |
| Strong AC reflection from sudden expansion | **Pass** | \(\mathcal T_E=4\epsilon/(1+\epsilon)^2=O(\epsilon)\) for \(\epsilon\ll1\). |
| Neumann reflection for generic scalar \(\psi\) | **Fail** | Pressure/potential-like scalar has \(R_p\to-1\), Dirichlet-like. |
| Neumann reflection for flow/displacement support variable | **Pass** | Dual variable has \(R_q\to+1\), \(q_w/q=ik\epsilon\to0\). |
| D/N ladder | **Conditional pass** | Valid if the finite-throat support field is the free-end flow/displacement coordinate or an impedance-transformed equivalent. |

---

## 6. Required wording patch

Replace the old cap language with:

> The finite throat is an open conduit. At \(w=L\), the exit radius is finite, \(R(L)>0\), and the throat opens into unconfined 4D spatial bulk. The half-shifted D/N support ladder is not produced by a hard geometric cap. It is produced by the low-impedance open-junction reflection of the free-end support coordinate. The conjugate pressure/potential variable sees the usual open-end Dirichlet-like phase, while the support displacement/flow coordinate sees a Neumann-like condition. The baseline/DC superfluid current exits through the same open junction and is tracked by the existing leakage/current bookkeeping.

This preserves both requirements:

\[
\text{AC support trapping}
\]

and

\[
\text{DC mass/current leakage}.
\]

---

## 7. Required geometrical/tapering condition if using a pressure-like scalar

If the program insists that the D/N field \(\psi\) is pressure-like or potential-like, then a sudden expansion does **not** generate the required Neumann condition.

To recover Neumann for that variable, one needs an effective high-impedance input:

\[
\frac{Z_L}{Z_t}\gg1.
\]

That requires one of the following:

1. a narrow high-impedance neck,
2. a quarter-wave impedance transformer,
3. a side-branch/stub resonator whose input impedance is high at the support frequencies,
4. or a variable redefinition in which the D/N support coordinate is the dual flow/displacement variable.

The simplest patch is option 4.

---

## 8. Volume 2 consequence

The rewritten V2-04 is not a cap-regularity problem anymore. It is an open-junction impedance theorem.

The next downstream stages should use:

\[
R(L)>0,
\]

not

\[
R(L)=0,
\]

and should distinguish:

- DC throughput and leakage, handled by continuity and \(S_{\rm leak}\),
- AC support modes, handled by impedance reflection,
- pressure/potential variables, which see Dirichlet-like open-end reflection,
- flow/displacement support variables, which see Neumann-like open-end reflection.

This is the clean way to preserve the D/N ladder without violating mass continuity.

stage_v2_06_07_poisson_newtonian_derivation.md
# Stage V2-06/V2-07 — Projected Continuity, Poisson Hook, and Newtonian Universality

## Purpose

This stage audits the path

\[
\text{bulk continuity}
\;\longrightarrow\;
\text{projected brane continuity}
\;\longrightarrow\;
\text{longitudinal identity}
\;\longrightarrow\;
\text{controlled Poisson hook}
\;\longrightarrow\;
\text{Newtonian source law}.
\]

The goal is not just to rederive the inverse-square profile.  The stronger question is whether the current moving-throat stack already proves a universal Newtonian gravity law, or whether additional branch data are still required.

## Source-status anchor

The 4D / PN stack treats projected continuity and the longitudinal identity as exact, while the Poisson equation appears only after a controlled near-zone reduction.  The 1PN and 2PN assemblies then import the Newtonian point-particle block inside a closure hierarchy.  This stage keeps that status firewall.

## 1. Exact projected continuity

Start from exact bulk continuity:

\[
\partial_t \rho+\partial_a j^a+\partial_w j^w=0,
\qquad a\in\{x,y,z\}.
\]

Define the projected brane density and current by a fixed measurement kernel \(W(w)\):

\[
\rho_{\rm br}(\mathbf x,t)
=
\int W(w)\rho(\mathbf x,w,t)\,dw,
\]

\[
J^a_{\rm br}(\mathbf x,t)
=
\int W(w)j^a(\mathbf x,w,t)\,dw.
\]

Then

\[
\partial_t\rho_{\rm br}
+
\nabla_3\cdot\mathbf J_{\rm br}
=
-\int W(w)\partial_w j^w\,dw.
\]

Using

\[
-W\partial_w j^w
=
-\partial_w(Wj^w)+W'j^w,
\]

the leakage term is

\[
\boxed{
S_{\rm leak}
=
-\left[Wj^w\right]_{w_1}^{w_2}
+
\int_{w_1}^{w_2}W'(w)j^w\,dw.
}
\]

On an infinite projection window with fast decay of \(Wj^w\),

\[
\boxed{
S_{\rm leak}
=
\int_{-\infty}^{+\infty}W'(w)j^w\,dw.
}
\]

The script verifies the product-rule identity exactly.

## 2. Exact longitudinal identity

Define the brane velocity by the ratio of projected quantities:

\[
\mathbf v_{\rm br}
=
\frac{\mathbf J_{\rm br}}{\rho_{\rm br}}.
\]

Use the Helmholtz decomposition

\[
\mathbf v_{\rm br}
=
\nabla_3\phi_{\rm br}+\mathbf v_T,
\qquad
\nabla_3\cdot\mathbf v_T=0.
\]

Projected continuity is

\[
\partial_t\rho_{\rm br}
+
\nabla_3\cdot(\rho_{\rm br}\mathbf v_{\rm br})
=
S_{\rm leak}.
\]

Expanding the divergence gives

\[
\nabla_3\cdot(\rho_{\rm br}\mathbf v_{\rm br})
=
(\nabla_3\rho_{\rm br})\cdot(\nabla_3\phi_{\rm br}+\mathbf v_T)
+
\rho_{\rm br}\nabla_3^2\phi_{\rm br}
+
\rho_{\rm br}\nabla_3\cdot\mathbf v_T.
\]

With \(\nabla_3\cdot\mathbf v_T=0\), the exact identity is

\[
\boxed{
\rho_{\rm br}\nabla_3^2\phi_{\rm br}
=
S_{\rm leak}
-
\partial_t\rho_{\rm br}
-
(\nabla_3\rho_{\rm br})\cdot(\nabla_3\phi_{\rm br}+\mathbf v_T).
}
\]

The script verifies the divergence expansion and the identity modulo continuity and \(\nabla_3\cdot\mathbf v_T=0\).

## 3. Controlled Poisson hook

The identity becomes a Poisson equation only in a controlled regime:

\[
\rho_{\rm br}\approx\rho_0,
\qquad
\partial_t\rho_{\rm br}\ \text{subleading},
\qquad
\nabla_3\rho_{\rm br}\ \text{subleading},
\qquad
\mathbf v_T\ \text{subleading or perturbative}.
\]

Then

\[
\boxed{
\nabla_3^2\phi_{\rm br}
\simeq
\frac{S_{\rm eff}}{\rho_0}.
}
\]

For a localized effective source

\[
S_{\rm eff}(\mathbf x)=S_0\delta^{(3)}(\mathbf x),
\]

the Green-function solution is

\[
\boxed{
\phi_{\rm br}(\mathbf x)
=
-\frac{S_0}{4\pi\rho_0 r}.
}
\]

The script verifies

\[
\nabla^2\left(-\frac{S_0}{4\pi\rho_0 r}\right)=0
\quad (r>0),
\]

and the sphere-flux normalization

\[
\oint \nabla\phi_{\rm br}\cdot d\mathbf S
=
\frac{S_0}{\rho_0},
\]

which gives the distributional identity

\[
\nabla^2\phi_{\rm br}
=
\frac{S_0}{\rho_0}\delta^{(3)}(\mathbf x).
\]

## 4. Newtonian normalization gate

Define the measured Newtonian potential by

\[
\Phi_N=\lambda_\Phi\phi_{\rm br}.
\]

For a compact defect of inertial mass \(m\), suppose the effective source amplitude is

\[
S_0=\eta m.
\]

Then

\[
\nabla^2\Phi_N
=
\frac{\lambda_\Phi\eta}{\rho_0}m\delta^{(3)}(\mathbf x).
\]

To match

\[
\nabla^2\Phi_N=4\pi Gm\delta^{(3)}(\mathbf x),
\]

one needs

\[
\boxed{
\eta=\frac{4\pi G\rho_0}{\lambda_\Phi}.
}
\]

Equivalently,

\[
\boxed{
\frac{\lambda_\Phi\eta}{4\pi\rho_0}=G.
}
\]

This is the first universality gate: the integrated source strength of every compact defect must be proportional to its inertial mass with the same coefficient \(\eta\).

## 5. Pair-counting gate

For two bodies \(A,B\), let their source coefficients be \(\eta_A,\eta_B\).  The potential at \(A\) from \(B\) is

\[
\Phi_B(A)
=
-\frac{\lambda_\Phi\eta_Bm_B}{4\pi\rho_0 r_{AB}}.
\]

The half-counted pair interaction is

\[
L_{\rm pair}
=
-\frac12\left[m_A\Phi_B(A)+m_B\Phi_A(B)\right].
\]

Therefore

\[
L_{\rm pair}
=
\frac{\lambda_\Phi(\eta_A+\eta_B)}{8\pi\rho_0}
\frac{m_Am_B}{r_{AB}}.
\]

So the effective pair coefficient is

\[
G_{AB}
=
\frac{\lambda_\Phi(\eta_A+\eta_B)}{8\pi\rho_0}.
\]

To get one universal \(G\) for arbitrary pairs, it is sufficient and necessary in this two-species test that

\[
\boxed{
\eta_A=\eta_B=\frac{4\pi G\rho_0}{\lambda_\Phi}.
}
\]

Then

\[
\boxed{
L_0
=
\frac12m_Av_A^2+\frac12m_Bv_B^2+\frac{Gm_Am_B}{r_{AB}}.
}
\]

The script verifies this algebra.

## 6. Response-mass / equivalence-principle gate

Let

\[
m_{g,A}=\kappa_A m_A,
\qquad
m_{g,B}=\kappa_B m_B.
\]

For a particle action

\[
L_A=\frac12m_Av_A^2-m_{g,A}\Phi_N,
\]

the acceleration is

\[
\ddot{\mathbf X}_A
=
-\kappa_A\nabla\Phi_N.
\]

Two species in the same field have acceleration difference

\[
\Delta \mathbf a
=
-(\kappa_A-\kappa_B)\nabla\Phi_N.
\]

Thus the model needs

\[
\boxed{
\kappa_A=\kappa_B=\kappa_\rho.
}
\]

The Newtonian ledger then chooses the normalization

\[
\boxed{\kappa_\rho=1.}
\]

This is not a consequence of projected continuity alone.  It is a compact-defect source/response theorem that must be proved from the actual branch, or treated as a declared closure datum.

## 7. Same-\(G\) ledger

The outgoing quadrupole target from the 2.5PN / 4PN bridge is

\[
\widehat m_0^2P_0
=
\frac{54Gc_s^5}{5a^5c^5}.
\]

If this sector used a different constant \(G_Q\) from the Newtonian source law \(G_N\), the mismatch would be

\[
\frac{54c_s^5}{5a^5c^5}(G_Q-G_N).
\]

So Volume 2 should carry a same-\(G\) constraint:

\[
\boxed{G_Q=G_N.}
\]

The script verifies the symbolic equality after imposing \(G_Q=G_N\).  This is a bookkeeping unification gate, not an independent derivation from the projection identity.

## 8. Verdict

The exact mathematical chain passes:

```text
EXACT_PROJECTED_CONTINUITY: PASS
EXACT_LONGITUDINAL_IDENTITY: PASS
CONTROLLED_POISSON_HOOK: PASS under stated quasi-static assumptions
INVERSE_SQUARE_GREEN_FUNCTION: PASS
```

The Newtonian source law receives a conditional pass:

```text
NEWTONIAN_SOURCE_LAW: CONDITIONAL PASS
```

Required branch-level facts:

\[
S_{\rm eff}^{(A)}
=
\eta m_A\delta^{(3)}(\mathbf x-\mathbf X_A)
+O(a_A^2\nabla^2\delta),
\]

with one universal \(\eta\), and

\[
m_{g,A}=m_A
\]

after one shared normalization.

The strict result is:

\[
\boxed{
\text{Projected continuity gives a Poisson hook, not yet a universal Newtonian source theorem.}
}
\]

The next repair theorem should compute the compact-defect source amplitude directly from the moving-throat branch and prove

\[
\boxed{
\frac{1}{m_A}\int S_{\rm eff}^{(A)}\,d^3x
=
\frac{4\pi G\rho_0}{\lambda_\Phi}
}
\]

for every admissible compact defect species in the intended point-particle class.

## 9. Script output summary

The accompanying SymPy audit reports:

```text
symbolic_checks: 11
symbolic_passes: 11
```

and isolates the remaining non-symbolic condition as branch-level source universality rather than algebraic inconsistency.

stage_v2_08_bdg_wall_schur_stability_derivation.md
# Stage V2-08 — BdG–Wall Schur Complement and Stability/Softening Audit

## Purpose

This stage audits the first matter-support closure in the Volume 2 roadmap.

The target is the reduced **stable BdG normal-mode** block coupled to a moving wall/worldtube mode. The question is not whether the full nonlinear moving-throat PDE is solved. The narrower question is:

> If the linearized GNLS/BdG sector is reduced to positive-energy stable normal coordinates, does integrating those modes out produce a consistent conservative wall kernel, and what stability constraints does that impose?

The answer is a conditional **yes**:

\[
\boxed{
\text{Stable positive-norm BdG support modes give a valid Schur-complement softening of the wall.}
}
\]

But the audit also finds three important restrictions:

1. the closure requires positive-Krein/positive-energy support modes;
2. near-softening must preserve \(D_0>0\);
3. the low-frequency BdG moments are constrained Stieltjes moments, not freely tunable coefficients.

---

## 1. One wall mode coupled to one stable support mode

Start with the reduced quadratic Lagrangian

\[
L
=
\frac12 M\dot q^2
-\frac12 Kq^2
+\frac12 \dot X^2
-\frac12 \varpi^2X^2
+gqX.
\]

The Euler–Lagrange equations are

\[
M\ddot q+Kq-gX=0,
\]

\[
\ddot X+\varpi^2X-gq=0.
\]

Using the frequency convention \(e^{-i\omega t}\),

\[
(K-M\omega^2)q-gX=0,
\]

\[
(\varpi^2-\omega^2)X-gq=0.
\]

Solving the second equation gives

\[
X=\frac{g}{\varpi^2-\omega^2}q.
\]

Substitution gives the exact effective wall operator

\[
\boxed{
D_{\rm eff}(\omega)
=
K-M\omega^2
-
\frac{g^2}{\varpi^2-\omega^2}.
}
\]

This is the one-mode Schur complement.

---

## 2. Low-frequency expansion

Expanding for \(|\omega|<\varpi\),

\[
\frac{1}{\varpi^2-\omega^2}
=
\frac1{\varpi^2}
+
\frac{\omega^2}{\varpi^4}
+
\frac{\omega^4}{\varpi^6}
+
O(\omega^6).
\]

Therefore

\[
D_{\rm eff}(\omega)
=
\left(K-\frac{g^2}{\varpi^2}\right)
-
\omega^2
\left(M+\frac{g^2}{\varpi^4}\right)
-
\omega^4
\left(\frac{g^2}{\varpi^6}\right)
+
O(\omega^6).
\]

Define

\[
K_{\rm eff}=K-\frac{g^2}{\varpi^2},
\]

\[
M_{\rm eff}=M+\frac{g^2}{\varpi^4},
\]

\[
N_{\rm eff}=\frac{g^2}{\varpi^6}.
\]

Then

\[
\boxed{
D_{\rm eff}(\omega)
=
K_{\rm eff}
-
M_{\rm eff}\omega^2
-
N_{\rm eff}\omega^4
+
O(\omega^6).
}
\]

The support mode therefore:

\[
\boxed{
\text{lowers static stiffness, raises inertia, and adds a positive higher even moment.}
}
\]

This is exactly the desired conservative support behavior.

---

## 3. Static positivity and softening bound

The static potential is

\[
V(q,X)
=
\frac12Kq^2
+
\frac12\varpi^2X^2
-
gqX.
\]

Completing the square,

\[
V(q,X)
=
\frac12\varpi^2
\left(
X-\frac{g}{\varpi^2}q
\right)^2
+
\frac12
\left(
K-\frac{g^2}{\varpi^2}
\right)q^2.
\]

So the static stability gate is

\[
\boxed{
K_{\rm eff}=K-\frac{g^2}{\varpi^2}>0.
}
\]

Equivalently,

\[
\boxed{
K\varpi^2-g^2>0.
}
\]

This is the exact Schur-complement positivity condition for the static Hessian

\[
H_{\rm stat}
=
\begin{pmatrix}
K & -g\\
-g & \varpi^2
\end{pmatrix}.
\]

The determinant is

\[
\det H_{\rm stat}=K\varpi^2-g^2.
\]

So conservative support softening is allowed only up to the stability boundary. A branch that uses softening to raise downstream normalization must keep

\[
D_0=K_{\rm eff}>0.
\]

For many support modes in a one-wall lane,

\[
D_0
=
K-\sum_\alpha\frac{g_\alpha^2}{\varpi_\alpha^2}.
\]

It is useful to define

\[
\epsilon_B
=
\frac{1}{K}
\sum_\alpha
\frac{g_\alpha^2}{\varpi_\alpha^2}.
\]

Then

\[
\boxed{
D_0=K(1-\epsilon_B),
\qquad
\epsilon_B<1.
}
\]

---

## 4. Dynamic pole stability

The exact dispersion relation is

\[
(K-M\omega^2)(\varpi^2-\omega^2)-g^2=0.
\]

Let

\[
\Omega_\eta^2=\frac{K}{M}.
\]

Then the two squared frequencies are

\[
\boxed{
\omega_\pm^2
=
\frac{
\Omega_\eta^2+\varpi^2
\pm
\sqrt{(\Omega_\eta^2-\varpi^2)^2+4g^2/M}
}{2}.
}
\]

Both poles are real. They are both positive exactly when

\[
\Omega_\eta^2\varpi^2-\frac{g^2}{M}>0,
\]

which is again

\[
\boxed{
K\varpi^2-g^2>0.
}
\]

So the same Schur-complement condition protects both static and dynamic stability.

If the support mode lies above the wall mode,

\[
\varpi^2=\Omega_\eta^2+\Delta,
\qquad
\Delta>0,
\]

and

\[
h=\frac{g^2}{M},
\]

then the wall-like pole shifts as

\[
\boxed{
\omega_-^2
=
\Omega_\eta^2
-
\frac{h}{\Delta}
+
\frac{h^2}{\Delta^3}
+
O(h^3).
}
\]

So an above-wall support mode softens the wall-like pole at weak coupling.

---

## 5. Matrix Schur complement

For wall coordinates \(Q^A\) and stable support coordinates \(X_\alpha\),

\[
L
=
\frac12\dot Q^TM\dot Q
-
\frac12 Q^TKQ
+
\frac12\dot X^T\dot X
-
\frac12 X^T\Omega^2X
+
Q^TCX.
\]

Eliminating \(X\) gives

\[
\boxed{
D_{\rm eff}(\omega)
=
K-\omega^2M
-
C(\Omega^2-\omega^2I)^{-1}C^T.
}
\]

The low-frequency expansion is

\[
\boxed{
D_{\rm eff}(\omega)
=
K_{\rm eff}
-\omega^2M_{\rm eff}
-\omega^4N_{\rm eff}
+O(\omega^6),
}
\]

with

\[
\boxed{
K_{\rm eff}=K-C\Omega^{-2}C^T,
}
\]

\[
\boxed{
M_{\rm eff}=M+C\Omega^{-4}C^T,
}
\]

\[
\boxed{
N_{\rm eff}=C\Omega^{-6}C^T.
}
\]

The positivity gate is the block-Hessian condition

\[
\begin{pmatrix}
K & -C\\
-C^T & \Omega^2
\end{pmatrix}>0.
\]

Since \(\Omega^2>0\), this is equivalent to

\[
\boxed{
K-C\Omega^{-2}C^T>0.
}
\]

That is the general matrix stability gate.

---

## 6. Negative-Krein / ghost warning

The reduced closure only works for positive-energy stable support modes.

If a support coordinate has negative kinetic/potential sign,

\[
L_X
=
-\frac12\dot X^2
+
\frac12\varpi^2X^2
+
gqX,
\]

then the static Hessian is effectively

\[
H_{\rm ghost}
=
\begin{pmatrix}
K & -g\\
-g & -\varpi^2
\end{pmatrix}.
\]

Its determinant is

\[
\det H_{\rm ghost}
=
-K\varpi^2-g^2<0.
\]

Thus the static energy is indefinite for \(K>0,\varpi^2>0\). The SymPy audit marks this as a failure mode.

So the Volume 2 statement must be:

\[
\boxed{
\text{V2-08 passes only for the positive-Krein stable normal-mode sector.}
}
\]

A later full BdG/Krein audit should verify that the modes retained in the reduced branch really have positive energy.

---

## 7. Moment constraints: \(B_0,B_2,B_4\) are not free

For scalar lane data, define

\[
B_0=\sum_\alpha\frac{g_\alpha^2}{\lambda_\alpha},
\]

\[
B_2=\sum_\alpha\frac{g_\alpha^2}{\lambda_\alpha^2},
\]

\[
B_4=\sum_\alpha\frac{g_\alpha^2}{\lambda_\alpha^3},
\]

where

\[
\lambda_\alpha=\varpi_\alpha^2>0.
\]

For two modes with positive weights \(w_i=g_i^2\),

\[
B_0B_4-B_2^2
=
\boxed{
\frac{
w_1w_2(\lambda_1-\lambda_2)^2
}{
\lambda_1^3\lambda_2^3
}
\ge0.
}
\]

This is a Stieltjes/Hankel moment constraint. It means the stable BdG front end cannot arbitrarily fit the conservative even coefficients. Any proposed branch data must satisfy these positivity/moment relations.

In particular, equality occurs only when the support effectively has one spectral scale.

---

## 8. Grouped real \(P_2\) isotropy

For each grouped real lane \(A\in\{20,21,22\}\),

\[
D_A(\omega)
=
D_{A0}+D_{A2}\omega^2+D_{A4}\omega^4+\cdots.
\]

The normalized response coefficient is

\[
u_2^{(A)}
=
-\frac{D_{A2}}{D_{A0}}.
\]

Grouped trace and anisotropy variables are

\[
\bar u_2
=
\frac{
u_2^{(20)}
+2u_2^{(21)}
+2u_2^{(22)}
}{5},
\]

\[
a_2
=
\frac{
2u_2^{(20)}
-u_2^{(21)}
-u_2^{(22)}
}{10},
\]

\[
b_2
=
\frac{
u_2^{(21)}
-u_2^{(22)}
}{2}.
\]

If all grouped lanes share identical stable support data,

\[
D_{20,n}=D_{21,n}=D_{22,n},
\]

then

\[
\boxed{
a_2=b_2=0.
}
\]

So the stable BdG Schur complement preserves grouped \(P_2\) isotropy if the support spectrum and overlap data are isotropic.

---

## 9. Verdict

The V2-08 audit passes as a **controlled reduced-sector gate**.

\[
\boxed{
\text{The Schur complement is algebraically exact and physically stable under positive-mode conditions.}
}
\]

The gate does **not** prove the full moving-throat branch. It establishes the conditions that any branch realization must satisfy:

\[
\boxed{
\Omega^2>0,\qquad
M>0,\qquad
K-C\Omega^{-2}C^T>0.
}
\]

It also establishes that downstream coefficient matching cannot treat the BdG moments as free:

\[
\boxed{
B_0B_4-B_2^2\ge0
}
\]

in the scalar two-mode prototype, with the corresponding positive-moment constraints in the multimode case.

## Carry-forward patch

Use this wording in Volume 2:

> The BdG-wall reduced closure is valid after projecting the linearized GNLS/BdG system onto positive-energy stable normal coordinates. Integrating those modes out gives an exact Schur-complement wall kernel. Stable support modes lower static stiffness, raise effective inertia, and generate positive higher even moments. The branch must satisfy \(K-C\Omega^{-2}C^T>0\); near-softened normalization enhancement is allowed only below that stability boundary. Negative-Krein modes are outside this closure and require a separate full BdG signature audit.

stage_v2_09_maxwell_mixed_kernel_derivation.md
# Stage V2-09 — Maxwell/mixed-sector kernel audit

## Status

**Verdict:** conditional pass inside the declared one-lane reduced Maxwell/mixed closure.

The audit checks the smallest reduced block that can carry an honest passive/outgoing channel:

\[
(Q,\ U,\ W),
\]

where \(Q\) is a wall/worldtube amplitude, \(U\) is a brane-like localized Maxwell coordinate, and \(W\) is the mixed-sector coordinate representing the \(A_w/F_{\mu w}/J^w\)-active block.

This stage does **not** prove that the full nonlinear moving-throat PDE realizes the needed branch. It proves that the reduced Maxwell/mixed mechanism has the right algebraic structure, stability gates, outgoing-transfer sign, and scalar-rescue condition.

---

## 1. Mixed-sector gauge invariants

Use the gauge convention

\[
A_0\mapsto A_0-\partial_t\chi,\qquad
A_a\mapsto A_a+\partial_a\chi,\qquad
A_w\mapsto A_w+\partial_w\chi.
\]

The mixed fields are

\[
E_w=-\partial_t A_w-\partial_w A_0,
\]

\[
C_a=\partial_a A_w-\partial_w A_a.
\]

The script verifies exactly that

\[
\delta_\chi E_w=0,
\qquad
\delta_\chi C_a=0.
\]

So the mixed variables used in this stage are gauge-invariant reduced observables, not gauge artifacts.

---

## 2. Reduced one-lane Maxwell/mixed model

The reduced quadratic Lagrangian is

\[
L
=
\frac12M\dot Q^2-\frac12KQ^2
+
\frac12\dot U^2-\frac12\Omega_U^2U^2
+
\frac12\dot W^2-\frac12\Omega_W^2W^2
+
R\,UW
+
g_UQU
+
g_WQW.
\]

Here:

- \(Q\) is the wall/worldtube amplitude,
- \(U\) is the localized brane-like Maxwell coordinate,
- \(W\) is the mixed \(A_w/F_{\mu w}/J^w\)-active coordinate,
- \(R\) mixes the two internal gauge-sector coordinates,
- \(g_U,g_W\) couple the wall to the internal gauge block.

With the convention \(e^{-i\omega t}\), define

\[
A(\omega)=\Omega_U^2-\omega^2,
\qquad
B(\omega)=\Omega_W^2-\omega^2,
\]

\[
\Delta(\omega)=A(\omega)B(\omega)-R^2.
\]

The internal equations are

\[
\begin{pmatrix}
A(\omega)&-R\\
-R&B(\omega)
\end{pmatrix}
\begin{pmatrix}
U\\W
\end{pmatrix}
=
\begin{pmatrix}
g_UQ\\g_WQ
\end{pmatrix}.
\]

The inverse internal block is

\[
\frac{1}{\Delta(\omega)}
\begin{pmatrix}
B(\omega)&R\\
R&A(\omega)
\end{pmatrix}.
\]

Therefore the conservative self-energy inherited by the wall is

\[
\boxed{
\Sigma_{\rm cons}(\omega)
=
\frac{
g_U^2B(\omega)+2g_Ug_WR+g_W^2A(\omega)
}{
\Delta(\omega)
}.
}
\]

The wall operator is

\[
\boxed{
D_{\rm cons}(\omega)=K-M\omega^2-\Sigma_{\rm cons}(\omega).
}
\]

The SymPy audit verifies the Schur-complement numerator and denominator exactly.

---

## 3. Conservative low-frequency expansion

Let

\[
\Delta_0=\Omega_U^2\Omega_W^2-R^2,
\]

\[
S_2=\Omega_U^2+\Omega_W^2,
\]

\[
Q_0=g_U^2\Omega_W^2+2g_Ug_WR+g_W^2\Omega_U^2,
\]

\[
G_2=g_U^2+g_W^2.
\]

Then

\[
\Sigma_{\rm cons}(\omega)
=
z_0+z_2\omega^2+z_4\omega^4+O(\omega^6),
\]

with

\[
\boxed{
z_0=\frac{Q_0}{\Delta_0},
}
\]

\[
\boxed{
z_2=\frac{Q_0S_2-G_2\Delta_0}{\Delta_0^2},
}
\]

\[
\boxed{
z_4=
\frac{
Q_0(S_2^2-\Delta_0)-S_2G_2\Delta_0
}{
\Delta_0^3
}.
}
\]

The script verifies this by multiplying the truncated series by the denominator and checking equality through \(O(\omega^4)\).

---

## 4. Stability and positivity gates

The static internal Hessian is

\[
H_{\rm int}
=
\begin{pmatrix}
\Omega_U^2&-R\\
-R&\Omega_W^2
\end{pmatrix}.
\]

The internal block is positive if

\[
\boxed{
\Omega_U^2>0,\qquad
\Omega_W^2>0,\qquad
\Delta_0>0.
}
\]

The full static Hessian for \((Q,U,W)\) is

\[
H_{\rm full}
=
\begin{pmatrix}
K&-g_U&-g_W\\
-g_U&\Omega_U^2&-R\\
-g_W&-R&\Omega_W^2
\end{pmatrix}.
\]

The script verifies

\[
\det H_{\rm full}
=
K\Delta_0-Q_0
=
\Delta_0\left(K-\frac{Q_0}{\Delta_0}\right).
\]

Thus the full reduced one-lane stability gate is

\[
\boxed{
\Delta_0>0,
\qquad
K-\Sigma_{\rm cons}(0)>0.
}
\]

Equivalently,

\[
\boxed{
K\Delta_0-Q_0>0.
}
\]

This is the Maxwell/mixed analogue of the BdG softening gate: conservative support can soften the wall, but the branch must stay on the positive side of the Schur boundary.

---

## 5. Outgoing-port dressing

Attach a passive outgoing channel to the mixed coordinate by replacing

\[
B(\omega)\mapsto B(\omega)-\Pi_{\rm out}(\omega).
\]

Then

\[
\Sigma_{\rm full}(\omega)
=
\frac{
g_U^2[B(\omega)-\Pi_{\rm out}]+2g_Ug_WR+g_W^2A(\omega)
}{
A(\omega)[B(\omega)-\Pi_{\rm out}]-R^2
}.
\]

Expanding to first order in \(\Pi_{\rm out}\),

\[
\Sigma_{\rm full}(\omega)
=
\Sigma_{\rm cons}(\omega)
+
\Pi_{\rm out}(\omega)\,N(\omega)
+
O(\Pi_{\rm out}^2),
\]

where the script verifies the exact transfer factor

\[
\boxed{
N(\omega)
=
\frac{[A(\omega)g_W+Rg_U]^2}{\Delta(\omega)^2}.
}
\]

At zero frequency,

\[
\boxed{
N_0
=
\frac{(\Omega_U^2g_W+Rg_U)^2}{\Delta_0^2}\ge0.
}
\]

So the reduced Maxwell/mixed block transfers the outgoing branch with a nonnegative static weight.

There is one important caveat:

\[
\Omega_U^2g_W+Rg_U=0
\]

is a **dark-port condition**. If this happens, then

\[
N_0=0,
\]

and the leading quadrupole outgoing channel is not transferred to the wall.

---

## 6. Transfer-factor low-frequency expansion

Let

\[
P=\Omega_U^2g_W+Rg_U.
\]

Then

\[
N(\omega)=N_0+N_2\omega^2+N_4\omega^4+O(\omega^6),
\]

with

\[
\boxed{
N_0=\frac{P^2}{\Delta_0^2},
}
\]

\[
\boxed{
N_2=
\frac{2P(PS_2-\Delta_0g_W)}{\Delta_0^3},
}
\]

\[
\boxed{
N_4=
\frac{
\Delta_0^2g_W^2
-2\Delta_0P^2
-4\Delta_0PS_2g_W
+3P^2S_2^2
}{
\Delta_0^4
}.
}
\]

These are the one-lane transfer moments that later grouped-\(P_2\) stages must generalize.

---

## 7. Outgoing \(l=2\) wall coefficient

The compact outgoing \(l=2\) branch has the normalized odd coefficient

\[
\Gamma_2^{\rm port}
=
\frac{a^5}{27c_s^5}.
\]

So

\[
\Pi_2^{\rm out}(\omega)
=
+i\frac{a^5}{27c_s^5}\omega^5+\cdots.
\]

Because

\[
D(\omega)=K-M\omega^2-\Sigma(\omega),
\]

the wall operator inherits

\[
\boxed{
\delta D_2^{\rm odd}(\omega)
=
-i\,N_0\,\frac{a^5}{27c_s^5}\omega^5
+
O(\omega^7).
}
\]

The negative imaginary sign is the passive damping sign in the \(e^{-i\omega t}\) wall-operator convention. Written in normalized response/admittance language, this corresponds to the positive outgoing sign used in the 2.5PN package.

---

## 8. Scalar derivative-coupling compatibility

A naive scalar outgoing port can generate an \(i\omega\) term. The reduced Maxwell/mixed block avoids that only if the scalar outlet is derivative-coupled.

Take

\[
g_U=0,\qquad
g_W(\omega)=\eta\omega.
\]

Then

\[
N_0^{\rm scalar}(\omega)
=
\frac{A(\omega)^2\eta^2\omega^2}{\Delta(\omega)^2}
=
\frac{\Omega_U^4\eta^2}{\Delta_0^2}\omega^2+O(\omega^4).
\]

If

\[
\Pi_0^{\rm out}(\omega)=i\gamma_1\omega+\cdots,
\]

then the wall-level correction is

\[
\delta D_0^{\rm odd}
=
-i\gamma_1
\frac{\Omega_U^4\eta^2}{\Delta_0^2}
\omega^3
+
O(\omega^5).
\]

So the scalar outlet is demoted from \(i\omega\) to \(i\omega^3\).

This is not a full scalar no-go theorem. It is a precise reduced compatibility condition:

\[
\boxed{
\text{direct non-derivative scalar port coupling is dangerous;}
}
\]

\[
\boxed{
\text{derivative scalar port coupling delays the leading scalar odd term.}
}
\]

---

## 9. Pass/fail ledger

| Gate | Result |
|---|---:|
| Mixed \(E_w,C_a\) gauge invariance | Pass |
| Conservative Schur self-energy | Pass |
| \(z_0,z_2,z_4\) low-frequency expansion | Pass |
| Outgoing transfer factor | Pass |
| \(N_0,N_2,N_4\) transfer expansion | Pass |
| Static stability determinant | Pass |
| \(l=2\) odd wall coefficient | Pass |
| Scalar derivative-coupling demotion | Pass |
| Full PDE branch realization | Open |

---

## 10. Carry-forward theorem gates

This stage leaves three concrete conditions for the full moving-throat branch:

1. **Stability gate**

   \[
   \Delta_0>0,
   \qquad
   K-\frac{Q_0}{\Delta_0}>0.
   \]

2. **Quadrupole transfer gate**

   \[
   N_0=
   \frac{(\Omega_U^2g_W+Rg_U)^2}{\Delta_0^2}
   \neq 0.
   \]

3. **Universal normalization gate**

   After grouped-\(P_2\) lifting,

   \[
   \widehat m_0^2\frac{N_0}{D_0}
   =
   \frac{54Gc_s^5}{5a^5c^5}.
   \]

This stage proves that the reduced Maxwell/mixed mechanism is algebraically capable of carrying the passive/outgoing \(l=2\) branch. It does not prove that the actual moving-throat PDE lands on the target value.

stage_v2_10_hamiltonian_stability_derivation.md
# Stage V2-10 — Hamiltonian / Stability Audit

## 0. Purpose

This stage combines the conservative wall, BdG, and localized Maxwell/mixed
reduced sectors into one stability ledger.

Earlier stages verified the individual Schur complements:

- stable BdG support modes soften the static wall stiffness and raise the
  effective inertia;
- the localized Maxwell/mixed block gives a conservative self-energy and a
  nonnegative outgoing-transfer factor;
- the outgoing `l=2` port supplies the odd `i omega^5` fingerprint.

The present stage asks the next necessary question:

> When all conservative pieces are active at once, what exact positivity gates
> keep the reduced branch from becoming a ghost, tachyon, or unstable
> near-softened denominator?

The result is a compact pass condition:

\[
M>0,\qquad \varpi^2>0,\qquad
\Omega_U^2>0,\qquad \Omega_W^2>0,
\]
\[
\Delta=\Omega_U^2\Omega_W^2-R^2>0,
\]
\[
D_0
=
K-\frac{c_B^2}{\varpi^2}
-
\frac{
g_U^2\Omega_W^2+2g_Ug_WR+g_W^2\Omega_U^2
}{
\Omega_U^2\Omega_W^2-R^2
}
>0.
\]

If these hold, the conservative reduced Hamiltonian is positive definite.

---

## 1. Conservative one-lane bundle

Use one representative grouped lane.  The reduced variables are

\[
q,\qquad X,\qquad U,\qquad W,
\]

where

- \(q\) is the wall/worldtube amplitude,
- \(X\) is a positive-energy BdG support coordinate,
- \(U\) is a brane-like localized Maxwell coordinate,
- \(W\) is the mixed \(A_w/F_{\mu w}/J^w\)-active coordinate.

The reduced Lagrangian is

\[
L
=
\frac12M\dot q^2-\frac12Kq^2
+
\frac12\dot X^2-\frac12\varpi^2X^2
+
c_BqX
\]
\[
+
\frac12\dot U^2-\frac12\Omega_U^2U^2
+
\frac12\dot W^2-\frac12\Omega_W^2W^2
+
RUW+g_UqU+g_WqW.
\]

The conservative potential Hessian for \((q,X,U,W)\) is

\[
P=
\begin{pmatrix}
K&-c_B&-g_U&-g_W\\
-c_B&\varpi^2&0&0\\
-g_U&0&\Omega_U^2&-R\\
-g_W&0&-R&\Omega_W^2
\end{pmatrix}.
\]

The kinetic matrix is

\[
T=\operatorname{diag}(M,1,1,1).
\]

So positive Hamiltonian kinetic energy requires \(M>0\).

---

## 2. Exact Schur-complement stability gate

The internal conservative block is

\[
A_{\rm int}
=
\begin{pmatrix}
\varpi^2&0&0\\
0&\Omega_U^2&-R\\
0&-R&\Omega_W^2
\end{pmatrix}.
\]

The mixed Maxwell sub-block is positive iff

\[
\Delta=\Omega_U^2\Omega_W^2-R^2>0
\]

with \(\Omega_U^2,\Omega_W^2>0\).

The wall-coupling vector is

\[
g_{\rm int}=(c_B,g_U,g_W)^T.
\]

The Schur complement of the internal block is

\[
D_0
=
K-g_{\rm int}^TA_{\rm int}^{-1}g_{\rm int}.
\]

The script verifies exactly that

\[
D_0
=
K-\frac{c_B^2}{\varpi^2}
-
\frac{
g_U^2\Omega_W^2+2g_Ug_WR+g_W^2\Omega_U^2
}{\Delta}.
\]

It also verifies

\[
\det P
=
\varpi^2\Delta D_0.
\]

Therefore the positive-energy conservative gate is

\[
\boxed{
M>0,\quad \varpi^2>0,\quad \Delta>0,\quad D_0>0.
}
\]

This is the precise version of “near-softening is allowed only below the
instability boundary.”

---

## 3. Low-frequency conservative moments

Define

\[
Q=
g_U^2\Omega_W^2+2g_Ug_WR+g_W^2\Omega_U^2,
\qquad
S=\Omega_U^2+\Omega_W^2,
\]
\[
H=g_U^2+g_W^2.
\]

The BdG moments are

\[
B_0=\frac{c_B^2}{\varpi^2},
\qquad
B_2=\frac{c_B^2}{\varpi^4},
\qquad
B_4=\frac{c_B^2}{\varpi^6}.
\]

The Maxwell/mixed conservative moments are

\[
Z_0=\frac{Q}{\Delta},
\]
\[
Z_2=\frac{QS-H\Delta}{\Delta^2},
\]
\[
Z_4=
\frac{Q(S^2-\Delta)-SH\Delta}{\Delta^3}.
\]

The script verifies the matrix identities

\[
Z_0=g^TA^{-1}g,\qquad
Z_2=g^TA^{-2}g,\qquad
Z_4=g^TA^{-3}g,
\]

where

\[
A=
\begin{pmatrix}
\Omega_U^2&-R\\
-R&\Omega_W^2
\end{pmatrix},
\qquad
g=(g_U,g_W)^T.
\]

Because \(A>0\), all these \(Z\)-moments are nonnegative quadratic forms.

The conservative wall operator is

\[
D(\omega)
=
D_0+D_2\omega^2+D_4\omega^4+O(\omega^6),
\]

with

\[
D_0=K-B_0-Z_0,
\]
\[
D_2=-(M+B_2+Z_2),
\]
\[
D_4=-(B_4+Z_4).
\]

The script verifies that this matches the direct low-frequency expansion of the
full conservative self-energy through \(O(\omega^4)\).

---

## 4. Outgoing-port passivity gate

For the Maxwell/mixed outgoing transfer factor,

\[
N_0
=
\frac{(\Omega_U^2g_W+Rg_U)^2}{\Delta^2}.
\]

This is a perfect square over \(\Delta^2\), so

\[
N_0\ge0
\]

provided the internal block is nonsingular.

For the compact outgoing \(l=2\) port,

\[
\Gamma_2^{\rm port}
=
\frac{a^5}{27c_s^5}.
\]

The wall-level odd coefficient is therefore

\[
\gamma_{\rm wall}
=
N_0\Gamma_2^{\rm port}
=
N_0\frac{a^5}{27c_s^5}
\ge0.
\]

With the convention \(e^{-i\omega t}\), the wall operator contribution

\[
\delta D_{\rm odd}(\omega)
=
-i\gamma_{\rm wall}\omega^5
\]

corresponds in time domain to

\[
+\gamma_{\rm wall} q^{(5)}.
\]

The script verifies the identity

\[
\dot q\,q^{(5)}
=
\frac{d}{dt}
\left(
\dot q\,q^{(4)}-\ddot q\,q^{(3)}
\right)
+
(q^{(3)})^2.
\]

Thus, after the usual Schott-energy reshuffling, the outgoing \(l=2\) term
dissipates

\[
\gamma_{\rm wall}(q^{(3)})^2\ge0.
\]

So the outgoing branch is passive whenever \(N_0\ge0\) and
\(\Gamma_2^{\rm port}>0\).

---

## 5. Failure modes

The audit isolates four concrete failure modes.

### 5.1 Internal mixed-block instability

If

\[
\Delta\le0,
\]

then the \(U/W\) internal block is indefinite.  The Maxwell/mixed support pair is
not a positive-energy conservative subsystem.

### 5.2 Static wall softening instability

If

\[
D_0\le0,
\]

then the Schur-complement stiffness of the wall is nonpositive.  The apparent
normalization enhancement from \(P_0=N_0/D_0\) has crossed into instability.

### 5.3 Ghost/Krein contamination

If the reduced BdG coordinates are not projected onto positive-energy stable
normal modes, the kinetic or symplectic signature assumption behind this
Hamiltonian ledger fails.  Such modes require a separate full BdG signature
audit.

### 5.4 Dark outgoing port

If

\[
\Omega_U^2g_W+Rg_U=0,
\]

then

\[
N_0=0.
\]

The branch can still be conservative and stable, but it cannot transfer the
leading outgoing \(l=2\) port to the wall at this order.

---

## 6. Numerical sanity checks

The script includes one stable rational example with

\[
D_0>0,\qquad \Delta>0,
\]

and all generalized eigenvalues \(\omega^2\) positive.

It then lowers \(K\) while holding the internal data fixed.  The resulting
example has

\[
D_0<0,
\]

and one generalized \(\omega^2\) eigenvalue becomes negative, confirming that
the Schur-complement gate is the static instability boundary.

It also includes a mixed-block example with \(\Delta<0\), showing that the
localized Maxwell/mixed internal block itself becomes indefinite.

---

## 7. Carry-forward statement

The combined reduced stability theorem is:

\[
\boxed{
T>0,\quad A_{\rm int}>0,\quad
D_0=K-g_{\rm int}^TA_{\rm int}^{-1}g_{\rm int}>0
\quad\Longrightarrow\quad
H_{\rm cons}>0.
}
\]

In scalar one-lane form,

\[
\boxed{
M>0,\quad
\varpi^2>0,\quad
\Omega_U^2\Omega_W^2-R^2>0,\quad
K-\frac{c_B^2}{\varpi^2}
-
\frac{
g_U^2\Omega_W^2+2g_Ug_WR+g_W^2\Omega_U^2
}{
\Omega_U^2\Omega_W^2-R^2
}>0.
}
\]

The outgoing \(l=2\) port is passive if

\[
\boxed{
\gamma_{\rm wall}
=
\frac{a^5}{27c_s^5}
\frac{(\Omega_U^2g_W+Rg_U)^2}
{(\Omega_U^2\Omega_W^2-R^2)^2}
\ge0.
}
\]

So the normalization strategy may use near-softening, but only with

\[
D_0>0
\]

kept as a hard branch gate.

stage_v2_11_grouped_p2_projectors_derivation.md
# Stage V2-11 — Grouped Real \(P_2\) Projector Calculus

## 0. Purpose

This stage freezes the algebra of the grouped real \(P_2\) bundle before the next normalization and source-map stages.

The goal is not to solve a radial/axial moving-throat branch. The goal is to verify the exact bookkeeping used whenever the five real quadrupole lanes

\[
(20,\;21c,\;21s,\;22c,\;22s)
\]

are compressed into the grouped triple

\[
x=(x_{20},x_{21},x_{22}).
\]

The \(21\) and \(22\) grouped lanes each contain two real components, so the grouped triple is not Euclidean with equal weights. Its natural metric is

\[
G_{\rm grp}=\operatorname{diag}(1,2,2).
\]

This stage verifies:

1. the grouped metric;
2. the trace/anomaly basis;
3. the exact projectors;
4. the inverse grouped-coordinate map;
5. the anisotropy norm;
6. the weak-axisymmetric splitting law;
7. first-order transport formulas for \(u_2=-D_2/D_0\) and \(P_0=N_0/D_0\).

The SymPy audit passes all checks.

---

## 1. Grouped metric from the five real lanes

Let

\[
y=(y_{20},y_{21c},y_{21s},y_{22c},y_{22s})^T
\]

and impose the grouped embedding

\[
y_{20}=x_{20},\qquad
y_{21c}=y_{21s}=x_{21},\qquad
y_{22c}=y_{22s}=x_{22}.
\]

The embedding matrix is

\[
E=
\begin{pmatrix}
1&0&0\\
0&1&0\\
0&1&0\\
0&0&1\\
0&0&1
\end{pmatrix}.
\]

Therefore the grouped inner product inherited from the five-lane Euclidean inner product is

\[
G_{\rm grp}=E^T E
=
\begin{pmatrix}
1&0&0\\
0&2&0\\
0&0&2
\end{pmatrix}.
\]

So every grouped trace, orthogonality condition, and norm must use \(G_{\rm grp}\), not the naive identity metric.

---

## 2. Trace/anomaly basis

The useful grouped basis is

\[
e_{\rm tr}=(1,1,1)^T,
\]

\[
e_a=(4,-1,-1)^T,
\]

\[
e_b=(0,1,-1)^T.
\]

The script verifies exact \(G_{\rm grp}\)-orthogonality:

\[
e_{\rm tr}^TG_{\rm grp}e_a=0,
\qquad
e_{\rm tr}^TG_{\rm grp}e_b=0,
\qquad
e_a^TG_{\rm grp}e_b=0.
\]

The exact squared norms are

\[
e_{\rm tr}^TG_{\rm grp}e_{\rm tr}=5,
\]

\[
e_a^TG_{\rm grp}e_a=20,
\]

\[
e_b^TG_{\rm grp}e_b=4.
\]

Thus every grouped triple decomposes uniquely as

\[
x=\bar x\,e_{\rm tr}+a_x\,e_a+b_x\,e_b.
\]

The exact coordinates are

\[
\boxed{
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}5,
}
\]

\[
\boxed{
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
}
\]

\[
\boxed{
b_x=\frac{x_{21}-x_{22}}2.
}
\]

The inverse map is

\[
\boxed{
x_{20}=\bar x+4a_x,
}
\]

\[
\boxed{
x_{21}=\bar x-a_x+b_x,
}
\]

\[
\boxed{
x_{22}=\bar x-a_x-b_x.
}
\]

---

## 3. Exact projectors

The \(G_{\rm grp}\)-orthogonal projector onto a basis vector \(e\) is

\[
P_e=\frac{e\,e^TG_{\rm grp}}{e^TG_{\rm grp}e}.
\]

The script verifies the three exact projectors:

\[
P_{\rm tr}
=
\frac15
\begin{pmatrix}
1&2&2\\
1&2&2\\
1&2&2
\end{pmatrix},
\]

\[
P_a
=
\frac1{20}
\begin{pmatrix}
16&-8&-8\\
-4&2&2\\
-4&2&2
\end{pmatrix},
\]

\[
P_b
=
\frac14
\begin{pmatrix}
0&0&0\\
0&2&-2\\
0&-2&2
\end{pmatrix}.
\]

They obey

\[
P_{\rm tr}+P_a+P_b=I_3,
\]

\[
P_i^2=P_i,
\]

\[
P_iP_j=0\qquad(i\neq j),
\]

and the \(G_{\rm grp}\)-self-adjointness condition

\[
P_i^TG_{\rm grp}=G_{\rm grp}P_i.
\]

So these are true \(G_{\rm grp}\)-orthogonal projectors.

---

## 4. Isotropy and anisotropy norm

The isotropic grouped branch is exactly

\[
a_x=0,\qquad b_x=0,
\]

which is equivalent to

\[
x_{20}=x_{21}=x_{22}.
\]

The grouped norm decomposes as

\[
x^TG_{\rm grp}x
=
5\bar x^2+20a_x^2+4b_x^2.
\]

The normalized anisotropy norm used downstream is the anisotropic part divided by the total grouped multiplicity \(5\):

\[
\boxed{
A_x^2
=
\frac{20a_x^2+4b_x^2}{5}
=
4a_x^2+\frac45 b_x^2.
}
\]

This is the exact norm behind the later grouped anisotropy measures.

---

## 5. Weak-axisymmetric splitting fingerprint

For a weak axisymmetric quadrupolar perturbation, the grouped splitting vector is

\[
\lambda=(1,\tfrac12,-1).
\]

The script verifies

\[
\bar\lambda=0,
\]

\[
a_\lambda=\frac14,
\]

\[
b_\lambda=\frac34.
\]

Therefore

\[
\boxed{
b_\lambda=3a_\lambda.
}
\]

So any first-order grouped coefficient of the form

\[
x_A=x^{(0)}+\epsilon\lambda_Ax^{(1)}
\]

has

\[
\bar x=x^{(0)},
\]

\[
a_x=\frac{\epsilon}{4}x^{(1)},
\]

\[
b_x=\frac{3\epsilon}{4}x^{(1)},
\]

and hence

\[
\boxed{b_x=3a_x.}
\]

This is the diagnostic line for pure weak-axisymmetric \(l=2\) symmetry breaking.

---

## 6. First-order transport for the normalized response \(u_2\)

Let

\[
u_2=-\frac{D_2}{D_0}.
\]

Perturb a grouped lane by

\[
D_{A0}=D_0+\epsilon\lambda_A\delta D_0,
\]

\[
D_{A2}=D_2+\epsilon\lambda_A\delta D_2.
\]

Then

\[
u_2^{(A)}
=
-\frac{D_{A2}}{D_{A0}}
=
u_2+\epsilon\lambda_A u_2^{(1)}+O(\epsilon^2),
\]

with

\[
\boxed{
u_2^{(1)}
=
-\frac{\delta D_2+u_2\,\delta D_0}{D_0}.
}
\]

Equivalently, using \(u_2=-D_2/D_0\),

\[
u_2^{(1)}
=
-\frac{D_0\delta D_2-D_2\delta D_0}{D_0^2}.
\]

On the weak-axisymmetric line, the \(u_2\) anisotropy obeys

\[
b_{u_2}=3a_{u_2}.
\]

---

## 7. First-order transport for the outgoing prefactor \(P_0\)

Let

\[
P_0=\frac{N_0}{D_0}.
\]

Perturb a grouped lane by

\[
N_{A0}=N_0+\epsilon\lambda_A\delta N_0,
\]

\[
D_{A0}=D_0+\epsilon\lambda_A\delta D_0.
\]

Then

\[
P_0^{(A)}
=
\frac{N_{A0}}{D_{A0}}
=
P_0+\epsilon\lambda_A P_0^{(1)}+O(\epsilon^2),
\]

with

\[
\boxed{
P_0^{(1)}
=
\frac{\delta N_0-P_0\delta D_0}{D_0}.
}
\]

Equivalently,

\[
P_0^{(1)}
=
\frac{D_0\delta N_0-N_0\delta D_0}{D_0^2}.
\]

On the weak-axisymmetric line, the outgoing-prefactor anisotropy also obeys

\[
b_{P_0}=3a_{P_0}.
\]

This is the exact \(P_1/P_0\) bookkeeping used later when \(\Xi_1\) is interpreted as the weak-axisymmetric outgoing-prefactor slope.

---

## 8. Result and carry-forward status

The script reports:

```text
checks_total: 41
checks_passed: 41
checks_failed: 0
```

So V2-11 gives a clean pass.

The frozen carry-forward packet is:

\[
G_{\rm grp}=\operatorname{diag}(1,2,2),
\]

\[
x=\bar x(1,1,1)+a_x(4,-1,-1)+b_x(0,1,-1),
\]

\[
A_x^2=4a_x^2+\frac45b_x^2,
\]

\[
\lambda_{\rm ax}=(1,\tfrac12,-1),
\qquad
b=3a,
\]

\[
u_2^{(1)}
=
-\frac{\delta D_2+u_2\delta D_0}{D_0},
\]

\[
P_0^{(1)}
=
\frac{\delta N_0-P_0\delta D_0}{D_0}.
\]

This stage is independent of the radial/axial open-throat boundary correction from V2-04. The open organ-pipe patch changes the radial/axial branch data inside \(D_n\) and \(N_n\), but not the grouped angular projector calculus.

stage_v2_12_stf_angular_source_map_derivation.md
# Stage V2-12 — STF Angular Source-Map Theorem

## 0. Purpose

This stage audits whether the surviving orbital/worldtube STF quadrupole branch has an **angular normalization defect** when mapped into the moving-throat grouped real \(P_2\) ports.

The result is:

\[
\boxed{\widehat m_{\rm ang}=1}
\]

in the canonical normalized real-STF basis. Therefore the remaining quadrupole-normalization problem is not angular. It is radial/axial, source-amplitude, and outgoing-port normalization data.

This stage is deliberately narrow. It does not compute the radial overlap integrals, the Maxwell/mixed transfer prefactor, or the final product \( \widehat m_0^2P_0 \). It proves that the angular map itself is an identity.

---

## 1. Canonical real STF basis

Let \(n^i\) be the unit direction on the throat/mouth sphere. Use five real symmetric trace-free tensors \(E_A^{ij}\), with

\[
A\in\{20,21c,21s,22c,22s\}.
\]

The script uses

\[
E_{20}
=
\frac{1}{\sqrt6}
\begin{pmatrix}
-1&0&0\\
0&-1&0\\
0&0&2
\end{pmatrix},
\]

\[
E_{21c}
=
\frac1{\sqrt2}
\begin{pmatrix}
0&0&1\\
0&0&0\\
1&0&0
\end{pmatrix},
\qquad
E_{21s}
=
\frac1{\sqrt2}
\begin{pmatrix}
0&0&0\\
0&0&1\\
0&1&0
\end{pmatrix},
\]

\[
E_{22c}
=
\frac1{\sqrt2}
\begin{pmatrix}
1&0&0\\
0&-1&0\\
0&0&0
\end{pmatrix},
\qquad
E_{22s}
=
\frac1{\sqrt2}
\begin{pmatrix}
0&1&0\\
1&0&0\\
0&0&0
\end{pmatrix}.
\]

The audit verifies

\[
{\rm Tr}(E_A)=0,
\]

and

\[
{\rm Tr}(E_AE_B)=\delta_{AB}.
\]

So these five tensors are an orthonormal real STF basis.

---

## 2. Normalized real \(l=2\) angular functions

Define

\[
Y_A(n)
=
\sqrt{\frac{15}{8\pi}}\,
E_A^{ij}n_in_j.
\]

The fourth-moment identity on the unit sphere is

\[
\int_{S^2}n_in_jn_kn_l\,d\Omega
=
\frac{4\pi}{15}
\left(
\delta_{ij}\delta_{kl}
+\delta_{ik}\delta_{jl}
+\delta_{il}\delta_{jk}
\right).
\]

Because \(E_A\) is trace-free,

\[
\int_{S^2}
(E_A^{ij}n_in_j)(E_B^{kl}n_kn_l)\,d\Omega
=
\frac{8\pi}{15}
{\rm Tr}(E_AE_B).
\]

Therefore

\[
\int_{S^2}Y_A(n)Y_B(n)\,d\Omega
=
{\rm Tr}(E_AE_B)
=
\delta_{AB}.
\]

The script verifies this twice:

1. by the fourth-moment identity;
2. by direct \((\theta,\phi)\) integration of the explicit real harmonics.

Both checks return the \(5\times5\) identity matrix.

---

## 3. Source-map theorem

Let the angular part of the orbital/worldtube STF source be

\[
S(n)=\sum_B S_B Y_B(n).
\]

The port projection onto the same real harmonic basis is

\[
S_A^{\rm port}
=
\int_{S^2}Y_A(n)S(n)\,d\Omega.
\]

Using orthonormality,

\[
S_A^{\rm port}
=
\sum_B S_B
\int_{S^2}Y_A(n)Y_B(n)\,d\Omega
=
S_A.
\]

Thus the angular source-map matrix is exactly

\[
\boxed{
M_{AB}^{\rm ang}=\delta_{AB}.
}
\]

There is no angular rescaling, no angular kernel, and no missing angular component.

---

## 4. Orbital/worldtube STF quadrupole reconstruction

For a point/worldtube quadrupole

\[
Q_{ij}
=
\mu
\left(
x_ix_j-\frac{r^2}{3}\delta_{ij}
\right),
\qquad
r^2=x^2+y^2+z^2,
\]

the real STF coefficients are

\[
Q_A={\rm Tr}(E_AQ).
\]

The script derives

\[
Q_{20}
=
-\frac{\sqrt6\mu}{6}(x^2+y^2-2z^2),
\]

\[
Q_{21c}=\sqrt2\,\mu xz,
\qquad
Q_{21s}=\sqrt2\,\mu yz,
\]

\[
Q_{22c}=\frac{\sqrt2\mu}{2}(x^2-y^2),
\qquad
Q_{22s}=\sqrt2\,\mu xy.
\]

It then checks the reconstruction

\[
Q_{ij}=\sum_A Q_AE_{Aij}.
\]

The residual is exactly zero.

It also checks the angular identity

\[
\sqrt{\frac{15}{8\pi}}Q_{ij}n^in^j
=
\sum_A Q_A Y_A(n),
\]

with exact zero residual.

---

## 5. Norm and grouped metric

For the full five-mode real STF source,

\[
\int_{S^2}S(n)^2\,d\Omega
=
S_{20}^2+S_{21c}^2+S_{21s}^2+S_{22c}^2+S_{22s}^2.
\]

When the two \(21\) real components are grouped together and the two \(22\) real components are grouped together, a representative grouped vector

\[
x_{\rm grp}=(x_{20},x_{21},x_{22})
\]

corresponds to the full vector

\[
(x_{20},x_{21},x_{21},x_{22},x_{22}).
\]

Therefore the full five-mode norm becomes

\[
x_{20}^2+2x_{21}^2+2x_{22}^2
=
x_{\rm grp}^T
\operatorname{diag}(1,2,2)
x_{\rm grp}.
\]

This is the grouped metric

\[
\boxed{
G_{\rm grp}=\operatorname{diag}(1,2,2).
}
\]

So the grouped metric used in the \(20/21/22\) projector calculus is exactly the norm inherited from the five real STF modes.

---

## 6. Convention bridge to the older \(\Pi\) notation

The 2.5PN package uses an older real-STF component convention in which the normalized port variables are related by

\[
q_{20}=\sqrt{\frac23}\Pi_{20},
\]

\[
q_{21c}=\Pi_{21c},
\qquad
q_{21s}=\Pi_{21s},
\]

\[
q_{22c}=2\Pi_{22c},
\qquad
q_{22s}=2\Pi_{22s}.
\]

The conversion matrix is

\[
B_{\Pi\to q}
=
\operatorname{diag}
\left(
\sqrt{\frac23},1,1,2,2
\right).
\]

The script verifies

\[
\det B_{\Pi\to q}
=
\frac{4\sqrt6}{3}
\neq0,
\]

and

\[
{\rm rank}(B_{\Pi\to q})=5.
\]

So even in the older \(\Pi\) convention, the map is invertible. It introduces only a convention conversion, not an angular obstruction.

---

## 7. Interpretation

The audit proves the angular part of the quadrupole map is complete and exactly normalized in the canonical real-STF basis:

\[
\boxed{
\widehat m_{\rm ang}=1.
}
\]

Therefore the remaining source normalization can be factored as

\[
\widehat m_0
=
\widehat m_{\rm ang}\widehat m_{\rm rad}
=
\widehat m_{\rm rad}.
\]

The unresolved product is still

\[
\widehat m_0^2 P_0
=
\frac{54Gc_s^5}{5a^5c^5},
\]

but V2-12 shows that any miss in this equation must come from

1. radial/axial overlap amplitudes,
2. mouth/worldtube source amplitude normalization,
3. the Maxwell/mixed transfer factor \(N_0\),
4. the conservative denominator \(D_0\),
5. or the final moving-throat branch selection.

It cannot be blamed on a missing angular STF component or an angular normalization mismatch.

---

## 8. SymPy result

The script reports:

```text
all_STF_traces_zero: PASS
STF_basis_orthonormal: PASS
angular_gram_from_moments_identity: PASS
direct_angular_gram_identity: PASS
source_projection_identity: PASS
orbital_STF_reconstructs_exactly: PASS
angular_function_reconstructs_exactly: PASS
source_norm_identity: PASS
grouped_metric_matches_full_pair_norm: PASS
Pi_to_q_map_invertible: PASS
```

Final verdict:

\[
\boxed{
\text{The STF angular source-map has no angular normalization defect.}
}
\]

\[
\boxed{
\text{The remaining quadrupole normalization gap is radial/axial and dynamical.}
}
\]

stage_v2_13_grouped_normalization_ratio_derivation.md
# Stage V2-13 — Grouped Normalization Ratio Audit

## Purpose

This stage freezes the exact algebraic bridge from grouped real `P2` conservative operator moments and outgoing-transfer moments to the normalized quantities used by the 2.5PN/4PN quadrupole-normalization package.

The stage treats one grouped lane first, then specializes to the isotropic full-bundle branch. It does **not** solve the moving-throat PDE. It identifies the exact target surface that the actual branch must hit.

The retained low-frequency inputs are

\[
D(\omega)=D_0+D_2\omega^2+D_4\omega^4+O(\omega^6),
\]

and

\[
N(\omega)=N_0+N_2\omega^2+N_4\omega^4+O(\omega^6).
\]

Here `D` is the conservative wall/worldtube operator and `N` is the outgoing-transfer numerator inherited from the Maxwell/mixed port.

A port-normalization factor \(\mathcal S_{\rm port}\) is carried explicitly. If `N` has already been gravitationally/port-normalized, set

\[
\mathcal S_{\rm port}=1.
\]

If `N` is the raw canonical mechanical transfer numerator, keep \(\mathcal S_{\rm port}\) until the comparison with the GR target.

---

## 1. Normalized conservative response

Define the normalized conservative response

\[
Y(\omega)=\frac{D_0}{D(\omega)}.
\]

Expanding,

\[
Y(\omega)=1+u_2\omega^2+u_4\omega^4+O(\omega^6).
\]

The exact conversion formulas are

\[
\boxed{u_2=-\frac{D_2}{D_0}},
\]

\[
\boxed{u_4=\frac{D_2^2-D_0D_4}{D_0^2}}.
\]

These formulas are purely algebraic and are independent of the microscopic origin of the coefficients.

---

## 2. Outgoing prefactor moments

The outgoing transfer appears through the internal prefactor

\[
{\rm Pref}(\omega)=\frac{D_0N(\omega)}{D(\omega)^2}.
\]

Write

\[
{\rm Pref}(\omega)=P_0+P_2\omega^2+P_4\omega^4+O(\omega^6).
\]

Then

\[
\boxed{P_0=\frac{N_0}{D_0}},
\]

\[
\boxed{P_2=\frac{D_0N_2-2D_2N_0}{D_0^2}},
\]

\[
\boxed{
P_4=
\frac{
D_0^2N_4
-2D_0(D_2N_2+D_4N_0)
+3D_2^2N_0
}{D_0^3}.
}
\]

Thus the static outgoing prefactor is the simple ratio

\[
\boxed{P_0=N_0/D_0}.
\]

This is the scalar that controls the leading universal quadrupole channel.

---

## 3. Constant-prefactor branch

The constant-prefactor branch through \(O(\omega^4)\) is defined by

\[
P_2=0,
\qquad
P_4=0.
\]

The exact branch conditions are

\[
\boxed{N_2=\frac{2D_2N_0}{D_0}},
\]

and

\[
\boxed{N_4=\frac{N_0(D_2^2+2D_0D_4)}{D_0^2}}.
\]

Equivalently, if the first condition is kept symbolic, the second may be written as

\[
N_4=
\frac{2D_0(D_2N_2+D_4N_0)-3D_2^2N_0}{D_0^2}.
\]

The important point is that the higher transfer moments do not need to vanish. They must be correlated with the conservative moments.

---

## 4. Multiplication by the compact outgoing \(l=2\) fingerprint

The compact outgoing quadrupole branch has normalized low-frequency form

\[
\widehat Y_2^{\rm out}(\omega)
=1+A\omega^2+B\omega^4+iG_5\omega^5+O(\omega^6),
\]

where, on the canonical compact branch,

\[
A=\frac{a^2}{9c_s^2},
\qquad
B=\frac{4a^4}{81c_s^4},
\qquad
G_5=\frac{a^5}{27c_s^5}.
\]

Multiplying by the internal prefactor gives

\[
\delta Y^{\rm out}(\omega)
={\rm Pref}(\omega)\widehat Y_2^{\rm out}(\omega).
\]

Thus

\[
K_0=P_0,
\]

\[
K_2=P_2+AP_0,
\]

\[
K_4=P_4+AP_2+BP_0,
\]

and

\[
\boxed{\Gamma_5=G_5P_0}.
\]

Therefore the leading odd 2.5PN coefficient depends only on the static prefactor \(P_0\), not on \(P_2\) or \(P_4\).

---

## 5. Isotropic full-bundle surface

On the isotropic grouped branch, the full coupled conservative operator is

\[
D_0=K-B_0-Z_0,
\]

\[
D_2=-(M+B_2+Z_2),
\]

\[
D_4=-(B_4+Z_4).
\]

Here:

- \(K,M\) are wall stiffness and inertia data,
- \(B_n\) are stable BdG support moments,
- \(Z_n\) are conservative Maxwell/mixed moments.

The one-pole conservative condition

\[
u_4=4u_2^2
\]

is equivalent to

\[
\boxed{
D_0(B_4+Z_4)=3(M+B_2+Z_2)^2.
}
\]

So the one-pole surface fixes the static denominator as

\[
D_0^{\rm one\ pole}
=
\frac{3(M+B_2+Z_2)^2}{B_4+Z_4}.
\]

The universal quadrupole target is

\[
\boxed{
\widehat m_0^{\,2}\,\mathcal S_{\rm port}\,P_0
=
\frac{54Gc_s^5}{5a^5c^5}.
}
\]

Since \(P_0=N_0/D_0\), this gives

\[
D_0^{\rm norm}
=
\frac{
\widehat m_0^{\,2}\mathcal S_{\rm port}N_0
}{T_{\rm GR}},
\]

with

\[
T_{\rm GR}=\frac{54Gc_s^5}{5a^5c^5}.
\]

The simultaneous one-pole and normalization compatibility surface is therefore

\[
\boxed{
\frac{
\widehat m_0^{\,2}\mathcal S_{\rm port}N_0
}{T_{\rm GR}}
=
\frac{3(M+B_2+Z_2)^2}{B_4+Z_4}.
}
\]

This is the exact isotropic scalar target surface for the reduced moving-throat branch.

---

## 6. GR coefficient check

Using

\[
\widehat m_0^{\,2}\mathcal S_{\rm port}P_0
=
\frac{54Gc_s^5}{5a^5c^5},
\]

and

\[
G_5=\frac{a^5}{27c_s^5},
\]

one obtains

\[
\widehat m_0^{\,2}\mathcal S_{\rm port}\Gamma_5
=
\frac{54Gc_s^5}{5a^5c^5}\frac{a^5}{27c_s^5}
=
\boxed{\frac{2G}{5c^5}}.
\]

So the normalization target is exactly equivalent to the standard quadrupole coefficient.

---

## 7. Weak-axisymmetric grouped transport

For a weak axisymmetric grouped perturbation with signature

\[
(20,21,22)\sim(1,\tfrac12,-1),
\]

let

\[
P_A=P_{\rm base}+\epsilon\lambda_A P_1.
\]

Then the grouped trace and anisotropy coordinates are

\[
\bar P=P_{\rm base},
\]

\[
a_P=\frac{\epsilon P_1}{4},
\]

\[
b_P=\frac{3\epsilon P_1}{4}.
\]

Therefore

\[
\boxed{b_P=3a_P}.
\]

This confirms that the weak-axisymmetric prefactor slope remains a one-scalar defect, not a generic three-lane deformation.

---

## 8. Stage result

The audit verifies 15 symbolic checks:

- normalized response formulas,
- outgoing prefactor formulas,
- constant-prefactor branch conditions,
- outgoing \(l=2\) multiplication,
- isotropic one-pole target surface,
- universal normalization to GR quadrupole coefficient,
- weak-axisymmetric grouped trace and \(b=3a\) law.

Final status:

\[
\boxed{\text{V2-13 passes as exact algebra.}}
\]

The remaining issue is not algebraic. It is branch realization: the actual moving-throat PDE must supply stable isotropic data satisfying the compatibility surface above.

stage_v2_14_outgoing_l2_fingerprint_derivation.md
# Stage V2-14 — Compact outgoing `l=2` fingerprint derivation

## Purpose

This stage rederives the compact outgoing quadrupole fingerprint

\[
\widehat Y_2^{\rm out}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6)
\]

from the outgoing spherical Hankel solution.  The main convention fixed here is:

\[
\boxed{
\widehat Y_2^{\rm out}
\text{ is the normalized inverse DtN compliance, not the raw DtN admittance.}
}
\]

If one uses the raw normalized Dirichlet-to-Neumann map, the signs of the
low-frequency even terms and the leading odd term are reversed before inversion.

---

## 1. Setup

Use the time convention

\[
e^{-i\omega t},
\]

so the outgoing radial solution in three spatial dimensions is the spherical
Hankel function

\[
h_l^{(1)}(kr).
\]

Let

\[
z=ka=\frac{a\omega}{c_s}.
\]

For \(l=2\), the exact outgoing function can be written as

\[
h_2^{(1)}(z)
=
i\,e^{iz}\frac{z^2+3iz-3}{z^3}.
\]

Its static singular behavior is

\[
h_2^{(1)}(z)\sim -\frac{3i}{z^3},
\]

which is the outgoing continuation of the exterior static quadrupole field
\(\sim r^{-3}\).

---

## 2. Raw normalized DtN map

For a boundary datum at \(r=a\), the outgoing exterior field is proportional to

\[
f(r)=h_2^{(1)}(kr).
\]

The raw DtN ratio is

\[
a\frac{f'(a)}{f(a)}
=
z\frac{h_2'(z)}{h_2(z)}.
\]

Since the static \(l=2\) exterior field scales as \(r^{-3}\), the normalized raw
DtN is

\[
\widehat D_2^{\rm raw}(z)
=
-\frac{z}{3}\frac{h_2'(z)}{h_2(z)}.
\]

SymPy gives

\[
\widehat D_2^{\rm raw}(z)
=
1-\frac{z^2}{9}
-\frac{z^4}{27}
-i\frac{z^5}{27}
+O(z^6).
\]

So the quoted bridge fingerprint is **not** this raw DtN object.

---

## 3. Normalized inverse DtN compliance

The compact outgoing bridge uses the inverse normalized compliance

\[
\widehat Y_2^{\rm out}(z)
=
\left(\widehat D_2^{\rm raw}(z)\right)^{-1}
=
-\,\frac{3h_2^{(1)}(z)}{z\,h_2^{(1)\prime}(z)}.
\]

The exact rational form is

\[
\widehat Y_2^{\rm out}(z)
=
\frac{
3(3-3iz-z^2)
}{
9-4z^2+i(z^3-9z)
}.
\]

Expanding at low frequency,

\[
\widehat Y_2^{\rm out}(z)
=
1+\frac{z^2}{9}
+\frac{4z^4}{81}
+i\frac{z^5}{27}
-\frac{11z^6}{729}
-i\frac{z^7}{243}
+O(z^8).
\]

Therefore, through the order needed by the 2.5PN and 4PN bridge,

\[
\boxed{
\widehat Y_2^{\rm out}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
}
\]

The leading odd coefficient is

\[
\boxed{
\Gamma_5^{\rm port}
=
\frac{a^5}{27c_s^5}.
}
\]

---

## 4. Channel interpretation

This result is a theorem for a **three-spatial-dimensional exterior outgoing
worldtube/STF quadrupole port**.

It is not automatically the result for an unrestricted \(4\)-spatial-dimensional
bulk radiation channel.  In \(d\) spatial dimensions, the small-\(z\) absorptive
power of a clean outgoing \(l\)-partial wave scales as

\[
z^{2l+d-2}.
\]

Thus

\[
d=3,\ l=2
\quad\Rightarrow\quad
z^5,
\]

while

\[
d=4,\ l=2
\quad\Rightarrow\quad
z^6
\]

(up to the usual integer-Bessel logarithmic subtleties).  Therefore the
\(i\omega^5\) branch is specifically the brane/worldtube \(3\)-space STF
quadrupole port.  A future full-PDE realization stage must show that the actual
moving-throat branch projects onto this port rather than an unrestricted
\(4\)-space bulk outgoing law.

---

## 5. Audit result

The SymPy script verifies:

1. the exact \(h_2^{(1)}\) closed form,
2. the static normalization,
3. the raw normalized DtN expansion,
4. the inverse-compliance expansion,
5. the coefficients \(1/9\), \(4/81\), and \(i/27\),
6. the leading coefficient \(\Gamma_5=a^5/(27c_s^5)\),
7. and the channel discriminator \(2l+d-2\).

The stage passes with one wording patch:

\[
\boxed{
\text{call the object a normalized outgoing compliance / inverse DtN response.}
}
\]

Calling the same expansion the raw DtN admittance would be incorrect.

stage_v2_15_25pn_4pn_interface_derivation.md
# Stage V2-15 — 2.5PN / 4PN Quadrupole-Normalization Interface Audit

## 0. Purpose

This stage checks whether the conservative 4PN hereditary/tail sector introduces a new quadrupole normalization problem, or whether it is controlled by the same canonically normalized outgoing STF quadrupole branch already isolated at 2.5PN.

The result is algebraic:

\[
C_{\rm tail}^{\rm GR}
=
\frac{GM}{2c^3}\gamma_{\rm GR},
\qquad
\gamma_{\rm GR}=\frac{2G}{5c^5},
\]

so

\[
C_{\rm tail}^{\rm GR}
=
\frac{G^2M}{5c^8}.
\]

Thus, once the model’s effective quadrupole coefficient equals the 2.5PN Burke--Thorne coefficient, the 4PN hereditary coefficient is fixed.

The only extra scalar that can remain in a toy-model tail bridge is a **tail-transport** factor, not a new grouped-\(P_2\) normalization datum.

---

## 1. Inputs from earlier Volume 2 stages

From V2-13 and V2-14, the canonically normalized moving-throat outgoing STF quadrupole coefficient is

\[
\gamma_{\rm quad}^{\rm eff}
=
P_{\rm eff}\frac{a^5}{27c_s^5},
\]

where

\[
P_{\rm eff}
=
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}.
\]

The universal 2.5PN target is

\[
\gamma_{\rm quad}^{\rm eff}
=
\gamma_{\rm GR}
=
\frac{2G}{5c^5}.
\]

Equivalently,

\[
P_{\rm eff}
=
\frac{54Gc_s^5}{5a^5c^5}.
\]

This is the same target previously written as

\[
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
=
\frac{54Gc_s^5}{5a^5c^5}.
\]

---

## 2. GR tail coefficient bridge

The conservative 4PN tail coefficient has the GR reference value

\[
C_{\rm tail}^{\rm GR}
=
\frac{G^2M}{5c^8}.
\]

Using the 2.5PN coefficient

\[
\gamma_{\rm GR}
=
\frac{2G}{5c^5},
\]

one obtains

\[
\frac{GM}{2c^3}\gamma_{\rm GR}
=
\frac{GM}{2c^3}\frac{2G}{5c^5}
=
\frac{G^2M}{5c^8}
=
C_{\rm tail}^{\rm GR}.
\]

So the exact interface identity is

\[
\boxed{
C_{\rm tail}^{\rm GR}
=
\frac{GM}{2c^3}\gamma_{\rm GR}.
}
\]

---

## 3. Moving-throat / toy-model tail bridge

A conservative symbolic bridge for the model is

\[
C_{\rm tail}^{\rm toy}
=
\Theta_{\rm tail}
\frac{GM}{2c_s^3}\gamma_{\rm quad}^{\rm eff}.
\]

Here:

- \(\gamma_{\rm quad}^{\rm eff}\) is the same STF quadrupole coefficient used at 2.5PN;
- \(c_s\) is the propagation/sound speed used in the outgoing port normalization;
- \(\Theta_{\rm tail}\) is a scalar monopole-scattering/tail-transport factor.

The ratio to the GR tail coefficient is

\[
\frac{C_{\rm tail}^{\rm toy}}{C_{\rm tail}^{\rm GR}}
=
\Theta_{\rm tail}
\left(\frac{c}{c_s}\right)^3
\frac{\gamma_{\rm quad}^{\rm eff}}{\gamma_{\rm GR}}.
\]

Therefore, if the 2.5PN target is already satisfied,

\[
\frac{\gamma_{\rm quad}^{\rm eff}}{\gamma_{\rm GR}}=1,
\]

the remaining tail condition is

\[
\boxed{
\Theta_{\rm tail}
\left(\frac{c}{c_s}\right)^3
=
1.
}
\]

On the \(c_s=c\) branch, this reduces to

\[
\boxed{
\Theta_{\rm tail}=1.
}
\]

So 4PN does not add another grouped-quadrupole normalization. It adds, at most, a scalar tail-transport gate.

---

## 4. Constant-prefactor check

The outgoing \(l=2\) fingerprint is

\[
\widehat Y_2^{\rm out}
=
1
+
\frac{a^2\omega^2}{9c_s^2}
+
\frac{4a^4\omega^4}{81c_s^4}
+
i\frac{a^5\omega^5}{27c_s^5}
+
O(\omega^6).
\]

Let the internal prefactor be

\[
P(\omega)=P_0+P_2\omega^2+P_4\omega^4+\cdots.
\]

Then

\[
P(\omega)\widehat Y_2^{\rm out}(\omega)
\]

has leading imaginary odd coefficient

\[
i\omega^5
\left[
P_0\frac{a^5}{27c_s^5}
\right].
\]

The \(P_2\) and \(P_4\) prefactor moments first affect higher odd orders:

\[
P_2\omega^2\cdot i\omega^5=O(i\omega^7),
\qquad
P_4\omega^4\cdot i\omega^5=O(i\omega^9).
\]

So the leading 2.5PN / 4PN interface uses only the static prefactor \(P_0=N_0/D_0\), not the higher even prefactor moments.

---

## 5. Residual bookkeeping

Define fractional deviations

\[
\frac{\gamma_{\rm quad}^{\rm eff}}{\gamma_{\rm GR}}=1+\delta_Q,
\]

and

\[
\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3=1+\delta_{\rm tail}.
\]

Then

\[
\frac{C_{\rm tail}^{\rm toy}}{C_{\rm tail}^{\rm GR}}-1
=
(1+\delta_Q)(1+\delta_{\rm tail})-1
=
\delta_Q+\delta_{\rm tail}+\delta_Q\delta_{\rm tail}.
\]

At linear order,

\[
\frac{\Delta C_{\rm tail}}{C_{\rm tail}^{\rm GR}}
\simeq
\delta_Q+\delta_{\rm tail}.
\]

This is a useful diagnostic split:

- \(\delta_Q\) is the quadrupole-normalization miss already seen at 2.5PN;
- \(\delta_{\rm tail}\) is the tail-transport/scattering miss.

They should not be conflated.

---

## 6. Dimension ledger

Use dimensions in \((M,L,T)\).

\[
[G]=M^{-1}L^3T^{-2},
\qquad
[c]=[c_s]=LT^{-1}.
\]

Then

\[
[\gamma_{\rm GR}]
=
\left[\frac{G}{c^5}\right]
=
M^{-1}L^{-2}T^3.
\]

The static outgoing prefactor target has

\[
[P_{\rm eff}]
=
\left[
\frac{Gc_s^5}{a^5c^5}
\right]
=
M^{-1}L^{-2}T^{-2}.
\]

Multiplying by

\[
\left[\frac{a^5}{c_s^5}\right]=T^5
\]

gives

\[
[\gamma_{\rm quad}^{\rm eff}]
=
M^{-1}L^{-2}T^3,
\]

matching \(\gamma_{\rm GR}\).

For the tail coefficient,

\[
[C_{\rm tail}^{\rm GR}]
=
\left[
\frac{G^2M}{c^8}
\right]
=
M^{-1}L^{-2}T^4.
\]

And

\[
\left[\frac{GM}{c^3}\right]=T,
\]

so

\[
\left[
\frac{GM}{2c^3}\gamma_{\rm GR}
\right]
=
M^{-1}L^{-2}T^4.
\]

The bridge is dimensionally consistent.

---

## 7. SymPy audit result

The audit script verifies:

1. \(P_{\rm eff}=54Gc_s^5/(5a^5c^5)\) converts \(\gamma_{\rm quad}^{\rm eff}\) into \(\gamma_{\rm GR}\).
2. \(C_{\rm tail}^{\rm GR}=(GM/2c^3)\gamma_{\rm GR}\).
3. The toy/GR tail ratio factorizes as
   \[
   \Theta_{\rm tail}(c/c_s)^3(\gamma_{\rm quad}^{\rm eff}/\gamma_{\rm GR}).
   \]
4. On \(c_s=c\), \(\Theta_{\rm tail}=1\), and the 2.5PN target, the toy tail equals the GR tail.
5. The leading \(i\omega^5\) coefficient ignores \(P_2\) and \(P_4\).
6. The tail interface ignores \(P_2\) and \(P_4\).
7. The fractional residual factorizes as
   \[
   \delta_Q+\delta_{\rm tail}+\delta_Q\delta_{\rm tail}.
   \]
8. All dimension checks pass.

The script reports:

```text
symbolic_checks_total: 9
symbolic_checks_passed: 9
dimension_checks_total: 3
dimension_checks_passed: 3
FINAL_STATUS: PASS with one explicit tail-transport gate Theta_tail*(c/c_s)^3 = 1
```

---

## 8. Carry-forward theorem statement

The clean V2-15 theorem statement is:

\[
\boxed{
C_{\rm tail}
=
\frac{GM}{2c^3}\gamma_{\rm quad}^{\rm eff}
}
\]

on the canonical \(c_s=c\), \(\Theta_{\rm tail}=1\) branch.

Therefore,

\[
\boxed{
\gamma_{\rm quad}^{\rm eff}
=
\frac{2G}{5c^5}
\quad\Longleftrightarrow\quad
C_{\rm tail}
=
\frac{G^2M}{5c^8}.
}
\]

So the conservative 4PN tail does not create a new grouped-\(P_2\) normalization target. It reuses the same 2.5PN STF quadrupole normalization.

The remaining PDE work is now split into two scalar gates:

\[
\boxed{
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
=
\frac{54Gc_s^5}{5a^5c^5}
}
\]

and

\[
\boxed{
\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3=1.
}
\]

On the \(c_s=c\) branch, the second gate is simply \(\Theta_{\rm tail}=1\).

stage_v2_16_branch_freeze_no_refit_derivation.md
# Stage V2-16 — Branch-freeze / no-refit protocol

## Purpose

This stage turns the anti-overfitting rule into an executable theorem gate.

The problem is not that the grouped `P2` algebra lacks enough variables. It has many. The problem is that, unless the branch is frozen before comparison, the same algebra contains enough freedom to make target matching look stronger than it is.

The rule for Volume 2 is therefore:

> Define the parent action, gauge convention, wall/interface action, open-exit boundary protocol, projection/source map, support-profile family, number of modes/ports, stability gates, and extraction formulas **before** evaluating any target residual.

Then, and only then, evaluate the residual packet.

---

## 1. Frozen protocol DAG

The stage encodes the protocol as a directed acyclic graph:

```text
parent_action
  -> gauge_convention
  -> wall_action_and_geometry
  -> open_boundary_protocol
  -> projection_and_source_map
  -> support_profile_family
  -> branch_solve
  -> coefficient_extraction
  -> target_residual_evaluation
  -> target_decision
```

The audit verifies three graph facts:

1. all allowed arrows point forward in the frozen order;
2. there is no path from target residuals or target decisions back into branch definitions;
3. a deliberately bad edge,
   ```text
   target_residual_evaluation -> support_profile_family
   ```
   is detected as a protocol violation.

This is the formal no-refit rule.

---

## 2. Algebraic residual packet

On the isotropic grouped `P2` branch define

\[
D_0=K-B_0-Z_0,
\]

\[
D_2=-(M+B_2+Z_2),
\]

\[
D_4=-(B_4+Z_4).
\]

The normalized conservative response moments are

\[
u_2=-\frac{D_2}{D_0},
\]

\[
u_4=\frac{D_2^2-D_0D_4}{D_0^2}.
\]

The outgoing prefactor moments are

\[
P_0=\frac{N_0}{D_0},
\]

\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]

\[
P_4=
\frac{
D_0^2N_4
-2D_0(D_2N_2+D_4N_0)
+3D_2^2N_0
}{D_0^3}.
\]

The residual packet used after freeze is:

### 2.1 One-pole conservative residual

\[
R_{\rm pole}
=
D_0(B_4+Z_4)
-
3(M+B_2+Z_2)^2.
\]

The script verifies

\[
(u_4-4u_2^2)D_0^2
=
R_{\rm pole}.
\]

So the one-pole condition is

\[
R_{\rm pole}=0.
\]

### 2.2 Constant-prefactor residuals

The constant-prefactor branch is

\[
P_2=0,
\qquad
P_4=0.
\]

The script verifies that these are equivalent to

\[
N_2=\frac{2D_2N_0}{D_0},
\]

and

\[
N_4=
\frac{
2D_0(D_2N_2+D_4N_0)-3D_2^2N_0
}{D_0^2}.
\]

### 2.3 Universal quadrupole normalization residual

With the port/source convention carried explicitly,

\[
P_{\rm eff}
=
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}.
\]

The target residual is

\[
R_{\rm norm}
=
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
-
\frac{54Gc_s^5}{5a^5c^5}.
\]

The script verifies that

\[
R_{\rm norm}=0
\]

is equivalent to

\[
\gamma_{\rm quad}^{\rm eff}
=
\frac{2G}{5c^5},
\]

because

\[
\gamma_{\rm quad}^{\rm eff}
=
P_{\rm eff}\frac{a^5}{27c_s^5}.
\]

### 2.4 4PN tail-transport residual

If the tail bridge is written as

\[
C_{\rm tail}^{\rm toy}
=
\Theta_{\rm tail}
\frac{GM}{2c_s^3}
\gamma_{\rm quad}^{\rm eff},
\]

then after the 2.5PN quadrupole normalization is met, the remaining tail residual is

\[
R_{\rm tail}
=
\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3-1.
\]

The script verifies that

\[
R_{\rm tail}=0
\]

makes the toy tail coefficient match

\[
C_{\rm tail}^{\rm GR}
=
\frac{GM}{2c^3}\gamma_{\rm GR}.
\]

---

## 3. Incidence matrix

The residual rows are

\[
(R_{\rm pole},R_{P2},R_{P4},R_{\rm norm},R_{\rm tail}).
\]

The branch/evaluation columns are

\[
(K,M,B_0,B_2,B_4,Z_0,Z_2,Z_4,N_0,N_2,N_4,\widehat m_0,\mathcal S_{\rm port},a,c_s,\Theta_{\rm tail}).
\]

The incidence matrix printed by the script is a dependency ledger, not a fit instruction. It shows which frozen branch data are used by which residual.

The key interpretation is:

- dependencies from branch data to residuals are allowed;
- dependencies from residuals or target decisions back to branch data are forbidden.

---

## 4. Why the no-refit rule is mandatory

The script computes the Jacobian of the five residuals with respect to the post-hoc knob set

\[
(K,N_2,N_4,N_0,\Theta_{\rm tail}).
\]

It finds

\[
\det
\frac{\partial(R_{\rm pole},R_{P2},R_{P4},R_{\rm norm},R_{\rm tail})}
{\partial(K,N_2,N_4,N_0,\Theta_{\rm tail})}
=
\frac{
\mathcal S_{\rm port}\widehat m_0^{\,2}(B_4+Z_4)
}{
D_0^3
}
\left(\frac{c}{c_s}\right)^3.
\]

Equivalently, as printed in the script’s native sign convention,

\[
-\frac{
\mathcal S_{\rm port}c^3\widehat m_0^{\,2}(B_4+Z_4)
}{
c_s^3(B_0-K+Z_0)^3
}.
\]

Under generic nonzero/stable conditions,

\[
D_0\neq0,
\qquad
B_4+Z_4\neq0,
\qquad
\widehat m_0\mathcal S_{\rm port}\neq0,
\qquad
c/c_s\neq0,
\]

this determinant is nonzero.

Therefore, if post-hoc fitting were allowed, five algebraic knobs could generically tune five residuals. That is the point of this stage: the algebra is powerful enough that target comparison is not trustworthy without a freeze certificate.

---

## 5. Freeze packet

The script creates a deterministic freeze packet and SHA256 certificate.

The packet requires freezing before target evaluation:

```text
parent action and current bookkeeping
gauge-fixing convention
wall/throat action or effective interface action
open-exit impedance boundary protocol
projection/source-map convention
support profile family and number of modes/ports
stability acceptance gates
coefficient extraction formulas
```

Only after these are frozen may the branch evaluate:

```text
one-pole residual
constant-prefactor residuals
universal quadrupole normalization
tail-transport scalar
weak-axisymmetric prefactor slope Xi_1
```

Forbidden after target evaluation:

```text
changing support-cardinality
changing boundary condition class
changing gauge convention
changing port/source normalization convention
dropping dark or unstable branches only after target miss unless this was predeclared
adding compensating modes because a target residual is nonzero
```

---

## 6. Stage result

The audit passes all symbolic and protocol checks.

The important warning is:

\[
\boxed{
\text{The target surface is algebraically tuneable unless the branch is frozen first.}
}
\]

So the correct Volume 2 wording is:

> A branch may be accepted only if the frozen pre-target packet produces stable coefficients satisfying the target residuals. A branch may not be rescued by changing mode support, boundary class, gauge convention, source normalization, or port structure after the residuals are known.

---

## 7. Carry-forward to the next stage

V2-17 can now use this freeze protocol for the weak-axisymmetric branch. The next audit should derive the axisymmetric splitting line

\[
(20,21,22)\sim(1,\tfrac12,-1),
\qquad
b=3a,
\]

and evaluate the transported prefactor slope

\[
\Xi_1=\frac{P_1}{P_0}
\]

without allowing post-target changes to the branch packet.

stage_v2_17_weak_axisymmetric_splitting_derivation.md
# Stage V2-17 — Weak-axisymmetric grouped-\(P_2\) splitting audit

## 0. Purpose and status

This stage audits the weak-axisymmetric grouped-\(P_2\) law that has been carried through the moving-throat program:

\[
(20,21,22)\sim \left(1,\frac12,-1\right),
\qquad
b=3a.
\]

The stage is **exact algebra inside the weak-axisymmetric grouped real \(P_2\) closure**. It does not solve the nonlinear moving-throat PDE. Its job is narrower:

1. derive the angular origin of the lane signature;
2. propagate that signature through the grouped trace/anomaly variables;
3. derive the first-order conservative response slope \(u_2^{(1)}\);
4. derive the outgoing-prefactor slope
   \[
   \Xi_1=\frac{P_1}{P_0};
   \]
5. identify the compensated branch where conservative even splitting is removed and only the prefactor-loading scalar remains.

The accompanying script is:

```text
stage_v2_17_weak_axisymmetric_splitting_sympy_audit.py
```

and its run output is:

```text
stage_v2_17_weak_axisymmetric_splitting_output.txt
```

The script reports:

```text
total_checks: 54
passed_checks: 54
failed_checks: 0
```

---

## 1. Angular source of the weak-axisymmetric signature

Use the normalized real \(l=2\) harmonics on the unit sphere:

\[
Y_{20}=\sqrt{\frac{5}{16\pi}}\,(2z^2-x^2-y^2),
\]

\[
Y_{21c}=\sqrt{\frac{15}{4\pi}}xz,
\qquad
Y_{21s}=\sqrt{\frac{15}{4\pi}}yz,
\]

\[
Y_{22c}=\sqrt{\frac{15}{16\pi}}(x^2-y^2),
\qquad
Y_{22s}=\sqrt{\frac{15}{4\pi}}xy.
\]

The script verifies orthonormality:

\[
\int_{S^2}Y_A Y_B\,d\Omega=\delta_{AB}.
\]

A weak axisymmetric quadrupole perturbation is proportional to \(Y_{20}\), so the lane-splitting matrix is controlled by

\[
M_{AA}^{(20)}=\int_{S^2}Y_A\,Y_{20}\,Y_A\,d\Omega.
\]

Using exact sphere monomial moments, the script obtains

\[
\int Y_{20}Y_{20}Y_{20}\,d\Omega
=
\frac{\sqrt5}{7\sqrt\pi},
\]

\[
\int Y_{21c}Y_{20}Y_{21c}\,d\Omega
=
\int Y_{21s}Y_{20}Y_{21s}\,d\Omega
=
\frac{\sqrt5}{14\sqrt\pi},
\]

\[
\int Y_{22c}Y_{20}Y_{22c}\,d\Omega
=
\int Y_{22s}Y_{20}Y_{22s}\,d\Omega
=
-\frac{\sqrt5}{7\sqrt\pi}.
\]

Therefore

\[
M^{(20)}
=
\kappa_\ast
\operatorname{diag}
\left(
1,\frac12,\frac12,-1,-1
\right),
\qquad
\kappa_\ast=\frac{\sqrt5}{7\sqrt\pi}.
\]

After grouping \(21c,21s\) into the \(21\) lane and \(22c,22s\) into the \(22\) lane, the grouped signature is

\[
\lambda_{20}=1,
\qquad
\lambda_{21}=\frac12,
\qquad
\lambda_{22}=-1.
\]

This is the angular origin of the weak-axisymmetric law.

---

## 2. Grouped trace/anomaly variables

The grouped metric is

\[
G_{\rm grp}=\operatorname{diag}(1,2,2).
\]

For any grouped vector

\[
x=(x_{20},x_{21},x_{22}),
\]

define

\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}5,
\]

\[
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\]

\[
b_x=\frac{x_{21}-x_{22}}2.
\]

The inverse map is

\[
x_{20}=\bar x+4a_x,
\]

\[
x_{21}=\bar x-a_x+b_x,
\]

\[
x_{22}=\bar x-a_x-b_x.
\]

For

\[
\lambda=\left(1,\frac12,-1\right),
\]

the script verifies

\[
\bar\lambda=0,
\qquad
a_\lambda=\frac14,
\qquad
b_\lambda=\frac34.
\]

Thus

\[
\boxed{b_\lambda=3a_\lambda.}
\]

So for any weak coefficient split

\[
x_A=x_0+\epsilon\lambda_Ax_1,
\]

one gets

\[
\bar x=x_0,
\]

\[
a_x=\frac{\epsilon x_1}{4},
\]

\[
b_x=\frac{3\epsilon x_1}{4},
\]

and therefore

\[
\boxed{b_x=3a_x.}
\]

This applies to \(D\)-moments, \(N\)-moments, response moments, prefactor moments, and any other grouped scalar whose first-order lane dependence is sourced by a pure axisymmetric \(l=2\) perturbation.

---

## 3. Conservative response transport

Let the conservative grouped operator moments split as

\[
D_{A0}=D_0+\epsilon\lambda_A D_{01},
\]

\[
D_{A2}=D_2+\epsilon\lambda_A D_{21},
\]

\[
D_{A4}=D_4+\epsilon\lambda_A D_{41}.
\]

The normalized conservative response is

\[
Y_A(\omega)
=
\frac{D_{A0}}{D_{A0}+D_{A2}\omega^2+D_{A4}\omega^4+\cdots}
=
1+u_{2,A}\omega^2+u_{4,A}\omega^4+\cdots.
\]

At zeroth order,

\[
u_2=-\frac{D_2}{D_0}.
\]

At first order,

\[
u_{2,A}
=
u_2+\epsilon\lambda_Au_2^{(1)}+O(\epsilon^2),
\]

with

\[
\boxed{
u_2^{(1)}
=
-\frac{D_{21}+u_2D_{01}}{D_0}.
}
\]

The script verifies the equivalent expanded form

\[
u_2^{(1)}
=
-\frac{D_0D_{21}-D_{01}D_2}{D_0^2}.
\]

Therefore the first conservative even-preserving condition is

\[
u_2^{(1)}=0
\quad\Longleftrightarrow\quad
D_{21}+u_2D_{01}=0.
\]

On the canonical outgoing branch where

\[
u_2=\frac19,
\]

this becomes

\[
\boxed{
D_{21}=-\frac{D_{01}}9.
}
\]

---

## 4. Hidden-even relation at \(O(\omega^4)\)

The next normalized response coefficient is

\[
u_4=\frac{D_2^2-D_0D_4}{D_0^2}.
\]

The canonical compact outgoing inverse-DtN branch has

\[
u_2=\frac19,
\qquad
u_4=\frac4{81}.
\]

Equivalently, in operator moments,

\[
D_2=-\frac{D_0}{9},
\qquad
D_4=-\frac{D_0}{27}.
\]

The one-pole/hidden-even first-order relation is the linearized form of

\[
u_4=4u_2^2.
\]

Since \(u_2=1/9\), the linearization is

\[
u_4^{(1)}=\frac89u_2^{(1)}.
\]

The script derives

\[
u_2^{(1)}
=
-\frac{D_{01}+9D_{21}}{9D_0},
\]

\[
u_4^{(1)}
=
-\frac{5D_{01}+18D_{21}+81D_{41}}{81D_0},
\]

and

\[
u_4^{(1)}-\frac89u_2^{(1)}
=
\frac{D_{01}+18D_{21}-27D_{41}}{27D_0}.
\]

Therefore the hidden-even gate is

\[
D_{01}+18D_{21}-27D_{41}=0,
\]

or

\[
\boxed{
D_{41}=\frac23D_{21}+\frac{D_{01}}{27}.
}
\]

On the even-preserving branch

\[
D_{21}=-\frac{D_{01}}9,
\]

the hidden-even condition reduces to

\[
\boxed{
D_{41}=-\frac{D_{01}}{27}.
}
\]

Thus the compensated conservative branch is

\[
\boxed{
D_{21}=-\frac{D_{01}}9,
\qquad
D_{41}=-\frac{D_{01}}{27}.
}
\]

---

## 5. Outgoing-prefactor transport and \(\Xi_1\)

Let the outgoing transfer numerator split as

\[
N_{A0}=N_0+\epsilon\lambda_A N_{01}.
\]

The static outgoing prefactor is

\[
P_A=\frac{N_{A0}}{D_{A0}}.
\]

Expanding,

\[
P_A=P_0+\epsilon\lambda_A P_1+O(\epsilon^2),
\]

with

\[
P_0=\frac{N_0}{D_0},
\]

\[
\boxed{
P_1=
\frac{D_0N_{01}-D_{01}N_0}{D_0^2}.
}
\]

The dimensionless transported prefactor slope is

\[
\boxed{
\Xi_1=\frac{P_1}{P_0}.
}
\]

Substituting \(P_0=N_0/D_0\),

\[
\boxed{
\Xi_1
=
\frac{N_{01}}{N_0}
-
\frac{D_{01}}{D_0}.
}
\]

This is the same load-mismatch form used later in the monomial/similarity-orbit package.

The prefactor lanes obey the same grouped weak-axisymmetric law:

\[
P_{20}=P_0+\epsilon P_1,
\]

\[
P_{21}=P_0+\frac{\epsilon}{2}P_1,
\]

\[
P_{22}=P_0-\epsilon P_1.
\]

Therefore

\[
\bar P=P_0,
\]

\[
a_P=\frac{\epsilon P_1}{4},
\]

\[
b_P=\frac{3\epsilon P_1}{4},
\]

and

\[
\boxed{
b_P=3a_P.
}
\]

The prefactor-isotropy condition is

\[
\boxed{
\Xi_1=0
\quad\Longleftrightarrow\quad
\frac{N_{01}}{N_0}=\frac{D_{01}}{D_0}.
}
\]

---

## 6. Final V2-17 gate packet

The weak-axisymmetric first-order gate packet is:

\[
\lambda=(1,\tfrac12,-1),
\]

\[
b=3a,
\]

\[
u_2^{(1)}
=
-\frac{D_{21}+u_2D_{01}}{D_0},
\]

\[
\Xi_1
=
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}
-
\frac{D_{01}}{D_0}.
\]

On the canonical compensated branch,

\[
D_{21}=-\frac{D_{01}}9,
\]

\[
D_{41}=-\frac{D_{01}}{27},
\]

and the only remaining first-order weak-axisymmetric normalization defect is

\[
\boxed{
\Xi_1
=
\frac{N_{01}}{N_0}
-
\frac{D_{01}}{D_0}.
}
\]

This is the scalar that future actual-branch calculations must compute.

---

## 7. Interpretation

V2-17 passes.

A pure weak axisymmetric \(l=2\) perturbation does not produce an arbitrary grouped anisotropy. It produces one fixed lane signature:

\[
(20,21,22)\sim(1,\tfrac12,-1).
\]

Therefore every first-order grouped defect produced by that angular mechanism must satisfy

\[
b=3a.
\]

After compensating the conservative even response, the remaining weak-axisymmetric first-order normalization problem is no longer a three-lane problem. It is one scalar loading mismatch:

\[
\Xi_1=\frac{P_1}{P_0}.
\]

That is the continuation target for the quotient/orbit stages.

stage_v2_18_monomial_quotient_similarity_orbit_derivation.md
# Stage V2-18 — Monomial Quotient and Similarity-Orbit Audit

## 0. Purpose

This stage turns the weak-axisymmetric grouped-`P2` defect into a quotient problem.

Earlier stages showed that the grouped weak-axisymmetric branch has a single leading
outgoing-prefactor slope

\[
\Xi_1=\frac{P_1}{P_0}
=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
\]

The purpose of this stage is to prove that this scalar is not an arbitrary
microscopic drift. It is one component of an exact quotient-residual packet

\[
q=(q_{\rm tr},q_{\rm nt},q_\eta)
=
\left(
\delta\ln\mathfrak C_{{\rm tr},*},
\delta\ln\mathfrak C_{{\rm nt},*},
\delta\ln\epsilon_\eta
\right),
\]

where the remaining microscopic drift directions are pure similarity-orbit
motion.

---

## 1. Original microscopic monomial map

Use the microscopic grouped drift vector

\[
\delta\mathbf x
=
(
\delta\ln\lambda_W,\,
\delta\ln c_{\eta U},\,
\delta\ln\gamma,\,
\delta\ln K_U,\,
\delta\ln K_\eta^{\rm eff},\,
\delta\ln K_W^{\rm eff},\,
\delta\ln\mu_W,\,
\delta\ln T_U
)^T.
\]

The exact monomials are

\[
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}},
\]

\[
\mathfrak C_{{\rm nt},*}
=
\left(
\frac{\lambda_W^2\mu_W}
{K_\eta^{\rm eff}(K_W^{\rm eff})^2}
\right)
\left(
\frac{\gamma^2\lambda_W^2\sigma}
{K_UK_W^{\rm eff}}
\right)^{E_*}
\left(
\frac{\pi^2T_U}{L^2K_U}
\right)^{-F_*},
\]

\[
\epsilon_\eta
=
\frac{c_{\eta U}^2}{K_UK_\eta^{\rm eff}}.
\]

Taking logarithmic differentials gives

\[
M_*\delta\mathbf x
=
\begin{pmatrix}
q_{\rm tr}\\
q_{\rm nt}\\
q_\eta
\end{pmatrix}.
\]

In the above variable order,

\[
M_*=
\begin{pmatrix}
0&1+\delta&1+\delta&-(2+\delta+\chi)&0&0&0&1+\chi\\
2+2E&0&2E&F-E&-1&-(2+E)&1&-F\\
0&2&0&-1&-1&0&0&0
\end{pmatrix},
\]

where \(\delta=\delta_{U,*}\), \(\chi=\chi_{0,*}\), \(E=E_*\), and \(F=F_*\).

The script verifies the quoted rank witness:

\[
\det M_*^{(T_U,K_\eta^{\rm eff},\mu_W)}
=
1+\chi.
\]

Therefore, on the physical branch \(1+\chi\ne0\),

\[
{\rm rank}(M_*)=3,
\qquad
\dim\ker M_*=5.
\]

Thus the microscopic drift space splits into a five-dimensional similarity
orbit plus three quotient directions.

---

## 2. Exact normal basis and similarity-orbit split

A convenient exact right inverse is

\[
n_{\rm tr}
=
\left(
0,0,0,0,0,0,
\frac{F_*}{1+\chi_{0,*}},
\frac{1}{1+\chi_{0,*}}
\right)^T,
\]

\[
n_{\rm nt}
=
(0,0,0,0,0,0,1,0)^T,
\]

\[
n_\eta
=
(0,0,0,0,-1,0,-1,0)^T.
\]

Let

\[
N=(n_{\rm tr},n_{\rm nt},n_\eta).
\]

The script verifies

\[
M_*N=I_3.
\]

Therefore any microscopic tangent has the exact split

\[
\delta\mathbf x
=
\delta\mathbf x_{\rm orbit}
+
Nq,
\]

with

\[
M_*\delta\mathbf x_{\rm orbit}=0.
\]

Equivalently,

\[
P_{\rm normal}=NM_*,
\qquad
P_{\rm orbit}=I-NM_*,
\]

are idempotent projectors, and

\[
M_*P_{\rm orbit}=0.
\]

This is the precise algebraic meaning of the similarity orbit: the orbit
directions preserve the three monomials, and only \(q\) moves the physical
quotient.

---

## 3. Normalized Stage-12/13 version

Using normalized variables

\[
(
d\ln G_W,\,
d\ln G_U,\,
d\ln R,\,
d\ln K,\,
d\ln M,\,
d\ln\Omega_U,\,
d\ln\Omega_W,\,
d\ln\delta_U
),
\]

the monomials become

\[
\mathfrak C_{{\rm tr},*}
=
\left(
\frac{RG_U}{\Omega_U^2G_W}
\right)^{1+\delta_U}
\delta_U^{1+\chi_0},
\]

\[
\mathfrak C_{{\rm nt},*}
=
\frac{MG_W^2}{K\Omega_W^4}
\left(
\frac{R^2\sigma}{\Omega_U^2\Omega_W^2}
\right)^{E_*}
\delta_U^{-F_*},
\]

\[
\epsilon_\eta
=
\frac{MG_U^2}{K\Omega_U^2}.
\]

The corresponding drift matrix is

\[
M_{\rm norm}=
\begin{pmatrix}
-(1+\delta)&1+\delta&1+\delta&0&0&-2(1+\delta)&0&1+\chi\\
2&0&2E&-1&1&-2E&-(4+2E)&-F\\
0&2&0&-1&1&-2&0&0
\end{pmatrix}.
\]

The script verifies a nonzero rank witness

\[
\det M_{\rm norm}^{(\delta_U,M,\Omega_W)}
=
2(E+2)(1+\chi).
\]

So the normalized quotient map also has rank three on the physical branch.

The zero-defect equations solve triangularly as

\[
d\ln\delta_U
=
-\frac{1+\delta_U}{1+\chi_0}
\left(
d\ln R+d\ln G_U-d\ln G_W-2d\ln\Omega_U
\right),
\]

\[
d\ln M
=
d\ln K-2d\ln G_U+2d\ln\Omega_U,
\]

\[
d\ln\Omega_W
=
\frac{
d\ln G_W-d\ln G_U+(1-E_*)d\ln\Omega_U
+E_*d\ln R-\frac{F_*}{2}d\ln\delta_U
}
{E_*+2}.
\]

The script substitutes these formulas back into \(M_{\rm norm}v\) and verifies
that all three monomial drifts vanish.

---

## 4. Support-blind extension

If explicit BdG-support primitive drift directions

\[
d\ln\lambda_B,\qquad d\ln\varpi
\]

are appended, the monomial matrix gains two exact zero columns.

The script verifies that these two columns are zero and that the extended
monomial kernel dimension becomes

\[
2+5=7.
\]

This is an important boundary on what the quotient theorem proves. It governs
the weak-axisymmetric normalization defect \(\Xi_1\), but it does not by itself
close every conservative even gate involving BdG moments.

---

## 5. Exact physical defect compiler

The quotient residuals map into the physical first-order defect triplet by

\[
\Theta_1
=
-C_{\rm tr}q_{\rm tr},
\]

\[
\Xi_1
=
A_{\rm tr}q_{\rm tr}+q_{\rm nt},
\]

\[
\mathcal R_1
=
-\Xi_1
-
\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta.
\]

In matrix form,

\[
\begin{pmatrix}
\Theta_1\\
\Xi_1\\
\mathcal R_1
\end{pmatrix}
=
\begin{pmatrix}
-C_{\rm tr}&0&0\\
A_{\rm tr}&1&0\\
-A_{\rm tr}&-1&-\dfrac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}
\end{pmatrix}
\begin{pmatrix}
q_{\rm tr}\\
q_{\rm nt}\\
q_\eta
\end{pmatrix}.
\]

The script verifies

\[
\det D_{\rm defect}
=
\frac{C_{\rm tr}\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}},
\]

so the compiler is invertible when

\[
C_{\rm tr}\ne0,\qquad
\epsilon_{\eta,*}\ne0,\qquad
\epsilon_{\eta,*}\ne1.
\]

The inverse is

\[
q_{\rm tr}
=
-\frac{\Theta_1}{C_{\rm tr}},
\]

\[
q_{\rm nt}
=
\Xi_1+\frac{A_{\rm tr}}{C_{\rm tr}}\Theta_1,
\]

\[
q_\eta
=
-\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}
(\mathcal R_1+\Xi_1).
\]

Thus

\[
\Theta_1=\Xi_1=\mathcal R_1=0
\quad\Longleftrightarrow\quad
q_{\rm tr}=q_{\rm nt}=q_\eta=0.
\]

And since \(q=M_*\delta\mathbf x\),

\[
\Theta_1=\Xi_1=\mathcal R_1=0
\quad\Longleftrightarrow\quad
\delta\mathbf x\in\ker M_*.
\]

That is the exact similarity-orbit zero-defect theorem.

---

## 6. Bridge back to the prefactor slope

The outgoing-prefactor slope is

\[
P_0=\frac{N_0}{D_0},
\]

\[
P_1
=
\frac{N_{01}D_0-N_0D_{01}}{D_0^2}.
\]

Therefore

\[
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}
-
\frac{D_{01}}{D_0}.
\]

The script verifies this identity exactly.

Combining this with the defect compiler gives

\[
\boxed{
\Xi_1
=
\frac{P_1}{P_0}
=
A_{\rm tr}q_{\rm tr}+q_{\rm nt}.
}
\]

So the weak-axisymmetric prefactor obstruction is not arbitrary. It is the
projection of the branch tangent onto the two quotient coordinates
\((q_{\rm tr},q_{\rm nt})\), while \(q_\eta\) enters the selected-branch residual
\(\mathcal R_1\).

---

## 7. Stage verdict

The audit passes.

The quotient variables

\[
(q_{\rm tr},q_{\rm nt},q_\eta)
\]

are exact physical coordinates on the defect quotient. The five remaining
microscopic directions are similarity-orbit directions. The physical
weak-axisymmetric first-order defects are exactly invertible functions of the
quotient residuals.

The immediate carry-forward theorem gate is therefore:

\[
\boxed{
\text{Extract the actual moving-throat branch tangent and test whether }
M_*\delta\mathbf x=0.
}
\]

Equivalently:

\[
\boxed{
\text{Compute }
\Xi_1=\frac{P_1}{P_0}
\text{ and the companion residuals from the quotient packet.}
}
\]

If \(M_*\delta\mathbf x=0\), the weak-axisymmetric prefactor problem is pure
similarity-gauge at first order. If not, the branch leaves the monomial-preserving
orbit, and the exact residual coordinates identify which physical defect
\((\Theta_1,\Xi_1,\mathcal R_1)\) failed.

stage_v2_19_isotropic_full_bundle_target_surface_derivation.md
# Stage V2-19 — Isotropic Full-Bundle Target Surface

## Purpose

This stage consolidates the isotropic moving-throat full-bundle target into a finite pass/fail packet for an actual branch extraction.

Earlier Volume 2 stages verified the ingredients separately:

- the conservative grouped response bridge,
- the outgoing `l=2` fingerprint,
- the 2.5PN / 4PN normalization interface,
- the branch-freeze protocol,
- the weak-axisymmetric `b=3a` line,
- and the quotient / similarity-orbit structure.

V2-19 turns those into one scalar target surface for the isotropic branch.

The stage is deliberately algebraic. It does **not** solve the nonlinear moving-throat PDE. It states what the actual PDE branch must output after the branch data have been frozen.

---

## 1. Isotropic full-bundle definitions

On the isotropic grouped real `P2` branch, all three grouped lanes share the same conservative operator moments and outgoing-transfer moments.

Define

\[
D(\omega)=D_0+D_2\omega^2+D_4\omega^4+O(\omega^6),
\]

with

\[
D_0=K-B_0-Z_0,
\]

\[
D_2=-(M+B_2+Z_2),
\]

\[
D_4=-(B_4+Z_4).
\]

It is useful to abbreviate

\[
A\equiv M+B_2+Z_2,
\]

\[
C\equiv B_4+Z_4.
\]

Then

\[
D_2=-A,
\qquad
D_4=-C.
\]

Here:

- `K,M` are the wall/worldtube static stiffness and inertia data;
- `B_0,B_2,B_4` are the stable BdG support moments;
- `Z_0,Z_2,Z_4` are the conservative Maxwell/mixed moments;
- `N_0,N_2,N_4` are the outgoing-transfer moments.

The outgoing prefactor is defined by

\[
P(\omega)=\frac{D_0N(\omega)}{D(\omega)^2},
\]

where

\[
N(\omega)=N_0+N_2\omega^2+N_4\omega^4+O(\omega^6).
\]

---

## 2. Conservative one-pole surface

The normalized conservative response is

\[
Y(\omega)=\frac{D_0}{D(\omega)}.
\]

Expanding gives

\[
Y(\omega)=1+u_2\omega^2+u_4\omega^4+O(\omega^6),
\]

with

\[
u_2=\frac{A}{D_0},
\]

\[
u_4=\frac{A^2+D_0C}{D_0^2}.
\]

The one-pole condition is

\[
u_4=4u_2^2.
\]

Substitution gives

\[
\frac{A^2+D_0C}{D_0^2}=4\frac{A^2}{D_0^2},
\]

so the exact one-pole surface is

\[
\boxed{
D_0C=3A^2.
}
\]

Equivalently,

\[
\boxed{
D_0=\frac{3A^2}{C}.
}
\]

Since

\[
D_0=K-B_0-Z_0,
\]

the required wall stiffness on the one-pole surface is

\[
\boxed{
K=B_0+Z_0+\frac{3(M+B_2+Z_2)^2}{B_4+Z_4}.
}
\]

This is a branch-output condition. It is not a license to refit `K` after target evaluation.

---

## 3. Universal quadrupole normalization

The universal 2.5PN / 4PN quadrupole target is

\[
\widehat m_0^2\mathcal S_{\rm port}P_0
=
\frac{54Gc_s^5}{5a^5c^5}.
\]

Define

\[
T_{\rm GR}
\equiv
\frac{54Gc_s^5}{5a^5c^5}.
\]

Since

\[
P_0=\frac{N_0}{D_0},
\]

the normalization gate is

\[
\boxed{
\widehat m_0^2\mathcal S_{\rm port}\frac{N_0}{D_0}=T_{\rm GR}.
}
\]

Equivalently, as a polynomial residual,

\[
\boxed{
R_{\rm norm}
=
\widehat m_0^2\mathcal S_{\rm port}N_0
-T_{\rm GR}D_0=0.
}
\]

Solving for `N_0` gives

\[
\boxed{
N_0=\frac{T_{\rm GR}D_0}{\widehat m_0^2\mathcal S_{\rm port}}.
}
\]

On the one-pole surface this becomes

\[
\boxed{
N_0
=
\frac{3T_{\rm GR}A^2}{\widehat m_0^2\mathcal S_{\rm port}C}.
}
\]

Writing out `T_GR`,

\[
\boxed{
N_0
=
\frac{162Gc_s^5(M+B_2+Z_2)^2}
{5\mathcal S_{\rm port}a^5c^5\widehat m_0^2(B_4+Z_4)}.
}
\]

The SymPy audit also verifies that this is exactly equivalent to

\[
\gamma_{\rm quad}^{\rm eff}
=
\frac{2G}{5c^5},
\]

where

\[
\gamma_{\rm quad}^{\rm eff}
=
\widehat m_0^2\mathcal S_{\rm port}P_0\frac{a^5}{27c_s^5}.
\]

---

## 4. Constant-prefactor outgoing branch

The prefactor expansion is

\[
P(\omega)=P_0+P_2\omega^2+P_4\omega^4+O(\omega^6).
\]

The exact coefficients are

\[
P_0=\frac{N_0}{D_0},
\]

\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]

\[
P_4=
\frac{
D_0^2N_4
-2D_0(D_2N_2+D_4N_0)
+3D_2^2N_0
}{D_0^3}.
\]

The constant-prefactor branch requires

\[
P_2=0,
\qquad
P_4=0.
\]

The exact moment constraints are therefore

\[
\boxed{
N_2=\frac{2D_2N_0}{D_0}.
}
\]

Since `D_2=-A`,

\[
\boxed{
N_2=-\frac{2AN_0}{D_0}.
}
\]

The next condition is

\[
\boxed{
N_4=\frac{N_0(D_2^2+2D_0D_4)}{D_0^2}.
}
\]

Since `D_2=-A` and `D_4=-C`,

\[
\boxed{
N_4=\frac{N_0(A^2-2D_0C)}{D_0^2}.
}
\]

On the one-pole surface `D_0C=3A^2`, this simplifies to

\[
\boxed{
N_4=-\frac{5A^2N_0}{D_0^2}.
}
\]

Using the normalization surface as well, the two outgoing-transfer target moments become

\[
\boxed{
N_2
=
-\frac{2AT_{\rm GR}}{\widehat m_0^2\mathcal S_{\rm port}},
}
\]

and

\[
\boxed{
N_4
=
-\frac{5CT_{\rm GR}}{3\widehat m_0^2\mathcal S_{\rm port}}.
}
\]

Equivalently, after writing out `T_GR`,

\[
\boxed{
N_2
=
-\frac{108Gc_s^5(M+B_2+Z_2)}
{5\mathcal S_{\rm port}a^5c^5\widehat m_0^2},
}
\]

\[
\boxed{
N_4
=
-\frac{18Gc_s^5(B_4+Z_4)}
{\mathcal S_{\rm port}a^5c^5\widehat m_0^2}.
}
\]

---

## 5. Optional tail-transport gate

If the 4PN tail transport is not already derived, carry the scalar gate

\[
\boxed{
R_{\rm tail}
=
\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3-1=0.
}
\]

Thus

\[
\boxed{
\Theta_{\rm tail}=\left(\frac{c_s}{c}\right)^3.
}
\]

On the branch `c_s=c`, this reduces to

\[
\boxed{
\Theta_{\rm tail}=1.
}
\]

This is separate from the grouped quadrupole normalization. It is not a new `P2` normalization slot.

---

## 6. Final V2-19 target packet

Let

\[
A=M+B_2+Z_2,
\qquad
C=B_4+Z_4.
\]

The isotropic full-bundle branch passes V2-19 only if the frozen branch output satisfies:

\[
\boxed{
D_0=K-B_0-Z_0=\frac{3A^2}{C},
}
\]

\[
\boxed{
\widehat m_0^2\mathcal S_{\rm port}\frac{N_0}{D_0}
=
\frac{54Gc_s^5}{5a^5c^5},
}
\]

\[
\boxed{
N_2=\frac{2D_2N_0}{D_0},
}
\]

\[
\boxed{
N_4=\frac{N_0(D_2^2+2D_0D_4)}{D_0^2},
}
\]

and, if the tail transport is not otherwise closed,

\[
\boxed{
\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3=1.
}
\]

---

## 7. Feasibility and stability implications

The one-pole surface gives

\[
D_0=\frac{3A^2}{C}.
\]

For a stable branch we require

\[
D_0>0.
\]

Therefore, for a nondegenerate branch with `A != 0`, we need

\[
\boxed{
C=B_4+Z_4>0.
}
\]

If `C <= 0`, the one-pole surface is incompatible with `D_0>0`, except for degenerate cases that would also collapse the target packet.

The normalized outgoing transfer is

\[
N_0=\frac{T_{\rm GR}D_0}{\widehat m_0^2\mathcal S_{\rm port}}.
\]

So if

\[
G>0,
\quad c_s>0,
\quad a>0,
\quad c>0,
\quad \widehat m_0^2>0,
\quad \mathcal S_{\rm port}>0,
\quad D_0>0,
\]

then the target branch has

\[
N_0>0.
\]

Thus it is not a dark port.

---

## 8. Local codimension check

The SymPy audit uses the residual vector

\[
R=
\begin{pmatrix}
R_{\rm pole}\\
R_{\rm norm}\\
R_{P2}\\
R_{P4}\\
R_{\rm tail}
\end{pmatrix}
\]

with output variables

\[
(K,N_0,N_2,N_4,\Theta_{\rm tail}).
\]

The Jacobian determinant is

\[
\boxed{
\det J
=
-\mathcal S_{\rm port}\,m_0^2
\left(\frac{c}{c_s}\right)^3
(B_4+Z_4)
(B_0-K+Z_0)^3.
}
\]

Since

\[
B_0-K+Z_0=-D_0,
\]

this determinant is nonzero whenever

\[
\mathcal S_{\rm port}\ne0,
\quad
m_0\ne0,
\quad
c_s\ne0,
\quad
B_4+Z_4\ne0,
\quad
D_0\ne0.
\]

On the target surface, the determinant becomes

\[
\boxed{
\det J\big|_{\rm surf}
=
\frac{
27\mathcal S_{\rm port}m_0^2c^3A^6
}{
c_s^3C^2
}.
}
\]

So the target packet is locally codimension five in those algebraic output slots. This is exactly why the branch-freeze rule matters: without freezing, the same slots could be tuned to hit the residuals.

---

## 9. SymPy audit result

The script performed 18 symbolic checks. All passed.

```text
checks_total: 18
checks_passed: 18
checks_failed: 0
```

The most important output formulas are:

```text
D0_surface = 3*(B2 + M + Z2)^2/(B4 + Z4)

K_surface = B0 + Z0 + 3*(B2 + M + Z2)^2/(B4 + Z4)

P0_surface = 54*G*c_s^5/(5*S_port*a_th^5*c^5*mhat0^2)

N0_surface = 162*G*c_s^5*(B2 + M + Z2)^2/
             (5*S_port*a_th^5*c^5*mhat0^2*(B4 + Z4))

N2_surface = -108*G*c_s^5*(B2 + M + Z2)/
             (5*S_port*a_th^5*c^5*mhat0^2)

N4_surface = -18*G*c_s^5*(B4 + Z4)/
             (S_port*a_th^5*c^5*mhat0^2)
```

---

## 10. Carry-forward statement

V2-19 gives the isotropic target sheet for the actual moving-throat branch:

> A frozen isotropic branch must output stable `D0`, non-dark `N0`, one-pole conservative moments, constant-prefactor outgoing moments, and the universal quadrupole normalization using one common `G`. If the tail scalar has not already been derived, the same branch must also satisfy the separate transport gate `Theta_tail(c/c_s)^3=1`.

This is the exact target packet that V2-20 should feed into a weak-form / numerical extraction plan.

stage_v2_20_weak_form_branch_extraction_derivation.md
# Stage V2-20 — Weak-form and numerical branch-extraction preparation

## Status

This stage is a **pre-solver theorem gate**, not a nonlinear branch solution.

The previous stages collapsed the isotropic grouped-\(P_2\) finish line to one scalar normalization test,

\[
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
=
\frac{54Gc_s^5}{5a^5c^5},
\]

together with conservative one-pole, constant-prefactor, tail-transport, stability, and freeze-order gates.

V2-20 answers the next practical question:

> When a future moving-throat solver outputs a branch, exactly what weak-form data must be extracted, and how do we turn those data into \(K,M,B_n,Z_n,N_n,\widehat m_0,\Theta_{\rm tail}\) without refitting the target?

## Architectural patch retained from V2-04

The throat is now treated as an **open finite-radius conduit**.

The exit data are:

\[
R(L)>0,
\]

not

\[
R(L)=0.
\]

The AC support modes may see an effective Neumann condition through an open-exit impedance mismatch, but this is not a hard geometric cap. The zero/DC superfluid current is allowed to exit and is recorded through the open-system leakage bookkeeping rather than through the finite-support AC mode equation.

In this stage, the far end is represented by a Robin/impedance condition,

\[
T_w q_w(L,\omega)+Y_L(\omega)q(L,\omega)=0.
\]

The Neumann support limit is

\[
Y_L(\omega)\to0.
\]

For a low-impedance open expansion, this is the support-coordinate version of the organ-pipe reflection result from V2-04.

## 1. Weak form for one wall harmonic

For a real harmonic sector \((l,m)\), write

\[
V_l(w)=K_\eta(w)+l(l+1)T_\Omega(w).
\]

The densitized wall equation is

\[
\mu_\eta q_{tt}
-\partial_w(T_w q_w)
+V_l q
=
S_l.
\]

Using \(q(w,t)=q(w,\omega)e^{-i\omega t}\), this becomes

\[
-\omega^2\mu_\eta q
-\partial_w(T_w q_w)
+V_lq
=
S_l.
\]

For a test function \(\varphi\) satisfying the mouth constraint, the weak form is

\[
\int_0^L
\left[
T_w q_w\varphi_w
+
(V_l-\omega^2\mu_\eta)q\varphi
\right]dw
+
Y_L(\omega)q(L)\varphi(L)
=
\int_0^L S_l\varphi\,dw.
\]

The SymPy script verifies the integration-by-parts identity

\[
Tq_w\varphi_w-\partial_w(Tq_w\varphi)
=
-\partial_w(Tq_w)\varphi.
\]

It also verifies the open-exit sign convention: the endpoint variation is

\[
\left[T_w q_w(L)+Y_Lq(L)\right]\varphi(L),
\]

so free exit variations impose

\[
T_w q_w(L)+Y_Lq(L)=0.
\]

## 2. Galerkin extraction matrices

Choose a frozen branch basis \(\{\chi_i(w)\}\) satisfying the mouth condition. In a numerical solver this basis may be finite elements, splines, spectral functions, or computed eigenfunctions. The extracted matrices are

\[
M_{ij}^{(l)}
=
\int_0^L\mu_\eta(w)\chi_i(w)\chi_j(w)\,dw,
\]

\[
K_{ij}^{(l)}
=
\int_0^L
\left[
T_w(w)\chi_i'(w)\chi_j'(w)
+
V_l(w)\chi_i(w)\chi_j(w)
\right]dw
+
Y_L(0)\chi_i(L)\chi_j(L).
\]

The script checks a two-function Dirichlet-mouth open-exit toy basis,

\[
\chi_1=w/L,\qquad \chi_2=(w/L)^2,
\]

and derives

\[
M
=
\begin{pmatrix}
L\mu/3 & L\mu/4\\
L\mu/4 & L\mu/5
\end{pmatrix},
\qquad
\det M=\frac{L^2\mu^2}{240}>0.
\]

The stiffness matrix is symmetric and is positive in the tested positive-coefficient sample.

For the AC support limit \(Y_L\to0\), the same script checks the D/N support seed,

\[
q(w)=\sin\left(\frac{\pi w}{2L}\right),
\]

with

\[
q(0)=0,\qquad q_w(L)=0.
\]

This keeps the support ladder while avoiding a capped geometry.

## 3. Wall/BdG/Maxwell extraction slots

For each grouped lane

\[
A\in\{20,21,22\},
\]

the solver must output or permit extraction of the following frozen primitives.

### 3.1 Wall data

For the selected wall/support coordinate,

\[
M_A,
\qquad
K_A.
\]

These are obtained from the weak-form matrices and the selected normalized wall eigenvector or response coordinate.

### 3.2 Stable BdG support data

For positive-energy stable BdG modes,

\[
c_{A\alpha},
\qquad
\varpi_{A\alpha}>0.
\]

The extracted moments are

\[
B_{A0}=\sum_\alpha\frac{c_{A\alpha}^2}{\varpi_{A\alpha}^2},
\]

\[
B_{A2}=\sum_\alpha\frac{c_{A\alpha}^2}{\varpi_{A\alpha}^4},
\]

\[
B_{A4}=\sum_\alpha\frac{c_{A\alpha}^2}{\varpi_{A\alpha}^6}.
\]

The script verifies the low-frequency Schur-complement expansion

\[
K-M\omega^2-\frac{g^2}{\varpi^2-\omega^2}
=
\left(K-\frac{g^2}{\varpi^2}\right)
-
\left(M+\frac{g^2}{\varpi^4}\right)\omega^2
-
\frac{g^2}{\varpi^6}\omega^4
+\cdots.
\]

It also verifies the two-mode Hankel/Stieltjes positivity identity

\[
B_0B_4-B_2^2
=
\frac{
w_1w_2(\lambda_1-\lambda_2)^2
}{
\lambda_1^3\lambda_2^3
}
\ge0.
\]

### 3.3 Conservative Maxwell/mixed data

For each mixed port \(r\), the numerical branch must output

\[
\Omega_{U,Ar},\quad
\Omega_{W,Ar},\quad
R_{Ar},\quad
g_{U,Ar},\quad
g_{W,Ar}.
\]

Define

\[
\Delta_{Ar}=\Omega_{U,Ar}^2\Omega_{W,Ar}^2-R_{Ar}^2,
\]

\[
S_{Ar}=\Omega_{U,Ar}^2+\Omega_{W,Ar}^2,
\]

\[
Q_{Ar}
=
g_{U,Ar}^2\Omega_{W,Ar}^2
+
2g_{U,Ar}g_{W,Ar}R_{Ar}
+
g_{W,Ar}^2\Omega_{U,Ar}^2,
\]

\[
H_{Ar}=g_{U,Ar}^2+g_{W,Ar}^2.
\]

Then

\[
Z_{A0}^{(r)}=\frac{Q_{Ar}}{\Delta_{Ar}},
\]

\[
Z_{A2}^{(r)}
=
\frac{Q_{Ar}S_{Ar}-H_{Ar}\Delta_{Ar}}{\Delta_{Ar}^2},
\]

\[
Z_{A4}^{(r)}
=
\frac{
Q_{Ar}(S_{Ar}^2-\Delta_{Ar})-S_{Ar}H_{Ar}\Delta_{Ar}
}{
\Delta_{Ar}^3
}.
\]

The full conservative mixed data are

\[
Z_{An}=\sum_r Z_{An}^{(r)}.
\]

The internal stability gate is

\[
\Delta_{Ar}>0
\]

for every active mixed port.

### 3.4 Outgoing transfer data

For each outgoing mixed port,

\[
P_{Ar}
=
\Omega_{U,Ar}^2g_{W,Ar}
+
R_{Ar}g_{U,Ar}.
\]

The outgoing-transfer moments are

\[
N_{A0}^{(r)}
=
\frac{P_{Ar}^2}{\Delta_{Ar}^2},
\]

\[
N_{A2}^{(r)}
=
\frac{
2P_{Ar}(P_{Ar}S_{Ar}-\Delta_{Ar}g_{W,Ar})
}{
\Delta_{Ar}^3
},
\]

\[
N_{A4}^{(r)}
=
\frac{
\Delta_{Ar}^2g_{W,Ar}^2
-2\Delta_{Ar}P_{Ar}^2
-4\Delta_{Ar}P_{Ar}S_{Ar}g_{W,Ar}
+3P_{Ar}^2S_{Ar}^2
}{
\Delta_{Ar}^4
}.
\]

The full transfer data are

\[
N_{An}=\sum_rN_{An}^{(r)}.
\]

The dark-port failure is

\[
P_{Ar}=0
\]

for all active ports, which gives

\[
N_{A0}=0.
\]

## 4. Total grouped operator and response extraction

The conservative lane operator is

\[
D_A(\omega)
=
D_{A0}+D_{A2}\omega^2+D_{A4}\omega^4+O(\omega^6),
\]

with

\[
D_{A0}=K_A-B_{A0}-Z_{A0},
\]

\[
D_{A2}=-(M_A+B_{A2}+Z_{A2}),
\]

\[
D_{A4}=-(B_{A4}+Z_{A4}).
\]

Then

\[
u_{2,A}=-\frac{D_{A2}}{D_{A0}},
\]

\[
u_{4,A}=\frac{D_{A2}^2-D_{A0}D_{A4}}{D_{A0}^2}.
\]

The outgoing prefactor is

\[
P_{0,A}=\frac{N_{A0}}{D_{A0}},
\]

\[
P_{2,A}=\frac{D_{A0}N_{A2}-2D_{A2}N_{A0}}{D_{A0}^2},
\]

\[
P_{4,A}
=
\frac{
D_{A0}^2N_{A4}
-2D_{A0}(D_{A2}N_{A2}+D_{A4}N_{A0})
+3D_{A2}^2N_{A0}
}{D_{A0}^3}.
\]

The script verifies these formulas algebraically and confirms that the constant-prefactor branch is

\[
N_2=\frac{2D_2N_0}{D_0},
\]

\[
N_4=\frac{N_0(D_2^2+2D_0D_4)}{D_0^2}.
\]

## 5. Grouped projectors and anisotropy extraction

The grouped metric is

\[
G_{\rm grp}=\mathrm{diag}(1,2,2).
\]

For any grouped triple \(x=(x_{20},x_{21},x_{22})\),

\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}5,
\]

\[
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\]

\[
b_x=\frac{x_{21}-x_{22}}2.
\]

The inverse map is

\[
x_{20}=\bar x+4a_x,
\]

\[
x_{21}=\bar x-a_x+b_x,
\]

\[
x_{22}=\bar x-a_x-b_x.
\]

The script verifies projector completeness and idempotence.

For a pure weak-axisymmetric perturbation,

\[
\lambda=(1,\tfrac12,-1),
\]

the script verifies

\[
\bar x=x_0,
\qquad
b_x=3a_x.
\]

So future branch output must be classified as:

- isotropic if \(a=b=0\),
- pure weak-axisymmetric if \(b=3a\),
- more general anisotropy otherwise.

## 6. Final target residual packet

The numerical branch extraction should report these residuals:

\[
R_{\rm pole}
=
D_0(B_4+Z_4)-3(M+B_2+Z_2)^2,
\]

\[
R_{\rm norm}
=
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
-
\frac{54Gc_s^5}{5a^5c^5},
\]

\[
R_{P2}=P_2,
\qquad
R_{P4}=P_4,
\]

\[
R_{\rm tail}
=
\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3-1.
\]

The script verifies that

\[
R_{\rm norm}=0
\]

is equivalent to

\[
\gamma_{\rm quad}^{\rm eff}
=
\frac{2G}{5c^5}.
\]

## 7. Mandatory stability gates

A branch cannot be accepted unless the following are true before target evaluation:

\[
\mu_\eta>0,
\qquad
T_w>0,
\]

\[
\Delta_{Ar}>0
\quad
\text{for every active mixed port},
\]

\[
D_0>0,
\]

\[
B_0B_4-B_2^2\ge0
\]

for the stable BdG moment sequence, and

\[
B_4+Z_4>0
\]

on a stable nondegenerate one-pole branch.

The branch also fails for the quadrupole route if the active outgoing ports are all dark:

\[
N_0=0.
\]

## 8. Freeze manifest

The script emits a freeze manifest requiring:

1. open finite-radius exit geometry;
2. no hard cap;
3. mouth boundary class;
4. exit impedance boundary class;
5. DC leakage policy;
6. weak-form convention;
7. densitized measure convention;
8. extraction slot list;
9. stability gates;
10. target residual packet;
11. no-refit rule.

The manifest hash produced by the script is

```text
067c626bb61456fad945f5b3f7fa4d10c19e38f9083bd835896e1d052261e390
```

This hash is not a physical invariant. It is a bookkeeping certificate for the V2-20 extraction schema. A future branch run should include its own branch-packet hash before target residuals are evaluated.

## 9. Solver protocol produced by this stage

A future numerical PDE branch test should run in this order:

1. **Freeze the packet**: parent action, gauge convention, open-exit boundary class, source projection, wall action, support basis, port list, and extraction formulas.
2. **Solve the stationary branch** with \(R(L)>0\).
3. **Linearize** the GNLS/BdG, Maxwell/mixed, and wall/interface equations.
4. **Select positive-energy stable modes** and record any excluded negative-Krein modes.
5. **Extract wall matrices** \(K_A,M_A\).
6. **Extract BdG moments** \(B_{A0},B_{A2},B_{A4}\).
7. **Extract mixed-sector moments** \(Z_{A0},Z_{A2},Z_{A4}\).
8. **Extract outgoing transfer moments** \(N_{A0},N_{A2},N_{A4}\).
9. **Project grouped data** into \(\bar x,a_x,b_x\).
10. **Compute residuals** \(R_{\rm pole},R_{\rm norm},R_{P2},R_{P4},R_{\rm tail}\).
11. **Report pass/fail** without changing the frozen branch packet.

## 10. SymPy result

The script reports:

```text
checks_total: 28
checks_passed: 28
checks_failed: 0
```

The algebraic scaffold passes.

The important limitation is unchanged:

> V2-20 prepares the branch extraction theorem gate. It does not prove that the actual moving-throat PDE branch lands on the target surface.

stage_v2_21_branch_extraction_fixture_derivation.md
# Stage V2-21 — Branch-Extraction Fixture and Observable Packet

## 0. Purpose

This stage turns the Volume 2 weak-form scaffold into an executable branch-extraction fixture.

The fixture does **not** solve the nonlinear moving-throat PDE. It does the next necessary thing: it accepts a frozen branch packet, extracts the grouped real `P2` observables, enforces open-throat and stability gates, and evaluates the target residuals without changing branch data after the residuals are known.

The implementation file is:

```text
stage_v2_21_branch_extraction_fixture.py
```

The fixture writes a machine-readable extraction packet:

```text
stage_v2_21_branch_extraction_packet.json
```

and includes a sample manifest:

```text
stage_v2_21_sample_branch_manifest.json
```

---

## 1. Input manifest

The branch manifest has the form

```json
{
  "schema": "stage_v2_21_branch_extraction_fixture/v1",
  "branches": [
    {
      "name": "branch name",
      "geometry": {
        "R_exit": 0.35,
        "boundary_class": "open_impedance",
        "Y_L_limit": 0.0
      },
      "target": {
        "constants": {
          "G": 1.0,
          "c_s": 1.0,
          "c": 1.0,
          "a": 1.0,
          "mhat0": 1.0,
          "S_port": 1.0,
          "theta_tail": 1.0
        }
      },
      "lanes": {
        "20": { "...": "lane data" },
        "21": { "...": "lane data" },
        "22": { "...": "lane data" }
      }
    }
  ]
}
```

Each grouped lane may be supplied in either of two ways.

### 1.1 Primitive reduced-mode input

```json
{
  "K": 1.6,
  "M": 0.9,
  "bdg_modes": [
    {"coupling": 0.42, "varpi": 2.7}
  ],
  "mixed_ports": [
    {"Omega_U": 3.0, "Omega_W": 4.0, "R": 0.65, "g_U": 0.28, "g_W": 0.52}
  ]
}
```

### 1.2 Direct coefficient input

```json
{
  "direct_coefficients": {
    "K": 1.0,
    "M": 1.0,
    "B0": 0.0,
    "B2": 0.0,
    "B4": 3.0,
    "Z0": 0.0,
    "Z2": 0.0,
    "Z4": 0.0,
    "N0": 10.8,
    "N2": -21.6,
    "N4": -54.0
  }
}
```

The direct-coefficient path is included only as an algebraic self-test and for importing future pre-extracted numerical coefficients. The primitive path is the intended route for real branch data.

---

## 2. Open-throat gate

The fixture enforces the open-organ-pipe patch by requiring

\[
R_{\rm exit}>0,
\]

and

\[
\texttt{boundary\_class}=\texttt{open\_impedance}.
\]

It rejects or fails hard-cap geometry:

\[
R(L)=0
\]

as a branch-realization condition for the support problem. In this fixture, the support-mode Neumann-like condition is represented by the open impedance limit

\[
T_w q_w(L,\omega)+Y_L(\omega)q(L,\omega)=0,
\qquad
Y_L\to0
\]

for AC support modes, while DC leakage is allowed through the finite exit.

---

## 3. Lane-level extraction formulas

For each grouped lane

\[
A\in\{20,21,22\},
\]

the fixture computes BdG moments, Maxwell/mixed moments, outgoing-transfer moments, conservative operator moments, response moments, and prefactor moments.

### 3.1 BdG support moments

For stable positive-energy support coordinates with couplings \(c_\alpha\) and frequencies \(\varpi_\alpha>0\),

\[
B_0=\sum_\alpha \frac{c_\alpha^2}{\varpi_\alpha^2},
\]

\[
B_2=\sum_\alpha \frac{c_\alpha^2}{\varpi_\alpha^4},
\]

\[
B_4=\sum_\alpha \frac{c_\alpha^2}{\varpi_\alpha^6}.
\]

These are the Schur-complement support moments from the stable BdG closure.

### 3.2 Conservative Maxwell/mixed moments

For a port with

\[
\Omega_U,
\qquad
\Omega_W,
\qquad
R,
\qquad
g_U,
\qquad
g_W,
\]

define

\[
\Delta=\Omega_U^2\Omega_W^2-R^2,
\]

\[
S=\Omega_U^2+\Omega_W^2,
\]

\[
Q=g_U^2\Omega_W^2+2g_Ug_WR+g_W^2\Omega_U^2,
\]

\[
H=g_U^2+g_W^2.
\]

Then

\[
Z_0=\frac{Q}{\Delta},
\]

\[
Z_2=\frac{QS-H\Delta}{\Delta^2},
\]

\[
Z_4=\frac{Q(S^2-\Delta)-SH\Delta}{\Delta^3}.
\]

For multiple mixed ports, the fixture sums these contributions over ports.

### 3.3 Outgoing-transfer moments

Define the outgoing-port numerator

\[
P=\Omega_U^2g_W+Rg_U.
\]

Then

\[
N_0=\frac{P^2}{\Delta^2},
\]

\[
N_2=\frac{2P(PS-\Delta g_W)}{\Delta^3},
\]

\[
N_4=
\frac{
\Delta^2 g_W^2
-2\Delta P^2
-4\Delta PSg_W
+3P^2S^2
}{\Delta^4}.
\]

Again, for multiple ports, the fixture sums over ports.

### 3.4 Conservative wall operator

The lane operator moments are

\[
D_0=K-B_0-Z_0,
\]

\[
D_2=-(M+B_2+Z_2),
\]

\[
D_4=-(B_4+Z_4).
\]

### 3.5 Normalized response

The normalized conservative response is

\[
Y(\omega)=\frac{D_0}{D_0+D_2\omega^2+D_4\omega^4+O(\omega^6)}.
\]

The fixture verifies symbolically that

\[
u_2=-\frac{D_2}{D_0},
\]

\[
u_4=\frac{D_2^2-D_0D_4}{D_0^2}.
\]

### 3.6 Outgoing prefactor

The outgoing prefactor is

\[
\mathrm{Pref}(\omega)
=
\frac{D_0\left(N_0+N_2\omega^2+N_4\omega^4\right)}{
\left(D_0+D_2\omega^2+D_4\omega^4\right)^2
}.
\]

The fixture verifies symbolically that

\[
P_0=\frac{N_0}{D_0},
\]

\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]

\[
P_4=
\frac{
D_0^2N_4
-2D_0(D_2N_2+D_4N_0)
+3D_2^2N_0
}{D_0^3}.
\]

---

## 4. Grouped real `P2` decomposition

For any grouped triple

\[
x=(x_{20},x_{21},x_{22}),
\]

the fixture computes

\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}{5},
\]

\[
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\]

\[
b_x=\frac{x_{21}-x_{22}}{2}.
\]

The inverse map is verified symbolically:

\[
x_{20}=\bar x+4a_x,
\]

\[
x_{21}=\bar x-a_x+b_x,
\]

\[
x_{22}=\bar x-a_x-b_x.
\]

The anisotropy norm is

\[
A_x^2=4a_x^2+\frac45b_x^2.
\]

The axisymmetric weak-splitting diagnostic is also reported:

\[
b_x-3a_x.
\]

---

## 5. Target residual packet

The fixture evaluates the isotropic trace residuals

\[
R_{\rm pole}
=
D_0(B_4+Z_4)-3(M+B_2+Z_2)^2,
\]

\[
R_{\rm norm}
=
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
-
\frac{54Gc_s^5}{5a^5c^5},
\]

\[
R_{P2}=P_2,
\]

\[
R_{P4}=P_4,
\]

\[
R_{\rm tail}=\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3-1.
\]

It also reports the equivalent quadrupole coefficient residual

\[
\gamma_{\rm eff}-\gamma_{\rm GR}
=
\widehat m_0^{\,2}\mathcal S_{\rm port}P_0\frac{a^5}{27c_s^5}
-
\frac{2G}{5c^5}.
\]

---

## 6. Stability gates

The fixture currently enforces the reduced stability gates:

\[
D_0>0,
\]

\[
B_4+Z_4>0,
\]

\[
N_0\ne0,
\]

\[
M>0,
\]

and for every primitive mixed port,

\[
\Delta=\Omega_U^2\Omega_W^2-R^2>0.
\]

These are reduced checks. A future PDE run should also supply the full positivity certificate for the BdG signature and full wall/mixed Hamiltonian.

---

## 7. Built-in fixture branches

The script includes two built-in branches.

### 7.1 Calibrated direct-coefficient branch

This is an algebraic self-test, not a physical primitive realization.

It uses normalized constants

\[
G=c=c_s=a=\widehat m_0=\mathcal S_{\rm port}=\Theta_{\rm tail}=1,
\]

so

\[
P_0^{\rm target}=\frac{54}{5}=10.8.
\]

The direct coefficients are chosen so that

\[
D_0=1,
\qquad
M+B_2+Z_2=1,
\qquad
B_4+Z_4=3,
\]

which gives

\[
D_0(B_4+Z_4)=3(M+B_2+Z_2)^2.
\]

Then

\[
N_0=10.8,
\]

and the constant-prefactor branch fixes

\[
N_2=\frac{2D_2N_0}{D_0}=-21.6,
\]

\[
N_4=\frac{N_0(D_2^2+2D_0D_4)}{D_0^2}=-54.
\]

The output shows this branch passing all gates up to floating-point roundoff.

### 7.2 Primitive open-throat demo branch

This branch uses actual reduced primitive data: two BdG support modes and one Maxwell/mixed port. It is stable and open-throat, but it is not tuned to the target surface.

The output is intentionally a falsifying packet:

- open gate passes,
- stability gate passes,
- grouped anisotropy vanishes because the lanes are isotropic,
- but the one-pole, normalization, and constant-prefactor residuals do not vanish.

This demonstrates the role of the fixture: it does not rescue a branch. It reports the residuals.

---

## 8. Output from this run

The script reported:

```text
symbolic_checks: 13/13 passed
```

For the calibrated branch:

```text
target_packet_pass: True
D0_bar: 1
N0_bar: 10.8
P0_bar: 10.8
P0_target: 10.8
R_pole: 0
R_norm: 0
R_P2: 0
R_P4: 1.4210854715202e-14
R_tail: 0
```

For the primitive open-throat demo branch:

```text
target_packet_pass: False
D0_bar: 1.546870256628791
N0_bar: 0.001146719334105218
P0_bar: 0.0007413157821033734
P0_target: 10.8
R_pole: -2.459877868885326
R_norm: -10.7992586842179
R_P2: 0.0009676829082644174
R_P4: 0.0008900366908790024
R_tail: 0
```

---

## 9. Interpretation

V2-21 changes the workflow from handwritten coefficient checks to a reproducible branch-extraction protocol.

A future PDE or numerical branch run should now do this:

1. freeze the branch manifest before target evaluation,
2. extract primitive wall, BdG, Maxwell/mixed, and outgoing-port data,
3. write them into the manifest schema,
4. run the fixture,
5. read the residual packet.

The important discipline is that the fixture has no rescue logic. It does not change support cardinality, boundary class, gauge convention, or port normalization after the residuals are known.

---

## 10. Carry-forward status

The fixture is ready for mock data and reduced numerical branch data.

The next technical stage should be either:

1. **V2-22A:** build a profile-to-coefficient adapter, where actual sampled profiles \(\chi_\eta(s),\phi_\alpha(s),u_r(s),w_r(s)\) are integrated numerically to produce \(B_n,Z_n,N_n\); or
2. **V2-22B:** build a weak-axisymmetric manifest adapter that adds first-order lane slopes and automatically reports \(\Xi_1=P_1/P_0\), \(K_1\), and the hidden-even residual.

The first option is the more direct continuation toward an actual PDE branch. The second is the more direct continuation toward the same-charge / weak-axisymmetric prefactor packet.

stage_v2_22a_profile_to_coefficient_adapter_derivation.md
# Stage V2-22A — Profile-to-Coefficient Adapter and Observable Packet

## 0. Purpose

This stage builds the missing interface between an actual moving-throat profile calculation and the V2-21 branch-extraction fixture.

The V2-21 fixture expects a frozen branch manifest containing lane-level reduced data:

\[
K_A,M_A,
\qquad
\{c_{A\alpha},\varpi_{A\alpha}\},
\qquad
\{\Omega_{U,Ar},\Omega_{W,Ar},R_{Ar},g_{U,Ar},g_{W,Ar}\}.
\]

A PDE or profile solver, however, naturally gives radial/axial profiles such as

\[
\chi_\eta(s),\qquad \phi_\alpha(s),\qquad u_r(s),\qquad w_r(s).
\]

This adapter computes the overlap integrals, converts them into the reduced couplings, builds the V2-21-compatible branch manifest, and then evaluates the same grouped `P2` observable packet as V2-21.

The stage does **not** solve the nonlinear moving-throat PDE. It makes the next numerical branch output machine-readable and falsifiable.

---

## 1. Open-throat geometry convention

The adapter assumes the V2-04 architectural correction:

\[
R_{\rm exit}>0,
\qquad
\texttt{boundary\_class}=\texttt{open\_impedance}.
\]

Thus the finite-throat support ladder is treated as an AC impedance-reflection / organ-pipe condition, not as a hard cap. The branch manifest generated by this stage therefore remains compatible with DC leakage through the open finite-radius exit.

---

## 2. Profile overlaps

Let \(s\in[0,L]\) be the axial coordinate and let \(\mu_s(s)\) be the one-dimensional measure/weight after densitization. The profile adapter computes

\[
I_{\eta\phi,\alpha}
=
\int_0^L\mu_s(s)\chi_\eta(s)\phi_\alpha(s)\,ds,
\]

\[
I_{\eta u,r}
=
\int_0^L\mu_s(s)\chi_\eta(s)u_r(s)\,ds,
\]

\[
I_{\eta w,r}
=
\int_0^L\mu_s(s)\chi_\eta(s)w_r(s)\,ds,
\]

\[
I_{uw,r}
=
\int_0^L\mu_s(s)u_r(s)w_r(s)\,ds.
\]

The built-in analytic D/N demonstration uses

\[
\chi_\eta(s)=\sqrt{\frac{2}{L}}\sin\frac{\pi s}{L},
\]

\[
\phi_{\rm DN}(s)=\sqrt{\frac{2}{L}}\sin\frac{\pi s}{2L},
\]

with

\[
u_r(s)=\chi_\eta(s),
\qquad
w_r(s)=\phi_{\rm DN}(s).
\]

SymPy verifies the exact overlap ledger

\[
\int_0^L\chi_\eta^2\,ds=1,
\]

\[
\int_0^L\phi_{\rm DN}^2\,ds=1,
\]

\[
\int_0^L\chi_\eta\phi_{\rm DN}\,ds
=
\frac{8}{3\pi}.
\]

Therefore the built-in branch has

\[
I_{\eta u}=1,
\qquad
I_{\eta\phi}=I_{\eta w}=I_{uw}=\frac{8}{3\pi}.
\]

---

## 3. Conversion from overlaps to reduced couplings

For each grouped lane \(A\in\{20,21,22\}\), the adapter maps overlap data into reduced coefficients by

\[
c_{A\alpha}=\lambda_{B,A\alpha}I_{\eta\phi,\alpha},
\]

\[
g_{U,Ar}=\lambda_{U,Ar}I_{\eta u,r},
\]

\[
g_{W,Ar}=\lambda_{W,Ar}I_{\eta w,r},
\]

\[
R_{Ar}=\lambda_{R,Ar}I_{uw,r}.
\]

These are exactly the primitive inputs used by the V2-21 branch-extraction fixture.

---

## 4. Lane-level coefficient extraction

The stable BdG support moments are

\[
B_{A0}=\sum_\alpha\frac{c_{A\alpha}^2}{\varpi_{A\alpha}^2},
\]

\[
B_{A2}=\sum_\alpha\frac{c_{A\alpha}^2}{\varpi_{A\alpha}^4},
\]

\[
B_{A4}=\sum_\alpha\frac{c_{A\alpha}^2}{\varpi_{A\alpha}^6}.
\]

For each Maxwell/mixed port,

\[
\Delta_{Ar}=\Omega_{U,Ar}^2\Omega_{W,Ar}^2-R_{Ar}^2,
\]

\[
S_{Ar}=\Omega_{U,Ar}^2+\Omega_{W,Ar}^2,
\]

\[
Q_{Ar}=g_{U,Ar}^2\Omega_{W,Ar}^2+2g_{U,Ar}g_{W,Ar}R_{Ar}+g_{W,Ar}^2\Omega_{U,Ar}^2,
\]

\[
H_{Ar}=g_{U,Ar}^2+g_{W,Ar}^2,
\]

\[
P_{Ar}=\Omega_{U,Ar}^2g_{W,Ar}+R_{Ar}g_{U,Ar}.
\]

Then

\[
Z_{A0}^{(r)}=\frac{Q_{Ar}}{\Delta_{Ar}},
\]

\[
Z_{A2}^{(r)}=\frac{Q_{Ar}S_{Ar}-H_{Ar}\Delta_{Ar}}{\Delta_{Ar}^2},
\]

\[
Z_{A4}^{(r)}=
\frac{Q_{Ar}(S_{Ar}^2-\Delta_{Ar})-S_{Ar}H_{Ar}\Delta_{Ar}}{\Delta_{Ar}^3},
\]

and

\[
N_{A0}^{(r)}=\frac{P_{Ar}^2}{\Delta_{Ar}^2},
\]

\[
N_{A2}^{(r)}=
\frac{2P_{Ar}(P_{Ar}S_{Ar}-\Delta_{Ar}g_{W,Ar})}{\Delta_{Ar}^3},
\]

\[
N_{A4}^{(r)}=
\frac{
\Delta_{Ar}^2g_{W,Ar}^2
-2\Delta_{Ar}P_{Ar}^2
-4\Delta_{Ar}P_{Ar}S_{Ar}g_{W,Ar}
+3P_{Ar}^2S_{Ar}^2
}{\Delta_{Ar}^4}.
\]

Summing over ports gives \(Z_{An}\) and \(N_{An}\).

The conservative operator moments are then

\[
D_{A0}=K_A-B_{A0}-Z_{A0},
\]

\[
D_{A2}=-(M_A+B_{A2}+Z_{A2}),
\]

\[
D_{A4}=-(B_{A4}+Z_{A4}).
\]

---

## 5. Grouped response and outgoing prefactor

The normalized response moments are

\[
u_{2,A}=-\frac{D_{A2}}{D_{A0}},
\]

\[
u_{4,A}=\frac{D_{A2}^2-D_{A0}D_{A4}}{D_{A0}^2}.
\]

The outgoing prefactor moments are

\[
P_{A0}=\frac{N_{A0}}{D_{A0}},
\]

\[
P_{A2}=\frac{D_{A0}N_{A2}-2D_{A2}N_{A0}}{D_{A0}^2},
\]

\[
P_{A4}=\frac{
D_{A0}^2N_{A4}
-2D_{A0}(D_{A2}N_{A2}+D_{A4}N_{A0})
+3D_{A2}^2N_{A0}
}{D_{A0}^3}.
\]

The adapter then performs the weighted grouped decomposition

\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}5,
\]

\[
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\]

\[
b_x=\frac{x_{21}-x_{22}}2.
\]

---

## 6. Weak-axisymmetric branch option

The adapter supports a weak-axisymmetric slope packet of the form

\[
x_A=x_0(1+\epsilon\lambda_As_x),
\]

with

\[
\lambda=(1,\tfrac12,-1).
\]

SymPy verifies the first-order grouped law

\[
\bar x=x_0,
\qquad
b_x=3a_x.
\]

The built-in weak-axisymmetric demo applies small slopes to wall, BdG, and Maxwell/mixed primitive parameters. Because the full extraction contains nonlinear ratios, the final observable packet obeys \(b\approx3a\) to first order, with higher-order residuals suppressed by the small \(\epsilon\).

---

## 7. Run results

The script reports

```text
symbolic_checks: 10/10 passed
```

The exact analytic overlaps are

```text
I_eta_phi = 8/(3*pi) = 0.8488263631567752
I_eta_u   = 1
I_eta_w   = 8/(3*pi) = 0.8488263631567752
I_u_w     = 8/(3*pi) = 0.8488263631567752
```

The sampled-grid copy reproduces the same overlaps within the expected trapezoidal error:

```text
I_eta_phi = 0.8488250214352183
I_eta_u   = 0.999997532602072
I_eta_w   = 0.8488250214352183
I_u_w     = 0.8488250214352183
```

The analytic isotropic branch is open and stable, but not target-calibrated:

```text
open_gate_pass: True
stability_gate_pass: True
target_packet_pass: False
D0_bar: 1.560684541862752
N0_bar: 0.0008248610026844818
P0_bar: 0.0005285251314785051
P0_target: 10.8
R_pole: -2.45246734045587
R_norm: -10.79947147486852
R_P2: 0.0006833304812588467
R_P4: 0.0006221211785235433
P0_anisotropy_norm_sq: 0
```

The weak-axisymmetric demo is also open and stable, and shows the expected small one-dimensional grouped-anisotropy signature:

```text
P0_anisotropy_norm_sq: 5.666182926653727e-14
P0_axisymmetric_b_minus_3a: -1.532472626228142e-10
```

The nonzero `b_minus_3a` is a second-order contamination from finite \(\epsilon\) and nonlinear ratios, not a failure of the first-order law.

---

## 8. Output artifacts

This stage produces:

```text
stage_v2_22a_profile_to_coefficient_adapter.py
stage_v2_22a_profile_adapter_output.txt
stage_v2_22a_profile_to_coefficient_adapter_derivation.md
stage_v2_22a_profile_input_manifest.json
stage_v2_22a_generated_v21_manifest.json
stage_v2_22a_observable_packet.json
```

The important handoff artifact is `stage_v2_22a_generated_v21_manifest.json`: it is a V2-21-compatible manifest built directly from profile overlap data.

---

## 9. Carry-forward status

V2-22A passes as an adapter / extraction bridge:

\[
\boxed{
\text{profiles}\to\text{overlaps}\to\text{reduced coefficients}\to\text{V2-21 observable packet}.
}
\]

It does not claim branch realization. Its purpose is to let the next actual PDE or numerical profile solve export data in a fixed format without changing the target definitions after residuals are known.

stage_v2_22b_solver_handoff_derivation.md
# Stage V2-22B — Numerical Branch Profile Schema and Solver Handoff

## 0. Purpose

Stage V2-22A created a profile-to-coefficient adapter:

\[
\text{profiles}\rightarrow \text{overlaps}\rightarrow B_n,Z_n,N_n\rightarrow D_n,u_n,P_n.
\]

Stage V2-22B freezes the **upstream solver export contract**.  Its job is to define exactly what a future moving-throat PDE solver must output before the V2-22A adapter and V2-21 extraction fixture are allowed to evaluate target residuals.

This is an anti-overfitting and branch-realization gate.  The solver must provide a **target-blind, pre-target-frozen, open-throat branch packet** with enough sampled data to compute the wall/BdG/Maxwell/mixed overlaps.  If the packet is a hard-cap branch, omits the mixed sector, has unstable primitive signs, or fails the grid/profile checks, it is rejected before any quadrupole target residual is computed.

## 1. Solver export schema

The frozen solver output schema is:

```text
stage_v2_22b_solver_handoff/v1
```

The required top-level fields are:

```text
schema
branch_id
freeze
geometry
constants
grid
wall
profiles
bdg_modes
mixed_ports
```

The required freeze block is:

```json
{
  "pre_target_freeze": true,
  "target_blind": true,
  "parent_action": "declared GNLS + localized Maxwell + wall/action convention",
  "gauge_convention": "declared gauge convention",
  "boundary_protocol": "open_impedance_AC_reflecting_DC_leaking"
}
```

The required geometry block is:

```json
{
  "L": "positive finite throat length",
  "R_mouth": "positive mouth radius",
  "R_exit": "strictly positive open exit radius",
  "boundary_class": "open_impedance",
  "Y_L_limit": "finite load-admittance parameter",
  "exit_model": "organ-pipe / impedance mismatch open exit model"
}
```

The validator explicitly forbids:

\[
R_{\rm exit}\le 0,
\qquad
\texttt{boundary\_class}=\texttt{hard\_cap}.
\]

So the old capped-tube interpretation cannot enter the branch-realization packet.

## 2. Mathematical checks inside the validator

The script verifies the D/N support profile convention:

\[
q(s)=\sin(ks),
\qquad
q(0)=0,
\qquad
q_s(L)=0
\quad\text{for}\quad
k=\frac{\pi}{2L}.
\]

It also checks the open-impedance endpoint law:

\[
T_w q_w(L,\omega)+Y_L(\omega)q(L,\omega)=0.
\]

In the AC support-mode reflection limit,

\[
Y_L\to0,
\]

this becomes

\[
q_w(L)=0,
\]

which is the Neumann end of the D/N ladder.  This is compatible with the V2-04 open-organ-pipe patch: AC support modes can reflect through impedance mismatch while the DC/background flow exits through the open finite-radius conduit.

The script also verifies the exact D/N overlap prototype:

\[
\chi_\eta(s)=\sqrt{\frac2L}\sin\frac{\pi s}{L},
\qquad
\phi_{DN}(s)=\sqrt{\frac2L}\sin\frac{\pi s}{2L},
\]

with

\[
\int_0^L \chi_\eta^2\,ds=1,
\qquad
\int_0^L \phi_{DN}^2\,ds=1,
\qquad
\int_0^L \chi_\eta\phi_{DN}\,ds=\frac{8}{3\pi}.
\]

For weak-axisymmetric grouped data, the script confirms the fixed signature:

\[
(20,21,22)\sim\left(1,\frac12,-1\right),
\qquad
b=3a.
\]

For the Maxwell/mixed block, it verifies that the leading outgoing transfer factor is a perfect square:

\[
N_0
=
\frac{(\Omega_U^2 g_W+Rg_U)^2}{(\Omega_U^2\Omega_W^2-R^2)^2}.
\]

So, once

\[
\Delta=\Omega_U^2\Omega_W^2-R^2>0,
\]

is enforced, the outgoing transfer has nonnegative static weight.

## 3. Runtime validation gates

For a solver output packet, the validator checks:

1. required schema keys and declared freeze metadata;
2. target-blind pre-target freeze status;
3. open-throat geometry;
4. monotone grid from \(0\) to \(L\);
5. finite sampled profiles with lengths matching the grid;
6. wall inertia/stiffness positivity;
7. BdG mode positivity \(\varpi>0\);
8. at least one mixed Maxwell/\(A_w\) port;
9. effective mixed-block positivity;
10. availability of the constants needed by V2-21/V2-22A.

The profile diagnostics compute, for example,

\[
\int_0^L \mu_s(s)\chi_\eta(s)^2\,ds
\]

and each BdG-mode norm.  Norm deviations generate warnings, not automatic failures, because a real PDE solver might output nonunit profiles plus explicit normalization conventions.  Structural issues such as hard-cap geometry or nonpositive \(\Delta\) are errors.

## 4. Conversion into the V2-22A profile manifest

If validation passes, the script converts the solver packet into the V2-22A profile manifest schema:

```text
stage_v2_22a_profile_adapter/v1
```

It maps:

```text
solver profiles.weight       -> profile "weight"
solver profiles.wall_chi_eta -> profile "wall"
solver bdg_modes[*]          -> sampled BdG profiles
solver mixed_ports[*].u/w    -> sampled U/W mixed profiles
```

and preserves the reduction coefficients:

\[
K,
\quad
M,
\quad
\lambda_B,
\quad
\varpi,
\quad
\lambda_U,
\quad
\Omega_U,
\quad
\lambda_W,
\quad
\Omega_W,
\quad
\lambda_R.
\]

The generated profile manifest is then immediately consumable by the V2-22A adapter.

## 5. Built-in valid sample result

The built-in sample solver export is an open D/N finite-throat prototype with:

```text
R_exit = 0.35
boundary_class = open_impedance
pre_target_freeze = true
target_blind = true
```

The V2-22B validation result was:

```text
symbolic_checks: 9/9 passed
validation_pass: True
error_count: 0
warning_count: 0
packet_hash: 959ba8b19d5c8afc683c006e0214c0fa2794f873ae7b611bf6d8a527df66b7d1
generated_profile_manifest_hash: 46da7fc420e753dd749e3e590f24614b82705399497fd79a3e227505dcba0a9c
```

A smoke test through V2-22A then produced:

```text
I_eta_phi::bdg_halfwave = 0.8488255450333207
I_eta_u::one_mixed_port = 1
I_eta_w::one_mixed_port = 0.8488255450333207
I_u_w::one_mixed_port = 0.8488255450333207
```

and the downstream observable packet remained open/stable but target-failing, as expected for an untuned prototype:

```text
open_gate_pass: True
stability_gate_pass: True
target_packet_pass: False
D0_bar: 1.560684600911219
N0_bar: 0.0008248594058984478
P0_bar: 0.0005285240883499759
P0_target: 10.8
R_norm: -10.79947147591165
```

This is the intended behavior.  V2-22B and V2-22A are not supposed to rescue a branch; they only validate, convert, and extract.

## 6. Built-in invalid hard-cap rejection

The invalid sample changes the exit to:

```text
R_exit = 0
boundary_class = hard_cap
boundary_protocol = hard_cap_DC_blocked
```

The validator rejected it:

```text
validation_pass: False
error_count: 4
```

The specific errors were:

```text
freeze.boundary_protocol: must declare open impedance AC/DC split
geometry.boundary_class: boundary_class must be open_impedance
geometry.R_exit: R_exit must be strictly positive; hard cap is forbidden
geometry.boundary_class: hard_cap geometry is forbidden by V2-04 patch
```

This confirms that the capped-tube mass-continuity bug cannot silently re-enter the PDE branch pipeline.

## 7. Status

V2-22B passes as a solver-handoff layer:

\[
\boxed{\text{PDE solver output}\rightarrow\text{validated frozen branch packet}\rightarrow\text{V2-22A profile manifest}.}
\]

It does not solve the moving-throat PDE.  It fixes the exact data contract required before an actual numerical branch can be judged.

## 8. Immediate continuation

The next natural stage is V2-22C:

```text
profile-to-fixture end-to-end smoke pipeline and tolerance budget
```

That stage should run the generated V2-22A manifest automatically through the V2-22A adapter and V2-21 extraction fixture, then emit a single end-to-end branch-realization report with:

```text
solver packet hash
profile manifest hash
V2-21 manifest hash
observable packet hash
open/stability gates
isotropy gates
target residuals
profile normalization diagnostics
```

That would turn the current validation stack into a one-command branch-realization pipeline.

stage_v2_22c_end_to_end_smoke_pipeline_derivation.md
# Stage V2-22C — End-to-end branch-realization smoke pipeline and tolerance budget

## Purpose

V2-22C closes the executable handoff chain:

\[
\text{solver export}
\rightarrow
\text{V2-22B validation}
\rightarrow
\text{V2-22A profile adapter}
\rightarrow
\text{V2-21 grouped observable extraction}
\rightarrow
\text{residual/tolerance packet}.
\]

This stage does **not** claim a solved moving-throat branch. Its purpose is to prove that a future PDE branch can be tested without changing the definitions after seeing the residuals.

The resulting executable stage distinguishes three classes of failure:

1. **handoff/validation failure**, such as a hard-cap geometry or nonpositive mixed-block determinant;
2. **open-throat/stability failure**, such as `D0 <= 0` or `Delta <= 0`;
3. **target-realization failure**, such as nonzero `R_norm`, `R_pole`, `R_P2`, or `R_P4`.

A stable branch is allowed to fail the target packet. That is a useful falsification, not a pipeline error.

---

## Inputs

The stage reuses the existing executable tools:

- `stage_v2_22b_solver_handoff_validator.py`
- `stage_v2_22a_profile_to_coefficient_adapter.py`
- `stage_v2_21_branch_extraction_fixture.py`

The built-in solver export is the open finite-throat D/N prototype:

\[
R_{\rm exit}>0,
\qquad
\texttt{boundary\_class}=\texttt{open\_impedance}.
\]

The invalid control is the forbidden hard-cap packet:

\[
R_{\rm exit}=0,
\qquad
\texttt{boundary\_class}=\texttt{hard\_cap}.
\]

---

## Formula checks performed in V2-22C

V2-22C runs a lightweight orchestration audit. It intentionally does not rerun the heavier component SymPy audits, because this stage is a handoff pipeline rather than another symbolic compiler.

The checked identities are:

### 1. Grouped inverse map

For

\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}5,
\]

\[
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\qquad
b_x=\frac{x_{21}-x_{22}}2,
\]

the inverse is

\[
x_{20}=\bar x+4a_x,
\]

\[
x_{21}=\bar x-a_x+b_x,
\]

\[
x_{22}=\bar x-a_x-b_x.
\]

### 2. Weak-axisymmetric grouped signature

For

\[
\lambda=(1,\tfrac12,-1),
\]

one has

\[
\bar x=x_0,
\qquad
b_x=3a_x.
\]

### 3. Constant-prefactor branch

The constant-prefactor conditions are

\[
N_2=\frac{2D_2N_0}{D_0},
\]

\[
N_4=\frac{N_0(D_2^2+2D_0D_4)}{D_0^2}.
\]

Substitution gives

\[
P_2=0,
\qquad
P_4=0.
\]

### 4. Quadrupole normalization equivalence

The target

\[
P_0^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}
\]

is equivalent to

\[
\gamma_{\rm quad}^{\rm eff}
=
\frac{2G}{5c^5}
\]

when

\[
\gamma_{\rm quad}^{\rm eff}
=\widehat m_0^{\,2}\mathcal S_{\rm port}P_0\frac{a^5}{27c_s^5}.
\]

The orchestration audit passed all `8/8` formula checks.

---

## Pipeline gates

### V2-22B validation gates

The solver export must satisfy:

\[
R_{\rm exit}>0,
\]

\[
\texttt{boundary\_class}=\texttt{open\_impedance},
\]

\[
\texttt{pre\_target\_freeze}=\texttt{true},
\qquad
\texttt{target\_blind}=\texttt{true}.
\]

The mixed block must satisfy

\[
\Delta=\Omega_U^2\Omega_W^2-R^2>0.
\]

The valid open-throat sample passed. The invalid hard-cap sample was rejected.

### V2-21 extraction gates

For each grouped lane, the extracted moments are

\[
D_0=K-B_0-Z_0,
\]

\[
D_2=-(M+B_2+Z_2),
\]

\[
D_4=-(B_4+Z_4).
\]

The normalized response and outgoing prefactor are

\[
u_2=-\frac{D_2}{D_0},
\]

\[
u_4=\frac{D_2^2-D_0D_4}{D_0^2},
\]

\[
P_0=\frac{N_0}{D_0},
\]

\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]

\[
P_4=
\frac{
D_0^2N_4
-2D_0(D_2N_2+D_4N_0)
+3D_2^2N_0
}{D_0^3}.
\]

The target residual packet is

\[
R_{\rm pole}=D_0(B_4+Z_4)-3(M+B_2+Z_2)^2,
\]

\[
R_{\rm norm}=\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
-
\frac{54Gc_s^5}{5a^5c^5},
\]

\[
R_{P2}=P_2,
\qquad
R_{P4}=P_4,
\]

\[
R_{\rm tail}=\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3-1.
\]

---

## Smoke-branch output

The open D/N smoke branch passed validation and stability, but it did not hit the target surface:

```text
open_gate_pass: True
stability_gate_pass: True
target_packet_pass: False
D0_bar: 1.560684600911219
N0_bar: 0.0008248594058984478
P0_bar: 0.0005285240883499759
P0_target: 10.8
R_pole: -2.452467306741504
R_norm: -10.79947147591165
R_P2: 0.0006833291043178326
R_P4: 0.0006221198975181069
R_tail: 0
```

This is the intended behavior. The sample branch is a stable open-throat profile fixture, not a calibrated realization of the universal quadrupole target.

---

## Calibrated direct-coefficient control

The calibrated direct-coefficient control passed the full target packet:

```text
target_packet_pass: True
R_pole: 0
R_norm: 0
R_P2: 0
R_P4: 1.4210854715202e-14
```

The tiny nonzero `R_P4` is floating-point roundoff.

This control proves that the gate can pass when the target surface is actually satisfied. Therefore the smoke branch failure is not caused by a broken extraction pipeline.

---

## Tolerance budget

The default extraction tolerance is

\[
\varepsilon_{\rm extract}=10^{-9}.
\]

The validator profile normalization tolerance is

\[
\varepsilon_{\rm norm}=5\times10^{-3}.
\]

The mixed-block determinant tolerance is

\[
\varepsilon_\Delta=10^{-12}.
\]

The smoke branch’s large normalized residual is dominated by `R_norm`, not by numerical roundoff:

\[
\left|\frac{R_{\rm norm}}{P_0^{\rm target}}\right|
\approx
0.999951.
\]

So the branch is very far below the required outgoing normalization.

---

## Final V2-22C conclusion

The executable branch-realization chain is now closed:

\[
\boxed{
\text{frozen solver packet}
\to
\text{validated profile manifest}
\to
\text{V2-21 observable packet}
\to
\text{residual/tolerance report}.
}
\]

The mechanical pipeline passes:

```text
mechanical_pipeline_pass: True
```

But no PDE branch realization is claimed:

```text
branch_target_realization_claimed: False
```

The next actual science step is to feed this pipeline a real moving-throat PDE branch export and see whether the resulting branch, without post-target refitting, lands on the isotropic full-bundle target surface.

stage_v2_23_minimal_branch_solver_derivation.md
# Stage V2-23 — Minimal open-throat branch solver and first real residual extraction

## 0. Purpose

This stage is the first reduced branch-realization prototype after the V2-22C
handoff pipeline.

It is **not** a full nonlinear GNLS/Maxwell/moving-wall PDE solution. It is a
target-blind, frozen, one-dimensional open-throat branch solve that moves beyond
mock profiles by solving:

1. a stationary finite-radius open throat profile,
2. a wall `l=2` Sturm-Liouville support profile,
3. a stable BdG support profile,
4. a brane-like localized gauge profile,
5. a mixed `A_w/F_{\mu w}/J^w` profile,

then computing overlap integrals and passing the resulting branch through the
same isotropic full-bundle residual packet used in V2-19 through V2-22C.

The branch is intentionally simple. Its job is to give Codex/local solvers a
concrete executable template for the next full PDE branch export.

## 1. Frozen branch definition

```text
branch_id: v2_23_minimal_open_throat_frozen_demo
branch_freeze_hash: 0feee821ffec9b23f79d80bb4e176139f01fbffb3f1358607dad6597032592ae
pre_target_freeze: true
target_blind: true
no_post_residual_refit: true
boundary_class: open_impedance
boundary_protocol: open_impedance_AC_reflecting_DC_leaking
```

The geometry is an open finite throat, not a capped tube:

\[
R(0)=a,\qquad R(L)>0.
\]

The stationary profile is obtained by minimizing

\[
E[R]
=
\frac12\int_0^L T_R (R')^2\,ds
+
\frac12\int_0^L K_R(R-R_{\rm pref})^2\,ds
+
\frac12Y_{\rm exit}(R(L)-R_{\rm exit,pref})^2,
\]

with fixed mouth radius \(R(0)=a\). The exit condition is a finite-radius open
penalty rather than a cap.

The solved branch gives

\[
R_{\rm mouth}=1,\qquad
R_{\rm exit}=0.452938042901,\qquad
R_{\rm min}=0.452938042901.
\]

The open gate passes because \(R_{\rm exit}>0\).

## 2. Axial measure

The solver uses the geometry-derived effective axial measure

\[
\mu_s(s)=R_0(s)^2\sqrt{1+R_0'(s)^2},
\]

renormalized so that

\[
\int_0^L\mu_s(s)\,ds=L.
\]

This remains a reduced one-dimensional model, but the overlap integrals now
depend on the solved open-throat geometry rather than on a flat hand-inserted
measure.

## 3. Solved Sturm-Liouville problems

Each profile is obtained from a finite-element Sturm-Liouville problem

\[
-\partial_s(T(s)\partial_s q)+V(s)q=\lambda\mu_s(s)q,
\]

with mouth Dirichlet condition

\[
q(0)=0,
\]

and open-end impedance condition

\[
T(L)q'(L)+Y_Lq(L)=0.
\]

For this first branch prototype the AC-reflecting organ-pipe limit is used:

\[
Y_L=0,
\]

so the exit is Neumann-like for AC support modes while remaining geometrically
open for DC/background flow.

The solved eigenvalues are:

\[
K=\lambda_{{\rm wall},l=2}=2.2393180779,
\]

\[
\varpi^2=3.16919632623,
\]

\[
\Omega_U^2=3.78509395487,\qquad
\Omega_W^2=4.05724378721.
\]

The FEM residuals are reported in `stage_v2_23_tolerance_report.json`.

## 4. Overlap integrals

The reduced overlap integrals are

\[
I_{\eta\phi}=\int_0^L\mu_s\chi_\eta\phi\,ds,
\]

\[
I_{\eta U}=\int_0^L\mu_s\chi_\eta u\,ds,
\]

\[
I_{\eta W}=\int_0^L\mu_s\chi_\eta w\,ds,
\]

\[
I_{UW}=\int_0^L\mu_s u w\,ds.
\]

For this frozen branch,

\[
I_{\eta\phi}=0.999946613757,
\qquad
I_{\eta U}=0.999180452256,
\]

\[
I_{\eta W}=0.994783990782,
\qquad
I_{UW}=0.998096444041.
\]

The reduced couplings are

\[
c_B=\lambda_B I_{\eta\phi},
\qquad
g_U=\lambda_U I_{\eta U},
\qquad
g_W=\lambda_W I_{\eta W},
\qquad
R=\lambda_R I_{UW}.
\]

## 5. Reduced coefficients

The stable BdG contribution is

\[
B_0=\frac{c_B^2}{\varpi^2},
\qquad
B_2=\frac{c_B^2}{\varpi^4},
\qquad
B_4=\frac{c_B^2}{\varpi^6}.
\]

The conservative Maxwell/mixed block uses

\[
\Delta=\Omega_U^2\Omega_W^2-R^2,
\]

\[
Q=g_U^2\Omega_W^2+2g_Ug_WR+g_W^2\Omega_U^2,
\]

\[
H=g_U^2+g_W^2,\qquad
S=\Omega_U^2+\Omega_W^2.
\]

Then

\[
Z_0=\frac Q\Delta,
\]

\[
Z_2=\frac{QS-H\Delta}{\Delta^2},
\]

\[
Z_4=\frac{Q(S^2-\Delta)-SH\Delta}{\Delta^3}.
\]

The outgoing-transfer moments are

\[
N_0=\frac{(\Omega_U^2g_W+Rg_U)^2}{\Delta^2},
\]

\[
N_2=
\frac{2P(P S-\Delta g_W)}{\Delta^3},
\qquad
P=\Omega_U^2g_W+Rg_U,
\]

\[
N_4=
\frac{\Delta^2g_W^2-2\Delta P^2-4\Delta PSg_W+3P^2S^2}{\Delta^4}.
\]

The total conservative operator is

\[
D_0=K-B_0-Z_0,
\]

\[
D_2=-(M+B_2+Z_2),
\]

\[
D_4=-(B_4+Z_4).
\]

The solved branch gives

\[
D_0=1.89448908712,\qquad
D_2=-1.09827027208,\qquad
D_4=-0.0282005090292.
\]

## 6. Observable packet

The normalized response and outgoing prefactor are

\[
u_2=-\frac{D_2}{D_0},
\]

\[
u_4=\frac{D_2^2-D_0D_4}{D_0^2},
\]

\[
P_0=\frac{N_0}{D_0},
\]

\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]

\[
P_4=
\frac{D_0^2N_4-2D_0(D_2N_2+D_4N_0)+3D_2^2N_0}{D_0^3}.
\]

For this branch,

\[
u_2=0.579718447339,\qquad
u_4=0.350959026601,
\]

\[
P_0=0.0197851167073,\qquad
P_2=0.0332589572512,\qquad
P_4=0.0365510452457.
\]

## 7. Target residuals

The residual packet is

\[
R_{\rm pole}=D_0(B_4+Z_4)-3(M+B_2+Z_2)^2,
\]

\[
R_{\rm norm}
=
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
-
\frac{54Gc_s^5}{5a^5c^5},
\]

\[
R_{P2}=P_2,\qquad R_{P4}=P_4,
\]

\[
R_{\rm tail}=\Theta_{\rm tail}\left(\frac c{c_s}\right)^3-1.
\]

The branch output is:

```text
R_pole = -3.56516721502
R_norm = -10.7802148833
R_P2   = 0.0332589572512
R_P4   = 0.0365510452457
R_tail = 0
```

The target packet therefore fails:

```text
target_packet_pass = False
```

The failure is expected. This first reduced branch is open and stable, but it was
not tuned to the isotropic full-bundle target surface.

The normalization ratio is

\[
\frac{P_0}{P_0^{\rm target}}=0.00183195525067.
\]

The one-pole ratio is

\[
\frac{D_0(B_4+Z_4)}{3(M+B_2+Z_2)^2}=0.0147641804366.
\]

So this branch undershoots the outgoing normalization and misses the conservative
one-pole surface.

## 8. Verdict

```text
open_gate_pass: True
stability_gate_pass: True
outgoing_transfer_gate_pass: True
target_packet_pass: False
```

This is the desired first-real-branch behavior:

- the open-throat geometry is valid,
- the reduced conservative branch is stable,
- the mixed/outgoing transfer is non-dark,
- but the branch fails the target residuals honestly.

That means the V2-22C pipeline is now usable for real reduced branch exports:
it does not merely process mock data, and it does not rescue a failing branch.

## 9. Next handoff to Codex

The next Codex/local continuation should replace the toy 1D operators with a
true moving-throat export while preserving this schema:

1. solve the stationary open-throat branch,
2. export `R0(s)`, wall coefficients, BdG profiles, gauge/mixed profiles,
3. compute `B_n,Z_n,N_n`,
4. run the same residual packet,
5. do not change the branch after reading residuals.

The most important data to improve are:

- `P0/P0_target`, currently far too small,
- `D0(B4+Z4)/(3A^2)`, currently far from the one-pole surface,
- the constant-prefactor conditions `P2=0`, `P4=0`.

stage_v2_25_actual_branch_protocol_derivation.md
# Stage V2-25 - Actual Branch Protocol and Notes Intake

## Purpose

This stage records the mathematical handoff created after the first
target-blind simulation miss and the intake of the unincorporated 5PN, barrier,
atom/lepton, and `P22` notes.

The point is to keep the next calculation from being lost in the notes:

- the reduced target algebra is already sharp;
- the current simulated branch families miss the target;
- the unincorporated notes add useful source and mouth physics;
- but those notes do not yet provide a calibrated physical exporter.

So the next load-bearing object is an actual moving-throat branch packet, not a
post-hoc retuning of the present simulation coefficients.

The executable companion for this note is
`simulation/diagnose_notes_intake.py`, which verifies the source anchors and
writes `simulation/output/notes_intake_report.json`.

---

## 1. Status after the target-blind simulation miss

The current simulation bundle freezes candidate packets before target
evaluation, then runs the existing V2-22B -> V2-22A -> V2-21 chain.

The current referee run gives:

- `0/192` target-passing reduced frozen candidates;
- `0/3` target-passing manufactured nonlinear candidates;
- reduced open-stable one-pole ratio
  `D0*C/(3*A^2)` between `0.0033775383274364888` and
  `0.1353664855760648`;
- reduced median required `C` or `D0` multiplier
  `16.30132163440465`;
- reduced median required `P0` multiplier `171.65261223353198`;
- projection stress shows that one-pole support alone and uniform outgoing
  amplitude scaling are both insufficient.

This is evidence against the current reduced and manufactured simulation
families.  It is not a theorem against every possible nonlinear moving-throat
branch.

---

## 2. Coefficient map

On the isotropic grouped real `P2` branch, write the conservative denominator
and outgoing numerator as

```text
D(omega) = D0 + D2 omega^2 + D4 omega^4 + O(omega^6),
N(omega) = N0 + N2 omega^2 + N4 omega^4 + O(omega^6).
```

The full-bundle aliases are

```text
D0 = K - B0 - Z0,
D2 = -(M + B2 + Z2),
D4 = -(B4 + Z4).
```

Equivalently,

```text
A = M + B2 + Z2 = -D2,
C = B4 + Z4 = -D4.
```

Here:

- `K,M` are wall/worldtube stiffness and inertia data;
- `B0,B2,B4` are stable BdG support moments;
- `Z0,Z2,Z4` are conservative Maxwell/mixed moments;
- `N0,N2,N4` are outgoing-transfer moments.

The normalized conservative response is

```text
Y(omega) = D0 / D(omega)
         = 1 + u2 omega^2 + u4 omega^4 + O(omega^6),
```

with

```text
u2 = A/D0,
u4 = (A^2 + D0*C)/D0^2.
```

The one-pole condition is `u4 = 4 u2^2`, so

```text
D0*C = 3*A^2,
```

or

```text
D0*C/(3*A^2) = 1.
```

That is the one-pole ratio reported by the simulation diagnostics.

---

## 3. Outgoing prefactor and moment-shape conditions

The outgoing prefactor is

```text
P(omega) = D0*N(omega) / D(omega)^2.
```

Expanding through `O(omega^4)` gives

```text
P0 = N0/D0,
P2 = (D0*N2 + 2*A*N0)/D0^2,
P4 = (D0^2*N4 + 2*D0*(A*N2 + C*N0) + 3*A^2*N0)/D0^3.
```

The constant-prefactor branch requires

```text
P2 = 0,
P4 = 0.
```

Equivalently,

```text
N2 = -2*A*N0/D0,
N4 = N0*(A^2 - 2*D0*C)/D0^2.
```

On the one-pole surface `D0*C = 3*A^2`, this reduces to

```text
N4 = -5*A^2*N0/D0^2.
```

This is why the projection-stress diagnostic matters.  Raising `C` or `D0` to
the one-pole surface is not enough.  Scaling the outgoing amplitude to fix `P0`
is also not enough.  The branch must realize the outgoing moment shape through
`N2` and `N4` on the same frozen branch.

---

## 4. Actual branch packets

The 5PN handoff reduces the remaining computation to the actual PDE-selected
branch.  The symbolic target surface is not the missing piece; the missing
piece is the realized branch data.

### Packet A

For each grouped lane `A in {20, 21, 22}`, export

```text
D_A0, D_A2, D_A4,
N_A0, N_A2, N_A4.
```

The same packet must also include

```text
mhat_0,
N_Q or chi_Q,
parent_action_status,
boundary_protocol,
stability_certificate,
source_hashes,
freeze_hash.
```

### Packet B

The same realized branch must export one equivalent orbit-lock representation:

```text
m_T, m_K, m_mu
```

or

```text
R_tr, R_nt, R_eta
```

or

```text
q_tr, q_nt, q_eta.
```

The four finish-line conditions are

```text
dln_R_tr = 0,
dln_R_target = 0,
dln_epsilon_eta = 0,
N_Q = 1.
```

The last condition may also be written as `chi_Q = 1` on the natural source-map
branch.

---

## 5. What the unincorporated notes add

The note intake found useful physical ingredients, but not a completed
exporter.

Promotable reduced ingredients:

- leakage/work lane: useful support-side source/work accounting;
- non-rigid `U/V` dressing: useful branch component for mouth/dressing and
  orbit-side response;
- microscopic export kernel: useful odd passive/export term in the active `V`
  equation;
- finite-throat atomic `P0/P2` source: useful replacement for point-source
  forcing;
- open/radiative scalar `P0` flux: possible outgoing-normalization hook;
- intrinsic `P22` bracing: useful mouth-shape source if the half-flux/mixed
  closure is realized.

But the missing calibrated maps remain:

```text
source physics -> D0, C,
source physics -> P0,
source physics -> N2, N4.
```

In particular:

- support/source enhancement is not the active bottleneck in the reduced 5PN
  stack;
- scalar `P0` hammering does not linearly drive the area-preserving `P22`
  mouth mode;
- outgoing transfer moments are support-blind in the explicit BdG pair, so
  support tuning alone cannot supply `N2/N4` control.

---

## 6. No-refit rules

The next exporter must preserve the V2-16 freeze discipline.

Before target residuals are evaluated, freeze:

- parent action status or effective closure declaration;
- gauge convention;
- open-exit boundary class;
- wall/interface action;
- support basis;
- mixed/outgoing port list;
- stability gates;
- extraction formulas;
- source terms admitted from the atom/lepton/barrier notes;
- source hashes and protocol hash.

After residuals are known, do not change:

- support cardinality;
- boundary class;
- gauge convention;
- source normalization;
- port normalization;
- outgoing branch;
- extraction formulas.

Do not project a realized `chi_Q != 1` branch back to `chi_Q = 1`.  Do not use
support enhancement to explain an orbit-lock miss.  Do not report an
algebraically projected zero-residual packet as a physical target-blind hit.

---

## 7. Next calculation

The next honest calculation is:

1. solve or continue a stationary open moving-throat branch;
2. linearize wall, BdG, Maxwell/mixed, and outgoing-port sectors;
3. extract `K,M`, `B_n`, `Z_n`, and `N_n`;
4. export Packet A and Packet B with freeze/source hashes;
5. run the unchanged target-blind guard and post-hoc residual chain;
6. report `R_pole`, `R_norm`, `R_P2`, and `R_P4` without refitting.

That is the point at which the unincorporated notes become testable branch
physics instead of promising but uncalibrated reduced ingredients.

---

## Current claim status

```text
Reduced target algebra: passed.
Target-blind reduced simulation: target miss.
Manufactured nonlinear readiness: target miss.
Notes intake: passed as source/protocol evidence.
Physical actual-branch exporter: not yet implemented.
Post-hoc retuning allowed: no.
Next required artifact: actual moving-throat branch packet.
```

stage_v2_26_program_status_after_audit.md
# Stage V2-26 - Program Status After the PDE Audit

## Purpose

This note records the honest program status after the PDE audit and the first
target-blind simulation layer.

The audit was useful because it did not merely check algebra.  It separated:

- what has been derived cleanly;
- what is conditional or reduced;
- what has executable support;
- what the current simulations fail to realize;
- and what physical derivations are still missing before the PDE can be
  claimed as a completed one-equation account.

This note is intended to be paper-facing.  It should prevent the audit result
from being summarized too strongly or too weakly.

For the physical ontology that should accompany this status statement, see
`stage_v2_28_physical_picture_and_ontology.md`.  That note records the
finite brane-bulk throat/puncture picture, the open-conduit condition, the
mouth/interior distinction, flux distinctions, and common misreadings to avoid.

For the superfluid material-closure gap, see
`stage_v2_29_superfluid_material_closure_gap.md`.  That note records the
difference between existing reduced EOS/sound-speed formulas and a completed
branch-level derivation of density, sound speed, effective light-speed
behavior, and flux feedback.

---

## 1. What can be said solidly

The current project contains real mathematical work.

The audit verifies many nontrivial pieces of the reduced/effective framework:

- parent-action versus effective-wall status;
- Maxwell localization and gauge normalization;
- dimension and port/source normalization ledgers;
- open finite-radius throat boundary logic;
- Poisson/Newtonian hook;
- stable BdG support Schur complement;
- Maxwell/mixed conservative and outgoing kernels;
- Hamiltonian positivity gates;
- grouped real `P2` projector calculus;
- STF angular source-map normalization;
- grouped normalization and constant-prefactor branch algebra;
- outgoing `l=2` fingerprint;
- 2.5PN/4PN interface and tail transport gate;
- branch-freeze/no-refit discipline;
- weak-axisymmetric splitting;
- similarity-orbit separation;
- isotropic full-bundle target surface;
- weak-form and branch-extraction scaffolding;
- fixture-backed solver handoff and negative controls;
- reduced FEM and manufactured nonlinear simulation checks.

These are not just narrative assertions.  Many of them are backed by executable
Python audits, Mathematica mirrors, fixtures, generated artifacts, and
target-blind guards.

So the proper positive claim is:

```text
The program has derived a growing reduced/effective framework whose algebraic
targets, adapters, stability gates, and no-refit protocol are internally
consistent and executable under stated assumptions.
```

---

## 2. What cannot be claimed yet

The audit does not prove that known physics has been derived from one completed
PDE.

The PDE is still incomplete as a physical branch-realization system.  In
particular, the current audited exporter does not yet derive one frozen
moving-throat branch that simultaneously outputs the required conservative,
outgoing, and orbit-lock data.

The following claim would be too strong:

```text
Known physics has been derived from a single PDE.
```

The current honest replacement is:

```text
Several reduced/effective pieces of known-physics structure have been derived
or targeted under explicit assumptions, and the audit now identifies the exact
branch-output conditions that a completed PDE must satisfy.
```

The distinction matters.  The target algebra can be correct while the actual
PDE branch remains unsolved.

---

## 3. Why the clean algebra and failed simulations are not contradictory

The algebra answers a conditional question:

```text
If a moving-throat branch is going to reproduce the desired GR-like packet,
what coefficient relations must it satisfy?
```

The simulations answer a different question:

```text
Do the currently exportable reduced/manufactured branch families land on that
target surface without target feedback?
```

The current answer to the second question is no.

The referee run gives:

```text
0/192 target-passing reduced frozen candidates.
0/3 target-passing manufactured nonlinear candidates.
```

Among reduced open-stable candidates, the one-pole ratio

```text
D0*C/(3*A^2)
```

has maximum `0.1353664855760648`, far below the target value `1`.

The required-deformation diagnostics show that the miss is not a small local
continuation:

```text
minimum reduced C-or-D0 multiplier: 7.387352901601946
median reduced C-or-D0 multiplier: 16.30132163440465
median reduced P0 multiplier: 171.65261223353198
```

The projection-stress diagnostic then shows that even fixing one-pole support
and uniform outgoing amplitude after the fact is insufficient.  A successful
branch needs outgoing moment-shape control through `N2` and `N4`, not just a
scalar multiplier.

So the simulation miss means:

```text
The current exported branch families are not the desired physical branch.
```

It does not mean:

```text
The target algebra is inconsistent.
```

It also does not mean:

```text
The desired branch exists.
```

That existence question remains the central unsolved problem.

---

## 4. What the audit says is missing

The audit has turned a broad uncertainty into a concrete missing-physics list.

The completed branch must realize, on one frozen branch:

```text
D0 = K - B0 - Z0,
A  = M + B2 + Z2,
C  = B4 + Z4,
D0*C/(3*A^2) = 1,
P0 = N0/D0,
P2 = 0,
P4 = 0.
```

Equivalently, the outgoing transfer moments must satisfy

```text
N2 = -2*A*N0/D0,
N4 = N0*(A^2 - 2*D0*C)/D0^2.
```

On the one-pole surface this becomes

```text
N4 = -5*A^2*N0/D0^2.
```

The 5PN/orbit-lock side additionally requires the actual branch to satisfy:

```text
dln_R_tr = 0,
dln_R_target = 0,
dln_epsilon_eta = 0,
N_Q = 1.
```

The current materials do not yet provide the physical exporter that determines
all of these quantities from the throat equations.

---

## 5. Physical mechanisms not yet represented

The five audit readouts `D0`, `C`, `P0`, `N2`, and `N4` are not a full physical
description of the throat.  They are low-order projected response moments.

They do not directly encode:

- the actual throat geometry `R0(s)`;
- wall shape and axial boundary-layer structure;
- superfluid intake and mass balance;
- DC leakage and open/radiative mouth flux;
- reservoir depletion or backreaction;
- nonlinear saturation and branch selection;
- full BdG/Krein mode signatures;
- complete Maxwell/mixed port spectra;
- finite-throat `P0/P2/P22` mouth geometry;
- intrinsic `P22` bracing and orientation;
- half-flux/mixed-sector closure;
- orbit-lock tangent data;
- thermodynamic or heat/export partition;
- superfluid material closure: EOS, stationary density, `c_s(rho)`, effective
  light-speed behavior if density dependent, and flux feedback;
- higher dissipative PN information;
- the physical origin of source/port normalization.

Those mechanisms may affect the five readouts, but the readouts alone do not
derive the mechanisms.

This is why the final PDE is harder than a coefficient-fitting problem.  It
must physically select the branch whose projected response moments satisfy the
audit packet.

---

## 6. What the unincorporated notes contribute

The notes-intake audit found real candidate ingredients:

- finite-throat `P0/P2` forcing from Hessian/tidal mouth loading;
- scalar open/radiative `P0` flux hooks;
- intrinsic `P22` bracing from mixed-sector/half-flux structure;
- non-rigid `U/V` dressing;
- leakage/work lanes;
- microscopic odd passive/export kernels;
- support/source branch data that is no longer the active bottleneck in the
  reduced 5PN stack.

These are useful.  They should be carried into the next derivation phase.

But they are not yet a completed exporter.  They do not yet provide a calibrated
map:

```text
throat/source physics -> D0, C, P0, N2, N4
```

and they do not yet prove that the same branch satisfies the orbit-lock packet.

So the notes are best understood as:

```text
source physics and protocol constraints for the next branch derivation,
not a solved target-passing PDE branch.
```

---

## 7. Publication framing

The audit should be published as a boundary-setting document.

It should claim:

- the algebraic and executable audit layers are internally consistent;
- the reduced/effective claims have been separated from strict parent-PDE
  claims;
- the current target-blind branch families fail the target;
- the miss is quantitatively diagnosed;
- post-hoc rescue/tuning is explicitly disallowed;
- the remaining PDE completion problem is now sharply specified.

It should not claim:

- a completed derivation of known physics from one PDE;
- an actual physical branch that already passes the target packet;
- that support tuning alone can explain the simulation miss;
- that scalar outgoing normalization alone can repair the branch;
- that algebraic projection after seeing residuals is physical evidence.

Recommended paper language:

```text
The present audit does not complete the moving-throat PDE program.  It
establishes a reproducible reduced/effective ledger, identifies the exact
coefficient packet required of any successful branch, and shows that the
currently exportable branch families do not realize that packet.  The remaining
problem is the derivation and target-blind export of an actual moving-throat
branch containing the missing intake, finite-mouth, outgoing, and orbit-lock
physics.
```

---

## 8. Remaining derivation program

The next phase should be derivation-first and gap-driven.

Required derivations:

1. Decide the strict parent/effective status of the wall/throat dynamics.
   Either promote `S_Sigma` or keep an explicit effective closure label.

2. Derive the actual stationary open-throat branch equation, including the
   physical terms admitted from the finite-mouth, scalar-flux, `P22`, `U/V`,
   leakage/work, and export-kernel notes.

3. Derive the linearized wall/BdG/Maxwell/mixed/outgoing-port exporter that
   maps the branch to `K,M,B_n,Z_n,N_n`.

4. Derive the source/port normalization law that produces `P0` as branch data
   rather than as a fitted scale.

5. Derive independent outgoing moment-shape control for `N2` and `N4`.

6. Derive the weak-axisymmetric tangent and orbit-lock packet on the same
   branch.

7. Freeze the complete protocol before target evaluation.

8. Run the unchanged V2-22B -> V2-22A -> V2-21 residual chain and report the
   result without refitting.

Higher PN work may still be useful, especially dissipative half-integer orders
such as 3.5PN, 4.5PN, and 5.5PN, because they may constrain outgoing
moment-shape physics.  But more PN algebra is secondary until the actual branch
exporter exists.

---

## Current status statement

```text
Real mathematical audit work: yes.
Real reduced/effective physics bridges: yes, under stated assumptions.
Completed one-PDE derivation of known physics: no.
Current physical exporter: incomplete.
Current target-blind simulated branch families: fail.
Main remaining problem: derive and export the actual moving-throat branch.
```

The audit therefore moves the program forward by changing the question from

```text
Can the algebra be arranged to look right?
```

to

```text
Can the completed PDE physically select a stable open branch whose frozen
projected response packet satisfies the audited target conditions?
```

stage_v2_27_executable_audit_implementation_summary.md
# Stage V2-27 - Executable Audit Implementation Summary

## Purpose

This note summarizes what the code in `research/pde_audit/` now accomplishes.
It is intended as a paper-writing handoff for the executable/reproducibility
side of the audit.

The math/status notes explain what the audit means.  This note explains what
was implemented, what artifacts are generated, and what the code does and does
not certify.

---

## 1. Top-level referee harness

The main entry point is:

```bash
bash research/pde_audit/run_all.sh
```

It runs the release bundle in this order:

1. fixture manifest verification;
2. Python audit scripts;
3. simulation bundle;
4. root JSON hygiene check;
5. serial Mathematica mirrors.

The combined summary is written to:

```text
research/pde_audit/output/_summary.txt
research/pde_audit/output/_summary.json
```

The last full run reported:

```text
PYTHON: 28/28 passed, 0 failed
SIMULATION: 16/16 passed, 0 failed
MATHEMATICA: 10/10 passed, 0 failed
ROOT_JSON_FILES: 0
REFEREE_PASS: True
```

This harness is a reproducibility check.  It does not prove a completed PDE
branch exists.  It proves that the stated audit scripts, fixtures, mirrors, and
simulation diagnostics execute consistently and produce the recorded verdicts.

---

## 2. Python audit layer

The Python audit layer lives in:

```text
research/pde_audit/scripts/
```

It is run by:

```bash
bash research/pde_audit/scripts/run_all_audits.sh
```

The current suite has 28 passing checks.  It covers:

- parent wall action versus effective linear wall closure;
- Maxwell gauge localization;
- dimension ledger and port/source normalization warnings;
- open finite-radius junction impedance;
- Poisson/Newtonian hook;
- stable BdG support Schur complement;
- Maxwell/mixed kernel;
- Hamiltonian positivity gates;
- grouped real `P2` projector calculus;
- STF angular source-map normalization;
- grouped normalization and constant-prefactor branch algebra;
- outgoing `l=2` fingerprint;
- 2.5PN/4PN interface;
- branch freeze and no-refit protocol;
- weak-axisymmetric splitting;
- monomial quotient and similarity-orbit separation;
- isotropic full-bundle target surface;
- weak-form branch-extraction scaffold;
- branch-extraction fixture;
- profile-to-coefficient adapter;
- solver handoff validator;
- end-to-end smoke pipeline;
- V2-23 formula, mesh convergence, and minimal reduced branch solver;
- V2-24 negative controls.

The layer writes text summaries to:

```text
research/pde_audit/scripts/output/*.txt
research/pde_audit/scripts/output/_summary.txt
research/pde_audit/scripts/output/_summary.json
```

and generated JSON artifacts to:

```text
research/pde_audit/scripts/output/artifacts/
```

Important generated artifacts include:

- V2-21 branch extraction packet;
- V2-22A generated V2-21 manifest and observable packet;
- V2-22B validation reports;
- V2-22C generated profile and V2-21 manifests;
- V2-23 reduced solver branch manifest, observable packet, tolerance report,
  and mesh convergence report;
- V2-24 negative-control report.

---

## 3. Fixture and negative-control layer

Fixtures live in:

```text
research/pde_audit/scripts/fixtures/
```

The fixture manifest is:

```text
research/pde_audit/scripts/fixtures/MANIFEST.json
```

The fixture verifier checks that fixture files are present and hash-stable.

The negative-control suite under

```text
research/pde_audit/scripts/fixtures/negative_controls/
```

tests malformed or invalid solver packets, including:

- bad boundary protocol;
- missing gauge convention;
- nonfinite solver residual;
- nonmonotone grid;
- nonpositive mixed denominator;
- missing pre-target freeze;
- profile length mismatch;
- target-blind false;
- target-output leakage.

The purpose is to show that the handoff pipeline rejects invalid branches
before target residuals are trusted.

---

## 4. Mathematica mirror layer

The Mathematica mirrors live in:

```text
research/pde_audit/mathematica/
```

They are run serially by:

```bash
bash research/pde_audit/mathematica/run_all_audits.sh
```

or through the top-level harness.  Mathematica is run one script at a time with
`math -script` because of licensing constraints.

The current mirror suite has 10 passing checks:

- V2-04 open junction impedance;
- V2-13 grouped normalization ratio;
- V2-16 branch freeze/no-refit;
- V2-17 weak-axisymmetric splitting;
- V2-19 isotropic full-bundle target surface;
- V2-21 branch extraction fixture;
- V2-22A profile-to-coefficient adapter;
- V2-22B solver handoff validator;
- V2-22C end-to-end smoke pipeline;
- V2-23 formula audit.

Outputs are written to:

```text
research/pde_audit/mathematica/output/
```

The mirrors are secondary execution coverage.  They are not claimed as
independent derivations of the entire program; they check load-bearing algebra
and fixture contracts through a separate CAS/runtime.

---

## 5. Simulation bundle

The simulation layer lives in:

```text
research/pde_audit/simulation/
```

It is run by:

```bash
bash research/pde_audit/simulation/run_all.sh
```

The simulation layer is target-blind during candidate generation.  Target
residuals are applied only after packets are frozen.

The current suite has 16 passing checks:

1. `verify_reduced_fem.py`
2. `verify_nonlinear_solver.py`
3. `verify_physical_model.py`
4. `diagnose_notes_intake.py`
5. `generate_nonlinear_packets.py`
6. nonlinear target-blind guard
7. nonlinear frozen sweep evaluator
8. nonlinear obstruction diagnostic
9. nonlinear required-deformation diagnostic
10. `generate_reduced_sweep.py`
11. reduced target-blind guard
12. reduced frozen sweep evaluator
13. reduced obstruction diagnostic
14. reduced required-deformation diagnostic
15. mechanism-gap diagnostic
16. projection-stress diagnostic

The simulation writes:

```text
research/pde_audit/simulation/output/_summary.txt
research/pde_audit/simulation/output/_summary.json
```

plus manifests, packets, candidate reports, and diagnostic reports.

---

## 6. Reduced FEM and frozen packet generation

The reduced FEM primitives are in:

```text
research/pde_audit/simulation/reduced_fem.py
```

`verify_reduced_fem.py` checks:

- matrix structure;
- manufactured D/N half-wave behavior;
- open-shape profile smoke behavior.

`generate_reduced_sweep.py` builds 192 frozen V2-22B-compatible reduced
open-throat packets from a predeclared `operator_v1` protocol.

Generated reduced packet artifacts are written under:

```text
research/pde_audit/simulation/output/packets/
research/pde_audit/simulation/output/manifest.json
```

The reduced packets are then classified post-hoc by `evaluate_frozen_sweep.py`.
The current result is:

```text
0/192 target-passing frozen candidates.
189/192 open and stable candidates.
```

This is an honest target-blind miss, not a pipeline error.

---

## 7. Nonlinear manufactured readiness lane

The nonlinear readiness protocol is implemented in:

```text
research/pde_audit/simulation/nonlinear_protocol.py
research/pde_audit/simulation/verify_nonlinear_solver.py
research/pde_audit/simulation/generate_nonlinear_packets.py
```

`verify_nonlinear_solver.py` checks:

- source import boundary;
- protocol and freeze hashes;
- Jacobian directional consistency;
- manufactured nonlinear mesh convergence;
- continuation sanity.

`generate_nonlinear_packets.py` emits a small target-blind manufactured
nonlinear packet set under:

```text
research/pde_audit/simulation/output/nonlinear_packets/
research/pde_audit/simulation/output/nonlinear_manifest.json
```

The current result is:

```text
0/3 target-passing nonlinear manufactured candidates.
3/3 open and stable candidates.
```

This lane verifies nonlinear mechanics and frozen export plumbing.  It is not
the final physical moving-throat exporter.

---

## 8. Physical-export boundary guard

The strict physical nonlinear exporter boundary is implemented in:

```text
research/pde_audit/simulation/physical_nonlinear_model.py
research/pde_audit/simulation/verify_physical_model.py
```

The current physical-model guard passes by confirming that physical export is
still blocked.

It verifies:

- no banned target-evaluation imports;
- status hash stability;
- strict parent dynamic wall status is not passed;
- effective wall closure is available;
- physical export is not permitted;
- no physical packets or physical manifest were emitted;
- cited ledger blockers and equation-inventory items have source hashes and
  evidence phrases.

The current summary is:

```text
physical_export_permitted: False
packets_emitted: False
blocker_count: 4
```

This guard prevents the manufactured nonlinear lane from being accidentally
described as a true physical branch realization.

---

## 9. Target-blind guards and frozen evaluation

`verify_target_blind.py` checks that packet generators do not import target
evaluation modules and that frozen packet manifests do not contain target
outputs.

The same guard is used for:

- reduced packet generation;
- nonlinear manufactured packet generation.

`evaluate_frozen_sweep.py` runs post-hoc classification through the existing
V2-22B -> V2-22A -> V2-21 chain after packets are frozen.

The target-blind separation is the key anti-overfitting feature of the
simulation code.

---

## 10. Post-hoc diagnostics

The simulation bundle includes post-hoc diagnostics that read already-frozen
candidate reports.  They do not generate, mutate, or refit candidates.

### Obstruction diagnostics

Implemented in:

```text
research/pde_audit/simulation/diagnose_obstruction.py
```

The diagnostic decomposes misses through the one-pole ratio:

```text
D0*C/(3*A^2).
```

Current reduced open-stable range:

```text
min 0.0033775383274364888
median 0.06134471930726503
max 0.1353664855760648
```

### Required-deformation diagnostics

Implemented in:

```text
research/pde_audit/simulation/diagnose_required_deformation.py
```

It reports the coefficient changes that would be required after the fact.
Current reduced open-stable values:

```text
C_multiplier_min 7.387352901601946
C_multiplier_median 16.30132163440465
P0_multiplier_median 171.65261223353198
local_continuation False
```

These are diagnostic numbers, not retuning instructions.

### Mechanism-gap diagnostic

Implemented in:

```text
research/pde_audit/simulation/diagnose_mechanism_gap.py
```

It combines deformation reports with the physical-model guard.  It classifies
the current miss as a large one-pole support deficit and records the physical
requirements for a future exporter.

### Projection-stress diagnostic

Implemented in:

```text
research/pde_audit/simulation/diagnose_projection_stress.py
```

It asks what would happen under post-hoc algebraic coefficient projections.
It records:

```text
one_pole_support_alone_is_insufficient: True
uniform_outgoing_amplitude_scale_is_insufficient: True
target_blind_hit_claimed: False
```

This is the code-level evidence that a future mechanism needs outgoing
moment-shape control, not only scalar support or scalar outgoing amplitude.

---

## 11. Notes-intake guard

The unincorporated-notes intake is implemented in:

```text
research/pde_audit/simulation/diagnose_notes_intake.py
```

It verifies 16 source anchors across:

- 5PN computational handoff;
- 5PN final packet shape;
- Family-1 support/source branch location;
- barrier audit summary;
- strict parent/effective wall status;
- no-refit status firewall;
- `U/V` dressing;
- microscopic export kernel;
- constant-prefactor `N2/N4` target equations;
- support-blind outgoing transfer split;
- leakage/work support-side lane;
- atomic finite-throat `P0/P2/P22` source;
- lepton scalar radiative `P0` flux hook;
- scalar hammer to `P22` veto;
- atom-work intrinsic `P22` bracing;
- simulation coefficient map.

Current output:

```text
pass: True
anchor_count: 16
failed_anchor_count: 0
actual_branch_packet_required: True
retune_current_candidates_allowed: False
outgoing_moment_shape_control_required: True
primary_next_artifact: actual_branch_protocol_v1
```

Generated artifacts:

```text
research/pde_audit/simulation/output/notes_intake_report.json
research/pde_audit/simulation/output/notes_intake_report.txt
```

This guard converts the unincorporated notes into evidence-backed source and
protocol constraints without claiming that they solve the branch.

---

## 12. Referee summary generator

The combined summary generator is:

```text
research/pde_audit/scripts/write_referee_summary.py
```

It reads:

```text
research/pde_audit/scripts/output/_summary.json
research/pde_audit/simulation/output/_summary.json
research/pde_audit/mathematica/output/_summary.json
```

and writes:

```text
research/pde_audit/output/_summary.json
research/pde_audit/output/_summary.txt
```

It also enforces the release hygiene condition that no stray root-level JSON
files remain under `research/pde_audit/`.

Recent fields added to the combined summary include:

- physical model export boundary status;
- notes-intake status;
- nonlinear export status;
- target-blind guard result;
- obstruction and required-deformation diagnostics;
- mechanism-gap classification;
- projection-stress result.

---

## 13. What the code establishes

The code establishes:

- the symbolic/reduced audit scripts execute successfully;
- fixtures and negative controls are stable and enforced;
- target-blind packet generation is separated from target evaluation;
- the current reduced and manufactured nonlinear branch families miss the
  target;
- the miss is quantitatively diagnosed;
- post-hoc coefficient projection is not claimed as evidence;
- strict physical nonlinear export remains blocked;
- unincorporated notes imply actual-branch extraction, not retuning.

The code does not establish:

- that a completed physical PDE branch exists;
- that known physics has been derived from one PDE;
- that the manufactured nonlinear lane is a physical branch;
- that the current miss can be repaired by tuning the existing candidates;
- that the unincorporated notes already provide calibrated `D0/C/P0/N2/N4`
  packets.

---

## 14. Paper-use summary

For paper writing, the implementation can be summarized as:

```text
The audit bundle contains a reproducible Python, Mathematica, fixture, and
simulation harness.  It verifies the stated reduced/effective algebra and
handoff contracts, rejects malformed branch packets, runs target-blind reduced
and manufactured nonlinear branch searches, and records the resulting miss
through obstruction, deformation, mechanism-gap, projection-stress, and
notes-intake diagnostics.  The final referee harness passes, but the current
simulation layer finds no target-passing frozen branch and explicitly blocks
reclassification of the manufactured nonlinear model as a completed physical
moving-throat exporter.
```

The short version:

```text
The code makes the audit reproducible and falsifiable.  It does not complete
the PDE.  It shows exactly where the current executable branch families fail
and what an actual physical exporter must supply next.
```

stage_v2_28_physical_picture_and_ontology.md
# Stage V2-28 - Physical Picture and Ontology Checklist

## Purpose

This note records the physical picture that the papers should state explicitly.
It exists because the model can be misread if the defect is described only by
reduced coefficients or mouth response variables.

The most important point:

```text
The particle is modeled as a finite brane-bulk throat defect: a puncture/open
conduit through the brane into the bulk, not a mere surface depression,
indentation, dimple, or capped pocket.
```

The code/audit variables `D0`, `C`, `P0`, `N2`, and `N4` are low-order
projected response readouts of this object.  They are not the full physical
description of the object.

---

## 1. Core ontology

The model treats a particle as a brane-bulk throat defect.

The brane sees a localized mouth field.  The interior carries throat/cavity
structure extending into the bulk.  Charge, mass, support moments, and outgoing
moments are reduction-layer quantities, not primitive labels attached to a
point particle.

The papers should say this early:

```text
The defect is a finite-radius opening of the brane into a bulk throat.  The
mouth is the brane-side cross-section of that throat.  The interior throat
supports bulk/superfluid, wall, BdG, Maxwell/mixed, and outgoing-port degrees
of freedom.  The reduced particle observables are extracted from this branch.
```

Do not describe the object as:

- a depression in the brane;
- a surface dimple;
- a shallow pocket;
- a closed bubble attached to the brane;
- a hard-capped tube, except when explicitly discussing obsolete toy models or
  negative controls.

---

## 2. Geometry

The effective moving-throat geometry is currently represented as a level-set or
shape-field lift:

```text
Sigma(X,t) = r - R(Omega,w,t),
```

where `w` is the transverse/bulk direction, `Omega` is the angular direction on
the brane-side sphere, and the finite throat surface is

```text
Sigma(X,t) = 0.
```

The convention carried by the notes is:

```text
exterior region: Sigma > 0,
interior/support region: Sigma < 0.
```

This geometry is an effective closure unless/until the parent throat action
`S_Sigma` is promoted and frozen.

The reference stationary branch has:

```text
mouth at w = 0,
finite depth 0 <= w <= L,
mouth radius R(0) = a,
open finite exit R(L) > 0.
```

The open-endpoint correction is load-bearing:

```text
R(0) = a,
R(L) > 0.
```

This is an open finite-radius conduit.  The old hard cap

```text
R(L) = 0
```

is an obsolete toy idealization unless explicitly declared as a negative
control or simplification.

---

## 3. Mouth versus interior

The mouth is the brane-side cross-section of the throat.  It is not the whole
defect.

Important distinctions:

- mouth geometry: radius, ellipticity, headless `P22` axis, area;
- interior geometry: axial throat profile, depth, support/cavity structure;
- open exit: finite-radius outlet into the bulk/reservoir side;
- projected branch observables: reduced coefficients extracted after solving
  or freezing the branch.

The papers should avoid language that collapses these into one scalar "size" or
"depth" unless the reduction being used has explicitly done that.

In particular, the same mouth radius can coexist with different internal
support, throughput, outgoing, and orbit-lock data.  A rigid-mouth reading is a
statement about the brane entrance geometry; it is not automatically a theorem
about every internal transfer or orbit-lock variable.

---

## 4. Open system and superfluid intake/output

The throat is an open-system object in the physical picture.

There are several flux-like quantities that must not be identified by fiat:

1. internal radial throughput entering the self-flow energy ledger;
2. projected brane-bulk exchange source built from transverse bulk current;
3. brane-side mouth-output flux measured through the mouth surface;
4. leakage/work terms through the open finite-radius junction;
5. outgoing/radiative port flux.

The superfluid material state is also not determined by these flux labels
alone.  A completed branch must separately specify or derive the density field
`rho`, EOS/internal energy, local sound speed `c_s(rho)`, and any effective
light-speed relation if the model makes the light cone density dependent.  See
`stage_v2_29_superfluid_material_closure_gap.md`.

The exact D/N trapped support mode is not itself a trans-brane current injector.
On that branch the support field can vanish at the mouth while its mouth
gradient is nonzero.  The first nontrivial mouth datum is therefore a quadratic
normal stress, not the mouth value.

The carried hammer-stress theorem gives a cycle-averaged normal mouth stress
after support-action normalization:

```text
T_nn_bar = pi*hbar*c_s/(2*L^2).
```

But the constitutive response coefficient that turns mouth stress into actual
mouth flux is not yet fully derived from the completed PDE.

A nonzero DC mouth source in the current notes appears only on an explicitly
open/radiative branch.  Closed conservative standing-wave support does not by
itself produce the required DC mouth flux.

---

## 5. Charge and circulation

The charge sign is tied to puncture orientation in the carried ontology:

```text
eta_Q = +/- 1.
```

The microscopic charge branch is

```text
q_* = eta_Q e_*.
```

The observable brane charge is thickness/localization dressed.

The fluxoid/circulation law is a topological law for loops surrounding the
mouth.  It quantizes the tangential vortical class.  It does not determine the
radial feed amplitude or the mouth-output strength by itself.

The papers should keep these separate:

- puncture orientation;
- microscopic charge branch;
- observable dressed charge;
- circulation/fluxoid sector;
- radial feed or throughput amplitude.

For the electromagnetic paper-facing status and claim boundary, see
`stage_v2_30_electromagnetic_ontology_and_status.md`.

---

## 6. Finite-mouth shape physics

Finite throat size is not a cosmetic correction.  It changes the first
shape-sensitive response problem.

For a finite throat, the first shape-sensitive external load is the partner
field Hessian across the mouth, not the raw scalar potential at a point.  This
load decomposes into:

```text
P0 trace/scalar channel,
P2 traceless quadrupole channel.
```

Under centering and area preservation, the first driven non-axisymmetric mouth
deformation is a real headless `P22` mouth ellipse.

The exact constant-area ellipse has semiaxes

```text
R1 = R_th exp(epsilon_ell),
R2 = R_th exp(-epsilon_ell),
R1*R2 = R_th^2,
Area = pi R_th^2.
```

Its first-order boundary perturbation is a real `P22` pair, not a dipole or
free monopole:

```text
delta r_m/R_th = epsilon_ell cos(2(phi - theta_mouth)).
```

The mouth quadrupole tensor is headless.  The axis satisfies

```text
theta_mouth ~ theta_mouth + pi.
```

This `P22` sector is essential for the atom/same-charge bridge notes, but it is
still reduced/conditional when continued to isolated same-charge structure.

---

## 7. Parent status and closure status

The physical picture must preserve the status firewall.

The strict parent action currently contains the GNLS/matter and localized
Maxwell sectors.  The moving throat is parent-level as a confinement-coupling
argument, but it is not yet parent-level as an autonomous dynamical field unless
`S_Sigma` is added to the total action.

Therefore the current wall/throat PDE should be described as:

```text
effective wall/throat closure
```

unless a paper explicitly derives and freezes the promoted parent action.

The reduced wall action and distributed wall equations are useful and
mathematically consistent inside the closure.  They should not be advertised as
already deriving the full parent moving-throat PDE.

---

## 8. What the response coefficients do not describe

The audit readouts

```text
D0, C, P0, N2, N4
```

are not physical ontology variables.  They are low-order projected response
moments extracted from a branch.

They do not directly describe:

- the throat surface `Sigma = 0`;
- the axial profile `R0(w)` or `R0(s)`;
- the mouth geometry;
- the open exit;
- superfluid intake;
- DC leakage;
- outgoing radiation;
- thermal/export partition;
- superfluid density, EOS, sound speed, or density-dependent effective
  light-speed behavior;
- internal support profiles;
- full BdG/Krein signatures;
- Maxwell/mixed mode shapes;
- finite-mouth `P22` orientation;
- branch selection;
- nonlinear saturation.

They are useful because they are the compressed packet that the GR-like audit
can test.  They are not a substitute for the physical mechanism that produces
the packet.

---

## 9. Paper checklist

Every paper in this bundle should make the following physical points explicit
when relevant:

1. The defect is a finite brane-bulk throat/puncture, not a surface depression.
2. The mouth is the brane-side cross-section, not the entire defect.
3. The throat is open in the branch-realization picture: `R(L)>0`.
4. Hard caps are obsolete toy idealizations or negative controls.
5. The wall/throat dynamics are currently an effective closure unless
   `S_Sigma` is promoted.
6. Charge sign, circulation, radial throughput, and mouth output are distinct
   physical quantities.
7. Electromagnetic language should preserve the charge/circulation firewall:
   puncture orientation carries `eta_Q`, while circulation belongs to the
   magnetic/vortical sector.
8. The exact D/N support mode gives a mouth gradient/stress, not automatically
   a mouth current.
9. Nonzero DC mouth output requires open/radiative or other dynamic
   rectification physics.
10. Finite-mouth response starts with Hessian/tidal loading, not point scalar
   depth.
11. The first non-axisymmetric finite-mouth deformation is headless real `P22`.
12. Rigid mouth is a statement about entrance geometry, not automatically about
   all internal transfer variables.
13. Reduced/effective/conditional/strict-parent claims must be labeled.

---

## 10. Recommended paper paragraph

A compact version suitable for an introduction or model section:

```text
In this program a particle is modeled as a finite brane-bulk throat defect: a
finite-radius puncture of the brane into an open bulk conduit, not a surface
depression or capped cavity.  The brane-side mouth is the entrance cross-section
of the throat, while the interior carries support, wall, Maxwell/mixed, and
outgoing-port degrees of freedom.  The current moving-throat equations are an
effective wall/throat closure unless the promoted parent action S_Sigma is
explicitly added.  Observable particle data are extracted only after a branch is
frozen; the audit coefficients D0, C, P0, N2, and N4 are therefore projected
response moments, not the full physical ontology of the defect.
```

---

## 11. Common misreadings to avoid

Misreading:

```text
The defect is just a depression in the brane.
```

Correction:

```text
The defect is a finite-radius throat/puncture through the brane into the bulk.
```

Misreading:

```text
The throat is a capped cavity.
```

Correction:

```text
The branch-realization geometry is an open finite-radius conduit with R(L)>0.
```

Misreading:

```text
The mouth value is the mouth source.
```

Correction:

```text
On the D/N support branch the mouth value can vanish; the first nontrivial
mouth datum is the normal stress from the boundary gradient.
```

Misreading:

```text
P0 scalar output can directly generate P22 shape control.
```

Correction:

```text
The scalar P0 hammer does not linearly drive the area-preserving P22 mouth
mode.  P22 requires genuine finite-mouth shape/bracing physics.
```

Misreading:

```text
D0, C, P0, N2, and N4 are the physical model.
```

Correction:

```text
They are the compressed audit readouts of a branch.  The physical model is the
throat, mouth, support, wall, mixed, flux, and outgoing-port system that must
produce those readouts.
```

stage_v2_29_superfluid_material_closure_gap.md
# Stage V2-29 - Superfluid Material Closure Gap

## Purpose

This note records a material-closure gap exposed by the PDE audit discussion.
It should be carried into the paper draft because it affects whether the final
PDE can be claimed as complete.

The short version:

```text
The notes contain partial superfluid EOS and sound-speed formulas, but the
audited branch machinery does not yet solve a self-consistent superfluid
material sector that determines rho, c_s(rho), and any density-dependent
effective light speed on the same frozen moving-throat branch.
```

This is not a small bookkeeping issue.  If the observable speed of light, the
tail transport speed, or the Maxwell localization speed depends on the
superfluid density, then density is a load-bearing branch variable rather than
a hidden constant.

---

## 1. What is already present

The project has not completely ignored the superfluid equation of state.
Several notes carry a frozen stiff-polytropic closure of the schematic form

```text
P(rho) = K_EOS rho^n
h(rho) = n K_EOS rho^(n-1)/(n-1)
c_s^2(rho) = (1/m) dP/drho.
```

For the commonly carried `n = 5` branch this gives

```text
c_s^2(rho) = (5 K_EOS/m) rho^4.
```

The notes also contain background-density formulas in reduced Thomas-Fermi or
stationary-branch settings, and some 2PN-side material uses normalized
sound-speed functions such as `C_s(u)=c_s^2(u)/c_s0^2`.

So the right statement is not:

```text
There is no EOS anywhere in the project.
```

The right statement is:

```text
The EOS appears as a reduced or frozen constitutive ingredient, but it is not
yet closed into the audited moving-throat branch exporter.
```

---

## 2. What is missing

The current audit and simulation packets mostly treat

```text
rho0, c_s, c, K_EOS, Z(w), W(w)
```

as symbolic constants, fixture inputs, normalization choices, or frozen branch
data.  That is acceptable for checking conditional algebra, but it is not a
completed physical PDE.

A completed branch realization needs equations that determine at least:

1. The superfluid density field:

   ```text
   rho = rho(X,t)
   ```

   including its stationary throat profile and perturbations.

2. The superfluid phase/current sector:

   ```text
   j_rho, v_s, continuity, intake, leakage, and outgoing/export terms.
   ```

3. The EOS or internal energy:

   ```text
   U(rho), P(rho), h(rho), chemical potential.
   ```

4. The local sound speed:

   ```text
   c_s^2(rho) = (1/m) dP/drho
   ```

   or its corrected equivalent if the final EOS is not the frozen
   stiff-polytropic one.

5. The effective light-speed relation, if the model assumes one:

   ```text
   c_eff = c_eff(rho, localization data, Maxwell/mixed branch data).
   ```

6. The way these material fields feed the audit readouts:

   ```text
   D0, A, C, P0, N0, N2, N4.
   ```

7. The way they feed the tail gate:

   ```text
   Theta_tail (c/c_s)^3 = 1
   ```

   or the branch-derived replacement for that gate.

Without those equations, `c_s` and `c` are parameters placed on the branch, not
outputs of the branch.

---

## 3. Why this can block the final PDE

This gap can prevent completion of the one-PDE program.

If density controls the propagation speeds, then a branch that matches the
target coefficients must also satisfy material constraints.  It is not enough
to find values of

```text
D0, C, P0, N2, N4
```

that pass the reduced audit.  The same branch must explain why the required
values of `c`, `c_s`, and `rho` are obtained and why they remain stable enough
to match observed physics.

In particular:

- If `c_eff` depends on `rho`, then density variations would generically change
  the effective light cone unless the branch has a stabilizing or screening
  mechanism.
- If `c_s` depends on `rho`, then the tail-transport gate is also a material
  equation, not just a scalar normalization condition.
- If `rho0` is left free, then source strength, port normalization, throat
  impedance, and outgoing transfer can hide untracked tuning.
- If the throat intake/output changes `rho`, then the open-system flux ledger
  must feed back into the coefficients rather than being appended after the
  coefficient extraction.

This gives a plausible reason why the reduced algebra can be clean while the
current simulations miss the target: the simulations are not solving the
material sector that would select or rule out the needed branch.

---

## 4. Referee-safe claim language

The paper should avoid saying:

```text
The present PDE derives the speed of light and sound speed from the superfluid
density.
```

unless the material sector is actually solved.

A safer statement is:

```text
The current audit treats c, c_s, and rho0 as frozen branch data or reduced
constitutive inputs.  Existing notes contain candidate EOS relations, including
a stiff-polytropic sound-speed law, but the moving-throat exporter has not yet
closed the superfluid material sector that would derive those quantities on the
same branch that produces the audit coefficients.
```

The strongest positive statement presently supported is:

```text
The audit identifies the additional material-closure equations that a completed
PDE must supply before the reduced target algebra can be promoted to a full
physical derivation.
```

---

## 5. Recommended next work package

The next derivation package should be a dedicated superfluid material closure.
It should produce a frozen branch packet only after the following are fixed
before target comparison:

```text
EOS choice or derivation
rho0(X) stationary throat profile
c_s(rho0) and perturbative c_s response
c_eff(rho0) or proof that c is density-insensitive in the relevant sector
intake/output continuity law
source and port normalization from rho and current data
tail-gate transport factor from branch data
effect of material response on D0, A, C, P0, N0, N2, N4
```

The important no-refit rule still applies: this material packet must be frozen
before comparing to the GR-like target surface.  Otherwise the material sector
would become a new tuning layer rather than a physical branch derivation.

---

## 6. Practical interpretation of the simulation miss

The current simulation miss should not be interpreted as proving the whole
program false.  It should also not be ignored.

The most accurate interpretation is:

```text
The reduced/manufactured branch families audited so far do not contain the
needed target branch.  One major missing ingredient is a self-consistent
superfluid material sector that determines density, sound speed, effective
light-speed behavior, and open-system flux feedback together with the throat
response coefficients.
```

This is exactly the kind of gap that can make a final PDE substantially harder
than the reduced algebra suggested.

stage_v2_30_electromagnetic_ontology_and_status.md
# Stage V2-30 - Electromagnetic Ontology and Status

## Purpose

This note records what the PDE/audit stack currently says about
electromagnetism, charge, circulation, and puncture orientation.

The short version:

```text
The program does contain an electromagnetic sector: a localized Maxwell field,
gauge-localization audit, mixed brane-bulk gauge invariants, charge variables,
and a circulation/fluxoid sector.  But electric charge and circulation are not
the same variable in the corrected ontology.
```

This matters for paper writing because older or informal language can make it
sound as if "charge = circulation" or "the puncture itself = electric charge."
The current stack is more precise.

---

## 1. What EM structure is present

The parent/reduced PDE materials include a Maxwell field

```text
A_M = (A_0, A_i, A_w),
F_MN = partial_M A_N - partial_N A_M,
```

with a localized kinetic weight `Z(w)` and source coupling to `J^M`.  The V2-02
audit checks the gauge-localization issue:

```text
S_EM ~ integral [-Z(w) F_MN F^MN/(4 mu0) - A_M J^M + gauge fixing].
```

The V2-09 audit then checks a reduced Maxwell/mixed block with variables

```text
Q, U, W.
```

Here:

- `Q` is the wall/worldtube amplitude;
- `U` is the localized brane-like Maxwell coordinate;
- `W` is the mixed `A_w/F_{mu w}/J^w` active coordinate.

The mixed gauge-invariant observables are

```text
E_w = -partial_t A_w - partial_w A_0,
C_a = partial_a A_w - partial_w A_a.
```

These are not gauge artifacts.  They are the tensor slots where brane-bulk
electromagnetic exchange, mixed-sector plumbing, and hidden-port response can
enter.

---

## 2. Electric charge variables

The corrected charge dictionary uses

```text
eta_Q = +/- 1,
q_* = eta_Q e_*,
q_eff = q_*/sqrt(Z_int).
```

Interpretation:

- `eta_Q` is the electric branch sign, tied to puncture orientation in the
  carried ontology;
- `q_*` is the microscopic branch charge;
- `q_eff` is the brane-observed/localization-dressed charge after zero-mode
  normalization;
- `Z_int = integral Z(w) dw` is the localization normalization entering the
  effective brane charge.

So the paper-safe phrasing is:

```text
The electric charge sign is carried by the oriented puncture branch eta_Q, and
the observed charge is the localization-dressed q_eff.
```

Do not write the stronger shorthand:

```text
The puncture is electric charge.
```

The puncture/throat is the physical defect.  Electric charge is one branch
label and coupling carried by that defect.

---

## 3. Circulation and magnetism

The corrected ontology keeps circulation separate from electric charge.

The circulation/fluxoid sector is the topological law for loops surrounding the
mouth.  In the charged superfluid notation it has the schematic form

```text
integral_C (partial_i theta - q_* A_i/hbar) dl^i = 2 pi n.
```

This law quantizes a tangential vortical/holonomy class.  It belongs to the
magnetic/vortical sector, not to the electric-charge dictionary.

Paper-safe phrasing:

```text
Circulation around or through the throat belongs to the magnetic/vortical
sector.  It is the natural place to encode magnetic/fluxoid topology, while
electric charge remains the separate puncture-orientation branch eta_Q with
microscopic charge q_*.
```

Avoid the stronger claim:

```text
Circulation is electric charge.
```

Also avoid claiming that the audit already derives all of magnetism from throat
intake circulation.  The files contain the Maxwell field, fluxoid/circulation
law, and mixed brane-bulk slots, but the full recirculation/plumbing law is not
closed yet.

---

## 4. What is reduced versus still open

The following pieces are present and audit-backed:

- localized Maxwell action and zero-mode normalization issues;
- finite gauge-localization warnings and patches;
- mixed gauge invariants `E_w` and `C_a`;
- reduced Maxwell/mixed kernel stability and transfer gates;
- conservative Maxwell/mixed moments `Z0,Z2,Z4`;
- charge dictionary `eta_Q`, `q_*`, `q_eff`;
- circulation/fluxoid sector as magnetic/vortical, not electric-charge
  defining.

The following pieces remain open:

- a complete moving-throat derivation of the source current `J^M`;
- a closed recirculation/plumbing law connecting exterior circulation, throat
  intake, mixed `A_w/F_{mu w}/J^w` transport, and brane magnetic fields;
- a derivation of the charge magnitude `e_*` rather than treating it as a
  branch parameter;
- a proof that the localized Maxwell sector itself emerges from the superfluid
  variables rather than being included as a parent sector;
- a completed same-charge internal twist/spin-like discretizer, if the lepton
  branch needs one.

---

## 5. Recommended paper paragraph

```text
The electromagnetic sector in the current PDE stack is represented by a
localized Maxwell field A_M with transverse localization profile Z(w), together
with mixed brane-bulk gauge-invariant components such as E_w and C_a.  Electric
charge is not identified with circulation.  The corrected charge dictionary
assigns the electric branch sign to the oriented puncture label eta_Q, with
microscopic charge q_* = eta_Q e_* and brane-observed charge
q_eff = q_*/sqrt(Z_int).  Circulation and fluxoid quantization around the mouth
belong to the magnetic/vortical sector and must be kept distinct from radial
throughput and mouth-output flux.  The present audit verifies the reduced
Maxwell/mixed algebra and gauge-localization consistency, but a full
moving-throat recirculation law deriving all electromagnetic response from the
superfluid branch remains future work.
```

---

## 6. Minimal claim boundary

Supported:

```text
The PDE program includes a localized Maxwell/mixed sector and a corrected
charge/circulation ontology.
```

Conditionally supported:

```text
Puncture orientation carries the electric charge sign eta_Q.
Circulation/fluxoid data belong to the magnetic/vortical sector.
```

Not yet supported as a completed theorem:

```text
The current PDE derives the whole of electromagnetism from superfluid
circulation entering the throat.
```

