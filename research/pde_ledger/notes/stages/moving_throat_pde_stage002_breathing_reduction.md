# Moving-Throat PDE — Phase 2A Breathing Reduction Back to the Old `(a,L)` Closure

## Purpose
This note is the first continuation after the Phase-1 scaffold and the finite-throat D/N benchmark.
Its job is to show that the new distributed wall field is not floating free from the old hierarchy.
At the conservative linear level, the axisymmetric part of the distributed wall theory reduces back to the same kind of matrix system that the parent 4D files already use for the effective throat variables `a(t)` and `L(t)`.

That reduction matters because it tells us we have not broken continuity with the existing program.
The distributed wall field is a **lift** of the old closure, not a replacement that forgets it.

---

## 1. Axisymmetric wall ansatz

Use the normalized real harmonic
\[
Y_{00}(\Omega)=\frac{1}{2\sqrt\pi},
\qquad
\int_{S^2}Y_{00}^2\,d\Omega=1,
\qquad
\frac{1}{4\pi}\int_{S^2}Y_{00}\,d\Omega=\frac{1}{2\sqrt\pi}.
\]
So if the physical mouth-average radius shift is `delta a(t)`, the corresponding monopole harmonic coefficient is
\[
q_{00}(0,t)=2\sqrt\pi\,\delta a(t).
\]
This is the normalization bridge between the distributed wall field and the old collective coordinate.

Take the axisymmetric wall perturbation to be the two-mode truncation
\[
\eta_{0}(w,t)=2\sqrt\pi\Big[\alpha_a(w)\,\delta a(t)+\alpha_L(w)\,\delta L(t)\Big],
\]
where:
- `alpha_a(w)` is the profile that changes the mouth radius,
- `alpha_L(w)` is the profile that changes the axial throat extent,
- higher axisymmetric wall modes are deferred.

This is the minimal truncation that can recover the old `(a,L)` closure.

---

## 2. Insert the two-mode ansatz into the quadratic wall action

On the axisymmetric branch, the quadratic wall action is
\[
S^{(2)}_{\eta,0}
=
\frac12\int dt\,dw\,d\Omega\,
\Big[
\mu_\eta(w)(\partial_t\eta_0)^2
-
T_w(w)(\partial_w\eta_0)^2
-
K_0(w)\eta_0^2
\Big],
\]
with
\[
K_0(w)=K_\eta(w).
\]
As in Stage 1, the remaining axial coefficients are written in the densitized
convention where the background surface measure has already been absorbed into
the effective one-dimensional fields and coefficients. So the reduced overlap
integrals below are taken with respect to `dw`, not with an additional explicit
`sqrt(gamma_0)` factor.
Using the two-mode ansatz and the normalized harmonic integral
\[
\int_{S^2}Y_{00}^2\,d\Omega=1,
\]
the reduced Lagrangian becomes
\[
L_{\rm red}^{(0)}
=
\frac12 M_{AB}\,\dot Q^A\dot Q^B
-
\frac12 K_{AB}\,Q^A Q^B,
\qquad
Q^A=(\delta a,\delta L),
\]
with the effective matrices
\[
M_{AB}=
4\pi\int dw\,\mu_\eta(w)\,\alpha_A(w)\alpha_B(w),
\]
\[
K_{AB}=
4\pi\int dw\,
\Big[
T_w(w)\,\alpha_A'(w)\alpha_B'(w)
+
K_0(w)\,\alpha_A(w)\alpha_B(w)
\Big].
\]
Indices `A,B` run over `a,L`.

So the distributed wall does exactly what it should do: it produces effective inertia and stiffness matrices by overlap integrals of the wall profiles.

---

## 3. Euler–Lagrange reduction

The reduced equations of motion are
\[
M_{AB}\,\ddot Q^B+K_{AB}\,Q^B=0.
\]
This is the conservative linearized version of the old geometry equations. In the parent 4D files the geometry sector is written schematically as
\[
M_{AB}\,\ddot Q^B+\Gamma_{AB}\,\dot Q^B=-\frac{\partial H_{\rm tot}}{\partial Q^A}.
\]
So the new distributed-wall reduction reproduces the old closure at the expected level:
- the conservative part comes directly from the wall action,
- any damping matrix `Gamma_AB` is an effective/open-system completion that would enter after coupling to matter/gauge/exterior channels.

This is exactly the relationship we wanted. The old `(a,L)` geometry equations are the lowest-mode truncation of the new distributed wall theory.

---

## 4. What happens to the grouped real `P2` sector at the same level

Now take one grouped real quadrupole component,
\[
\eta_{2m}(\Omega,w,t)=\beta_2(w)\,q_{2m}(t)\,Y_{2m}^{\rm real}(\Omega),
\]
with
\[
-\Delta_{S^2}Y_{2m}^{\rm real}=6Y_{2m}^{\rm real}.
\]
Using the same quadratic wall action, the reduced one-mode Lagrangian is
\[
L_{2m}
=
\frac12 M_2\,\dot q_{2m}^2
-
\frac12 K_2\,q_{2m}^2,
\]
with
\[
M_2=
\int dw\,\mu_\eta(w)\,\beta_2(w)^2,
\]
\[
K_2=
\int dw\,
\Big[
T_w(w)\,\beta_2'(w)^2
+
\bigl(K_\eta(w)+6T_\Omega(w)\bigr)\beta_2(w)^2
\Big].
\]
So before any symmetry breaking or matter/gauge coupling, every real `P2` component has the same uncoupled operator
\[
M_2\,\ddot q_{2m}+K_2\,q_{2m}=0.
\]
That is the microscopic reason the grouped real `P2` block is degenerate on the isotropic reference throat.

---

## 5. Why this reduction matters for the roadmap

This reduction gives three concrete answers.

First, it proves that the distributed wall lift is compatible with the old geometry sector. The old `a,L` equations are not being discarded; they are being reinterpreted as the lowest axisymmetric collective truncation.

Second, it shows that the grouped real `P2` bundle is on exactly the same footing as the old scalar geometry variables. It is not an artificial add-on. It is the next harmonic family of the same wall PDE.

Third, it tells us what the next coupled calculation has to do. Once the BdG matter sector, the localized Maxwell sector, and the mixed channels are turned back on, they must:
- renormalize the effective matrices,
- split or preserve the `P2` degeneracy,
- shift the pole data,
- and produce the passive/outgoing odd parts that the uncoupled wall theory cannot generate on its own.

---

## 6. Script-backed status

The accompanying SymPy audit verifies the concrete algebraic claims used here:
- the real-harmonic normalization and zero-average facts,
- the axisymmetric mouth-average extraction rule,
- the reduced two-mode Euler–Lagrange matrix form,
- and the one-mode grouped-`P2` reduction.

Supporting file:
- `scripts/moving_throat_pde_stage002_breathing_reduction_sympy_audit.py`
- `mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl`
  mirrors the same symbolic reduction in a second CAS so the `Y_00` bridge, the
  `4\pi` overlap factor, the conservative `(a,L)` matrix system, and the
  grouped-real `P2` degeneracy are all dual-CAS checked.
