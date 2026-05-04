# Parent Throat Action — Minimal Nonlinear Candidate and Quadratic Recovery

## Purpose

This note continues the parent-action branch after the wall-action audit.
The goal is not to claim a unique final throat action. The goal is narrower:
write the **smallest parent-complete throat action candidate** that

1. turns the throat into an autonomous field at parent level,
2. reduces exactly to the audited quadratic wall action around a stationary isotropic branch,
3. preserves the old `(a,L)` reduction as the lowest axisymmetric truncation,
4. keeps the grouped real `P2` lane literal.

The project files already fixed the status split:

- the present parent action with `V_conf(X;Sigma)` gives a **wall force**, not an autonomous wall PDE;
- the quadratic wall action `S_eta^(2)` is a mathematically consistent **effective linear wall closure**;
- the moving-throat program wants a promoted throat action `S_Sigma` whose quadratic limit reproduces that closure.

This note supplies the minimal gauge-fixed candidate that does exactly that.

---

## 1. Candidate parent-complete throat action

Use the same shape-field language already adopted by the moving-throat program:

\[
\Sigma(X,t)=r-R(\Omega,w,t),
\qquad
r=\sqrt{x^2+y^2+z^2}.
\]

Take the throat field itself to be the graph variable `R(Ω,w,t)` and add a
new autonomous action

\[
S_\Sigma[R]
=
\int dt\,dw\,d\Omega\;
\mathcal L_\Sigma,
\]

with the minimal gauge-fixed density

\[
\boxed{
\mathcal L_\Sigma
=
\frac12\,\mu_\Sigma(R,w)\,(\partial_t R)^2
-
\frac12\,T_{w,\Sigma}(R,w)\,(\partial_w R)^2
-
\frac12\,T_{\Omega,\Sigma}(R,w)\,|\nabla_\Omega R|^2
-
U_\Sigma(R,w).
}
\]

This is the smallest parent-complete promotion in the current throat
coordinates. It is **not** the most geometric possible choice. A later,
stronger upgrade could replace it by a reparameterization-invariant worldvolume
form built from the induced throat metric and extrinsic curvature invariants.
But that extra structure is not needed to recover the audited linear wall PDE.

The total parent-level action becomes

\[
\boxed{
S_{\rm total}=S_\psi[\psi,A;\Sigma]+S_{\rm EM}[A]+S_\Sigma[R].
}
\]

---

## 2. Exact Euler–Lagrange equation of the promoted throat field

Varying `R` gives the autonomous nonlinear throat equation

\[
\partial_t\!\bigl(\mu_\Sigma R_t\bigr)
-
\partial_w\!\bigl(T_{w,\Sigma} R_w\bigr)
-
\nabla_\Omega\!\cdot\!\bigl(T_{\Omega,\Sigma}\nabla_\Omega R\bigr)
-
\frac12\,\partial_R\mu_\Sigma\;R_t^2
+
\frac12\,\partial_R T_{w,\Sigma}\;R_w^2
+
\frac12\,\partial_R T_{\Omega,\Sigma}\;|\nabla_\Omega R|^2
+
\partial_R U_\Sigma
=
\mathcal S_\Sigma[\psi,A;R].
\]

Here `\mathcal S_\Sigma` is the source obtained by varying the other sectors
with respect to `R`. In particular, the promoted confinement coupling

\[
V_{\rm conf}(X;\Sigma)=V_{\rm wall}(\Sigma/\ell_c)
\]

already produces the linearized source used in the Stage-1 wall/BdG reduction,
so the new action is compatible with the existing coupling route.

The updated SymPy script now derives this Euler-Lagrange operator directly from
the nonlinear density by `sympy.calculus.euler.euler_equations` in a local
orthonormal angular chart and then matches it back to the boxed coordinate-free
form written above.

This is the missing parent-level feature that the audit said the current
confinement-only action does not supply: the throat now has its **own** kinetic,
axial-stiffness, angular-stiffness, and restoring terms before any matter or
gauge reduction is performed.

---

## 3. Static isotropic background equation

Take a stationary isotropic branch

\[
R(\Omega,w,t)=R_0(w).
\]

Then the background equation from `S_\Sigma` alone is

\[
\boxed{
\frac{d}{dw}\Bigl(T_{w,\Sigma}(R_0,w)R_0'\Bigr)
-
\frac12\,\partial_R T_{w,\Sigma}(R_0,w)(R_0')^2
-
\partial_R U_\Sigma(R_0,w)
=0,
}
\]

before the matter/gauge source correction is reinserted.

In the audited SymPy derivation this comes from the explicit linear
integration-by-parts identity
\[
-T_w(R_0,w)R_0'\eta_w
=
\partial_w\!\bigl[-T_w(R_0,w)R_0'\eta\bigr]
+
\partial_w\!\bigl(T_w(R_0,w)R_0'\bigr)\eta,
\]
so the boundary term \([ -T_w(R_0,w)R_0'\eta ]_{\partial w}\) is now carried
explicitly rather than being silently discarded. The script also evaluates that
boundary term on a concrete decaying Gaussian test profile by explicit
`sp.limit`, verifying that the discharge really happens in the example rather
than only in formal notation.

So the background throat profile is no longer a purely external constitutive
input. It is the stationary branch of a genuine throat action.

---

## 4. Quadratic expansion around the stationary branch

Now write

\[
R(\Omega,w,t)=R_0(w)+\eta(\Omega,w,t).
\]

Expanding the candidate action to quadratic order and using the background
stationarity equation gives the fluctuation density

\[
\mathcal L_\Sigma^{(2)}
=
\frac12\,\mu_\eta(w)\,\eta_t^2
-
\frac12\,T_w(w)\,\eta_w^2
-
\frac12\,T_\Omega(w)\,|\nabla_\Omega \eta|^2
-
\frac12\,K_\eta(w)\,\eta^2,
\]

with the exact identifications

\[
\mu_\eta(w)=\mu_\Sigma(R_0,w),
\qquad
T_w(w)=T_{w,\Sigma}(R_0,w),
\qquad
T_\Omega(w)=T_{\Omega,\Sigma}(R_0,w),
\]

and

\[
\boxed{
K_\eta(w)
=
\partial_R^2 U_\Sigma(R_0,w)
-
\frac{d}{dw}\Bigl(\partial_R T_{w,\Sigma}(R_0,w)\,R_0'\Bigr)
+
\frac12\,\partial_R^2 T_{w,\Sigma}(R_0,w)\,(R_0')^2.
}
\]

This is exactly the form that the wall-action audit needed.

The cross term is handled by the explicit identity
\[
-T_{w,\Sigma,R}(R_0,w)R_0'\eta\eta_w
=
\partial_w\!\Bigl[-\frac12\,T_{w,\Sigma,R}(R_0,w)R_0'\eta^2\Bigr]
+
\frac12\,\partial_w\!\bigl(T_{w,\Sigma,R}(R_0,w)R_0'\bigr)\eta^2,
\]
so the quadratic boundary term
\[
\Bigl[-\frac12\,T_{w,\Sigma,R}(R_0,w)R_0'\eta^2\Bigr]_{\partial w}
\]
is also tracked explicitly before the bulk \(K_\eta\) formula is read off, and
the script again checks its vanishing on a concrete decaying profile by
explicit `sp.limit`. The boundary-value helper is separately checked on
\(\arctan w\), where the endpoint discharge is nonzero, so these vanishing
boundary claims are tied to the chosen decaying profiles rather than to a
degenerate boundary operator. The script also repeats the linear and quadratic
boundary-discharge checks on a second non-Gaussian decaying profile
\(\eta(w)=1/(1+w^2)\). For the linear check it uses a separate
\(B(w)=e^{-w^2}\) coefficient so the Lorentzian denominator is not cancelled;
for the quadratic check the rational decay is retained against the audited
\(A(w)\) coefficient. The script also tests a finite-endpoint Lorentzian probe
with boundary discharge \(-2\), so the boundary operator distinguishes
nonzero finite limits instead of only returning zero on decaying profiles. The
boundary discharge is therefore not supported only by one Gaussian test
function.

So the audited linear closure

\[
S_\eta^{(2)}
=
\frac12\int dt\,dw\,d\Omega
\Bigl[
\mu_\eta\eta_t^2
-
T_w\eta_w^2
-
T_\Omega\eta(-\Delta_{S^2})\eta
-
K_\eta\eta^2
\Bigr]
\]

is the quadratic limit of the promoted nonlinear `S_\Sigma`, not a bolt-on
spring law.

---

## 5. Modal split and recovery of the audited wall PDE

Projecting onto spherical harmonics gives

\[
\mu_\eta q_{\ell m,tt}
-
\partial_w\bigl(T_w q_{\ell m,w}\bigr)
+
\bigl[K_\eta+\ell(\ell+1)T_\Omega\bigr]q_{\ell m}
=
S_{\ell m}.
\]

Therefore:

- the scalar lane `l=0` has restoring operator `K_\eta`,
- the grouped real `P2` lane `l=2` has restoring operator `K_\eta+6T_\Omega`.

So the parent promotion automatically preserves the wall audit’s `l=0` /
`l=2` split.

The audit script now verifies this in a fused way. It still checks the genuine
\(Y_{20}\) spherical-Laplacian eigenvalue, but it also inserts
\(R(t,w,\Omega)=q(t,w)Y_{20}(\Omega)\) directly into the angular part of the
quadratic action, performs the \(S^2\) integral, and derives the reduced modal
Euler-Lagrange equation
\[
\mu_\eta q_{tt} - \partial_w(T_w q_w) + (K_\eta + 6T_\Omega)q = 0
\]
from that projected Lagrangian itself.

---

## 6. Recovery of the old `(a,L)` closure

Using the same axisymmetric two-profile truncation already used in the moving-throat notes,

\[
\eta_0(w,t)=2\sqrt\pi\bigl[\alpha_a(w)\,\delta a(t)+\alpha_L(w)\,\delta L(t)\bigr],
\]

the reduced matrices are

\[
M_{AB}=4\pi\int dw\,\mu_\eta\,\alpha_A\alpha_B,
\qquad
K_{AB}=4\pi\int dw\,[T_w\alpha_A'\alpha_B'+K_\eta\alpha_A\alpha_B].
\]

So the old geometry sector remains what the program wanted it to be:
not the fundamental throat law, but the lowest axisymmetric collective
truncation of the promoted parent throat action.

---

## 7. Grouped real `P2` branch data from the promoted action

For one grouped real quadrupole profile

\[
\eta_{2m}(\Omega,w,t)=\beta_2(w)\,q_{2m}(t)\,Y^{\rm real}_{2m}(\Omega),
\]

the parent-complete wall contribution becomes

\[
M_2=\int dw\,\mu_\eta\,\beta_2^2,
\qquad
K_2=\int dw\,[T_w(\beta_2')^2+(K_\eta+6T_\Omega)\beta_2^2].
\]

This matters because the isotropic full-bundle target surface later uses

\[
D_0=K-B_0-Z_0,
\qquad
D_2=-(M+B_2+Z_2),
\qquad
D_4=-(B_4+Z_4),
\]

and the promotion of `S_\Sigma` turns `K` and `M` from pure closure knobs into
actual branch integrals of a parent-level throat field.

The target equations themselves do **not** change. What changes is the meaning
of the wall entries inside them.

---

## 8. Positivity / stability gates

At quadratic order the local wall positivity gates remain

\[
\mu_\eta>0,
\qquad
T_w>0,
\qquad
K_\eta+\ell(\ell+1)T_\Omega\ge 0.
\]

So the promotion preserves the audit’s linear stability requirements while
removing the status gap about where those coefficients come from.

---

## 9. Honest status statement

This candidate does **not** solve the full PDE branch-realization problem.
It does something more specific:

1. it closes the strict parent/effective status gap identified in the audit,
2. it shows that the linear wall action can be the quadratic limit of a real parent throat field,
3. it keeps the existing moving-throat reduction stack intact.

What it does **not** do by itself is fix the support, Maxwell/mixed, or outgoing
branch data. Those still have to come from the coupled branch realization.

So the right reading is:

> a promoted `S_\Sigma` is feasible right now at the reduced derivation level,
> and the smallest gauge-fixed nonlinear candidate already reproduces the audited
> wall PDE exactly at quadratic order.
