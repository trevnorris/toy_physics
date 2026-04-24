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
