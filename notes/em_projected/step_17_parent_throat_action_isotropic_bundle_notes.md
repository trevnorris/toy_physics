# Parent Throat Action — Isotropic Full-Bundle Bridge

## Purpose

This note splices the promoted parent throat action directly into the isotropic grouped-`P2` bundle formulas.
The aim is to rewrite the existing one-pole and outgoing-normalization targets so that the wall inputs come from the parent throat field rather than from abstract wall coefficients.

---

## 1. Parent wall inputs from the promoted throat action

Start from the promoted nonlinear throat action

\[
S_{\rm total} = S_\psi[\psi,A;\Sigma] + S_{\rm EM}[A] + S_\Sigma[R],
\qquad
\Sigma = r - R(\Omega,w,t),
\]

with quadratic wall limit around an isotropic background `R_0(w)`

\[
\mathcal L_\Sigma^{(2)}
=
\frac12\mu_\eta(w)\eta_t^2
-
\frac12 T_w(w)\eta_w^2
-
\frac12 T_\Omega(w)|\nabla_\Omega\eta|^2
-
\frac12 K_\eta(w)\eta^2.
\]

The grouped real `P2` wall mode is written as

\[
\eta_{2m}(\Omega,w,t) = \beta_2(w)\,q_{2m}(t)\,Y^{\rm real}_{2m}(\Omega).
\]

The parent wall block therefore contributes the two exact branch integrals

\[
\boxed{
M_\Sigma = \int dw\,\mu_\eta(w)\,\beta_2(w)^2,
}
\]

\[
\boxed{
K_\Sigma = \int dw\,\Big[T_w(w)\,(\beta_2'(w))^2 + \bigl(K_\eta(w)+6T_\Omega(w)\bigr)\beta_2(w)^2\Big].
}
\]

Using the parent-candidate formula

\[
K_\eta
=
U_{\Sigma,RR}(R_0,w)
-
\partial_w\!\bigl(T_{w,\Sigma,R}(R_0,w)R_0'\bigr)
+
\frac12 T_{w,\Sigma,RR}(R_0,w)(R_0')^2,
\]

this becomes the explicit parent-complete stiffness integral

\[
K_\Sigma
=
\int dw\,\Big[
T_{w,\Sigma}(R_0,w)(\beta_2')^2
+
\Big(U_{\Sigma,RR}-\partial_w(T_{w,\Sigma,R}R_0')+\frac12 T_{w,\Sigma,RR}(R_0')^2+6T_{\Omega,\Sigma}\Big)\beta_2^2
\Big].
\]

So the wall input to the bundle is fully determined by the chosen stationary branch `R_0` and wall profile `\beta_2`.

---

## 2. Full isotropic bundle with parent wall inputs

Carry the already-frozen support and conservative mixed moments

\[
B_0,B_2,B_4,
\qquad
Z_0,Z_2,Z_4,
\qquad
N_0,N_2,N_4.
\]

Then the isotropic denominator moments become

\[
\boxed{
D_0 = K_\Sigma - B_0 - Z_0,
\qquad
D_2 = -(M_\Sigma + B_2 + Z_2),
\qquad
D_4 = -(B_4 + Z_4).
}
\]

The normalized conservative response is

\[
Y(\omega) = \frac{D_0}{D_0 + D_2\omega^2 + D_4\omega^4 + O(\omega^6)}
= 1 + u_2\omega^2 + u_4\omega^4 + O(\omega^6),
\]

with

\[
\boxed{u_2 = -\frac{D_2}{D_0}},
\qquad
\boxed{u_4 = \frac{D_2^2 - D_0D_4}{D_0^2}}.
\]

The outgoing prefactor is

\[
P(\omega) = \frac{D_0\,(N_0+N_2\omega^2+N_4\omega^4+\cdots)}{(D_0+D_2\omega^2+D_4\omega^4+\cdots)^2}
= P_0 + P_2\omega^2 + P_4\omega^4 + O(\omega^6),
\]

with

\[
\boxed{P_0 = \frac{N_0}{D_0}},
\]

\[
\boxed{P_2 = \frac{D_0N_2 - 2D_2N_0}{D_0^2}},
\]

\[
\boxed{P_4 = \frac{D_0^2N_4 - 2D_0(D_2N_2 + D_4N_0) + 3D_2^2N_0}{D_0^3}}.
\]

---

## 3. Exact one-pole condition in parent-complete form

The conservative one-pole condition is

\[
u_4 = 4u_2^2.
\]

Substituting the exact response moments gives

\[
\boxed{
D_0\,(B_4+Z_4) = 3\,(M_\Sigma + B_2 + Z_2)^2.
}
\]

Equivalently,

\[
\boxed{
K_\Sigma
=
B_0+Z_0 + \frac{3\,(M_\Sigma+B_2+Z_2)^2}{B_4+Z_4}.
}
\]

So the one-pole surface fixes the required parent wall stiffness once the parent wall inertia and the support/mixed bundle are known.

One may also solve for the parent wall inertia,

\[
\boxed{
M_\Sigma
=
\sqrt{\frac{(K_\Sigma-B_0-Z_0)(B_4+Z_4)}{3}} - (B_2+Z_2),
}
\]

with the opposite algebraic root also solving the one-pole equation. The
positive square-root branch is therefore not fixed by one-pole algebra alone.
The updated script now makes the additional sign criterion explicit:
on the one-pole surface
\[
u_2=\frac{M_\Sigma+B_2+Z_2}{D_0}
\]
takes the values
\[
u_2^{(+)}=+\frac{1}{D_0}\sqrt{\frac{D_0(B_4+Z_4)}{3}},
\qquad
u_2^{(-)}=-\frac{1}{D_0}\sqrt{\frac{D_0(B_4+Z_4)}{3}}.
\]
So under the stable-pole sign convention \(D_0>0\) and \(B_4+Z_4>0\), the
positive branch is exactly the one with \(u_2>0\), while the negative branch
gives \(u_2<0\). The script still verifies the full factorization and Vieta
identities for the two algebraic roots, but the branch choice is now tied to a
stated response-sign criterion rather than left as an unexamined convention.
It also evaluates the two branches at a pinned stable-pole sample,
\(D_0=10\) and \(B_4+Z_4=12\), giving
\[
u_2^{(+)}\approx 0.6324555320336759,
\qquad
u_2^{(-)}\approx -0.6324555320336759.
\]
The self-test repeats the same sign check on three stable-pole samples,
including small-\(D_0\)/large-tail and large-\(D_0\)/small-tail regimes. That
numerical sign check is the load-bearing branch-sign guard; the symbolic
factorization and Vieta identities remain supporting algebra.

---

## 4. Exact outgoing normalization condition

The universal grouped-`P2` outgoing target is

\[
\widehat m_0^2 P_0 = \frac{54Gc_s^5}{5a^5c^5}.
\]

Using `P_0 = N_0/D_0`, the isotropic normalization condition becomes

\[
\boxed{
\widehat m_0^2\,\frac{N_0}{K_\Sigma-B_0-Z_0}
=
\frac{54Gc_s^5}{5a^5c^5}.
}
\]

Equivalently,

\[
\boxed{
K_\Sigma
=
B_0+Z_0 + \frac{N_0}{P_{0,\rm target}},
\qquad
P_{0,\rm target} = \frac{54Gc_s^5}{5a^5c^5\,\widehat m_0^2}.
}
\]

So the isotropic branch fixes the same wall stiffness in a second way.

---

## 5. Exact compatibility equation

Because the one-pole surface and the outgoing normalization target both fix the same quantity `K_\Sigma`, the parent-complete branch must satisfy the compatibility identity

\[
\boxed{
\frac{N_0}{P_{0,\rm target}}
=
\frac{3\,(M_\Sigma+B_2+Z_2)^2}{B_4+Z_4}.
}
\]

Substituting the wall inertia integral gives the most explicit parent-complete form:

\[
\boxed{
\frac{N_0}{P_{0,\rm target}}
=
\frac{3\,\Big(\int dw\,\mu_\eta\beta_2^2 + B_2 + Z_2\Big)^2}{B_4+Z_4}.
}
\]

So the isotropic parent-throat problem is no longer “pick `K` and `M` until the target works.”
It is “does the actual throat branch produce wall integrals that satisfy this equation together with the support and mixed moments on the same frozen branch?”

---

## 6. Constant-prefactor branch conditions

If the outgoing branch is constant-prefactor through `O(\omega^4)`, then

\[
P_2 = 0,
\qquad
P_4 = 0.
\]

The exact conditions are

\[
\boxed{
N_2 = 2\,\frac{D_2N_0}{D_0} = -2\,\frac{(M_\Sigma+B_2+Z_2)N_0}{K_\Sigma-B_0-Z_0},
}
\]

\[
\boxed{
N_4 = \frac{2D_0(D_2N_2 + D_4N_0) - 3D_2^2N_0}{D_0^2}.
}
\]

On the isotropic one-pole surface, this simplifies to

\[
\boxed{
N_4 = -5\,\frac{(M_\Sigma+B_2+Z_2)^2}{(K_\Sigma-B_0-Z_0)^2}\,N_0.
}
\]

Equivalently, using
\(D_0(B_4+Z_4)=3(M_\Sigma+B_2+Z_2)^2\), the same condition can be written as

\[
N_4
=
-5\,\frac{(B_4+Z_4)^2}{9(M_\Sigma+B_2+Z_2)^2}\,N_0.
\]

So the parent-complete branch must not only hit the one-pole surface and the static prefactor target. It must also land the higher transfer moments on the correct shape relation if the constant-prefactor outgoing branch is the one realized.

The script no longer checks these conditions by solving for \(N_2,N_4\) and
then substituting that same solution back into \(P_2,P_4\). Instead it shows:

- the cleared \((P_2,P_4)\) system is triangular with determinant \(D_0^3\),
- \(P_2=(N_2-N_2^{\rm branch})/D_0\),
- \(P_4|_{N_2=N_2^{\rm branch}}=(N_4-N_4^{\rm branch})/D_0\),

so the branch formulas are now structural factorizations of the outgoing
conditions rather than a solve-and-reinsert tautology.

The script also mutates the constant-prefactor closures by replacing
\(N_2^{\rm branch}\) and \(N_4^{\rm branch}\) with
\(N_2^{\rm branch}+\epsilon\) and \(N_4^{\rm branch}+\epsilon\). Both mutated
factorizations are required to become nonzero, so the zero checks are
load-bearing rather than decorative.

The audit script also includes a concrete wall-integral sanity check with
\(\beta_2(w)=e^{-w^2/2}\), \(\mu_\eta=1\), \(T_w=1\), \(T_\Omega=1/6\), and
\(K_\eta=0\), giving

\[
M_\Sigma=\sqrt{\pi},
\qquad
K_\Sigma=\frac{3\sqrt{\pi}}{2}.
\]

---

## 7. What to export from an actual branch

A single frozen isotropic branch should export:

1. `R_0(w)` and the grouped wall profile `\beta_2(w)`,
2. the promoted-action coefficient functions `\mu_\eta,T_w,T_\Omega,K_\eta`,
3. the two wall integrals `M_\Sigma,K_\Sigma`,
4. the support moments `B_0,B_2,B_4`,
5. the conservative mixed moments `Z_0,Z_2,Z_4`,
6. the outgoing moments `N_0,N_2,N_4`,
7. the source factor `\widehat m_0`.

After that, the entire isotropic branch test is algebraic.
