# Parent Throat Action — Weak-Axisymmetric Packet Bridge

## Purpose

This note pushes the promoted throat action into the weak-axisymmetric packet.
The goal is to replace the abstract wall slopes by explicit parent-action wall integrals and then solve the live first-order gates.

---

## 1. Parent wall slopes from a `Y_{20}` deformation

On the promoted throat-action branch, let a weak axisymmetric perturbation induce first-order variations

\[
\mu_\eta \to \mu_\eta + \varepsilon\,\delta\mu_\eta(w)Y_{20}(\Omega),
\]

\[
T_w \to T_w + \varepsilon\,\delta T_w(w)Y_{20}(\Omega),
\qquad
T_\Omega \to T_\Omega + \varepsilon\,\delta T_\Omega(w)Y_{20}(\Omega),
\qquad
K_\eta \to K_\eta + \varepsilon\,\delta K_\eta(w)Y_{20}(\Omega).
\]

For the grouped `P2` wall profile `\beta_2(w)`, define the wall slope integrals

\[
\boxed{
\delta M_\Sigma = \int dw\,\delta\mu_\eta(w)\,\beta_2(w)^2,
}
\]

\[
\boxed{
\delta K_\Sigma = \int dw\,\Big[\delta T_w(w)(\beta_2'(w))^2 + (\delta K_\eta(w)+6\delta T_\Omega(w))\beta_2(w)^2\Big].
}
\]

A pure `Y_{20}` angular deformation still carries the grouped signature

\[
\lambda = \left(1,\frac12,-1\right),
\]

so every first-order grouped wall quantity lies on the weak-axisymmetric line `b = 3a`.
The accompanying SymPy audit derives this signature from the Wigner/Gaunt
triple-overlap coefficients, including the real-harmonic \((-1)^m\) factor
that converts the complex-\(Y_{2m}\) overlap into the squared real-harmonic
lane weight, and then verifies the compensated packet algebra.

---

## 2. Full weak-axisymmetric denominator slopes

Carry the already-derived support and conservative mixed slopes

\[
B_{01},B_{21},B_{41},
\qquad
Z_{01},Z_{21},Z_{41},
\qquad
N_{01}.
\]

Then the full first-order denominator moments are

\[
\boxed{
D_{01} = \delta K_\Sigma - B_{01} - Z_{01},
}
\]

\[
\boxed{
D_{21} = -\,(\delta M_\Sigma + B_{21} + Z_{21}),
}
\]

\[
\boxed{
D_{41} = -\,(B_{41}+Z_{41}).
}
\]

The transported weak-axisymmetric packet is

\[
\boxed{
K_1 = D_{21} + \frac{D_{01}}{9},
}
\]

\[
\boxed{
H_{\rm even} = D_{41} - \frac{2}{3}D_{21} - \frac{D_{01}}{27},
}
\]

\[
\boxed{
\Xi_1 = \frac{N_{01}}{N_0} - \frac{D_{01}}{D_0},
\qquad
D_0 = K_\Sigma - B_0 - Z_0.
}
\]

So the parent promotion does not change the gate definitions.  It changes only the meaning of the wall-side slots.

---

## 3. Exact wall solve of the even gates

Impose the canonical even-compensation conditions

\[
K_1 = 0,
\qquad
H_{\rm even} = 0.
\]

Solving for the parent wall slopes gives

\[
\boxed{
\delta K_\Sigma = B_{01}+Z_{01}+27(B_{41}+Z_{41}),
}
\]

\[
\boxed{
\delta M_\Sigma = -\,(B_{21}+Z_{21}) + 3(B_{41}+Z_{41}).
}
\]

This is the sharpest new result of the continuation.

It says that on the canonical compensated branch, the parent wall block has no first-order freedom left once the support and conservative mixed anisotropy slopes are fixed.
The wall promotion therefore removes an important source of post-target retuning: one can no longer independently choose arbitrary `\delta K` and `\delta M` after the support/mixed branch is known.

---

## 4. Canonical compensated branch relations

Substituting the exact wall-slope solve back into the denominator moments gives

\[
\boxed{
D_{01} = 27(B_{41}+Z_{41}),
}
\]

\[
\boxed{
D_{21} = -3(B_{41}+Z_{41}),
}
\]

\[
\boxed{
D_{41} = -(B_{41}+Z_{41}).
}
\]

Hence the canonical compensated relations follow automatically,

\[
\boxed{D_{21} = -\frac{D_{01}}{9}},
\qquad
\boxed{D_{41} = -\frac{D_{01}}{27}}.
\]

So the familiar weak-axisymmetric denominator pattern reappears, but now as a consequence of the parent wall solve rather than as a purely abstract bookkeeping condition.

---

## 5. Residual prefactor slope on the parent-complete branch

After even compensation, the only remaining first-order normalization defect is

\[
\boxed{
\Xi_1
=
\frac{N_{01}}{N_0}
-
\frac{27(B_{41}+Z_{41})}{K_\Sigma-B_0-Z_0}.
}
\]

So the residual weak-axisymmetric theorem gate depends on three inputs only:

1. the isotropic denominator `D_0 = K_\Sigma-B_0-Z_0`,
2. the combined conservative even anisotropy `B_{41}+Z_{41}`,
3. the outgoing slope `N_{01}/N_0`.

The parent wall slopes `\delta K_\Sigma,\delta M_\Sigma` no longer appear explicitly after compensation because they have already been fixed by the even-gate solve.

---

## 6. Parent-complete interpretation

This continuation clarifies the role of the promoted throat action.

### What the parent promotion now controls directly

- The isotropic denominator `D_0` through `K_\Sigma`.
- The wall contributions to the weak-axisymmetric packet through `\delta K_\Sigma,\delta M_\Sigma`.
- The exact wall-side part of the even-gate compensation.

### What still must come from the realized branch

- The BdG slope packet `B_{01},B_{21},B_{41}`.
- The conservative mixed slope packet `Z_{01},Z_{21},Z_{41}`.
- The outgoing slope `N_{01}`.

So the parent promotion tightens the branch protocol, but it does not eliminate the need to compute the support and mixed sectors on the same frozen branch.

---

## 7. Actual-branch export packet

A frozen weak-axisymmetric branch should export:

1. the isotropic wall data `K_\Sigma,M_\Sigma`,
2. the wall slope integrals `\delta K_\Sigma,\delta M_\Sigma`,
3. the support slope packet `B_{01},B_{21},B_{41}`,
4. the conservative mixed slope packet `Z_{01},Z_{21},Z_{41}`,
5. the outgoing slope `N_{01}`.

The algebraic checks are then:

\[
K_1 = 0,
\qquad
H_{\rm even} = 0,
\qquad
\Xi_1 = 0 \;\text{or at least admissibly small}.
\]

That is the complete first-order parent-complete packet.
