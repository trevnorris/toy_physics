
# 2PN geometry-Hessian monopole breathing closure: current result

## What this step adds

The previous Family-1 wall work compressed the full static even sector to

\[
\text{[local support/source constitutive law]} + \frac{109}{280}\,\mathbb P_{00},
\]

so the only unresolved static piece was the **global monopole breathing auxiliary**.

This step shows that the missing monopole projector can be generated directly from a reduced **geometry Hessian** in the throat variables \((a,L)\), with a natural volume-work coupling. So the monopole add-on does **not** need to be inserted as an independent ad hoc channel.

---

## 1. Minimal geometry-side model

Use the 4D cylinder-like throat geometry

\[
V(a,L)=\frac{4\pi}{3}a^3L,
\qquad
A(a,L)=4\pi a^2L+\frac{8\pi}{3}a^3,
\]

and the minimal curvature-completed geometry energy

\[
E_{\rm geom}(a,L)
=
P_{\rm vac}V(a,L)+\sigma A(a,L)+\kappa_b\frac{a^2}{L}.
\]

The new ingredient is the last term. It is the smallest explicit curvature/bending completion that can repair the known \(P_{\rm vac}V+\sigma A\) limitation.

---

## 2. Natural monopole coupling and exact projector coefficient

Treat the uniform monopole port as a pressure-like source that couples through volume work,

\[
\delta W = -p\,\delta V,
\qquad
g=\nabla_{(a,L)}V.
\]

Let \(H_0\) be the Hessian of \(E_{\rm geom}\) at the reference point
\[
(a_0,L_0)=(a_0,\Lambda a_0).
\]

Then integrating out \((\delta a,\delta L)\) gives the exact reduced geometry-side monopole coefficient

\[
\boxed{
\Delta K_{00}^{\rm geom}
=
\frac{g^T H_0^{-1} g}{V_0^2}.
}
\]

So the old global monopole projector is reinterpreted as a genuine low-frequency **geometry compressibility**.

---

## 3. Baseline no-go from the actual geometry Hessian

Write

\[
\sigma=\frac{\Sigma}{a_0^3},
\qquad
P_{\rm vac}=\rho\,\frac{\Sigma}{a_0^4},
\qquad
\kappa_b=\beta\,\frac{\Sigma}{a_0}.
\]

Then the exact reduced coefficient is

\[
\boxed{
\Delta K_{00}^{\rm geom}
=
\frac{2\pi\Lambda^2\rho+5\pi\Lambda^2-2\pi\Lambda-4\beta}
{2\pi\Sigma\left(\pi\Lambda^3\rho^2+4\pi\Lambda^3\rho+4\pi\Lambda^3-2\Lambda\beta\rho-3\Lambda\beta-2\beta\right)}.
}
\]

For the baseline model with no curvature completion \((\beta=0)\), the reduced Hessian is

\[
H_{\rm base}=
\begin{pmatrix}
8\pi(\Lambda\rho+\Lambda+2) & 4\pi(\rho+2) \\
4\pi(\rho+2) & 0
\end{pmatrix},
\]

with determinant

\[
\det H_{\rm base} = -16\pi^2(\rho+2)^2 < 0.
\]

So the literal \(P_{\rm vac}V+\sigma A\) geometry energy is **not** a passive 2DOF geometry Hessian and cannot by itself close the monopole channel. That exactly matches the limitation flagged in the throat notes.

---

## 4. Exact positivity conditions for the minimal curvature completion

The curvature-completed Hessian is passive / positive-definite when

\[
\beta>0
\]
and
\[
\boxed{
\beta >
\beta_{\rm stab}(\Lambda,\rho)
=
\frac{\pi\Lambda^3(\rho+2)^2}{2\Lambda\rho+3\Lambda+2}.
}
\]

For the resulting monopole coefficient to be positive as well, one also needs

\[
\boxed{
\beta >
\beta_{\Delta}(\Lambda,\rho)
=
\frac{\pi\Lambda(2\Lambda\rho+5\Lambda-2)}{4}.
}
\]

So the geometry-side closure exists whenever

\[
\beta > \max(\beta_{\rm stab},\beta_\Delta).
\]

---

## 5. EM-branch worked example

Take the EM-cavity aspect ratio
\[
\Lambda_{\rm EM}=\frac{\sqrt{2}\pi}{x_{01}},
\qquad
x_{01}=2.40482555769577\ldots,
\]
so
\[
\Lambda_{\rm EM}\approx 1.847486577120128.
\]

A simple positive example is

\[
\rho=\frac{1}{10},
\qquad
\beta=12,
\]
for which

\[
\beta_{\rm stab}(\Lambda_{\rm EM},\rho)\approx 11.0420171,
\qquad
\beta_\Delta(\Lambda_{\rm EM},\rho)\approx 11.0377513.
\]

So this point lies safely inside the passive/positive-support region.

Matching the exact target
\[
\Delta K_{00}^{\rm geom}=\frac{109}{280}
\]
then fixes the overall geometry stiffness scale to

\[
\boxed{
\Sigma_* \approx 0.2076143291835488854.
}
\]

At this point the exact reduced geometry coefficient is

\[
\Delta K_{00}^{\rm geom} = \frac{109}{280}
\]
to machine precision.

---

## 6. Why the earlier single breathing auxiliary was already the right idea

At the worked point above, the 2DOF geometry Hessian eigenvalues are

\[
\lambda_1 \approx 0.10664211,
\qquad
\lambda_2 \approx 24.42044437.
\]

Using the natural normalized coupling vector

\[
\hat g = \frac{\nabla V}{V_0} = \left(3,\frac{1}{\Lambda_{\rm EM}}\right),
\]

the mode-resolved contributions to \(\Delta K_{00}^{\rm geom}\) are

\[
(0.00878310,\ 0.38050262),
\]
which sum to

\[
0.389285714285714\ldots = \frac{109}{280}.
\]

So the dominant mode contributes

\[
97.743791\%
\]
of the total static monopole response.

That is the key interpretation:

> the earlier single global monopole auxiliary is a **controlled reduction** of the full \((a,L)\) geometry sector, not an extra arbitrary assumption.

At this stage the monopole closure is structurally:

\[
\boxed{
\text{local Family-1 wall law}
\;+\;
\text{dominant-coupling reduction of the reduced geometry Hessian}.
}
\]

---

## 7. What remains

This does **not** finish the PDE derivation of the monopole channel, but it narrows it a lot.

The remaining task is now:

1. derive the curvature completion \(\kappa_b a^2/L\) (or its more accurate soft-wall analog) from the actual Family-1 soft-wall / traction physics,
2. and derive the corresponding breathing inertia if we want the full pole
   \[
   Y_{\rm mono}(\omega)=\frac{109/280}{1-\omega^2/\Omega_{\rm mono}^2}.
   \]

So the static monopole story is no longer “mysterious projector needed.”
It is:

- baseline \(PV+\sigma A\) is insufficient,
- minimal curvature completion repairs the geometry Hessian,
- and the passed \(109/280\) coefficient is realizable directly from that repaired geometry sector.
