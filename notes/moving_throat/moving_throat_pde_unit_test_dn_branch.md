# Moving-Throat PDE — Unit-Test D/N Benchmark

## Purpose
This note records the **first actual calculation** to accompany the Phase-1 linearized scaffold.
It is the frozen-wall, finite-throat benchmark that the full moving-throat PDE must reduce to in the appropriate limit.

---

## 1. Benchmark problem

Take a straight finite throat of length \(L_0\) with longitudinal coordinate
\[
s\in[0,L_0].
\]
Assume a scalar support field \(\phi(s,t)\) obeying the linear wave equation
\[
\partial_t^2\phi-c_s^2\partial_s^2\phi=0.
\]
Pass to frequency space,
\[
\phi(s,t)=\Re\bigl[\hat\phi(s,\omega)e^{-i\omega t}\bigr],
\qquad
k\equiv \omega/c_s,
\]
so that
\[
\hat\phi''+k^2\hat\phi=0.
\]

Use the finite-throat boundary conditions
- prescribed mouth datum at \(s=0\): \(\hat\phi(0,\omega)=\hat\phi_m(\omega)\),
- Neumann cap at \(s=L_0\): \(\hat\phi'(L_0,\omega)=0\).

---

## 2. Exact solution

The general solution is
\[
\hat\phi(s)=A\cos(ks)+B\sin(ks).
\]
Imposing \(\hat\phi'(L_0)=0\) gives
\[
-kA\sin(kL_0)+kB\cos(kL_0)=0
\quad\Longrightarrow\quad
B=A\tan(kL_0).
\]
Then the mouth datum \(\hat\phi(0)=A=\hat\phi_m\) yields
\[
\boxed{
\hat\phi(s,\omega)
=
\hat\phi_m(\omega)
\frac{\cos(k(L_0-s))}{\cos(kL_0)}.
}
\]

---

## 3. Exact mouth operator

Define the scalar mouth operator as outgoing normal derivative over mouth value,
\[
Z_{00}^{\rm DN}(\omega)
\equiv
-\frac{\partial_s\hat\phi(0,\omega)}{\hat\phi(0,\omega)}.
\]
Using the exact solution,
\[
\partial_s\hat\phi(0,\omega)=k\tan(kL_0)\,\hat\phi_m(\omega),
\]
so
\[
\boxed{
Z_{00}^{\rm DN}(\omega)
=
-\frac{\omega}{c_s}\tan\!\left(\frac{\omega L_0}{c_s}\right).
}
\]

This is the exact finite-throat benchmark mouth operator.

---

## 4. Pole ladder

Poles occur when
\[
\cos(kL_0)=0,
\]
so the benchmark throat resonances are
\[
\boxed{
\omega_n^{\rm pole}
=
\frac{\pi c_s}{L_0}\left(n+\frac12\right),
\qquad n=0,1,2,\dots
}
\]
This is the characteristic half-shifted finite-throat ladder already identified in the earlier finite-throat reinterpretation.

---

## 5. Low-frequency expansion

Using
\[
\tan x=x+\frac{x^3}{3}+\frac{2x^5}{15}+O(x^7),
\]
we get
\[
Z_{00}^{\rm DN}(\omega)
=
-\frac{\omega}{c_s}
\left[
\frac{\omega L_0}{c_s}
+\frac13\left(\frac{\omega L_0}{c_s}\right)^3
+\frac{2}{15}\left(\frac{\omega L_0}{c_s}\right)^5
+O(\omega^7)
\right].
\]
Therefore
\[
\boxed{
Z_{00}^{\rm DN}(\omega)
=
-\frac{L_0}{c_s^2}\,\omega^2
-\frac{L_0^3}{3c_s^4}\,\omega^4
-\frac{2L_0^5}{15c_s^6}\,\omega^6
+O(\omega^8).
}
\]
So the first benchmark coefficients are
\[
\boxed{Z_2=-\frac{L_0}{c_s^2},
\qquad
Z_4=-\frac{L_0^3}{3c_s^4}.}
\]

---

## 6. Why this benchmark matters

This benchmark is the correct first unit test for the full linearized moving-throat PDE because:

1. it is genuinely **finite-throat**, not a noncompact soft-wall surrogate;
2. it reproduces the exact half-shifted pole ladder expected on the finite branch;
3. it gives a clean analytic low-frequency series against which any numerically extracted response operator can be checked;
4. it is the simplest branch on which the future wall-field lift must collapse when wall motion is frozen.

If the full coupled linearized PDE does not reduce to this branch in the geometry-frozen limit, then the interface/geometry bookkeeping is wrong before any grouped-\(P_2\) or quadrupole normalization work begins.

---

## 7. Immediate next benchmark extensions

The next benchmark extensions should be tackled in this order:

1. allow a breathing wall mode \(\eta_{00}(s,t)\) and derive the first correction to \(Z_{00}^{\rm DN}(\omega)\);
2. switch on one grouped real \(P_2\) wall mode at a time and define the corresponding benchmark operators \(Z_{2m}(\omega)\);
3. test isotropy on the reference branch by checking whether the grouped \(P_2\) low-frequency coefficients collapse to one common value.

That is the cleanest way to start turning the scaffold into an actual response calculation.
