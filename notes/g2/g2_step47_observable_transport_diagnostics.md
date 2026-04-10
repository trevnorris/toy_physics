# Step 47 — The adiabatic anomaly defines a one-parameter flow on the outgoing target surface and can be read from any even bundle coefficient

## Goal

Step 46 compressed the minimal isotropic outgoing branch to an exact algebraic
surface in the observable bundle \((K_0,K_2,K_4,\Gamma_5)\):

\[
\Gamma_5=\frac{2G}{5c^5},
\qquad
K_2^2=\frac14 K_0K_4.
\]

But that still does not tell us how the **adiabatic anomaly parameter**
\[
\ell=\ln(1+\Lambda_1 f)
\]
moves the bundle *along* this surface.

The next honest move is therefore:

> derive the exact one-parameter transport law on that target surface, and then
> invert it so that a future PDE output can be tested directly against the anomaly
> parameter without first reconstructing \((a,c_s)\).

This step does that.

---

## Step 47A — Exact finite flow on the target surface

From Step 46,
\[
K_0 \propto \left(\frac{c_s}{a}\right)^5,\qquad
K_2 \propto \left(\frac{c_s}{a}\right)^3,\qquad
K_4 \propto \left(\frac{c_s}{a}\right),\qquad
\Gamma_5=\text{const}.
\]

On the anomaly-only adiabatic branch we already fixed
\[
\delta\ln a = 0,
\qquad
\delta\ln c_s = \frac{\ell}{5},
\]
so
\[
\frac{c_s}{a}\mapsto e^{\ell/5}\frac{c_s}{a}.
\]

Therefore the outgoing bundle moves exactly as

\[
\boxed{
K_0 \mapsto e^\ell K_0,
\qquad
K_2 \mapsto e^{3\ell/5} K_2,
\qquad
K_4 \mapsto e^{\ell/5} K_4,
\qquad
\Gamma_5 \mapsto \Gamma_5.
}
\]

So the even bundle coefficients do not drift independently. They are locked to a
single one-parameter flow, while the odd quadrupole coefficient stays universal.

---

## Step 47B — Direct inversion formulas for the anomaly parameter

Because the bundle flow is one-parameter, the anomaly parameter can be read from
**any one** of the even coefficients:

\[
\boxed{
\ell
=
\ln\!\frac{K_0}{K_{0,*}}
=
\frac{5}{3}\ln\!\frac{K_2}{K_{2,*}}
=
5\ln\!\frac{K_4}{K_{4,*}},
}
\]
where \((K_{0,*},K_{2,*},K_{4,*})\) are the reference coefficients on the same
constant-prefactor branch.

That is a very practical theorem gate.

A future PDE calculation does **not** need to reconstruct \(a\) and \(c_s\)
first. It can compute any one of the even bundle coefficients, compare it to the
reference value, and infer \(\ell\) immediately.

Then the other two even coefficients are predicted automatically.

---

## Step 47C — Exact differential transport laws

Linearizing the finite flow gives the exact logarithmic transport laws

\[
\boxed{
\delta\ln K_0 = \ell,
\qquad
\delta\ln K_2 = \frac35\ell,
\qquad
\delta\ln K_4 = \frac15\ell,
\qquad
\delta\ln\Gamma_5 = 0.
}
\]

Eliminating \(\ell\) gives purely observable slope identities:

\[
\boxed{
\delta\ln K_0 = 5\,\delta\ln K_4,
}
\]

\[
\boxed{
\delta\ln K_2 = 3\,\delta\ln K_4,
}
\]

\[
\boxed{
\delta\ln K_0 = \frac53\,\delta\ln K_2.
}
\]

So the adiabatic anomaly track is now an extremely sharp bundle diagnostic:

- if a computed branch does **not** satisfy these slope laws, it is not the same
  adiabatic constant-prefactor track;
- if it **does**, then the whole even outgoing bundle is already locked.

---

## Step 47D — Electron-point numbers

Using the carried electron value
\[
\ell=\ln(1+\Lambda_1 f),
\qquad
\Lambda_1\approx 0.279605891931464,
\qquad
f\approx 0.001161409732093,
\]
gives

\[
\ell \approx 3.24684288391064\times 10^{-4}.
\]

So the exact bundle ratios are

\[
\boxed{
\frac{K_0}{K_{0,*}} = e^\ell \approx 1.00032473700404,
}
\]

\[
\boxed{
\frac{K_2}{K_{2,*}} = e^{3\ell/5} \approx 1.00019482954985,
}
\]

\[
\boxed{
\frac{K_4}{K_{4,*}} = e^{\ell/5} \approx 1.00006493896612,
}
\]

\[
\boxed{
\frac{\Gamma_5}{\Gamma_{5,*}} = 1.
}
\]

So the electron anomaly predicts an exact ppm pattern in the outgoing bundle:

- \(K_0\): `324.737` ppm,
- \(K_2\): `194.830` ppm,
- \(K_4\): `64.939` ppm,
- \(\Gamma_5\): `0` ppm.

That is a much cleaner and more testable signature than the original staggered
charge/inertia picture.

---

## Step 47E — What is established, and what is still conditional

### Established here

Once we keep the already-chosen ingredients,

1. adiabatic wall,
2. pure anomaly increment on top of the fixed ground state,
3. constant isotropic outgoing prefactor branch,

the outgoing bundle is not just on an algebraic surface. It follows a single exact
flow with invariant odd coefficient.

### Still conditional

What is still **not** yet proven is that the true moving-throat PDE selects this
constant-prefactor branch rather than a more general outgoing branch with
nonzero \(P_2\) or \(P_4\).

So this step should be read as:

> the exact observable transport law of the minimal outgoing branch.

That is still exactly the right next theorem gate, because it reduces the PDE test
to a handful of coefficient comparisons.

---

## Main result of the step

The adiabatic anomaly defines a one-parameter flow on the minimal outgoing target
surface:

\[
\boxed{
K_0 \mapsto e^\ell K_0,
\qquad
K_2 \mapsto e^{3\ell/5} K_2,
\qquad
K_4 \mapsto e^{\ell/5} K_4,
\qquad
\Gamma_5 \mapsto \Gamma_5.
}
\]

Equivalently,

\[
\boxed{
\ell
=
\ln\!\frac{K_0}{K_{0,*}}
=
\frac53\ln\!\frac{K_2}{K_{2,*}}
=
5\ln\!\frac{K_4}{K_{4,*}},
}
\]

with differential form

\[
\boxed{
\delta\ln K_0 = 5\,\delta\ln K_4,
\qquad
\delta\ln K_2 = 3\,\delta\ln K_4,
\qquad
\delta\ln\Gamma_5 = 0.
}
\]

So the next PDE-facing test is now extremely concrete:

> compute the outgoing bundle and check whether it lands on this flow with a fixed universal \(\Gamma_5\).
