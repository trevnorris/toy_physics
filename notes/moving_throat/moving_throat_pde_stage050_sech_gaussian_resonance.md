
# Moving-Throat PDE — Stage 50: Exact Sech–Gaussian Coherence Resonance Benchmark

## Purpose

The memo you shared suggests a concrete transverse profile benchmark for the parent source/support coherence:

- fluid/compressional source profile:
  `chi_sigma(y) = sech(y / w_f)`,
- geometric/support profile:
  `chi_phi(y) = exp(-y^2 / w_g^2)`.

This stage checks that benchmark carefully and places it in the right role inside the reduced theorem program.

The result is useful, but more limited than the memo’s strongest interpretation.

The useful theorem is:

1. the benchmark has an **exact self-duality**,
2. that self-duality forces an exact stationary point at
   `w_g / w_f = sqrt(pi)`,
3. the corresponding maximal coherence is
   `C_res^2 = 0.994418836451529...`,
4. so the best independent sech–Gaussian mismatch branch gets within about `0.56%` of the exact matched-layer ideal,
5. but this does **not** by itself prove threshold survival, because the survival theorem still depends on the wall/source figure of merit from Stages 44–49, not on coherence alone.

So the memo is genuinely helpful as an explicit benchmark family, but it is not itself the final theorem.

---

## 1. Exact norms

Let

`chi_sigma(y) = sech(y / w_f),`

`chi_phi(y) = exp(-y^2 / w_g^2),`

with positive widths `w_f`, `w_g`, and define the ratio

`r := w_g / w_f.`

The exact norms are

`N_(sigma sigma) = int dy chi_sigma^2 = 2 w_f,`

`N_(phi phi)     = int dy chi_phi^2 = w_g sqrt(pi/2).`

The overlap is

`O_(sigma phi) = int dy chi_sigma chi_phi`
`             = w_f I(r),`

with the dimensionless integral

`I(r) := int_{-inf}^{inf} dx sech(x) exp(-x^2 / r^2).`

Therefore the coherence factor is

`C^2(r) = O_(sigma phi)^2 / [ N_(sigma sigma) N_(phi phi) ]`
`       = I(r)^2 / [ r sqrt(2 pi) ].`

---

## 2. Exact self-duality

Using the Fourier-transform identity

`FT[ sech(x) ](k) = pi sech(pi k / 2),`

together with the Gaussian transform

`FT[ exp(-x^2 / r^2) ](k) = sqrt(pi) r exp(-r^2 k^2 / 4),`

Parseval gives

`I(r) = (sqrt(pi) r / pi) int dt sech(t) exp(-r^2 t^2 / pi^2)`
`     = (r / sqrt(pi)) I(pi / r).`

So the overlap obeys the exact duality law

`I(r) = (r / sqrt(pi)) I(pi / r).`

Substituting into the coherence formula gives

`C^2(r) = C^2(pi / r).`

That is the key exact structural fact behind the resonance.

---

## 3. Exact stationary resonance point

Because `C^2(r) = C^2(pi / r)`, the self-dual point is

`r_* = sqrt(pi).`

Differentiating the duality relation gives

`dC^2/dr |_(r = sqrt(pi)) = 0.`

So the ratio

`w_g / w_f = sqrt(pi)`

is not just numerically close to optimal. It is the **exact self-dual stationary point** of the explicit sech–Gaussian benchmark family.

A numerical monotonicity audit then shows that this stationary point is the unique global maximum on the constructive branch.

---

## 4. Numerical maximum

Evaluating the exact dimensionless overlap at the self-dual point gives

`C_res^2 := C^2(sqrt(pi))`
`        = 0.9944188364515293487...`

So the best independent sech–Gaussian profile mismatch falls short of perfect support/source matching by only

`1 - C_res^2 = 0.0055811635484706513...`

or about `0.56%`.

Equivalently, the resonance penalty factor is

`P_res := 1 / C_res^2`
`      = 1.0056124877605762169...`

That is the exact multiplicative penalty by which the explicit independent-profile family misses the ideal matched-layer branch.

---

## 5. Why this does not supersede Stage 47

This benchmark is important, but it does **not** supersede the earlier parent equilibrium result.

Stage 47 already showed that on the parent equilibrium-aligned source/support branch,

`chi_sigma(y) = g_phi chi_phi(y) / H(y),`

and in the thin active layer where `H(y)` is nearly constant across the support region, the coherence becomes

`C_(sigma phi)^2 = 1.`

So:

- the present sech–Gaussian family is a strong **independent-profile benchmark**,
- the Stage-47 branch is the **equilibrium-matched** benchmark.

The new resonance result therefore tells us that even a fairly natural profile mismatch can get extremely close to the exact matched limit. But it does not replace the matched-limit theorem.

---

## 6. What Stage 50 changes

Before this stage, the mismatch between an ideal matched source/support branch and a concrete profile family was still diffuse.

After this stage, the first explicit off-equilibrium profile family is sharply characterized:

- exact stationary resonance at `w_g / w_f = sqrt(pi)`,
- exact duality `C^2(r) = C^2(pi/r)`,
- maximal coherence `C_res^2 = 0.994418836451529...`,
- penalty factor `P_res = 1.005612487760576...`.

So the user memo is mathematically useful — but the strongest justified conclusion is not “automatic threshold survival.” The justified conclusion is:

> the first explicit independent-profile family almost saturates the ideal matched source/support branch.

The actual survival theorem still has to be stated in terms of the Stage-44/49 gain and wall-figure thresholds.
