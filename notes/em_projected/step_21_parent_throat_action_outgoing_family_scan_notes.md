# Parent Throat Action — Outgoing-Family Extension Scan

## Purpose

`step_20` showed that the three-coordinate reduced branch family

\[
(\log A,\ \log \alpha_R,\ \log \alpha_\beta)
\]

cannot cancel all four isotropic residuals at first order near the frozen
branch. This step tests the next natural enlargement:

add one upstream **outgoing-profile curvature coordinate** while keeping the
same reduced parent-background wall operator and the same target-blind branch
discipline.

So this is not a post-target retune of `(N_0,N_2,N_4)`. It is a declared
reduced family extension before the residual check.

---

## 1. Frozen baseline

The script declares

- `branch_id = v2_local_parent_background_outgoing_family_scan`,
- `pre_target_freeze = true`,
- `target_blind = true`,
- `no_post_residual_refit = true`,
- `boundary_class = open_impedance_demo`,
- `branch_freeze_hash = 0da71b0e29736d98`.

At the baseline point it reproduces the same 4-mode residual packet as
`step_19` and `step_20`:

\[
(R_{\rm pole},R_{\rm norm},R_{P2},R_{P4})
\approx
(-13.134593938872376,\,-10.33719584868593,\,
0.37009844569768474,\,
0.8889149882257381).
\]

Its norm is
\[
\|R\|_2 \approx 16.742231591665213.
\]

---

## 2. Extended reduced family

The outgoing source is promoted from the fixed baseline form

\[
\phi_N(w)=\left(1+\frac{w^2}{2}\right)\beta_2(w)
\]

to the one-parameter family

\[
\phi_N(w)=\bigl(1+(\tfrac12+\delta_N)w^2\bigr)\beta_2(w),
\]

with family coordinates

\[
(\log A,\ \log \alpha_R,\ \log \alpha_\beta,\ \delta_N).
\]

The wall block and the conservative support/mixed sources remain otherwise
unchanged:

\[
R_0(w)=A\,e^{-\alpha_R w^2/2},
\qquad
\beta_2(w)=A\,e^{-\alpha_\beta w^2/2},
\]
\[
\phi_B=\beta_2,
\qquad
\phi_Z=R_0\beta_2.
\]

So the added coordinate is a reduced outgoing-profile family coordinate, not a
retroactive moment repair.

---

## 3. Jacobian and linear rank

With central-difference step `h = 0.02`, the residual Jacobian becomes

\[
J
\approx
\begin{pmatrix}
 -70.025142145431 & -2.246449416019 & 11.693618046138 & 0 \\
 -0.568305619312 & -0.117170794621 & -0.315922071187 & 0.438756463756 \\
 -0.499694612692 & 0.003100955096 & -0.440318331254 & 0.371693302795 \\
 -1.403341014894 & -0.360050582391 & -0.952565606312 & 0.854814803381
\end{pmatrix}.
\]

Its singular values are

\[
\sigma(J)\approx
(71.044213096589,\ 1.715640711598,\ 0.141855721241,\ 0.108888020712),
\]

so the numerical rank is now `4`.

That is the first important change from `step_20`: after adding the outgoing
curvature coordinate, the local family has enough first-order directions to
span the whole four-component isotropic residual packet.

The least-squares linearized correction is

\[
\delta_{\rm LS}
\approx
(8.677800811209,\ 2.872115918746,\ 53.640446430441,\ 74.190514652389),
\]

and the linearized residual drops to numerical zero:

\[
\|R + J\delta_{\rm LS}\|_2
\approx 4.4642745655749085\times 10^{-13}.
\]

So there is **no first-order rank obstruction** in this enlarged reduced
family.

---

## 4. Actual finite steps

The linear picture is not the whole story, because the required correction is
large. The script therefore evaluates actual branch points along that direction.

At scale `0.01`:

\[
\|R\|_2 \approx 16.601833564436067,
\]
\[
R \approx
(-12.93103029,\,-10.39720022,\,
0.19154916,\,
0.52045143).
\]

At scale `0.02`:

\[
\|R\|_2 \approx 16.447455615657052,
\]
\[
R \approx
(-12.59169074,\,-10.58009881,\,
0.01933555,\,
0.17104359).
\]

At scale `0.03`:

\[
\|R\|_2 \approx 16.123181731091783,
\]
\[
R \approx
(-12.05718987,\,-10.70411312,\,
-0.01928138,\,
0.05245994).
\]

At scale `0.05`:

\[
\|R\|_2 \approx 14.73799806685835,
\]
\[
R \approx
(-10.0437788,\,-10.7856852,\,
-0.00693044,\,
0.00639645).
\]

At scale `0.08`:

\[
\|R\|_2 \approx 11.40374534724855,
\]
\[
R \approx
(-3.66355306,\,-10.7992494,\,
-4.46772936\times10^{-4},\,
3.28092630\times10^{-4}).
\]

At scale `0.10`:

\[
\|R\|_2 \approx 11.328790887276268,
\]
\[
R \approx
(3.42109084,\,-10.7998908,\,
-6.61056196\times10^{-5},\,
4.64262836\times10^{-5}).
\]

Among the tested scales, `0.10` gives the smallest residual norm.

---

## 5. Interpretation

This step changes the diagnosis in a useful way.

1. The earlier three-parameter family had a genuine local rank obstruction.
2. After adding one outgoing-family coordinate, that first-order obstruction is
   gone.
3. But the actual finite-step scan shows a new practical obstruction:
   `R_{P2}` and `R_{P4}` can be driven essentially to zero, and `R_{\rm pole}`
   can be reduced substantially, while
   \[
   R_{\rm norm}\approx -10.8
   \]
   remains dominant.

So in this enlarged reduced family the main unresolved isotropic bottleneck is
no longer the constant-prefactor packet. It is the **static normalization**
surface.

That is the theorem-facing takeaway of `step_21`.
