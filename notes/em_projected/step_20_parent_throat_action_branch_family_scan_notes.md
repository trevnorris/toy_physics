# Parent Throat Action — Reduced Branch-Family Scan

## Purpose

This step does not replace the frozen `step_19` branch audit. It sits one layer
above it.

- `step_19` is the fixed reduced branch packet and target-blind negative control.
- `step_20` asks whether a nearby reduced branch family can move the isotropic
  residual packet toward zero, or whether the miss already looks structural in
  that local family.

So this note is a **local branch-sensitivity diagnostic**, not a theorem that
the full moving-throat PDE branch exists or fails.

---

## 1. Frozen baseline

The script imports the reduced parent-background residual packet from the live
4-mode checkpoint owner in `step_19` and declares the V2-style audit metadata

- `branch_id = v2_local_parent_background_reduced_family_scan`,
- `pre_target_freeze = true`,
- `target_blind = true`,
- `no_post_residual_refit = true`,
- `boundary_class = open_impedance_demo`,
- `branch_freeze_hash = 74152b0cef38e0ec`.

The frozen baseline packet agrees numerically with the `step_19` 4-mode export:

\[
M_\Sigma \approx 3.340256569115096,
\qquad
K_\Sigma \approx 4.124157928219886,
\]
\[
(B_0,B_2,B_4)
\approx
(0.7816273402373898,\,-0.6449331985854893,\,0.5321531025367580),
\]
\[
(Z_0,Z_2,Z_4)
\approx
(0.5083076399936368,\,-0.4143227047792867,\,0.3408544970748320),
\]
\[
(N_0,N_2,N_4)
\approx
(1.3116901460788468,\,-1.0623716467074260,\,0.8725297055800808).
\]

The corresponding isotropic residual packet is

\[
(R_{\rm pole},R_{\rm norm},R_{P2},R_{P4})
\approx
(-13.134593938872369,\,-10.33719584868593,\,
0.3700984456976848,\,
0.8889149882257383),
\]

with residual norm
\[
\|R\|_2 \approx 16.742231591665213.
\]

---

## 2. Local reduced branch family

The scan varies three upstream reduced branch coordinates:

\[
(\log A,\ \log \alpha_R,\ \log \alpha_\beta),
\]

where

\[
R_0(w)=A\,e^{-\alpha_R w^2/2},
\qquad
\beta_2(w)=A\,e^{-\alpha_\beta w^2/2}.
\]

The parent wall block remains

\[
\mu_\eta(w)=1+\Bigl(1+\frac{w^2}{4}\Bigr)R_0(w),
\]
\[
T_w(w)=1+\Bigl(1+\frac{w^2}{6}\Bigr)R_0(w),
\qquad
T_\Omega(w)=\frac{1+\left(1+\frac{w^2}{8}\right)R_0(w)}{6},
\]

and \(K_\eta(w)\) is recomputed from the background equation for each family
point. The conservative support/mixed/outgoing sources remain

\[
\phi_B=\beta_2,
\qquad
\phi_Z=R_0\beta_2,
\qquad
\phi_N=\left(1+\frac{w^2}{2}\right)\beta_2,
\]

so the scan changes the branch family **before** evaluating targets rather than
retuning the target moments afterward.

---

## 3. Jacobian at the frozen branch

Using a central finite-difference step `h = 0.02`, the script computes the
Jacobian of the isotropic residual packet with respect to
\((\log A,\log \alpha_R,\log \alpha_\beta)\):

\[
J
\approx
\begin{pmatrix}
-70.025142145431 & -2.246449416019 & 11.693618046138 \\
-0.568305619312 & -0.117170794621 & -0.315922071187 \\
-0.499694612692 & 0.003100955096 & -0.440318331254 \\
-1.403341014894 & -0.360050582391 & -0.952565606312
\end{pmatrix}.
\]

Its singular values are

\[
\sigma(J)\approx(71.044210206226,\ 1.376398462469,\ 0.141312291580),
\]

so the numerical rank is `3`.

This already says something useful: in this local reduced family there are only
three first-order control directions for a four-component isotropic residual
packet.

---

## 4. Best linearized correction and its limit

The least-squares first-order correction is

\[
\delta_{\rm LS}
\approx
(-0.122264570277,\,-4.364517656908,\,-0.446335906159).
\]

But the best linearized residual is still

\[
R + J\,\delta_{\rm LS}
\approx
(0.012386508789964168,\,-9.615310840436756,\,
0.6141891008979157,\,
3.057105231659007),
\]

with irreducible linearized norm
\[
\|R+J\delta_{\rm LS}\|_2 \approx 10.108287522271972.
\]

So even after the best first-order correction inside this 3-parameter family,
the residual packet stays far from zero. That is a real local obstruction for
this family, not just a lack of motion.

---

## 5. Actual trial steps

The script then tests finite steps along the least-squares direction.

For scale `0.05`:

\[
\|R\|_2 \approx 16.279110841839092,
\]
\[
R \approx
(-12.55527084058592,\,-10.311239147614597,\,
0.36476518818401554,\,
0.9591233072995973).
\]

For scale `0.10`:

\[
\|R\|_2 \approx 15.965938937899807,
\]
\[
R \approx
(-12.153998893881553,\,-10.302369208407521,\,
0.3356940669736331,\,
0.9695438325968013).
\]

For scale `0.20`:

\[
\|R\|_2 \approx 15.828401363791917,
\]
\[
R \approx
(-11.96804731342165,\,-10.315781715587908,\,
0.2629782100552601,\,
0.9053304890309445).
\]

So the family is not completely flat: it can reduce the residual norm somewhat.
But it does **not** come close to the isotropic target surface in the tested
local window.

---

## 6. Sign-Flip Mutation Guard

The script also mutates the parent wall measure sign,

\[
\mu_\eta(w)=1+\left(1+\frac{w^2}{4}\right)R_0(w)
\quad\longrightarrow\quad
\mu_\eta^{\rm flip}(w)=1-\left(1+\frac{w^2}{4}\right)R_0(w).
\]

This is not treated as a physical branch. It is a sign-convention guard for the
mass-weighted reduced operator. The flipped packet has

\[
\|R_{\rm flip}\|_2 \approx 10.337242025167106,
\]

which is smaller than the baseline scalar norm, so the test deliberately does
not assert that the wrong sign must make the norm worse. Instead it checks the
packet movement itself:

\[
\|R_{\rm flip}-R_{\rm base}\|_2
\approx 13.139644483637701.
\]

That residual-packet separation is the relevant guard: flipping the
\(\mu_\eta\) contribution is not inert, and future edits cannot silently erase
that sign dependence.

---

## 7. Interpretation

The strongest honest reading is:

1. the frozen `step_19` branch packet is a real target-failing negative control;
2. nearby reduced branch motion can improve the packet quantitatively;
3. but this particular three-parameter family has a nonzero first-order
   obstruction to killing all four isotropic residuals at once.

So the next move should not be “retune the same packet harder.” It should be at
least one of:

- enlarge the reduced family by adding another physically motivated upstream
  branch coordinate,
- derive a more faithful outgoing/mixed branch export from the PDE side,
- or move to a genuinely different realized branch family.

That is the useful theorem-facing output of this step.
