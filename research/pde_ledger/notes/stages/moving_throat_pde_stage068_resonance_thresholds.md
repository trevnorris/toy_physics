
# Moving-Throat PDE — Stage 68: Resonance-Corrected Thresholds for the Sech–Gaussian Benchmark Family

## Purpose

Stage 067 verified that the explicit independent sech–Gaussian profile family has an exact self-dual stationary point at

`w_g / w_f = sqrt(pi),`

with maximal coherence

`C_res^2 = 0.994418836451529...`

and penalty factor

`P_res = 1 / C_res^2 = 1.005612487760576...`

The next question is how that benchmark changes the actual support/source survival thresholds derived earlier.

This stage gives the exact answer:

- on the Stage-064 equilibrium-matched branch, the ideal coherence remains `C^2 = 1`,
- on the independent sech–Gaussian benchmark family, the best possible coherence is `C_res^2`,
- so the explicit independent-profile family modifies the Stage-066 thresholds only by the tiny factor `P_res`.

That means the resonance benchmark is useful, but it does **not** rewrite the theorem structure. It only sharpens it.

---

## 1. Resonance-corrected gain

The general parent gain from Stages 062–063 has the form

`G_micro = [ rho_* g_phi^2 N_(phi phi) / (m c_(s,*)^2 K_X) ] C_(sigma phi)^2.`

On the Stage-064 matched equilibrium branch,

`C_(sigma phi)^2 = 1,`

so the matched-branch gain is

`G_match = rho_* g_phi^2 N_(phi phi) / (m c_(s,*)^2 K_X).`

Stage 066 repackaged this on the first explicit thin-wall branch into the dimensionless wall figure of merit

`W_wall = kappa G_match`
`       = 4 pi a^2 L^2 J_1 V0^2 / (T_X ell).`

For the explicit independent sech–Gaussian family, the actual gain is therefore

`G_res(r) = C^2(r) G_match,`

and the corresponding wall figure is

`W_res(r) = C^2(r) W_wall.`

At the self-dual resonance,

`W_res,* = C_res^2 W_wall.`

So the memo’s profile family simply inserts one multiplicative coherence factor into the already-frozen theorem window.

---

## 2. Exact resonance-family thresholds

Stage 066 gave the universal matched-branch window

`W_wall <= Pe_req / Delta_inf`  -> matched-branch fail,

`W_wall >= Pe_req / Delta_0`    -> matched-branch succeed.

For the explicit independent sech–Gaussian family, replace `W_wall` by `W_res(r)=C^2(r)W_wall`.
That gives the exact profile-family conditions

`W_wall <= Pe_req / [ C^2(r) Delta_inf ]`  -> fail on the profile family,

`W_wall >= Pe_req / [ C^2(r) Delta_0 ]`    -> succeed on the profile family.

At the resonance point `r = sqrt(pi)` this becomes

`W_wall <= P_res Pe_req / Delta_inf`  -> resonance-family fail,

`W_wall >= P_res Pe_req / Delta_0`    -> resonance-family succeed,

with

`P_res = 1 / C_res^2 = 1.005612487760576...`

So the first explicit independent-profile family changes the wall thresholds by only about `0.56%`.

---

## 3. What the memo does and does not prove

This is the critical interpretation point.

The resonance family does **not** by itself prove support-threshold survival.

The true threshold test still depends on:

- `Pe_req`,
- the axial functions `Delta_0`, `Delta_inf`,
- and the actual wall/source figure of merit `W_wall`.

The memo only sharpens the source/support coherence factor inside that test.

So the strongest justified statement is:

> the explicit independent sech–Gaussian family comes within `0.56%` of the exact matched-layer ideal.

That is strong evidence that profile mismatch is probably not the dominant reduced-theorem bottleneck.
But it is **not** the same thing as proving that the support/source branch clears the threshold.

---

## 4. Small profile-sensitivity band

Because the resonance family differs from the matched branch only by the factor `P_res = 1.005612...`, the only region where profile choice can change the reduced verdict is a very narrow band.

On the success side:

- matched-branch guaranteed success:
  `W_wall >= Pe_req / Delta_0,`
- resonance-family guaranteed success:
  `W_wall >= P_res Pe_req / Delta_0.`

So the success-side profile-sensitive band is only

`Pe_req / Delta_0 <= W_wall < P_res Pe_req / Delta_0,`

a width of about `0.56%`.

On the failure side:

- matched-branch guaranteed failure:
  `W_wall <= Pe_req / Delta_inf,`
- resonance-family guaranteed failure:
  `W_wall <= P_res Pe_req / Delta_inf.`

So again the profile-sensitive band is only about `0.56%`.

That is the sharpest practical consequence of the memo.

---

## 5. What Stage 068 changes

Before this stage, the independent-profile benchmark from the memo could have been interpreted as potentially introducing a qualitatively new source/support theorem.

After this stage, that is no longer the right reading.

The memo’s explicit sech–Gaussian resonance family:

- is mathematically real,
- provides a strong benchmark,
- and nearly saturates the matched branch,

but it changes the Stage-066 theorem window only by the small factor

`P_res = 1.005612487760576...`

So the dominant unresolved reduced-theorem question is still not transverse profile mismatch.
It is the actual wall/axial branch data entering `W_wall`, `Delta_0`, and `Delta_inf`.
