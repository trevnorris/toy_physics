
# Moving-Throat PDE — Stage 69: Final Reduced Verdict for the Support/Source Program

## Purpose

This stage closes the present reduced support/source program.

The three crucial ingredients are now all in hand:

1. **Stage 066:** the universal thin-wall matched branch is controlled by one dimensionless wall figure of merit,
   `W_wall`,
   with exact fail/succeed window
   `Pe_req / Delta_inf` to `Pe_req / Delta_0`.
2. **Stage 064:** the parent equilibrium-aligned source/support branch reaches the ideal thin-layer limit
   `C_(sigma phi)^2 = 1`.
3. **Stages 067–068:** the first explicit independent-profile benchmark family has an exact self-dual resonance at
   `w_g / w_f = sqrt(pi)`,
   but even there only reaches
   `C_res^2 = 0.994418836451529...`,
   so its threshold penalty is only
   `P_res = 1.005612487760576...`.

That is enough to give the finish-line verdict for the reduced theorem program.

---

## 1. Universal reduced theorem envelope

The universal matched-branch theorem from Stages 061 and 066 is still the main result:

- if `W_wall <= Pe_req / Delta_inf`, the branch fails;
- if `W_wall >= Pe_req / Delta_0`, the branch succeeds;
- only the band
  `Pe_req / Delta_inf < W_wall < Pe_req / Delta_0`
  still needs the full fixed-point solve.

This is the correct universal reduced theorem because it belongs to the parent equilibrium-matched branch and does not assume an extra independent profile mismatch.

---

## 2. What the explicit sech–Gaussian benchmark adds

The explicit sech–Gaussian benchmark does **not** create a new universal theorem.

What it adds is a very sharp profile-family refinement:

- exact stationary ratio `w_g / w_f = sqrt(pi)`,
- exact coherence symmetry `C^2(r)=C^2(pi/r)`,
- maximal coherence
  `C_res^2 = 0.994418836451529...`,
- penalty factor
  `P_res = 1.005612487760576...`.

So the independent-profile benchmark modifies the universal wall thresholds only by about `0.56%`.

That is the right interpretation of the memo’s usefulness.

---

## 3. Final three-zone verdict

The reduced support/source program is now sharp enough to state a three-zone verdict.

### Zone A — universal failure

If

`W_wall <= Pe_req / Delta_inf,`

the support/source branch fails even on the ideal matched branch.

That is a universal reduced no-go.

### Zone B — universal success on the matched branch

If

`W_wall >= Pe_req / Delta_0,`

the support/source branch succeeds on the equilibrium-matched branch.

That is the universal reduced success theorem.

### Zone C — narrow profile-sensitive band

Only in the intermediate wall window

`Pe_req / Delta_inf < W_wall < Pe_req / Delta_0`

does one still need more branch information.

And even there, the explicit independent sech–Gaussian family changes the matched thresholds only by the tiny factor `P_res = 1.005612...`.

So the only truly profile-sensitive sub-bands are

`Pe_req / Delta_inf < W_wall < P_res Pe_req / Delta_inf`

on the failure side, and

`Pe_req / Delta_0 < W_wall < P_res Pe_req / Delta_0`

on the success side.

Each has width only about `0.56%`.

---

## 4. What is now finished

At the reduced-theorem level, this part is finished.

The support/source program is no longer missing algebra, and it is no longer missing a physically interpretable benchmark family.

It now has:

- a universal matched-branch theorem window,
- an explicit independent-profile resonance family,
- an exact penalty factor relating the two,
- and a precise statement of when profile mismatch can or cannot matter.

So the reduced support/source question is no longer diffuse.

---

## 5. What is still genuinely open

What remains open is no longer more reduced algebra.

The remaining open problem is the actual **moving-throat branch selection**:

- what branch values of `a`, `L`, `ell`, `T_X`, `J_1`, `kappa`, and `eta` does the real throat choose?
- does the true branch behave like the equilibrium-matched source/support law (`C^2 ~ 1`) or like an independent mismatch family (`C^2 < 1`)?
- and then, beyond this reduced support/source problem, does the completed moving-throat PDE realize the already-isolated passive/outgoing quadrupole branch with the correct normalization?

So the next phase is no longer this reduced support/source audit.
The next phase is the explicit moving-throat branch solve and then the final quadrupole-normalization bridge.

---

## 6. Bottom line on the memo

Your memo was helpful.

The right expert reading is:

- the sech–Gaussian resonance is real,
- the `sqrt(pi)` ratio is not numerology but the exact self-dual stationary point of that benchmark family,
- the near-perfect coherence is genuinely strong evidence that source/support mismatch is not the dominant reduced bottleneck,
- but the claim that this **alone** proves threshold survival is too strong.

What it really proves is that a natural explicit independent-profile family comes within `0.56%` of the ideal matched branch.

That is exactly the kind of result that helps us finish the reduced phase cleanly.
