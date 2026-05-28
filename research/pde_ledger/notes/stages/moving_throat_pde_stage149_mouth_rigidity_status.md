
# Moving-Throat PDE — Stage 149: Non-Exponential Mouth-Rigidity Status

The mouth-side branch problem is now significantly tighter than it was at Stage 145.

## What is now explicit

1. Any positive normalized mouth deformation near the canonical exponential branch can be
   represented as
   \[
   \Sigma_\epsilon=(1-\epsilon)\Sigma_*+\epsilon\varsigma,
   \qquad
   \varsigma\ge0,
   \qquad
   \int_0^1\varsigma=1.
   \]

2. At first order, the entire deformation enters only through the two source moments
   \[
   \bar g_\varsigma=\int_0^1\varsigma\,c,
   \qquad
   \bar S_\varsigma=\int_0^1\varsigma\,K_q.
   \]

3. Retuning the electrochemical bias to stay on the canonical lower branch gives
   \[
   \delta\Pi
   =
   -\epsilon\frac{\bar g_\varsigma-\mathfrak g_*}{\mathfrak g'_*},
   \]
   and the physical traction shift reduces to
   \[
   \delta \widehat T_m
   =
   \epsilon\Big[
   A_T(\bar g_\varsigma-\mathfrak g_*)
   +
   B_T(\bar S_\varsigma-\mathcal S_*)
   \Big],
   \]
   with
   \[
   A_T\approx -4.27263956256927,
   \qquad
   B_T\approx 0.134875005736706.
   \]

4. So the canonical branch is controlled primarily by the overlap channel, not by the
   mixed-kernel channel:
   \[
   |A_T|/B_T\approx 31.68.
   \]

5. For representative positive non-exponential families at this same first-order
   deformation level:
   - uniform broadening raises the canonical point,
   - derivative matching lowers it,
   - and the zero-shift mixture coincides almost exactly with the earlier Stage 126
     exact compensation broadening fraction.

## Updated interpretation

Inside the explicit Family-1 positive mouth-layer closure, and at the same
first-order deformation level used in Stages 146-148, the mouth-side ambiguity is now
no longer a branch-selection problem and no longer a large shape-space uncertainty.

It has been reduced to:

1. a single regular canonical branch inside the explicit Family-1 closure,
2. one explicit rigidity kernel \(\mathcal W_*(x)\),
3. and perturbative finite shifts under smooth positive non-exponential deformations.

The stage remains a perturbative rigidity summary, not yet a full nonlinear
mouth theorem for arbitrary finite deformations.

So the remaining PDE-facing question is not “which mouth branch?” but rather:

\[
\boxed{
\text{what small non-exponential correction does the real moving-throat mouth layer induce around }(\Pi_*,\widehat T_{m,*})?
}
\]

That is a much smaller target than the earlier mouth-source ambiguity.
