# Parent Throat Action — Branch-B Parameter Burden

## Purpose

`step_29` identified a very promising moderate branch-B corridor:

\[
\sigma=-19 \quad\text{to}\quad \sigma=-49,
\]

with large normalization relief and almost no movement in the isotropic defect.

The next question is whether that corridor is still outrageously expensive when
translated into actual outlet parameters.

This step answers that by mapping the corridor into the exact Stage-95
branch-B parameters and comparing them against the forbidden natural-source
Robin burdens from `step_26`.

The moderate corridor records are imported from `step_29`, so this comparison
does not keep a second local copy of the \(\lambda_{\rm out}=20,50\) patch
points.

---

## 1. Exact branch-B parameter map

On the canonical Stage-95 branch-B family one has

\[
\rho=4\sigma,
\qquad
\kappa=\frac13,
\qquad
\gamma=\frac19,
\qquad
\chi_Q=1.
\]

So the moderate corridor points from `step_29` map to

\[
\sigma=-19 \quad\Rightarrow\quad \rho=-76,
\]
\[
\sigma=-49 \quad\Rightarrow\quad \rho=-196.
\]

These are still significant deformations, but they are not remotely comparable
to the natural-source burdens that appeared when `lambda_out` was interpreted as
\(N_Q\).

The script also checks that the sign of the \(\rho=4\sigma\) map is
load-bearing. Flipping it would move the two moderate-corridor \(\rho\) values
by residuals \(152\) and \(392\), so that sign convention cannot silently
regress.

---

## 2. Comparison against the forbidden natural-source burdens

From `step_26`, the forbidden natural-source interpretation required minimal
Robin burdens of approximately

\[
\rho_{\rm nat}(0.09)\approx -113614.01985490319,
\]
\[
\rho_{\rm nat}(0.092)\approx -137733.1209385474.
\]

So the moderate branch-B points are cheaper by factors of about

\[
\frac{|\rho_{\rm nat}(0.09)|}{76}
\approx 1494.9213138803052,
\]

and

\[
\frac{|\rho_{\rm nat}(0.092)|}{196}
\approx 702.7200047885071.
\]

That is the real quantitative distinction between the viable branch-B story and
the rejected \(N_Q\) story.

---

## 3. Interpretation

This step closes another important ambiguity.

1. The branch-B escape hatch is not only algebraically different from the
   forbidden natural-source interpretation.
2. It is also vastly cheaper in the underlying outlet parameters.
3. The moderate corridor `sigma = -19` to `sigma = -49` now looks like a real
   targetable regime rather than a disguised extreme-load corner.

So the next realization question is no longer whether moderate branch-B points
are absurd on parameter grounds. They are not. The next question is whether a
realized branch can actually generate that moderate outlet family while
remaining compatible with the rest of the PDE-side transport structure.
