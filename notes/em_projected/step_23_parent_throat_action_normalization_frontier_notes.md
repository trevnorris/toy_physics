# Parent Throat Action — Target-Blind Normalization Frontier

## Purpose

`step_22` was deliberately non-target-blind: it imposed a static-normalized
slice and then isolated the remaining isotropic packet.

This step returns to the V2 branch discipline:

- `pre_target_freeze = true`,
- `target_blind = true`,
- `no_post_residual_refit = true`.

Its job is to answer a narrower but important question:

> Along the best frozen outgoing-family ray found in `step_21`, how rapidly
> does the isotropic packet improve as the required static normalization
> increases?

So this note is a target-blind **normalization frontier**, not a repaired
branch.

---

## 1. Frozen setup

The script declares

- `branch_id = v2_local_parent_background_normalization_frontier_scan`,
- `branch_freeze_hash = 2243249329a9d5ab`,
- `boundary_class = open_impedance_demo`.

It reuses the `step_21` outgoing-family least-squares direction

\[
\delta
\approx
(8.677800811209,\ 2.872115918746,\ 53.640446430441,\ 74.190514652389),
\]

in the coordinates

\[
(\log R_0\text{ amplitude},\ \log R_0\text{ inverse-width},\ \log \beta\text{
 inverse-width},\ \delta_N).
\]

The sampled ray is frozen in advance:

\[
\text{scale}\in\{0,\ 0.03,\ 0.05,\ 0.08,\ 0.088,\ 0.09,\ 0.092\}.
\]

At each sample point the script records:

- the isotropic packet
  \[
  (R_{\rm pole},R_{\rm norm},R_{P2},R_{P4}),
  \]
- the reduced isotropic defect
  \[
  Q_{\rm iso}:=\sqrt{R_{\rm pole}^2+R_{P2}^2+R_{P4}^2},
  \]
- the full residual norm
  \[
  \|R\|_2:=\sqrt{R_{\rm pole}^2+R_{\rm norm}^2+R_{P2}^2+R_{P4}^2},
  \]
- and the required static normalization
  \[
  \widehat m_0^{\rm req}:=\sqrt{\frac{54/5}{P_0}},
  \qquad P_0>0.
  \]

---

## 2. Sampled frontier

The baseline point is

\[
\widehat m_0^{\rm req}\approx 4.8307354536946585,
\qquad
Q_{\rm iso}\approx 13.169840570593973,
\]
\[
\|R\|_2\approx 16.742231591665213.
\]

The frozen sampled ray then gives:

### Scale `0.03`

\[
\widehat m_0^{\rm req}\approx 10.612856360546363,
\qquad
Q_{\rm iso}\approx 12.057319412929672,
\qquad
\|R\|_2\approx 16.123181731091783.
\]

### Scale `0.05`

\[
\widehat m_0^{\rm req}\approx 27.467490672562466,
\qquad
Q_{\rm iso}\approx 10.043783258330908,
\qquad
\|R\|_2\approx 14.73799806685835.
\]

### Scale `0.08`

\[
\widehat m_0^{\rm req}\approx 119.94850420735051,
\qquad
Q_{\rm iso}\approx 3.6635530969600296,
\qquad
\|R\|_2\approx 11.40374534724855.
\]

### Scale `0.088`

\[
\widehat m_0^{\rm req}\approx 176.7200985570496,
\qquad
Q_{\rm iso}\approx 1.1434078630513964,
\qquad
\|R\|_2\approx 10.860014360884977.
\]

### Scale `0.09`

\[
\widehat m_0^{\rm req}\approx 194.6081703105869,
\qquad
Q_{\rm iso}\approx 0.4513177752288337,
\qquad
\|R\|_2\approx 10.809140954536216.
\]

### Scale `0.092`

\[
\widehat m_0^{\rm req}\approx 214.2709506975902,
\qquad
Q_{\rm iso}\approx 0.26705543084121786,
\qquad
\|R\|_2\approx 10.803066122095556.
\]

Every sampled point lies on the Pareto frontier of the frozen ray: moving
forward decreases `Q_iso` but increases the required normalization.

---

## 3. Key diagnostics

Two numbers matter most.

First, the **first** sampled point where

\[
Q_{\rm iso}<1
\]

is

\[
\text{scale}=0.09,
\qquad
\widehat m_0^{\rm req}\approx 194.6081703105869,
\qquad
Q_{\rm iso}\approx 0.4513177752288337.
\]

Second, the **best** sampled `Q_iso` value is

\[
\text{scale}=0.092,
\qquad
\widehat m_0^{\rm req}\approx 214.2709506975902,
\qquad
Q_{\rm iso}\approx 0.26705543084121786.
\]

So on this target-blind outgoing-family ray:

- the isotropic packet can be made much smaller,
- but only by paying a very steep normalization cost.

---

## 4. Interpretation

This sharpens the branch-realization picture again.

1. The enlarged outgoing family is not failing because `R_{P2}` and `R_{P4}`
   refuse to cooperate; they become tiny.
2. The main target-blind tradeoff is now explicit:
   reducing `Q_iso` drives the required normalization rapidly upward.
3. On the sampled ray,
   \[
   Q_{\rm iso}<1
   \]
   appears only once
   \[
   \widehat m_0^{\rm req}\sim 2\times 10^2.
   \]

So the next serious question is no longer “can this reduced outgoing family
approach the isotropic packet?” It can.

The real next question is:

> can one add a physically motivated branch coordinate that lowers
> \(\widehat m_0^{\rm req}\) substantially without reopening the
> `R_{P2}` / `R_{P4}` packet?

That is the correct next target.  
