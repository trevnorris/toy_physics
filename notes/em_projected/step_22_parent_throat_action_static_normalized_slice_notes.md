# Parent Throat Action — Static-Normalized Slice Diagnostic

## Purpose

This step is deliberately **not** another target-blind branch audit.

`step_21` showed that after adding one outgoing-family coordinate, the reduced
family can make

\[
R_{P2},\ R_{P4}
\]

tiny and reduce `R_pole` substantially, but `R_norm` remains the dominant
obstruction. So the next useful question is narrower:

> If one explicitly imposes the static normalization slice `R_norm = 0` inside
> that reduced outgoing-family direction, what obstruction remains?

That is a residual-isolation diagnostic, not a realized-branch theorem.

---

## 1. Status label

The script declares

- `branch_id = v2_local_parent_background_static_normalized_slice`,
- `branch_freeze_hash = 63fe3a626ecf00c7`,
- `pre_target_freeze = false`,
- `target_blind = false`,
- `no_post_residual_refit = false`,
- `declared_slice = R_norm == 0`.

So this artifact should be read as:

- useful for local obstruction analysis,
- not admissible as the main V2 target-blind branch packet.

---

## 2. Underlying outgoing-family direction

The starting point is the exact least-squares direction from `step_21`:

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

The underlying target-blind baseline packet is still

\[
(R_{\rm pole},R_{\rm norm},R_{P2},R_{P4})
\approx
(-13.134593938872369,\,-10.33719584868593,\,
0.3700984456976848,\,
0.8889149882257383).
\]

---

## 3. Static-normalized root

The script follows the `step_21` outgoing-family direction and brackets the
one-pole sign change on the declared static-normalized slice in

\[
\text{scale}\in(0.09,\ 0.092).
\]

After bisection, it finds

\[
\text{root scale} \approx 0.091265380859375.
\]

At that point the required static-normalizing source factor is

\[
\widehat m_0^{\rm static}\approx 206.8317356734255.
\]

The resulting isotropic packet on this declared slice is

\[
(R_{\rm pole},R_{\rm norm},R_{P2},R_{P4})
\approx
(8.45325692733212\times10^{-5},\ 0,\ -1.5242626961102732\times10^{-4},\ 1.0887180997006705\times10^{-4}),
\]

with total packet norm

\[
\|R\|_2 \approx 2.055057029418912\times10^{-4}.
\]

So on this declared slice the isotropic packet can be made numerically tiny.

---

## 4. Interpretation

This sharpens the reduced-branch story considerably.

1. The target-blind 3-coordinate family in `step_20` had a real local rank
   obstruction.
2. The target-blind 4-coordinate outgoing family in `step_21` removed that
   first-order obstruction and made `R_{P2},R_{P4}` nearly vanish on actual
   finite steps, but `R_norm` stayed large.
3. The non-target-blind static-normalized slice in this step shows that once
   `R_norm` is forced to zero, the remaining isotropic packet can be made
   extremely small.

So the main remaining issue in this reduced neighborhood is no longer algebraic
compatibility of the isotropic packet. It is the price of the required source
normalization:

\[
\widehat m_0 \sim 2\times10^2.
\]

That is why this step is diagnostic rather than celebratory. It says:

- the reduced family is close to isotropic closure in residual space,
- but only on a slice that currently asks for an implausibly large normalization.

That is exactly the kind of distinction the PDE program needs.
