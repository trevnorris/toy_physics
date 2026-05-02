# Parent Throat Action — Outgoing-Amplitude Normalization Frontier

## Purpose

`step_23` showed that on the frozen `step_21` outgoing-family ray, the isotropic
defect

\[
Q_{\rm iso}:=\sqrt{R_{\rm pole}^2+R_{P2}^2+R_{P4}^2}
\]

can be reduced only at the price of rapidly increasing the required static
normalization \(\widehat m_0^{\rm req}\).

This step tests one more **target-blind upstream coordinate**:

an outgoing-transfer amplitude scale \(\lambda_{\rm out}\) multiplying the whole
outgoing packet.

So the question is:

> if the branch is allowed a target-blind outgoing amplitude coordinate, does
> the normalization frontier remain bad?

---

## 1. Frozen setup

The script declares

- `branch_id = v2_local_parent_background_outgoing_amplitude_frontier_scan`,
- `branch_freeze_hash = b7026d6357b399bf`,
- `pre_target_freeze = true`,
- `target_blind = true`,
- `no_post_residual_refit = true`.

It reuses the frozen `step_21` least-squares ray

\[
\delta
\approx
(8.677800811209,\ 2.872115918746,\ 53.640446430441,\ 74.190514652389),
\]

and keeps the same sample scales

\[
\{0,\ 0.03,\ 0.05,\ 0.08,\ 0.088,\ 0.09,\ 0.092\}.
\]

The new outgoing amplitude grid is

\[
\lambda_{\rm out}\in
\{1,\ 5,\ 20,\ 50,\ 100,\ 200,\ 500,\ 1000,\ 2000\}.
\]

At fixed scale, this coordinate rescales

\[
(P_0,P_2,P_4)\mapsto
(\lambda_{\rm out}P_0,\ \lambda_{\rm out}P_2,\ \lambda_{\rm out}P_4),
\]

while \(R_{\rm pole}\) is unchanged.

So it is an upstream outgoing-family coordinate, not a post-target repair.

---

## 2. Best sampled isotropic defect at normalization budgets

On the frozen sampled grid, the best `Q_iso` values under explicit
normalization budgets are:

\[
\widehat m_0^{\rm req}\le 50
\quad\Rightarrow\quad
(\text{scale},\lambda_{\rm out})
=(0.092,\ 20),
\]
\[
\widehat m_0^{\rm req}\approx 47.912441136331765,
\qquad
Q_{\rm iso}\approx 0.2670781822626949.
\]

\[
\widehat m_0^{\rm req}\le 20
\quad\Rightarrow\quad
(\text{scale},\lambda_{\rm out})
=(0.092,\ 200),
\]
\[
\widehat m_0^{\rm req}\approx 15.151244224955443,
\qquad
Q_{\rm iso}\approx 0.2693266571833183.
\]

\[
\widehat m_0^{\rm req}\le 10
\quad\Rightarrow\quad
(\text{scale},\lambda_{\rm out})
=(0.092,\ 500),
\]
\[
\widehat m_0^{\rm req}\approx 9.582488227266353,
\qquad
Q_{\rm iso}\approx 0.28094980884298526.
\]

\[
\widehat m_0^{\rm req}\le 5
\quad\Rightarrow\quad
(\text{scale},\lambda_{\rm out})
=(0.092,\ 2000),
\]
\[
\widehat m_0^{\rm req}\approx 4.7912441136331765,
\qquad
Q_{\rm iso}\approx 0.4394839373049669.
\]

So even with a very small required normalization, the isotropic defect on this
sampled grid stays comfortably below `1`.

---

## 3. Best sampled normalization at defect thresholds

Looking the other way around:

\[
Q_{\rm iso}\le 1
\quad\Rightarrow\quad
(\text{scale},\lambda_{\rm out})
=(0.09,\ 2000),
\]
\[
\widehat m_0^{\rm req}\approx 4.351570977913287,
\qquad
Q_{\rm iso}\approx 0.618690285150578.
\]

\[
Q_{\rm iso}\le 0.5
\quad\Rightarrow\quad
(\text{scale},\lambda_{\rm out})
=(0.092,\ 2000),
\]
\[
\widehat m_0^{\rm req}\approx 4.7912441136331765,
\qquad
Q_{\rm iso}\approx 0.4394839373049669.
\]

This is a major change relative to the one-coordinate frontier in `step_23`,
where `Q_iso < 1` first appeared only once
\[
\widehat m_0^{\rm req}\sim 1.95\times 10^2.
\]

---

## 4. Interpretation

This step materially changes the diagnosis.

1. On the frozen outgoing-family ray alone, reducing `Q_iso` demanded a very
   large required normalization.
2. Once a target-blind outgoing amplitude coordinate is added, that frontier
   changes dramatically.
3. On the sampled grid, `Q_iso <= 1` is compatible with
   \[
   \widehat m_0^{\rm req}\approx 4.35,
   \]
   and `Q_iso <= 0.5` with
   \[
   \widehat m_0^{\rm req}\approx 4.79.
   \]

So static normalization no longer looks structurally fatal **inside this
enlarged reduced outgoing family**.

That does **not** yet prove realized-branch success. The next question is now
sharper:

> is such a large outgoing-amplitude enhancement physically admissible as a
> genuine branch-derived export, or is it only a reduced-family convenience?

That is the correct next frontier.
