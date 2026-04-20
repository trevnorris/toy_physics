
# Moving-Throat PDE — Stage 119: Mouth-Layer Fixed-Point Status After the Coupled Solve

## What is now fixed

After Stages 116–118, the mouth-source selection problem has narrowed again.

1. The actual GNLS + localized-Maxwell mouth layer is no longer described by an
   effective slope alone. It is governed by the exact coupled fixed-point law
   \[
   \Pi
   =
   M_+\mathcal S(\Pi,\kappa_+)+M_-\mathcal S(\Pi,\kappa_-).
   \]

2. On the first explicit Family-1 branch with one static shell lane and one mixed
   D/N half-wave lane,
   \[
   \Pi = M_s + M_q\,\mathcal S_q(\Pi),
   \qquad
   \kappa_s=0,
   \qquad
   \kappa_q=\frac{\pi}{2}.
   \]

3. The canonical source-shape compensation point is now an exact gain line:
   \[
   M_s \approx 1.50882951349316 - 0.658075937605429\,M_q.
   \]

4. If the mouth-layer gains inherit the same \(4:-1\) shell-to-mixed weighting as
   the nontrivial compensated throat-core outlet, the selection law becomes
   \[
   \Pi
   =
   \Sigma_m\left[4-\mathcal S_q(\Pi)\right],
   \]
   with canonical Family-1 value
   \[
   \Sigma_m^* \approx 0.451485277739090.
   \]

## Meaning

The source-shape problem is no longer:

- which positive family to use,
- which sign branch is physical,
- or what effective combination
  \( \partial_z\delta V_{{\rm conf}}-q_*\partial_zA_0 \)
  should be guessed.

It is now simply:

\[
\boxed{
\text{what dimensionless mouth gain pair }(M_s,M_q)\text{ — or, under the
outlet-consistent closure, what single gain }\Sigma_m\text{ — does the real
moving-throat mouth layer select?}
}
\]

So the open microscopic bias problem has collapsed from a free profile question to
a small, explicit fixed-point law.
