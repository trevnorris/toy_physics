# 5PN / moving-throat continuation — Stage 161–170 bundle

This bundle fills the previously missing numbered continuation after Stage 160.

## What is in this bundle

The scripts are split one stage at a time so the executable numbering now matches the note chain:

- Stage 161 — outgoing-port co-loading theorem
- Stage 162 — wall-normalized transfer-shape theorem
- Stage 163 — effective transfer-shape collapse
- Stage 164 — coherent tracking-branch defect law and support-blindness
- Stage 165 — microscopic coherent-kernel slippage decomposition
- Stage 166 — exact triangular normal form
- Stage 167 — branch-invariant coordinates
- Stage 168 — direct microscopic monomials and compatibility ledger
- Stage 169 — exact microscopic similarity orbit
- Stage 170 — exact orbit–quotient closure

Each script has a paired `_output.txt` file captured from a clean run in this environment.

## Structural summary

The chain now reads:

1. `Xi_1` is the mismatch between the outgoing-weighted static transfer slope and the wall-baseline slope.
2. That mismatch is exactly twice the weighted transfer-shape slope.
3. The whole many-port problem collapses to one effective transfer shape `T_eff^2 = N_0 / K`.
4. On the coherent branch, the support ratio drops out of the weak-axisymmetric defect exactly.
5. The remaining defect depends only on microscopic mixed/outgoing placement drifts.
6. Those drifts collapse to the three branch-adapted coordinates
   `Sigma_tr`, `Sigma_nt`, `Sigma_eta`.
7. Those coordinates are the logarithmic drifts of three exact direct microscopic monomials
   `C_tr,*`, `C_nt,*`, and `epsilon_eta`.
8. Their zero-defect set is the tangent space of an exact five-parameter multiplicative similarity orbit `G_*`.
9. At the finite level, the level sets of `(C_tr,*, C_nt,*, epsilon_eta)` are exactly the `G_*` orbits.
10. So the remaining microscopic question is purely branch-selective: whether the true moving-throat branch preserves those three quotient coordinates.

## Practical continuation point

The smallest honest next theorem gate after Stage 170 is:

- compute the actual branch drift of the three quotient coordinates
  `(C_tr,*, C_nt,*, epsilon_eta)`
  from the real moving-throat PDE branch;
- equivalently, determine whether the real branch stays on a single `G_*` similarity orbit.

If it does, the coherent first-order grouped weak-axisymmetric defect vanishes automatically.