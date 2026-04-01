# Moving-Throat PDE — Stage 86: Conditional Reduced 2.5PN Closure

## Statement

Inside the present reduced hierarchy, the following are now fixed:

- the actual isotropic grouped-`P2` conservative branch is geometry-clean through `O(omega^4)`,
- the conservative quadrupole split is exactly `3/4 + 1/4`,
- the selected support/source ratio is `rho_alpha = 4/3`,
- the explicit Family-1 support/source branch is automatically sufficient,
- the natural source-map branch gives `mhat_0 = 1 + O(a^2/r^2)`,
- and all higher odd retarded data beginning at `O(omega^7)` are irrelevant to the point-particle 2.5PN theorem.

Therefore the reduced 2.5PN theorem is conditionally closed by one and only one remaining branch datum:

`chi_Q`.

More precisely:

- if the actual passive/outgoing grouped-`P2` quadrupole branch satisfies `chi_Q = 1`,
  then the reduced GR-like point-particle 2.5PN theorem is closed;
- if `chi_Q != 1`, then the entire remaining reduced failure is measured by

`Delta_Q := chi_Q - 1`.

So the remaining PDE-facing problem is no longer “derive 2.5PN somehow.”
It is:

> compute the leading outgoing `omega^5` normalization of the actual passive/outgoing grouped-`P2` quadrupole branch and determine whether it equals the canonical compact outgoing value.
