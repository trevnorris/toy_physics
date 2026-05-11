# Moving-Throat PDE — Stage 093: Updated Status After the Direct Grouped-`P2` + Geometry Derivation

## Purpose

Stages 74–75 close the exact next step that was left open after Stage 73.

They do two things:

1. derive the minimal isotropic `3/4 + 1/4` conservative quadrupole module directly from the grouped-`P2` + geometry split;
2. isolate the exact obstruction formula if the geometry lane carries dynamic `omega^2` or `omega^4` moments.

This note records the status update.

---

## 1. What is now settled

Inside the reduced hierarchy, the following implication is now exact:

- if the actual isotropic grouped-`P2` conservative branch is one pole,
- and if the geometry lane is static through `O(omega^4)`,

then the conservative quadrupole module is forced to be

`Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`,

and therefore

`rho_alpha = 4/3`,

`zeta_req = 1/3`.

So the Stage-71/72 Family-1 verdict is now a direct corollary of the grouped-`P2` + static-geometry realization rather than an additional standalone assumption.

---

## 2. What remains genuinely open

The remaining reduced theorem gap is now extremely narrow.

It is no longer support-source sufficiency.
It is no longer the existence of the minimal isotropic module in the abstract.
It is now exactly this:

> does the real moving-throat branch realize
>
> - one isotropic grouped-`P2` pole,
> - and a geometry lane that is static through `O(omega^4)`?

If yes, the `3/4 + 1/4` split and the `rho_alpha = 4/3` verdict follow immediately.
If not, the exact deviation is still controlled by the Stage-75 obstruction formula

`c_pole = (1 + eps_4) / [ 4 (1 + eps_2)^2 ]`.

So the next real PDE-side task is no longer broad. It is to determine the two geometry-contamination numbers `eps_2` and `eps_4`, or prove that they vanish on the natural branch.

---

## 3. Best current expert verdict after Stage 76

At this point the reduced program is in its sharpest state yet.

- The explicit Family-1 support/source side is finished.
- The minimal isotropic quadrupole module is no longer just imported from the outgoing-moment analysis; it is directly recovered from the grouped-`P2` + static-geometry split.
- The only remaining reduced ambiguity is whether the real geometry lane is dynamically inert through `O(omega^4)` on the natural branch.

So the next derivation phase should be aimed squarely at the geometry lane, not at reopening the already-solved support/source side.
