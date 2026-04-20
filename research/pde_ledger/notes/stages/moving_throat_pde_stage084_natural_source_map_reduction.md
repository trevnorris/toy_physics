# Moving-Throat PDE — Stage 84: On the Natural Source-Map Branch the Last Reduced 2.5PN Obstruction is Purely Outgoing

## Purpose

Stage 83 showed that the full odd normalization condition factorizes exactly as

`mhat_0^2 chi_Q N_Q = 1`.

The next honest step is to combine that with the already-carried source-map result of the natural orbital/worldtube STF branch:

`mhat_0 = 1 + O(a^2/r^2)`.

This stage shows that in the strict point-particle limit the remaining reduced 2.5PN obstruction is no longer a mixed source/outgoing product. It is purely the outgoing-normalization factor `chi_Q`.

---

## 1. Point-particle natural branch

On the natural source-map branch isolated by the 2.5PN audit, one has

`mhat_0 = 1 + O(a^2/r^2)`.

So in the strict point-particle limit,

`mhat_0 -> 1`.

Then the exact Stage-83 factorization becomes

`N_Q = 1/chi_Q`.

So:

- if `chi_Q = 1`, then `N_Q = 1`,
- if `chi_Q != 1`, then the entire reduced 2.5PN mismatch is just the inverse outgoing renormalization.

---

## 2. Canonical compact outgoing branch

Stage 80 already fixed the unique minimal passive/outgoing grouped-`P2` one-pole completion that matches the exact compact outgoing `l=2` fingerprint:

`Yhat_Q^ret(omega)`
`= 3/4 + (1/4)/(1 - omega^2/Omega_Q^2 - i sigma_Q^can omega^5) + O(omega^6),`

with

`sigma_Q^can = 9/(8 Omega_Q^5) = 4 a^5/(27 c_s^5)`.

By definition, that canonical compact outgoing branch has

`chi_Q = 1`.

Therefore on the natural source-map branch,

`N_Q = 1`.

So the reduced 2.5PN theorem closes automatically **if** the actual passive/outgoing moving-throat branch is exactly the canonical compact one-pole `l=2` completion.

---

## 3. Small deviation parameter

Define the outgoing-normalization defect by

`Delta_Q := chi_Q - 1`.

Then on the natural source-map branch,

`N_Q = 1/(1 + Delta_Q)`

and therefore

`N_Q - 1 = -Delta_Q/(1 + Delta_Q)`.

For a small outgoing deviation,

`N_Q - 1 = -Delta_Q + O(Delta_Q^2).`

So the last reduced theorem gap is linearly controlled by the first outgoing-normalization defect.

---

## 4. Meaning

At this point the reduced hierarchy has done everything it can without solving the actual passive/outgoing DtN problem.

What remains is no longer a broad question about conservative structure, support sufficiency, or source-map ambiguity. It is just this:

> does the actual passive/outgoing moving-throat quadrupole branch realize `chi_Q = 1`, or does it carry a nontrivial outgoing-normalization defect `Delta_Q`?

That is the cleanest reduced 2.5PN finish line reached so far.
