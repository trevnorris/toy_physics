# Moving-Throat PDE — Stage 102: At 2.5PN the Only Live Retarded Obstruction is the Leading `omega^5` Outgoing Normalization

## Purpose

Stage 84 reduced the last finite reduced 2.5PN question to one outgoing-normalization factor `chi_Q`.

A natural objection remains:

> what if the true moving-throat PDE contains extra retarded structure beyond the canonical one-pole `l=2` denominator?

This stage shows that, at the level of the 2.5PN theorem, any extra retarded structure that first enters at `O(omega^7)` or higher is irrelevant.

So the only live retarded obstruction at 2.5PN is the leading `omega^5` outgoing-normalization factor `chi_Q`.

---

## 1. Generalized one-pole denominator with a higher odd tail

Take the most general one-pole grouped-`P2` denominator consistent with the already-fixed conservative branch but allowing one extra higher odd term:

`Yhat_Q^ret(omega)`
`= 3/4 + (1/4)/(1 - omega^2/Omega_Q^2 - i chi_Q sigma_Q^can omega^5 - i tau_Q omega^7) + O(omega^8).`

Expanding through `O(omega^5)` gives

`Yhat_Q^ret(omega)`
`= 1 + omega^2/(4 Omega_Q^2) + omega^4/(4 Omega_Q^4)`
`  + i chi_Q 9/(32 Omega_Q^5) omega^5 + O(omega^6).`

The extra `tau_Q` term does not appear until `O(omega^7)`.

So at 2.5PN order, the only retarded parameter that matters is `chi_Q`.

---

## 2. Consequence for the theorem ledger

Combining Stages 83–84 with the `O(omega^7)` irrelevance statement gives the exact reduced 2.5PN hierarchy:

1. `eps_2 = eps_4 = 0` kills the conservative geometry contamination,
2. `rho_alpha = 4/3` and the Family-1 analysis make support/source automatic,
3. `mhat_0 = 1 + O(a^2/r^2)` removes the source-map ambiguity on the natural branch,
4. every higher odd retarded coefficient starting at `omega^7` is invisible to the 2.5PN theorem,
5. so the only remaining reduced obstruction is

`Delta_Q = chi_Q - 1`.

If `Delta_Q = 0`, the reduced point-particle 2.5PN theorem is closed.

---

## 3. Best current finish-line statement

Inside the present reduced hierarchy, the moving-throat PDE program is now reduced to one sharply isolated PDE-facing question:

> does the actual passive/outgoing grouped-`P2` quadrupole branch have the canonical compact outgoing `omega^5` coefficient, i.e. `chi_Q = 1`?

Everything else that could have obstructed the reduced 2.5PN theorem has either been fixed or pushed above the relevant order.
