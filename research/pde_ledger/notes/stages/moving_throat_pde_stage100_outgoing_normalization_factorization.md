# Moving-Throat PDE — Stage 100: Exact Factorization of the Last 2.5PN Defect into Conservative and Outgoing Pieces

## Purpose

Stage 82 reduced the full isotropic point-particle 2.5PN problem to one scalar defect

`N_Q := Kbar_0 / Kbar_0^target`.

That statement was made on the *minimal passive/outgoing grouped-`P2` branch*, where the outgoing `l=2` odd coefficient already matches the canonical compact value. The next honest step is to relax that one assumption slightly and ask:

> if the actual passive/outgoing moving-throat branch is still a one-pole grouped-`P2` branch, but its leading outgoing odd coefficient is renormalized, how does the final defect factorize?

This stage shows that the answer is exact.

Introduce one dimensionless outgoing-normalization factor `chi_Q` by writing the actual normalized retarded branch as

`Yhat_Q^ret(omega)`
`= 3/4 + (1/4)/(1 - omega^2/Omega_Q^2 - i chi_Q sigma_Q^can omega^5) + O(omega^6),`

where the canonical compact outgoing coefficient is

`sigma_Q^can = 9/(8 Omega_Q^5) = 4 a^5/(27 c_s^5)`.

Then the invariant low-frequency tuple is

`Kbar_2 = Kbar_0/(4 Omega_Q^2),`

`Kbar_4 = Kbar_0/(4 Omega_Q^4),`

`Gammabar_5 = chi_Q * 9 Kbar_0/(32 Omega_Q^5).`

So the even conservative defect and the odd retarded defect separate cleanly:

- the even moments depend only on `N_Q`,
- the odd Burke–Thorne coefficient depends on `chi_Q N_Q`.

That is the exact factorization of the last 2.5PN obstruction.

---

## 1. Conservative target ratios remain scalar

Define the canonical invariant GR targets on the isotropic branch by

`Kbar_0^target = 64 G Omega_Q^5/(45 c^5),`

`Kbar_2^target = Kbar_0^target/(4 Omega_Q^2),`

`Kbar_4^target = Kbar_0^target/(4 Omega_Q^4).`

Then with

`N_Q := Kbar_0/Kbar_0^target`,

the exact even ratios are still

`Kbar_2/Kbar_2^target = N_Q,`

`Kbar_4/Kbar_4^target = N_Q.`

So the conservative/even side remains one-number clean.

---

## 2. Odd branch ratio picks up exactly one extra factor

The 2.5PN target odd coefficient is

`Gammabar_5^target = 2 G/(5 c^5)`.

With the renormalized outgoing branch,

`Gammabar_5 = chi_Q * 9 Kbar_0/(32 Omega_Q^5)`

and therefore

`Gammabar_5 / Gammabar_5^target = chi_Q N_Q.`

So all deviations from the exact compact outgoing `l=2` fingerprint are captured by the single multiplier `chi_Q`.

---

## 3. Observable point-particle normalization condition

The actual point-particle observable branch includes the source-map factor `mhat_0`, so the full odd normalization condition is

`mhat_0^2 Gammabar_5 = 2 G/(5 c^5)`.

Substituting the renormalized one-pole branch gives

`mhat_0^2 chi_Q N_Q = 1.`

Equivalently,

`N_Q = 1/(mhat_0^2 chi_Q).`

This is the exact factorized form of the last reduced 2.5PN defect.

---

## 4. Meaning

Stage 82 said the remaining reduced theorem gap was one scalar normalization defect `N_Q - 1`.

This stage refines that statement:

- if the branch is only *conservatively* known, the even defect is `N_Q - 1`;
- if the branch is *retarded* but not yet proven canonical, the full 2.5PN odd defect is `mhat_0^2 chi_Q N_Q - 1`;
- and all genuinely new retarded uncertainty sits in one number only:

`chi_Q`.

So the problem is no longer “some unknown outgoing structure.” It is exactly the leading outgoing-normalization factor of the actual grouped-`P2` one-pole branch.
