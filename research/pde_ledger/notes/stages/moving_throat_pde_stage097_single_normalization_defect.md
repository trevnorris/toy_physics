# Moving-Throat PDE — Stage 097: The Actual Isotropic Passive/Outgoing Branch Collapses to a Single Normalization Defect

## Purpose

Stage 79 finished the last serious **conservative** ambiguity on the actual isotropic branch:

- the grouped real `P2` carrier is the only dynamic quadrupole lane,
- the geometry lane is static through `O(omega^4)`,
- and the conservative module is therefore

`Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`.

The next honest step is to combine that actual conservative result with the already-frozen 2.5PN outgoing theorem target.

This stage shows that, once the actual branch is both

1. isotropic,
2. grouped-`P2` one-pole,
3. and passively/outgoingly completed in the natural `l=2` way,

then the entire reduced 2.5PN normalization problem collapses to **one scalar defect**.

---

## 1. Actual branch conservative module

Write the actual isotropic conservative quadrupole response in canonical invariant form as

`Kbar_Q^cons(omega)`
`= Kbar_0 [ 3/4 + (1/4)/(1 - omega^2/Omega_Q^2) ]`.

Then its low-frequency coefficients are exactly

`Kbar_2 = Kbar_0 / (4 Omega_Q^2),`

`Kbar_4 = Kbar_0 / (4 Omega_Q^4).`

So once `Kbar_0` and `Omega_Q` are known, the entire even conservative ledger is fixed.

---

## 2. Outgoing odd coefficient on the same branch

On the minimal isotropic outgoing `l=2` branch, the 2.5PN audit already fixed the odd coefficient algebraically as

`Gammabar_5 = 9 Kbar_2^(5/2) / Kbar_0^(3/2)`.

Substituting the actual one-pole conservative relations gives

`Gammabar_5 = 9 Kbar_0 / (32 Omega_Q^5).`

So the odd Burke–Thorne coefficient is not independent either.
It is already determined by the same two quantities `(Kbar_0, Omega_Q)`.

---

## 3. The GR target branch

The GR target on the same isotropic outgoing branch is

`Kbar_0^target = 64 G Omega_Q^5 / (45 c^5),`

equivalently, after using the already-carried geometric pole

`Omega_Q = 3 c_s / (2 a)`,

it becomes

`Kbar_0^target = 54 G c_s^5 / (5 a^5 c^5).`

Then automatically

`Kbar_2^target = Kbar_0^target / (4 Omega_Q^2),`

`Kbar_4^target = Kbar_0^target / (4 Omega_Q^4),`

`Gammabar_5^target = 2 G / (5 c^5).`

So the actual branch and the GR target have exactly the same algebraic structure.
The only possible mismatch is the overall normalization of `Kbar_0`.

---

## 4. Single normalization defect

Define the actual branch normalization defect by

`N_Q := Kbar_0 / Kbar_0^target`.

Then all three low-frequency targets scale by the same factor:

`Kbar_2 = N_Q Kbar_2^target,`

`Kbar_4 = N_Q Kbar_4^target,`

`Gammabar_5 = N_Q Gammabar_5^target = N_Q * 2 G / (5 c^5).`

Therefore the actual isotropic passive/outgoing branch satisfies the full reduced GR-like point-particle 2.5PN theorem **iff**

`N_Q = 1`.

Equivalently, the following defects are all the same number:

`R_0 := Kbar_0 / Kbar_0^target - 1,`

`R_2 := Kbar_2 / Kbar_2^target - 1,`

`R_4 := Kbar_4 / Kbar_4^target - 1,`

`R_5 := Gammabar_5 / (2G/(5c^5)) - 1,`

with

`R_0 = R_2 = R_4 = R_5 = N_Q - 1`.

So the final reduced theorem gap is now a **single scalar normalization defect**.

---

## 5. What this means physically

At this point, the reduced program has separated into two sharply different questions:

1. **support/source sufficiency** of the explicit branch,
2. **radiative normalization** of the actual outgoing quadrupole module.

The present stage shows that the second question is no longer a multi-parameter branch-selection problem.
Once the actual isotropic grouped-`P2` one-pole structure is accepted, it is only the one-number defect `N_Q`.

That is as narrow a reduced theorem gate as one can reasonably ask for before solving the full moving-throat PDE.
