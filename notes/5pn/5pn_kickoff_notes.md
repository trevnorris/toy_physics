# 5PN kickoff from the moving-throat PDE stack

## What is now fixed

The first conservative 5PN gate is not a generic coefficient hunt.
Using the grouped real `P2` handoff variables,

- `D_A^(cons)(omega) = D_{A0} + D_{A2} omega^2 + D_{A4} omega^4 + O(omega^6)`,
- `Y_A(omega) = D_{A0} / D_A^(cons)(omega) = 1 + u_2^(A) omega^2 + u_4^(A) omega^4 + ...`,

with

- `u_2^(A) = - D_{A2}/D_{A0}`,
- `u_4^(A) = (D_{A2}^2 - D_{A0} D_{A4}) / D_{A0}^2`,

the exact grouped trace/anomaly variables are

- `xbar = (x20 + 2 x21 + 2 x22)/5`,
- `a_x  = (2 x20 - x21 - x22)/10`,
- `b_x  = (x21 - x22)/2`.

So the 5PN conservative test is:

1. compute `(ubar2, a2, b2)` and `(ubar4, a4, b4)`,
2. test isotropy: `a2 = b2 = a4 = b4 = 0`,
3. on the isotropic branch, test the single-pole identity `u4 = 4 u2^2`.

## What the later moving-throat stages sharpened

If the actual conservative grouped-`P2` branch is carried by one isotropic pole and the
geometry lane is static through `O(omega^4)`, then the conservative quadrupole module is forced to be

`Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`.

Equivalently,

- `K_pole = K0/4`,
- `K_geom = 3 K0/4`,
- `u2 = 1/(4 Omega_Q^2)`,
- `u4 = 1/(4 Omega_Q^4) = 4 u2^2`.

If geometry instead carries dynamic even contamination, the exact obstruction is

`c_pole = (1 + eps4) / (4 (1 + eps2)^2)`

with

- `eps2 = Omega_Q^2 K_(g,2) / K_pole`,
- `eps4 = Omega_Q^4 K_(g,4) / K_pole`.

So any nonzero `eps2` or `eps4` is the clean reduced way that the 5PN conservative check can fail.

## Actual reduced-branch verdict carried forward

On the actual natural isotropic branch now frozen in the reduced hierarchy,

- `eps2 = eps4 = 0`,
- the geometry lane is dynamically inert through `O(omega^4)`,
- the conservative quadrupole carrier is exactly the `3/4 + 1/4` module.

That means the still-live reduced ambiguity has already moved away from the conservative 5PN identity and into the single normalization datum

`N_Q = Kbar_0 / Kbar_0^target`.

With outgoing renormalization included,

`mhat_0^2 chi_Q N_Q = 1`.

On the natural source-map branch `mhat_0 -> 1`, so the remaining point-particle obstruction is purely outgoing:

`N_Q = 1/chi_Q`.

## Diagnostic for non-isotropy

For a weak axisymmetric `Y_20` contamination,

- `x20 = x^(0) + eps x^(1)`,
- `x21 = x^(0) + (eps/2) x^(1)`,
- `x22 = x^(0) - eps x^(1)`,

which forces

- `b = 3 a`.

So if a future PDE extraction produces grouped anisotropy and it does **not** satisfy `b = 3a`, the contamination is not a simple weak axisymmetric branch distortion.

## Immediate next exact step

The next computation should be the overlap-extraction layer for the actual grouped bundle:

1. compute `D_{A0}, D_{A2}, D_{A4}` for `A in {20,21,22}` from the coupled wall/BdG/Maxwell/mixed branch,
2. project those lane values to `(ubar2,a2,b2)` and `(ubar4,a4,b4)`,
3. test isotropy and `u4 = 4 u2^2`,
4. then move directly to the normalization datum `N_Q` rather than reopening the whole conservative branch.
