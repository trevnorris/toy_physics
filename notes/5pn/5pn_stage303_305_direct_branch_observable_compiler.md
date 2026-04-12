# 5PN continuation notes — Stages 303–305

These stages take the Stage 300–302 continuum-ratio orbit and push it one step closer to the actual moving-throat branch observables.

The conceptual simplification is stronger than I expected.
After using the exact coherent-branch identities already established in the moving-throat program,
the first-order weak-axisymmetric orbit-lock test no longer needs the full five-variable reduced ratio state

a) `(chi_0, delta_U, epsilon_eta, Z_W, Lambda)`,

and it does not even need the finite quotient packet `(q_tr, q_nt, q_eta)` as a primary object.

It can be written directly in the three exact branch observables

a) `R_tr`,

b) `R_target`,

c) `epsilon_eta`.

That is now the cleanest actual-branch chart for the coherent weak-axisymmetric endgame.

## Stage 303 — exact finite quotient coordinates from direct branch observables

Using the exact branch-invariant definitions

a) `T_* = R_tr^(-C_*)`,

b) `N_* = Lambda_0 (1 - epsilon_eta) R_tr^(B_*) / R_target`,

with

a) `C_* = ((1+chi_{0,*})(1+delta_{U,*})(1+chi_{0,*}+delta_{U,*}))/(chi_{0,*} delta_{U,*})`,

b) `B_* = 2(1+chi_{0,*}+delta_{U,*})/delta_{U,*}`,

the finite quotient coordinates relative to a reference branch become exactly

`q_tr = - C_* ln(R_tr / R_tr,ref)`,

`q_nt = B_* ln(R_tr / R_tr,ref)
        + ln[(1-epsilon_eta)/(1-epsilon_eta,ref)]
        - ln(R_target / R_target,ref)`,

`q_eta = ln(epsilon_eta / epsilon_eta,ref)`.

So the Stage 300 finite quotient packet has an exact direct branch-observable form.
The inverse map is also exact:

`R_tr = R_tr,ref exp(-q_tr / C_*)`,

`epsilon_eta = epsilon_eta,ref exp(q_eta)`,

`R_target = R_target,ref
            exp[-q_nt - (B_*/C_*) q_tr]
            * (1-epsilon_eta)/(1-epsilon_eta,ref)`.

This proves that the coherent branch can be charted equivalently by either

- `(R_tr, R_target, epsilon_eta)`, or
- `(q_tr, q_nt, q_eta)`.

## Stage 304 — exact first-order defect compiler in direct branch language

Linearizing the Stage 303 finite quotient map gives

`q_tr = - C_* delta ln R_tr`,

`q_nt = B_* delta ln R_tr
        - [epsilon_eta,* / (1-epsilon_eta,*)] delta ln epsilon_eta
        - delta ln R_target`,

`q_eta = delta ln epsilon_eta`.

Composing this with the exact Stage 302 quotient-to-defect compiler yields a completely triangular first-order map:

`Theta_1 = delta ln R_tr`,

`Xi_1 = - delta ln R_target
        - [epsilon_eta,* / (1-epsilon_eta,*)] delta ln epsilon_eta`,

`R_1 = delta ln R_target`.

The inverse map is equally simple:

`delta ln R_tr = Theta_1`,

`delta ln R_target = R_1`,

`delta ln epsilon_eta = -[(1-epsilon_eta,*)/epsilon_eta,*] (R_1 + Xi_1)`.

So the physical first-order defect packet on the actual coherent branch is now exactly equivalent to the three direct observable drifts

`(delta ln R_tr, delta ln R_target, delta ln epsilon_eta)`.

## Stage 305 — exact three-observable and support-blind theorem

The coherent support-compensation side from Stages 30–35 only rescales the total baseline through

`M_tr = M_mix S(zeta; epsilon)`,

`S(zeta; epsilon) = 1 + zeta(1-epsilon)/(1-zeta epsilon)`.

But Stage 303 shows that the finite quotient packet depends only on

- `R_tr`,
- `R_target`,
- `epsilon_eta`,

and Stage 304 shows the same is true for the first-order defect packet.
Therefore both packets are exactly blind to the total baseline and to the coherent support-enhancement ratio:

`partial_(M_tr) (q_tr, q_nt, q_eta) = 0`,

`partial_(M_tr) (Theta_1, Xi_1, R_1) = 0`,

`partial_(zeta) (q_tr, q_nt, q_eta) = 0`,

`partial_(zeta) (Theta_1, Xi_1, R_1) = 0`.

So the support-compensation branch and the weak-axisymmetric orbit-lock branch are exact but distinct pieces of the 5PN endgame.

### Final structural consequence of Stages 303–305

The actual coherent weak-axisymmetric zero-defect theorem now has its sharpest direct form so far:

`Theta_1 = Xi_1 = R_1 = 0`

if and only if

`delta ln R_tr = 0`,

`delta ln R_target = 0`,

`delta ln epsilon_eta = 0`.

So the actual moving-throat branch no longer needs to be tested against the older eight-variable microscopic ledger, or even explicitly against the five continuum ratios, to answer the first-order orbit-lock question.
It only has to answer three exact branch-observable questions:

1. is `R_tr` invariant?
2. is `R_target` invariant?
3. is `epsilon_eta` invariant?

And importantly, coherent support enhancement does not enter that test.
It matters for whether the selected quadrupole branch can hit the normalization target, but not for whether the branch stays on the exact first-order similarity orbit.
