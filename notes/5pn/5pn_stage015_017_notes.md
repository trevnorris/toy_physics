# 5PN continuation notes — stages 15 through 17

These stages finally splice the Stage-5 primitive obstruction triplet into the later Stage-10/11 monomial-orbit package.

The main point is that the three objects are **not** on the same footing:

- `Xi_load = N01/N0 - D01/D0 = P1/P0` is exactly the weak-axisymmetric normalization defect governed by the monomial/similarity theorem;
- `K1 = D21 + D01/9` and
- `H_even = D41 - (2/3) D21 - D01/27`

are separate conservative even gates that survive after the normalization defect is killed.

## Stage 15 — the exact `Xi_load` bridge

The exact Stage-5 prefactor identity is
\[
P_0 = \frac{N_0}{D_0},
\qquad
P_1 = \frac{N_{01}D_0 - N_0 D_{01}}{D_0^2},
\qquad
\Xi_{\rm load}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}
=
\frac{P_1}{P_0}.
\]
So the Stage-5 load defect is **exactly** the same weak-axisymmetric prefactor slope already isolated later as
\[
\Xi_1 = \frac{P_1}{P_0}.
\]

Rewriting Stages 10–14 in one place, the direct monomial drifts are
\[
\Sigma_{\rm tr}=\delta\ln \mathfrak C_{{\rm tr},*},
\qquad
\Sigma_{\rm nt}=\delta\ln \mathfrak C_{{\rm nt},*},
\qquad
\Sigma_\eta=\delta\ln \epsilon_\eta,
\]
and the weak-axisymmetric normalization defect has the exact triangular form
\[
\Xi = A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt}.
\]
The Stage-13 normalized kernel annihilates all three monomial drifts, and the Stage-14 extension shows the explicit BdG directions `dln lambda_B` and `dln varpi` are zero columns of the monomial matrix. Therefore
\[
\Xi_{\rm load}=0
\]
on:

- the full normalized similarity kernel,
- its injected support-blind extension,
- the explicit BdG amplitude/frequency directions,
- the matched wall-only co-scaling direction.

So the monomial/similarity theorem really is a theorem about `Xi_load`.

## Stage 16 — why that still does **not** solve 5PN

The next check is to compare `Xi` with the conservative even gates.

On the matched wall-only co-scaling direction,
\[
d\ln K = d\ln M = 1,
\]
we still get
\[
K1_{\rm wall}=\frac{K}{9}-M,
\qquad
H_{{\rm even,wall}}=-\frac{K}{27}+\frac{2M}{3},
\]
which are generically nonzero even though `Xi = 0`.

On the explicit support-blind BdG directions,
\[
d\ln \lambda_B = 1,
\qquad d\ln \varpi = 0,
\]
or
\[
d\ln \lambda_B = 0,
\qquad d\ln \varpi = 1,
\]
we again have `Xi = 0`, but both `K1_B` and `H_even,B` are nonzero.

Even on the pure BdG self-similar branch
\[
d\ln \lambda_B = d\ln \varpi = \delta,
\]
which is the branch that kills the BdG load defect, one still gets
\[
K1_B = \frac{2B_0\delta}{\varpi^2},
\qquad
H_{{\rm even},B} = \frac{4B_0\delta(3-\varpi^2)}{3\varpi^4},
\]
so the even gates remain open unless `\delta = 0` (or an extra tuning is imposed).

That is the cleanest executable statement yet that:
\[
\text{similarity-orbit rigidity} \neq \text{full 5PN closure}.
\]
The orbit theorem kills `Xi`; it does **not** automatically close `K1` or `H_even`.

## Stage 17 — the first lower-bound intersection calculation

Once the Stage-13/14 monomial-rigid kernel is parameterized by
\[
(\alpha_K,\alpha_W,\alpha_U,\alpha_R,\alpha_{\Omega_U},\beta_B,\beta_{\varpi}),
\]
we can impose the **lower-bound** conservative even gates obtained by keeping only the explicit wall-only and pure-BdG pieces of `K1` and `H_even`.

This gives an exact `2 x 7` gate matrix of generic rank `2`, so the lower-bound even-gate intersection has dimension `5`.

A convenient exact solve is to determine
\[
\alpha_K,
\qquad
\beta_B,
\]
in terms of the remaining five free directions.

The corresponding null basis has five directions:

1. pure `\alpha_W`,
2. pure `\alpha_R`,
3. matched `\alpha_U/\alpha_{\Omega_U}`,
4. BdG-amplitude deformation `\beta_B` with compensating `\alpha_K,\alpha_U`,
5. BdG-frequency deformation `\beta_{\varpi}` with compensating `\alpha_K,\alpha_U`.

The important caution is that this is only a **lower-bound** solve, because the conservative mixed-sector `Z_2,Z_4` moments have not yet been reinstated. In particular, the fact that `\alpha_W` and `\alpha_R` survive untouched here is telling us exactly where the omitted mixed-sector moments still have to act.

## Net result after Stage 17

The 5PN continuation point is now much sharper:

1. `Xi_load = P1/P0` is fully absorbed into the similarity-orbit / monomial-rigidity theorem;
2. the real conservative 5PN bottleneck is the pair of even gates `K1` and `H_even`;
3. the explicit wall-only and BdG-only pieces show those even gates survive after `Xi` is killed;
4. and the next honest task is therefore:

> reinstate the conservative mixed-sector `Z_2,Z_4` moments and determine how they cut the remaining lower-bound tangent family.

That is the cleanest next theorem gate we have reached so far.
