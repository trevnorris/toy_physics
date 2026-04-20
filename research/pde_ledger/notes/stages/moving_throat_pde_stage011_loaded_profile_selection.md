# Moving-Throat PDE — Stage 11: Loaded Axial Profile Selection, Exact 2x2 Eigenproblem, and Why the Blind-Angle Branch Is Dynamically Disfavored

## Purpose

Stage 10 replaced the constant axial wall/brane branch by the first nonconstant N/N family and showed that the whole normalization problem depends on one profile angle `theta` through

`kappa(theta) = kappa_0 cos(theta) + kappa_1 sin(theta)`.

That clarified the overlap algebra, but it still left one unsatisfactory feature:

> the profile angle `theta` was being chosen by hand rather than produced by a reduced moving-throat operator.

The next honest step is therefore to make `theta` an **output** of the first loaded axial eigenproblem.

The main result of this stage is that the minimal loaded wall operator already selects the profile direction in a very rigid way.

Using the first two N/N wall modes,

- the bare wall geometry is diagonal,
- the support/mixed branch enters as a rank-1 attractive load proportional to the overlap vector `(kappa_0, kappa_1)`,
- the profile angle satisfies one exact formula,
- and for positive loading it is driven toward the **max-coupling angle**, not toward the Stage-10 blind-angle no-go branch.

So after this stage the theorem gap narrows again:

- the blind-angle branch is still an exact algebraic no-go,
- but in the minimal positive-loading model it is also dynamically disfavored,
- and the real PDE now needs to tell us only how large the effective loading is and whether it stays in the same sign/class as this minimal model.

---

## 1. Minimal loaded wall basis

Work in the first two exact N/N wall modes,

`u_0(s) = 1 / sqrt(L)`,

`u_1(s) = sqrt(2/L) cos(pi s / L)`.

Write the wall quadrupole profile as

`chi(s) = q_0 u_0(s) + q_1 u_1(s)`,

with normalized coefficient vector

`q = (q_0, q_1)^T`,

`q_0^2 + q_1^2 = 1`.

The natural angle parameterization is

`q_0 = cos(theta)`,

`q_1 = sin(theta)`.

The bare wall operator carried from the earlier stages is

`G_eta = -T_w d_s^2 + K_eta + 6 T_Omega`.

So its matrix in the `{u_0,u_1}` basis is exactly diagonal,

`K_bare = [[K_0, 0], [0, K_1]]`,

with

`K_0 = K_eta + 6 T_Omega`,

`K_1 = K_eta + 6 T_Omega + Delta K_ax`,

`Delta K_ax = T_w pi^2 / L^2`.

The D/N half-wave coupling vector from Stage 10 is

`v = (kappa_0, kappa_1)^T`,

with

`kappa_0 = 2 sqrt(2) / pi`,

`kappa_1 = -4 / (3 pi)`.

So the Stage-10 profile overlap is simply

`kappa(theta) = v . q`.

---

## 2. Minimal effective loading from the support/mixed branch

At low frequency, the support/mixed branch lowers the wall energy in proportion to the square of the overlap with the D/N half-wave. The minimal reduced model is therefore

`E_eff[q] = (1/2) q^T K_bare q - (alpha/2) (v . q)^2`,

with one positive loading susceptibility

`alpha > 0`.

In matrix form,

`K_eff = K_bare - alpha v v^T`.

So explicitly,

`K_eff = [[K_0 - alpha kappa_0^2,      -alpha kappa_0 kappa_1],`
`         [    -alpha kappa_0 kappa_1,  K_1 - alpha kappa_1^2 ]].`

This is the first actual loaded moving-throat axial operator in the program.

It is still reduced and minimal, but it already converts the free profile angle `theta` into a dynamical eigenvector problem.

---

## 3. Exact eigenvalues and exact profile-angle equation

The effective trace and determinant are

`tr(K_eff) = K_0 + K_1 - alpha (kappa_0^2 + kappa_1^2)`,

`det(K_eff) = K_0 K_1 - alpha ( K_1 kappa_0^2 + K_0 kappa_1^2 )`.

So the exact eigenvalues are

`lambda_(+-)`
` = (1/2) [ K_0 + K_1 - alpha (kappa_0^2 + kappa_1^2)`
`          +- sqrt( (Delta K_ax + alpha (kappa_0^2 - kappa_1^2))^2`
`                   + 4 alpha^2 kappa_0^2 kappa_1^2 ) ]`.

The lower branch `lambda_-` is the physically relevant loaded wall stiffness.

Parameterizing the eigenvector by `theta`, the exact stationarity condition is

`tan(2 theta)`
` = 2 alpha kappa_0 kappa_1 / ( Delta K_ax + alpha (kappa_0^2 - kappa_1^2) )`.

Since

`kappa_0^2 - kappa_1^2 = 56 / (9 pi^2) > 0`,

and

`2 kappa_0 kappa_1 = -16 sqrt(2) / (3 pi^2) < 0`,

the sign structure is rigid.

For any positive loading `alpha > 0`, the denominator is positive and the numerator is negative, so

`theta < 0`.

This already proves an important selection statement:

> the first loaded wall eigenprofile is driven away from the flat Stage-9 branch in the **negative** `u_1` direction.

That is exactly the direction of the Stage-10 max-coupling angle, not the positive blind-angle branch.

---

## 4. Small-loading and strong-loading limits

### 4.1 Weak loading

For small loading,

`theta = alpha kappa_0 kappa_1 / Delta K_ax + O(alpha^2)`.

Because `kappa_0 kappa_1 < 0`, the first correction to the constant branch is a small negative angle.

So the Stage-9 flat branch is perturbatively stable, but it is not the exact loaded eigenprofile once the support/mixed load is turned on.

### 4.2 Strong loading

For large positive loading,

`tan(2 theta) -> 2 kappa_0 kappa_1 / (kappa_0^2 - kappa_1^2)`.

That is precisely the two-angle relation satisfied by the normalized coupling vector `v / |v|`.

So in the strong-loading limit the selected profile approaches

`q -> v / |v|`,

which means

`theta -> theta_max`,

with

`tan(theta_max) = kappa_1 / kappa_0 = - sqrt(2) / 3`.

This is a major simplification.

The Stage-10 max-coupling angle is **not** an arbitrary hand-picked member of the family. It is the exact strong-loading eigenvector of the first reduced loaded wall operator.

So the minimal loaded model already explains why the physically preferred profile should move toward the max-coupling branch rather than stay flat forever.

---

## 5. Why the blind-angle no-go branch is dynamically disfavored

Stage 10 showed that the blind angle

`tan(theta_blind) = 3 sqrt(2) / 2`

forces

`kappa(theta_blind)=0`,

and therefore kills the outgoing quadrupole normalization bridge.

Stage 11 sharpens that conclusion.

Because the loaded wall angle satisfies

`tan(2 theta)`
` = 2 alpha kappa_0 kappa_1 / ( Delta K_ax + alpha (kappa_0^2 - kappa_1^2) )`,

and the right-hand side is always **negative** for positive loading,

the selected `theta` always lies in the negative-coupling direction.

But the blind angle is positive,

`theta_blind > 0`.

So in the minimal attractive-loading model,

> the blind-angle branch is not only an algebraic no-go for normalization; it is dynamically disfavored by the profile-selection eigenproblem itself.

That is the strongest theorem-like statement we have so far about profile selection.

---

## 6. Exact softening threshold

The lower eigenvalue crosses zero when

`det(K_eff)=0`.

So the exact loading threshold is

`alpha_crit = K_0 K_1 / ( K_1 kappa_0^2 + K_0 kappa_1^2 )`.

This is the first concrete softening threshold in the moving-throat PDE program.

Interpretation:

- for `alpha < alpha_crit`, the loaded quadrupole wall mode remains stable,
- for `alpha = alpha_crit`, the branch goes marginal,
- for `alpha > alpha_crit`, the minimal reduced wall operator predicts a quadrupole instability/condensation.

So the future PDE does not just need to say which angle is chosen.
It must also place the physical branch on the stable or near-softened side of this exact threshold.

---

## 7. What this means for the normalization bridge

The outgoing-normalization equation from Stage 10 depends on

- the geometry-side wall stiffness,
- and the profile-dressed overlap `kappa(theta)`.

Stage 11 shows that these are not independent in the minimal loaded model.
They are linked because the same support/mixed loading that wants a larger `|kappa|` also rotates the wall toward `theta_max` and softens the lower wall eigenvalue.

So the next theorem question is no longer

> can I pick a `theta` that makes the algebra work?

The sharper question is

> does the actual loaded wall eigenmode selected by the throat operator carry the right pair `(lambda_-, kappa(theta))` to satisfy the outgoing quadrupole normalization target while remaining on the stable branch?

That is a much more physical and much more falsifiable question.

---

## 8. Best current summary after Stage 11

The first reduced loaded wall eigenproblem already resolves a major ambiguity.

It shows that:

- the Stage-10 profile angle is not a free nuisance parameter,
- the support/mixed load selects a definite direction in the N/N wall basis,
- the selected profile rotates away from the constant branch toward the max-coupling branch,
- the blind-angle no-go branch is dynamically disfavored rather than merely algebraically bad,
- and there is an exact softening threshold `alpha_crit` for the loaded wall quadrupole mode.

So the next honest derivation step is now very sharply defined:

> derive the actual effective loading strength `alpha` and the associated loaded wall eigenmode from the coupled wall/BdG/Maxwell/mixed operator, then insert that eigenpair into the Stage-10 normalization equation.

At this point the roadmap is no longer “invent a PDE.”
It is starting to look like a specific spectral problem.
