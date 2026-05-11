# Moving-Throat PDE — Stage 048: Exact Support-Compensation Theorem on the Physical Tracking Branch

## Purpose

Stage 30 compressed the coherent local D/N kernel to one exact support-enhancement factor,

`M_tr = M_mix S(zeta;eps),`

with

`S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps).`

So the reduced problem is no longer “what support closure do we use?”
It is now a direct kernel-strength question:

> can the coherent support lane overcome the exact tracking-branch normalization deficit before the selected branch reaches its softening limit?

This stage answers that exactly.

The main result is that there is **no reduced-level support no-go** on the physical tracking branch.

For any finite target ratio `R_target > 1`, any stable geometry `(delta,R_tr)`, and any mixed-only baseline below the tracking critical load, there is a unique coherent support ratio `zeta_req` that hits the target before softening.

So the remaining theorem gap is not a reduced closure obstruction. It is whether the completed moving-throat PDE supplies a physical `zeta` large enough to reach that exact required value.

---

## 1. Universal tracking branch and its critical load

On the coherent branch the selected mode obeys the exact tracking law

`M_tr = G_tr(xi,delta;R_tr),`

with

`G_tr(xi,delta;R)`
`= 9 xi (xi + delta) / [ 9 delta + (9 + 2 R^2) xi ].`

The stable selected branch lives on

`0 < xi < 1.`

The exact derivative is

`dG_tr/dxi`
`= 9 [ 2 R^2 xi^2 + 9 delta^2 + 18 delta xi + 9 xi^2 ]`
`  / [ 2 R^2 xi + 9 delta + 9 xi ]^2 > 0.`

So `G_tr` is strictly increasing along the stable branch.

Its exact softening-limit value is therefore

`M_crit(delta,R)`
`:= G_tr(1,delta;R)`
`= 9 (1 + delta) / [ 9 delta + 9 + 2 R^2 ].`

Equivalently,

`M_crit - G_tr(xi,delta;R)`
`= 9 (1-xi) [ 2 R^2 xi + 9 delta^2 + 9 delta xi + 9 delta + 9 xi ]`
`  / [ (2 R^2 + 9 delta + 9) (2 R^2 xi + 9 delta + 9 xi) ] > 0`

for `0 < xi < 1`.

So every stable selected-branch point lies strictly below the critical load.

---

## 2. Exact support-enhancement factor and its inverse

From Stage 30,

`M_tr = M_mix S(zeta;eps),`

`S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps),`

with the stability window

`0 < eps < 1,`

`0 <= zeta < 1/eps.`

The exact derivative is

`dS/dzeta = (1-eps)/(1-zeta eps)^2 > 0,`

so `S` is strictly increasing.

The endpoint values are

`S(0;eps) = 1,`

`lim_(zeta -> 1/eps^-) S(zeta;eps) = +infinity.`

So the support lane can realize any finite enhancement factor above `1` while staying below the blocking pole.

The inverse map is exact. If `S_req > 1`, then the unique support ratio producing that enhancement is

`zeta_req = (S_req - 1) / [ 1 + eps (S_req - 2) ].`

It also obeys the exact stability margin

`1/eps - zeta_req`
`= (1-eps) / [ eps ( 1 + eps (S_req - 2) ) ] > 0.`

So every finite required enhancement lies strictly inside the stability window.

---

## 3. Existence of a stable-side target point

The exact tracking normalization function is

`F_tr(xi,delta;R)`
`= [ 9 delta + (9 + 2 R^2) xi ]^2 [ 9 delta + (9 + 2 R) xi ]^2`
`  / [ 81 (1 - xi) ( 9 delta^2 + 18 delta xi + (9 + 2 R^2) xi^2 )^2 ].`

The second bracket is intentionally linear in `R`, not `R^2`: on the coherent
tracking branch it comes from the source-side direction factor, while the first
bracket and the denominator carry the conservative `R^2` dependence.

Two endpoint facts are exact:

`F_tr(0,delta;R) = 1,`

`lim_(xi -> 1^-) F_tr(xi,delta;R) = +infinity.`

Therefore, by continuity, for every finite target ratio

`R_target > 1`

there exists at least one stable-side root

`xi_req in (0,1)`

such that

`F_tr(xi_req,delta;R_tr) = R_target.`

Define the corresponding required total load by

`M_req := G_tr(xi_req,delta;R_tr).`

Because `0 < xi_req < 1` and `G_tr` is strictly increasing,

`0 < M_req < M_crit(delta,R_tr).`

So the normalization target, if finite, is always associated with a finite stable-side required load.

---

## 4. Exact support-compensation theorem

Suppose the mixed-only coherent branch is stable,

`0 < M_mix < M_crit(delta,R_tr).`

There are then two cases.

### Case A — mixed-only branch already strong enough

If

`M_mix >= M_req,`

then the target is already reachable with

`zeta_req = 0.`

### Case B — support enhancement is needed

If

`M_mix < M_req,`

define the exact required enhancement factor

`S_req := M_req / M_mix > 1.`

Then the unique coherent support ratio that hits the target is

`zeta_req = (S_req - 1) / [ 1 + eps (S_req - 2) ].`

And because `M_req < M_crit`, one has

`S_req < S_crit := M_crit / M_mix,`

so the support ratio lies strictly below the critical support ratio

`zeta_crit = (S_crit - 1) / [ 1 + eps (S_crit - 2) ].`

Indeed,

`zeta_crit - zeta_req`
`= (S_crit - S_req)(1-eps)`
`  / [ (1 + eps(S_crit - 2))(1 + eps(S_req - 2)) ] > 0.`

So the target is reached **before** the selected branch softens out.

---

## 5. Exact monotone wall-softening response to support enhancement

Combining Stages 30 and 31 gives the exact implicit selected-branch law

`M_mix S(zeta;eps) = G_tr(xi_phys,delta;R_tr).`

Differentiating implicitly gives

`dxi_phys/dzeta`
`= M_mix (dS/dzeta) / (dG_tr/dxi)|_(xi_phys)`
`> 0.`

So coherent support enhancement always drives the selected branch to larger softening depth.

This is the exact reduced statement behind the compensation theorem:

- support does not alter `R_target`,
- it only increases the available baseline,
- and that increase moves the physical branch monotonically deeper into the tracking family.

---

## 6. Best current theorem statement after Stage 31

The coherent local D/N branch now has an exact reduced support-feasibility theorem.

### What is exact now

- the tracking branch has a finite critical load `M_crit(delta,R_tr)`,
- the coherent support-enhancement factor `S(zeta;eps)` is strictly increasing and invertible,
- every finite target ratio `R_target > 1` corresponds to at least one stable-side branch point `xi_req in (0,1)`,
- the corresponding required load `M_req` is always strictly below `M_crit`,
- and whenever `M_mix < M_req`, there is a unique coherent support ratio `zeta_req < zeta_crit < 1/eps` that reaches the target before softening.

### What this means physically

The exact tracking-branch deficit from Stages 28–29 is real, but it is **not** a reduced-level no-go.

The first coherent local D/N support lane can compensate that deficit exactly.
The remaining question is no longer whether compensation is possible in principle.
It is narrower:

> **does the actual moving-throat PDE produce a physical support ratio `zeta` large enough to reach the exact required value `zeta_req` on the real branch?**
