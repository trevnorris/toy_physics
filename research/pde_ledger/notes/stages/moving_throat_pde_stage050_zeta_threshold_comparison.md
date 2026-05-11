# Moving-Throat PDE — Stage 050: Exact Comparison of the Physical Coherent `zeta` Against `zeta_req`

## Purpose

Stage 31 proved that the support-compensation problem is exactly

`zeta_phys >= zeta_req`,

while Stage 32 derived the physical coherent D/N support ratio on the first explicit finite-throat kernel,

`zeta_n^(phys)`
`= (K_W^(eff) / K_(phi,n)^(eff)) / (2n+1)^2,`

and, on the same-operator twin family,

`zeta_n^(twin)`
`= 1 / [ (2n+1)^2 ( 1 + x n(n+1) ) ].`

So the comparison that was still abstract at Stage 31 can now be carried out exactly.

The main results are:

1. the symmetric lowest twin lane `n=0` has
   
   `zeta_0^(twin)=1`,
   
   so its support-enhancement factor is exactly
   
   `S_0 = 2`,
   
   independent of `eps`;
2. therefore the lowest symmetric twin branch succeeds **iff**
   
   `S_req <= 2`,
   
   equivalently
   
   `zeta_req <= 1`;
3. every higher D/N support harmonic is bounded by
   
   `zeta_n^(twin) < 1/(2n+1)^2`,
   
   so if
   
   `zeta_req > 1/(2n+1)^2`,
   
   that harmonic is ruled out immediately;
4. when it is not ruled out immediately, the exact stiffness threshold is
   
   `x <= x_max(n; zeta_req)`
   `   = [ 1 / ( (2n+1)^2 zeta_req ) - 1 ] / [ n(n+1) ].`

So the exact support-threshold comparison is now sharp enough to distinguish the physically viable twin lane from the strongly suppressed higher-harmonic branches.

---

## 1. Exact support-threshold inequality in Stage-048 variables

Stage 31 defines the required enhancement factor

`S_req := M_req / M_mix`

and the exact required support ratio

`zeta_req = (S_req - 1) / [ 1 + eps (S_req - 2) ],`

with `0 < eps < 1` and `S_req > 1` whenever support is really needed.

Because the coherent support-enhancement factor is

`S(zeta;eps) = 1 + zeta (1-eps) / (1 - eps zeta),`

the support branch succeeds exactly when

`zeta_phys >= zeta_req.`

This is the comparison formula every explicit kernel has to meet.

---

## 2. Lowest symmetric twin lane: exact doubling theorem

For the same-operator twin branch, Stage 32 gave

`zeta_0^(twin) = 1.`

Substituting into the exact support-enhancement factor gives

`S_0 = S(1;eps)`
`    = 1 + (1-eps)/(1-eps)`
`    = 2.`

So the symmetric lowest twin lane doubles the mixed baseline exactly, independent of `eps`.

This immediately yields the exact comparison theorem.

### Exact equivalence

`zeta_req <= 1`
`<=> S_req <= 2.`

So the lowest symmetric twin lane succeeds exactly when the Stage-048 required enhancement does not exceed a factor of two.

This is the cleanest physical interpretation yet of the coherent support problem.

The first explicit finite-throat twin lane is not just “supportive.” It is a universal **doubling** branch.

---

## 3. Higher D/N support harmonics: exact impossibility bound

For `n>=1`, Stage 32 gives

`zeta_n^(twin)`
`= 1 / [ (2n+1)^2 ( 1 + x n(n+1) ) ],`

with `x>0`.

Therefore

`zeta_n^(twin) < 1 / (2n+1)^2.`

So an exact necessary condition for the `n`th support harmonic to succeed is

`zeta_req <= 1 / (2n+1)^2.`

If this fails, the `n`th D/N support harmonic is ruled out before any finer stiffness modeling is done.

This gives the first explicit microscopic no-go tower:

- `n=1` impossible if `zeta_req > 1/9`,
- `n=2` impossible if `zeta_req > 1/25`,
- `n=3` impossible if `zeta_req > 1/49`,
- and so on.

So the physical support problem is already strongly biased toward the lowest coherent twin lane.

---

## 4. Exact stiffness threshold when a higher harmonic is not yet ruled out

When

`zeta_req <= 1 / (2n+1)^2`,

one can solve the exact inequality

`zeta_n^(twin) >= zeta_req`

for the stiffness parameter `x`.

For `n>=1`, the result is

`x <= x_max(n; zeta_req)`
`   = [ 1 / ( (2n+1)^2 zeta_req ) - 1 ] / [ n(n+1) ].`

So higher support harmonics are viable only when the corresponding support lane is sufficiently soft.

This is the first exact microscopic support-feasibility threshold written directly in terms of the Stage-048 support demand.

---

## 5. Exact enhancement bounds for the higher harmonics

The physical enhancement produced by the `n`th twin harmonic is

`S_n^(twin)`
`= 1 + zeta_n^(twin) (1-eps) / (1 - eps zeta_n^(twin)).`

Because `S` is strictly increasing in `zeta`, the upper bound on `zeta_n^(twin)` gives the exact enhancement ceiling

`S_n^(twin) < S_n^(max)`
`:= 1 + (1-eps) / [ (2n+1)^2 - eps ].`

So, for example,

`S_1^(twin) < 1 + (1-eps)/(9-eps) = (10 - 2 eps)/(9 - eps),`

which is only a modest enhancement above `1`.

By contrast,

`S_0^(twin) = 2`

exactly.

So the lowest twin lane is not just the strongest support branch. It is qualitatively different from the higher harmonics in the size of the enhancement it can deliver.

---

## 6. Best current theorem statement after Stage 33

The support-compensation problem is now resolved down to an explicit harmonic-selection test on the first finite-throat coherent kernel.

### What is now exact

- the lowest symmetric twin lane has `zeta=1` and therefore `S=2` exactly,
- `zeta_req <= 1` is equivalent to `S_req <= 2`,
- every higher harmonic obeys the exact impossibility bound `zeta_req <= (2n+1)^(-2)`,
- when that bound is satisfied, the exact softness threshold is `x <= x_max(n;zeta_req)`,
- and the higher-harmonic support enhancement is strictly bounded above by `S_n^(max)`.

### What this means physically

The first explicit moving-throat kernel does more than produce a physical `zeta`.
It also orders the entire support tower.

- If the exact Stage-048 demand satisfies `zeta_req <= 1`, the symmetric lowest twin D/N lane is already sufficient.
- If `zeta_req > 1`, the lowest symmetric twin lane fails and stronger-than-twin asymmetry or a different kernel family is required.
- If `zeta_req > 1/9`, the first excited D/N support harmonic is already impossible.
- And the higher harmonics become rapidly less viable.

So the next honest PDE question is now extremely narrow:

> **does the completed moving-throat operator place the physical coherent support lane on the lowest twin branch with `zeta_req <= 1`, or does the real branch need non-twin asymmetry to overcome the exact support threshold?**
