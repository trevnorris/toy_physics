# Moving-Throat PDE — Stage 052: Exact Non-Twin Asymmetry Requirement Once the Symmetric Lowest Twin Fails

## Purpose

Stage 051 compressed the lowest-twin question to one exact product test:

`Pi_tr(xi_req,delta;R_tr) <= 16 Lambda (1-eps) / pi^2.`

That tells us **when** the symmetric lowest twin lane is enough.
But the next physical question is just as important:

> if the symmetric lowest twin fails, what exact non-twin deformation is required to beat the support threshold?

This stage answers that.

The main result is that the support problem now splits into three exact regimes controlled by the same branch product `Pi_tr` and mixed product scale

`C_mix := 8 Lambda (1-eps) / pi^2.`

They are:

1. **mixed-only enough**

   `Pi_tr <= C_mix;`
2. **symmetric lowest twin enough**

   `C_mix < Pi_tr <= 2 C_mix;`
3. **non-twin asymmetry required**

   `Pi_tr > 2 C_mix.`

And once the third regime is entered, the exact required lowest-lane support ratio is

`zeta_req = (Pi_tr - C_mix) / [ C_mix - eps (2 C_mix - Pi_tr) ].`

So the residual theorem gap is no longer just “do we need asymmetry?”
It is the much sharper question:

> **what exact overlap boost and/or support softening is required once the symmetric twin doubling branch is no longer enough?**

---

## 1. Exact regime classification in branch-product form

Define again

`C_mix := 8 Lambda (1-eps) / pi^2,`

`Pi_tr := F_tr G_tr.`

Then the exact support-demand factor is

`S_req = Pi_tr / C_mix,`

and the exact required coherent support ratio is

`zeta_req = (S_req - 1) / [ 1 + eps (S_req - 2) ]`
`         = (Pi_tr - C_mix) / [ C_mix - eps (2 C_mix - Pi_tr) ].`

This formula is relevant once support is actually needed, i.e. once `Pi_tr > C_mix`.

The exact derivative is

`dzeta_req / dPi_tr`
`= C_mix (1-eps) / [ C_mix - eps (2 C_mix - Pi_tr) ]^2 > 0,`

so `zeta_req` grows strictly with the branch product.

The key anchor values are

`Pi_tr = C_mix   => zeta_req = 0,`

`Pi_tr = 2 C_mix => zeta_req = 1.`

So the support regimes are exact:

- `Pi_tr <= C_mix`: mixed-only branch already enough;
- `C_mix < Pi_tr <= 2 C_mix`: support is needed, but the symmetric lowest twin lane still suffices;
- `Pi_tr > 2 C_mix`: the symmetric lowest twin fails, and non-twin asymmetry is required.

This is the cleanest classification yet of the support problem.

---

## 2. Exact excess beyond the symmetric twin

When `Pi_tr > 2 C_mix`, the extra support demand beyond the symmetric twin is

`Delta_zeta := zeta_req - 1.`

Substituting the exact formula above gives

`Delta_zeta`
`= (1-eps) (Pi_tr - 2 C_mix)`
`  / [ C_mix - eps (2 C_mix - Pi_tr) ].`

So `Delta_zeta` is strictly positive exactly when the symmetric twin fails.

This means the “distance beyond the twin branch” is not heuristic.
It is an explicit rational function of the same branch product and the same blocking parameter `eps`.

---

## 3. General lowest-lane support ratio and the exact asymmetry threshold

The lowest support lane need not remain perfectly twin-symmetric.
For a general lowest support lane, define the relative overlap boost

`Omega_0 := I_(phi,0) / I_W,`

and keep the effective stiffness ratio explicit. Then the exact physical lowest-lane support ratio is

`zeta_0^(phys) = (K_W^(eff) / K_(phi,0)^(eff)) * Omega_0^2.`

So the exact success condition is

`(K_W^(eff) / K_(phi,0)^(eff)) * Omega_0^2 >= zeta_req.`

This yields two equivalent threshold forms.

### Overlap-asymmetry threshold at fixed stiffness

If the stiffness ratio is held fixed, the exact required overlap boost is

`Omega_0^2 >= Omega_(0,req)^2`
`:= zeta_req * K_(phi,0)^(eff) / K_W^(eff).`

### Support-softening threshold at fixed overlap

If the overlap ratio is held fixed, the exact required lowest-lane softening is

`K_(phi,0)^(eff) <= K_(phi,0)^(req)`
`:= K_W^(eff) * Omega_0^2 / zeta_req.`

So once `zeta_req > 1`, the lowest support lane must become non-twin in one of two exact ways:

- enhanced coupling to the lowest support mode (`Omega_0 > 1`),
- softer effective lowest-lane support stiffness (`K_(phi,0)^(eff) < K_W^(eff)`),
- or a combination of both.

---

## 4. Exact twin-failure diagnostics

The perfectly symmetric same-operator twin lane has

`Omega_0 = 1,`

`K_(phi,0)^(eff) = K_W^(eff),`

so

`zeta_0^(twin) = 1.`

If `zeta_req > 1`, this branch fails identically.

Two immediate exact diagnostics follow.

### Pure-overlap rescue at equal stiffness

If `K_(phi,0)^(eff) = K_W^(eff)`, then the exact overlap threshold is

`Omega_0 >= sqrt(zeta_req).`

So any `zeta_req > 1` forces a true overlap asymmetry above the symmetric twin value.

### Pure-softening rescue at equal overlap

If `Omega_0 = 1`, then the exact stiffness threshold is

`K_(phi,0)^(eff) <= K_W^(eff) / zeta_req.`

So any `zeta_req > 1` forces the lowest support lane to be physically softer than the mixed lane.

The exact softening fraction in that equal-overlap case is

`1 - K_(phi,0)^(req) / K_W^(eff)`
`= 1 - 1 / zeta_req`
`= (1-eps) (Pi_tr - 2 C_mix) / (Pi_tr - C_mix).`

So the size of the required non-twin softening is also explicit.

---

## 5. Best current theorem statement after Stage 052

The support bottleneck is no longer vague and no longer merely qualitative.

### What is exact now

- the required coherent support ratio is

  `zeta_req = (Pi_tr - C_mix) / [ C_mix - eps (2 C_mix - Pi_tr) ],`
- `zeta_req` is strictly increasing in the branch product `Pi_tr`,
- the exact regime split is

  `Pi_tr <= C_mix`  : mixed-only enough,

  `C_mix < Pi_tr <= 2 C_mix` : symmetric lowest twin enough,

  `Pi_tr > 2 C_mix` : non-twin asymmetry required,
- the exact excess beyond the symmetric twin is `Delta_zeta = zeta_req - 1`,
- and the exact lowest-lane rescue thresholds are

  `Omega_0^2 >= zeta_req K_(phi,0)^(eff) / K_W^(eff),`

  `K_(phi,0)^(eff) <= K_W^(eff) Omega_0^2 / zeta_req.`

### What this means physically

The first explicit coherent kernel has now answered the question posed at the end of Stage 050.

If the physical branch lands in the window `Pi_tr <= 2 C_mix`, then the symmetric lowest twin lane is enough.

If it lands above that, then the missing theorem gap is no longer “find some better support somehow.”
It is:

> **derive from the completed moving-throat operator whether the lowest support lane acquires the exact overlap boost and/or stiffness softening required by the non-twin threshold formulas above.**
