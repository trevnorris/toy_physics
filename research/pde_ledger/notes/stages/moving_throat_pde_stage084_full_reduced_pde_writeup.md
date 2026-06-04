# Moving-Throat PDE — Stage 084: Full Reduced Moving-Throat PDE Write-Up Skeleton and the Final Remaining Theorem Gap

## Purpose

Stages 001–083 have now done two different kinds of work:

1. they built the reduced moving-throat PDE hierarchy layer by layer, and
2. they compressed the remaining support/source branch question to one explicit quadrupole residual.

So this stage is not another isolated calculation. It is the first honest **write-up skeleton** of the full reduced moving-throat PDE program as it currently stands.

The point is to say clearly what we now have, what is already exact or SymPy-backed, and what the actual remaining theorem gap is.

---

## 1. The reduced moving-throat PDE hierarchy now in hand

### 1.1 Exact parent bulk system

The parent theory already carries:

- the gauged 4D GNLS matter equation,
- the localized `4+1` Maxwell equation,
- the geometry variables `(a,L)` or their distributed wall lift,
- and the open-system projection/leakage machinery.

That is the exact parent backbone inherited from the action-based 4D program.

### 1.2 Linearized support/source throat operator

The explicit reduced throat-support sector is now an actual PDE/ODE operator, not just a symbolic closure:

`partial_t sigma + partial_s J = 0,`

`J = -D_sigma partial_s sigma + v_sigma sigma,`

`-T_X partial_s^2 phi + K_X phi = Lambda_phi sigma,`

with Robin/Neumann boundaries and the operator-selected fixed-point law

`Pe_* = Xi Delta(Pe_*;kappa,eta).`

This is the exact reduced support/source operator carried from Stage 058 and
reassembled into the master residual form at Stage 082.

### 1.3 Exact explicit lowest-lane support family

The explicit lowest support lane from Stages 056-057, specialized in Stage 079 and
used operationally in Stages 082-083, is now fully reduced to physical variables

`(Pe, eta, kappa)`

through

`zeta_phys(Pe,eta;kappa)`
`= Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2),`

`Omega_Pe`
`= pi Pe (2 Pe exp(Pe) + pi)`
`  / [ (4 Pe^2 + pi^2) (exp(Pe)-1) ],`

`Omega_0 = 1,`

`y tan y = eta.`

This family is no longer heuristic. Its monotonicity, closure window, overlap ceiling, Robin softening law, and exact operator-selected bounds are all SymPy-backed.

### 1.4 Exact selected-branch quadrupole demand map

On the outgoing quadrupole side, the selected-branch demand from Stages 052 and 082
is now compressed to

`zeta_req(Pi_tr,C_mix,eps_blk)`
`= (Pi_tr - C_mix) / [ C_mix - eps_blk (2 C_mix - Pi_tr) ]`

with exact inversion

`Pi_tr = C_mix Q(zeta_req;eps_blk),`

`Q(zeta;eps_blk) = [1 + (1 - 2 eps_blk) zeta]/[1 - eps_blk zeta].`

So the outgoing branch has also been reduced to one exact scalar demand function.

---

## 2. The full reduced theorem gate

Putting those two sides together, the reduced moving-throat PDE is now governed by one scalar residual:

`R_quad(Pi_tr,C_mix,eps_blk ; Xi,eta,kappa)`
`= zeta_req(Pi_tr,C_mix,eps_blk)`
` - zeta_phys(Pe_*(Xi,eta,kappa),eta;kappa).`

This is the current write-up version of the full reduced theorem gate from
Stage 082.

- `R_quad < 0`  ->  explicit support/source supply exceeds selected quadrupole demand;
- `R_quad = 0`  ->  exact saturation of the explicit branch;
- `R_quad > 0`  ->  the explicit branch fails.

And because the support/source fixed-point already obeys the exact Stage-058
brackets,

`Xi Delta_0 <= Pe_* <= Xi Delta_inf,`

the reduced PDE also has exact success/failure bounds before any full root solve.

---

## 3. The explicit Family-1 specialization is now closed

The first real branch is now no longer abstract.
On Family-1 we have

`kappa_F1 = 12321/5,`

`eta_F1 = 37,`

`Xi_F1 = 1369 Upsilon_w = 136900 Theta_w.`

Using the exact Stage-077 wall-depth extraction from the frozen `n=5` branch, the
Stage-073-075 Family-1 data, and the Stage-083 direct operator evaluation, the
natural shell-weighted and conservative-floor branches produce direct
operator-selected windows in:

- transport bias `Pe`,
- support ratio `zeta`,
- and selected branch product `Pi_tr/C_mix`.

Numerically, the natural shell-weighted branch at `lambda_mu = 1` already pushes the explicit support window to

`zeta_-^(chi) ≈ 2.46622291347846,`

`zeta_+^(chi) ≈ 2.46752913273870,`

while the conservative floor gives

`zeta_-^(J) ≈ 2.44257571477179,`

`zeta_+^(J) ≈ 2.46752736855058,`

against the hard Stage-079 ceiling

`zeta_max^(F1) ≈ 2.46752922945601.`

So the explicit Family-1 support/source branch is effectively saturated.

That is why the support/source side is no longer the serious unresolved issue.

---

## 4. What is still not solved

Even after all this reduction, the project still does **not** yet have a first-principles theorem of the full moving-throat PDE.

What remains open is narrower:

1. the actual passive/outgoing quadrupole branch still has to produce its physical selected product `Pi_tr`,
2. the physical mixed baseline `C_mix` and blocking ratio `eps_blk` still have to be fixed on the true moving-throat branch,
3. and that final branch point must then be inserted into `R_quad`.

So the unresolved part is no longer the support/source operator, no longer the wall-depth supply, and no longer the explicit Family-1 reduction.

It is the **outgoing quadrupole-normalization branch itself**.

---

## 5. Honest current status

So we are indeed very close to the point where a full moving-throat PDE write-up becomes natural — but in the precise sense below.

### What is now ready to write up

A full **reduced moving-throat PDE write-up** is now ready:

- exact parent system,
- exact reduced support/source operator,
- exact explicit support family,
- exact selected-branch demand map,
- exact master residual `R_quad`,
- exact Family-1 specialization,
- and exact success/failure bounds.

### What is not yet honest to write up as finished

A full **first-principles theorem of the complete moving-throat quadrupole branch** is not yet ready, because the final outgoing-normalization branch is still the one unresolved ingredient.

So the correct statement is:

> the reduced moving-throat PDE program is now fully organized and almost fully written; the one remaining theorem gap is the actual passive/outgoing quadrupole-normalization branch that fixes `Pi_tr` on the true moving-throat solution.

---

## 6. The clean next move from here

The next derivation is now as narrow as it can reasonably get:

> derive the physical outgoing quadrupole product `Pi_tr` (and its accompanying `C_mix`, `eps_blk`) from the actual passive/outgoing moving-throat branch and evaluate the sign of `R_quad`.

That is the remaining finish line for the present reduced moving-throat PDE program.
