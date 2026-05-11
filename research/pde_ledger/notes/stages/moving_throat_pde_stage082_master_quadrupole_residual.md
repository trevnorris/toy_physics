# Moving-Throat PDE — Stage 082: Master Quadrupole Residual of the Full Reduced Moving-Throat PDE

## Purpose

The explicit Family-1 branch phase ended with a hard support ceiling in the selected-branch product variable `Pi_tr`, but that result was still being carried through several intermediate variables (`Pe_req`, `zeta_req`, `Pi_tr/C_mix`).

The next honest step toward a full moving-throat PDE write-up is to state the entire remaining quadrupole-normalization problem as **one reduced PDE residual** and then show how the earlier exact threshold theorems sit inside it.

That is what this stage does.

The main result is that the whole reduced moving-throat PDE, on the surviving passive/outgoing quadrupole branch, is now summarized by one scalar residual

`R_quad(Pi_tr,C_mix,eps_blk ; Xi,eta,kappa)`
`:= zeta_req(Pi_tr,C_mix,eps_blk) - zeta_phys(Pe_*(Xi,eta,kappa),eta;kappa),`

where

`zeta_req(Pi_tr,C_mix,eps_blk)`
`= (Pi_tr - C_mix) / [ C_mix - eps_blk (2 C_mix - Pi_tr) ]`

is the exact selected-branch demand from Stages 34–35,

and

`zeta_phys(Pe,eta;kappa)`
`= Omega_Pe^2 (kappa + pi^2/4) / (kappa + y(eta)^2),`

`y tan y = eta,`

is the exact explicit lowest-lane support ratio from Stages 39–40,
with `Pe_*` the operator-selected transport bias solving the Stage-41 fixed-point equation

`Pe_* = Xi Delta(Pe_*;kappa,eta).`

So the remaining reduced theorem question is no longer diffuse:

> **does the completed moving-throat quadrupole branch make `R_quad` nonpositive on the physical branch?**

---

## 1. Full reduced system entering the quadrupole residual

The reduced moving-throat PDE hierarchy now has five live layers.

### 1.1 Parent bulk matter/gauge sector

The exact parent theory still carries

`i hbar D_t psi`
`= [ -(hbar^2/2m) D_A D_A + V_conf(X;a,L) + h(|psi|^2) ] psi,`

and the localized `4+1` Maxwell sector

`partial_M ( Z(w) F^(MN) ) + (1/xi) partial^N(partial·A) = mu_0 J^N,`

with the geometry variables `(a,L)` or their distributed lift already embedded in the confinement sector.

### 1.2 Lowest support/source throat operator

On the reduced axial support/source branch, the same moving-throat program already isolated the coupled operator

`partial_t sigma + partial_s J = 0,`

`J = -D_sigma partial_s sigma + v_sigma sigma,`

`-T_X partial_s^2 phi + K_X phi = Lambda_phi sigma,`

with Robin/Neumann support boundaries

`T_X phi_s(0) = K_m phi(0),`

`phi_s(L) = 0.`

### 1.3 Stationary operator-selected transport bias

On the stationary zero-flux branch, the source density is

`Sigma_Pe(x) = Pe exp(Pe x)/(exp(Pe)-1),`

with `x=s/L in [0,1]`, and the support/source operator selects the transport bias through the exact fixed-point law

`Pe_* = Xi Delta(Pe_*;kappa,eta),`

where

`kappa = K_X L^2 / T_X,`

`eta = K_m L / T_X,`

`Xi = mu_sigma Lambda_phi^2 L^2 / (D_sigma T_X).`

### 1.4 Exact explicit lowest-lane support ratio

Once the support/source branch is selected, the explicit lowest support lane carries the exact ratio

`zeta_phys(Pe,eta;kappa)`
`= Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2),`

with the transport overlap factor

`Omega_Pe`
`= pi Pe (2 Pe exp(Pe) + pi)`
`  / [ (4 Pe^2 + pi^2) (exp(Pe)-1) ]`

and `y(eta)` the unique Robin root `y tan y = eta`, `0<y<pi/2`.

### 1.5 Exact selected-branch quadrupole demand

On the outgoing quadrupole side, the selected-branch support demand remains

`zeta_req(Pi_tr,C_mix,eps_blk)`
`= (Pi_tr - C_mix) / [ C_mix - eps_blk (2 C_mix - Pi_tr) ],`

which is just Stage 35 rewritten as the exact demand variable of the branch-product formulation.

So the reduced moving-throat PDE closes to one scalar comparison:

`zeta_req ? zeta_phys`.

---

## 2. The master quadrupole residual

Define the exact residual

`R_quad(Pi_tr,C_mix,eps_blk ; Xi,eta,kappa)`
`:= zeta_req(Pi_tr,C_mix,eps_blk)`
` - zeta_phys(Pe_*(Xi,eta,kappa),eta;kappa).`

Then the selected outgoing quadrupole branch satisfies:

- `R_quad < 0`  ->  support supply exceeds demand;
- `R_quad = 0`  ->  exact saturation of the surviving quadrupole branch;
- `R_quad > 0`  ->  the explicit lowest-lane support branch fails.

This is the cleanest full reduced-PDE statement reached so far.

It absorbs all of the earlier stage variables into one exact scalar diagnostic.

---

## 3. Exact bounded version using the support/source operator brackets

Stage 41 already proved the exact operator-selected bracket

`Xi Delta_0(kappa,eta) <= Pe_* <= Xi Delta_inf(kappa,eta),`

with

`Delta_0(kappa,eta)`
`= eta (cosh(alpha)-1) / [ alpha^2 ( alpha sinh(alpha) + eta cosh(alpha) ) ],`

`Delta_inf(kappa,eta)`
`= [ cosh(alpha) + (eta/alpha) sinh(alpha) - 1 ]`
`  / [ alpha sinh(alpha) + eta cosh(alpha) ],`

`alpha = sqrt(kappa).`

Since `zeta_phys(Pe,eta;kappa)` is strictly increasing in `Pe`, the physical support ratio obeys the exact bounds

`zeta_-(Xi,eta;kappa)`
`:= zeta_phys( Xi Delta_0(kappa,eta), eta; kappa )`
`<= zeta_phys(Pe_*,eta;kappa)`
`<= zeta_phys( Xi Delta_inf(kappa,eta), eta; kappa )`
`: = zeta_+(Xi,eta;kappa).`

Therefore the exact residual is trapped between

`R_-(...) := zeta_req - zeta_+ <= R_quad <= zeta_req - zeta_- =: R_+(...).`

So the reduced PDE already has theorem-level success/failure bounds:

- if `R_+ <= 0`, success is guaranteed;
- if `R_- > 0`, failure is guaranteed;
- only the strip `R_- <= 0 < R_+` requires the full fixed-point solve.

---

## 4. Direct thresholds in the branch-product variable

Because `zeta_req` is strictly increasing in `Pi_tr`, the bounded residual can be translated back into exact product thresholds.

Write

`Q(zeta;eps_blk)`
`:= [ 1 + (1 - 2 eps_blk) zeta ] / [ 1 - eps_blk zeta ].`

Then the inverse map is

`Pi_tr = C_mix Q(zeta_req;eps_blk).`

So the bounded reduced-PDE theorem is:

`Pi_tr <= Pi_suff(Xi,eta,kappa ; C_mix,eps_blk)`  -> guaranteed success,

`Pi_tr >= Pi_fail(Xi,eta,kappa ; C_mix,eps_blk)`  -> guaranteed failure,

with

`Pi_suff := C_mix Q( zeta_-(Xi,eta;kappa) ; eps_blk ),`

`Pi_fail := C_mix Q( zeta_+(Xi,eta;kappa) ; eps_blk ).`

This is the first exact product-window theorem derived directly from the coupled support/source operator rather than through the intermediate `Pe_req` bookkeeping.

---

## 5. Family-1 specialization of the master residual

The explicit Family-1 branch now drops out immediately by inserting the fixed geometry/support data

`kappa_F1 = 12321/5,`

`eta_F1 = 37,`

and the wall/source strength relation

`Xi_F1 = W_wall = Upsilon_w Lambda_ell^2 = 1369 Upsilon_w = 136900 Theta_w.`

So the Family-1 residual is simply

`R_quad^(F1)`
`= zeta_req(Pi_tr,C_mix,eps_blk)`
` - zeta_phys(Pe_*(Xi_F1,37,12321/5), 37; 12321/5).`

That is exactly the residual the later explicit branch stages had been approaching indirectly.

---

## 6. Why Stage 65 matters

This stage is the first place where the whole reduced moving-throat PDE can be written in one line without hiding the actual theorem gap.

The remaining unresolved problem is not another support/source reduction and not another wall-depth estimate. It is the sign of one scalar residual:

`R_quad(Pi_tr,C_mix,eps_blk ; Xi,eta,kappa).`

Everything else in the reduced program now feeds this one object.
