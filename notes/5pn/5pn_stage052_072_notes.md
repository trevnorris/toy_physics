# 5PN stages 52–72 notes

This bundle continues the moving-throat/support-source chain from the resonance-threshold point
through the explicit Family-1 branch closure and the first direct comparison with the minimal
isotropic passive/outgoing quadrupole module.

## Stage 52 — final reduced verdict for the support/source program

The matched-branch theorem remains

- `W_wall <= Pe_req / Delta_inf`  -> fail,
- `W_wall >= Pe_req / Delta_0`    -> succeed,

and the explicit sech–Gaussian resonance family only perturbs those thresholds by the tiny factor
`P_res = 1.005612487760576...`, so the genuinely profile-sensitive sub-bands are only about `0.56%`
wide.

## Stage 53 — explicit GNLS wall-shell reduction

Projecting the parent GNLS quadratic shell energy onto the wall-support mode `q(s) chi_phi(y)` gives

- `T_X = hbar^2 N_(phi phi) / (4 m rho_w)`,
- `K_X = H_w N_(phi phi) + hbar^2 G_(phi phi)/(4 m rho_w)`,

and on the matched thin-wall branch the support/source fixed-point coupling is exactly the wall
figure of merit:
`Xi = W_wall`.

## Stage 54 — canonical tanh-wall branch

For the canonical wall
`f(xi) = (1 + tanh xi)/2`, `chi_phi = f'(xi) = (1/2) sech^2 xi`,
the exact moments are

- `I_f = 1/3`,
- `I_g = 4/15`,
- `I_g / I_f = 4/5`.

With the natural local Robin closure `K_m = T_X / ell`, the explicit branch is

- `eta = L/ell`,
- `kappa = 4 (m c_(s,w) L / hbar)^2 + (4/5) (L/ell)^2`,
- `W_wall = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2)`.

## Stage 55 — explicit branch thresholds

Writing
`chi_s = m c_(s,w) L / hbar`,
`Lambda_ell = L/ell`,
`Upsilon_w = 4 rho_w^2 V0^2 /(hbar^2 c_(s,w)^2)`,
the branch map is

- `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2`,
- `eta = Lambda_ell`,
- `W_wall = Upsilon_w Lambda_ell^2`.

So the explicit thresholds are

- `Upsilon_fail = Pe_req / [Lambda_ell^2 Delta_inf(kappa,eta)]`,
- `Upsilon_suff = Pe_req / [Lambda_ell^2 Delta_0(kappa,eta)]`.

The shell-gradient and compression-dominated asymptotics from the notes are verified directly.

## Stages 56–57 — Family-1 geometry map and healing lock

Using the carried reference values

- `L/a = 37/20`,
- `ell/a = 1/20`,

the Family-1 branch fixes

- `Lambda_ell = 37`,
- `eta = 37`.

With the healing-width closure `ell = hbar/(2 m c_(s,w))`, the same branch fixes

- `chi_s = 37/2`,
- `kappa = 12321/5`,
- `alpha = sqrt(kappa) = 111/sqrt(5) ≈ 49.6407091`.

So the only remaining explicit branch input is the wall-depth amplitude.

## Stage 58 — explicit Family-1 threshold window

At `(kappa,eta) = (12321/5, 37)`, the exact support/source scales are

- `Delta_0 ≈ 1.73302079021525e-4`,
- `Delta_inf ≈ 2.01447565540522e-2`.

Hence

- `Upsilon_fail ≈ 0.0362605617972939 Pe_req`,
- `Upsilon_suff ≈ 4.21495341569977 Pe_req`,
- `Theta_fail ≈ 3.62605617972939e-4 Pe_req`,
- `Theta_suff ≈ 4.21495341569977e-2 Pe_req`

after `Upsilon_w = 100 Theta_w`.

## Stages 59–60 — exact `n=5` wall-depth lock and shell-weighted extraction

For the frozen `n=5` EOS,
`h(rho) = m c_s(rho)^2 / 4`, so with
`mu_* = lambda_mu h_w`
and the healing-width lock one gets

- `Theta_w = lambda_mu^2 rho_w^2 /(16 ell^2)`.

On the normalized Family-1 wall this becomes
`Theta_w = 25 lambda_mu^2 rho_w^2`.

Using the explicit Family-1 wall profile and canonical support weight gives

- `<rho_r>_chi ≈ 0.192619005556493`,
- `<rho_r^2>_chi ≈ 0.162745294003265`,
- `Theta_w^(chi) ≈ 4.06863235008162 lambda_mu^2`,
- `Theta_w^(J)   ≈ 0.927552032539308 lambda_mu^2`.

## Stage 61 — explicit Family-1 wall-depth verdict

Comparing the extracted `Theta_w` values to the Stage-58 window gives

- `Pe_suff^(chi) ≈ 96.5285247264386 lambda_mu^2`,
- `Pe_fail^(chi) ≈ 11220.5441626259 lambda_mu^2`,
- `Pe_suff^(J)   ≈ 22.0062226330754 lambda_mu^2`,
- `Pe_fail^(J)   ≈ 2558.01892349205 lambda_mu^2`.

So the explicit Family-1 wall-depth supply is not naturally starved; the remaining question becomes the
quadrupole-side demand `Pe_req`.

## Stages 62–64 — Family-1 demand map in `Pe`, `zeta_req`, and `Pi_tr/C_mix`

For the Family-1 branch,

- `y_F1 tan y_F1 = 37`,
- `A_F1 = (kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2) ≈ 1.00005192880220`.

The exact constructive support map is
`zeta_F1(Pe) = A_F1 Omega_Pe^2`,
with hard ceiling
`zeta_max^(F1) = A_F1 pi^2/4 ≈ 2.46752922945601`.

Converting the Stage-61 `Pe_req` windows through this map gives, at `lambda_mu = 1`,

- `zeta_suff^(chi) ≈ 2.46622291347846`,
- `zeta_fail^(chi) ≈ 2.46752913273870`,
- `zeta_suff^(J)   ≈ 2.44257571477179`,
- `zeta_fail^(J)   ≈ 2.46752736855058`.

Using
`Pi_tr = C_mix Q(zeta;eps_blk)`,
`Q(zeta;eps_blk) = [1 + (1 - 2 eps_blk) zeta]/[1 - eps_blk zeta]`,
the unblocked (`eps_blk = 0`) explicit product window becomes

- `Pi_tr/C_mix <= 3.46622291347846`  -> guaranteed success,
- `Pi_tr/C_mix >= 3.46752913273870`  -> guaranteed failure,
- `Pi_tr/C_mix <  3.46752922945601`  -> hard explicit ceiling.

## Stages 65–70 — master reduced quadrupole residual and pure loading-ratio collapse

The whole reduced moving-throat PDE is compressed to one scalar residual
`R_quad = zeta_req - zeta_phys(Pe_*)`,
with `Pe_*` selected by the support/source fixed-point law.

After the exact demand-cancellation step, the normalized support theorem depends only on the loading ratio

`rho_alpha = alpha_req / alpha_mix`,

with support demand
`zeta_req = (rho_alpha - 1) / [1 - eps_blk (2 - rho_alpha)]`.

So the explicit Family-1 theorem is finally

- `rho_alpha <= rho_suff^(chi)(lambda_mu;eps_blk)`  -> guaranteed success,
- `rho_alpha >= rho_fail^(chi)(lambda_mu;eps_blk)`  -> guaranteed failure,
- `rho_alpha < rho_max^(F1)(eps_blk)`               -> hard ceiling,

and in the natural unblocked case

- `rho_alpha <= 3.46622291347846`  -> guaranteed success,
- `rho_alpha >= 3.46752913273870`  -> guaranteed failure,
- `rho_alpha <  3.46752922945601`  -> hard ceiling.

So by Stage 70 the explicit support/source side is completely reduced to one number:
the outgoing-branch loading ratio `rho_alpha`.

## Stages 71–72 — loading ratio from the minimal isotropic conservative quadrupole module

The natural contact-plus-pole reading of the minimal isotropic conservative precursor

`Y_Q^cons(omega) = c0 + c1/(1 - omega^2/Omega_Q^2)`

gives

- `c0 = 1/rho_alpha`,
- `c1 = (rho_alpha - 1)/rho_alpha`,
- `zeta_req = c1/c0`.

For the minimal isotropic module
`c0 = 3/4`, `c1 = 1/4`,
the exact loading data are

- `rho_alpha = 4/3`,
- `zeta_req = 1/3`,
- `Pi_tr = (4/3) C_mix`.

This lies in the symmetric-lowest-twin regime and, on the explicit Family-1 branch, it already satisfies

- `zeta_req = 1/3 < A_F1 ≈ 1.00005192880220`.

So the explicit Family-1 branch succeeds at **zero transport bias** on the natural minimal isotropic passive/outgoing quadrupole module.

## Bottom line from Stage 72

By the end of this bundle, the explicit Family-1 support/source side is no longer the active reduced bottleneck. On the natural minimal isotropic contact-plus-pole branch it is already comfortably satisfied. The remaining reduced theorem question becomes whether the actual grouped-`P2` / geometry moving-throat branch really realizes that minimal isotropic passive/outgoing quadrupole module, which is the right bridge back into the 5PN grouped-`P2` program.
