# 5PN continuation notes — Stages 21–23 and 27–29

This note fills the earlier numbering gap and records the next continuation past the Stage 26 placement map.

## Missing stages 21–23

### Stage 21 — Dimensionless continuum placement map
The continuum-selected quadrupole placement problem compresses to the exact dimensionless ratios

a) `eps_eta = c_(etaU)^2 / (K_U K_eta^(eff))`

b) `eps_W = c_(UW)^2 sigma / (K_U K_W^(eff))`

c) `rho = c_(UW) c_(etaU) / (K_U c_(etaW))`

d) `Z_W = c_(etaW)^2 / (K_eta^(eff) K_W^(eff))`

e) `delta_0 = pi^2 T_w / (L^2 K_eta^(eff))`

with radiative scale

`Lambda = 27 pi^2 G c_s^5 K_W^(eff) / (20 a^5 c^5 mu_W)`.

The exact placement coordinates are

`delta = delta_0 / (1 - eps_eta)`,

`M_mix = 8 Z_W (1 + rho)^2 / [pi^2 (1 - eps_eta)(1 - eps_W)]`,

`R_target = Lambda (1 - eps_eta)(1 - eps_W)^2 / [Z_W (1 + rho)^2]`.

The key exact product law is

`R_target M_mix = 8 Lambda (1 - eps_W) / pi^2`

`= 54 G c_s^5 K_W^(eff) (1 - eps_W) / (5 a^5 c^5 mu_W)`.

So Stage 21 factorizes the continuum map into a geometry lane, a mixed-stability/product lane, and a redistribution lane.

### Stage 22 — First non-flat U doublet
Turning on the first axial structure of the internal `U` sector preserves the scalar placement map but breaks the old one-direction geometry.

The new exact quantities are

`delta_split = [delta_0 + eps_eta delta_U/(1+delta_U)] / (1-eps_eta)`,

`eps_W_split = eps_W [1 - (2/11) delta_U/(1+delta_U)]`,

`R_U = [1 + rho_0/(1+delta_U)] / (1 + rho_0)`.

The exact direction-splitting invariant is

`D_dir = - kappa_0 kappa_1 g_W rho_0 delta_U / (1+delta_U)`.

So collinearity survives iff `delta_U = 0` or `rho_0 = 0`.

### Stage 23 — Generalized selected-branch normalization
With distinct loading vector `z` and source vector `s`, the selected branch survives exactly, but the old flat-`U` function pair `(F,G)` deforms to

`F_(q,eta)(xi,delta)` and `G_q(xi,delta)`.

For the split-`U` continuum, the whole deformation collapses to the single ratio `R_U`, giving the exact one-parameter family

`F_U(xi,delta;R_U)`, `G_U(xi,delta;R_U)`.

Setting `R_U = 1` recovers the Stage 18–19 flat branch exactly.

The exact first deformation about the flat branch is

`F_U/F_flat = 1 + eps H_F + O(eps^2)`,

`G_U/G_flat = 1 + eps H_G + O(eps^2)`,

with the closed-form `H_F` and `H_G` verified in the companion SymPy audit.

## Continuation past Stage 26: Stages 27–29

### Stage 27 — Continuum-selected rank-2 closure
Once the actual support-direction data from Stage 26 are inserted, the physical selected wall branch is pinned by an exact quadratic for the softening depth `xi`:

`xi^2 + B_cont xi + C_cont = 0`,

with

`B_cont = delta - M_mix (1+t^2 R_U^2) - M_supp (1+t^2 R_phi^2)`,

`C_cont = - delta (M_mix + M_supp) + t^2 M_mix M_supp (R_U - R_phi)^2`.

So the continuum-selected normalization theorem gate becomes simply

`R_target = F_cont(xi_phys)`.

The exact special surfaces are:

- minimal-kernel source-tied surface: `sigma_0 = 0`,
- interference-matched tracking surface: `sigma_0 = rho_0`.

### Stage 28 — Coherent local D/N support kernel
The first coherent local D/N kernel forces

`g_B g_R = g_W g_S`,

so it lands exactly on the tracking surface. Therefore

`R_phi = R_U = R_tr`,

and the full Stage 27 rank-2 closure collapses to the one-parameter tracking laws

`M_tr = G_tr(xi,delta;R_tr)`,

`R_target = F_tr(xi,delta;R_tr)`.

The coherent branch also has the exact constructive bound

`1/(1+delta_U) < R_tr < 1`.

### Stage 29 — Tracking-branch bounds and residual comparison
On the coherent tracking branch, at fixed `(xi,delta)`,

`dG_tr/dR < 0`,

`dF_tr/dR > 0`.

So for the constructive split-`U` branch, where `R_tr < 1`, one has the exact inequalities

`G_tr > G_flat`,

`F_tr < F_flat`.

The endpoint formulas are

`G_tr(R=0) = xi`,

`F_tr(R=0) = 1/(1-xi)`.

The exact residual comparison theorem is

`E_tr - E_flat = F_flat - F_tr > 0`.

So the first coherent local D/N kernel does **not** rescue the Stage 18–19 normalization target. It makes the target harder: more total loading is required, while the normalized selected-branch response is lower than on the flat branch.

## Bottom line after this batch

The selected-branch problem is now narrower than before:

1. the physical softening depth is set by an exact quadratic;
2. the first coherent local D/N kernel forces exact tracking rather than a generic intermediate rank-2 closure;
3. and the resulting split-`U` deformation worsens the normalization residual relative to the old flat branch.

So the next theorem gate is no longer “which closure should we use?” It is whether the later support-enhancement / compensation mechanisms in the notes are strong enough to overcome this exact tracking-branch deficit.

## Further continuation — Stages 30–31

### Stage 30 — Coherent-kernel dimensionless map
The first coherent local D/N kernel compresses to a very small exact parameter set:

- `eps_eta`
- `delta_U`
- `chi_0`
- `eps_W`
- `Z_W`
- `zeta`
- `Lambda`

The exact coherent tracking factor is

`R_tr = [1 + chi_0/(1+delta_U)] / (1 + chi_0)`.

The total baseline splits as

`M_mix = 8 Z_W (1 + chi_0)^2 / [pi^2 (1 - eps_eta)(1 - eps)]`,

`M_supp = 8 zeta Z_W (1 + chi_0)^2 / [pi^2 (1 - eps_eta)(1 - zeta eps)]`,

with

`eps = eps_W [1 - (2/11) delta_U/(1+delta_U)]`.

So the support lane enters only through the exact enhancement factor

`S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps)`,

and

`M_tr = M_mix S(zeta;eps)`.

The target ratio remains independent of `zeta`:

`R_target = Lambda (1 - eps_eta)(1 - eps)^2 / [Z_W (1 + chi_0)^2]`.

### Stage 31 — Support-compensation theorem
On the coherent tracking branch,

`M_tr = G_tr(xi,delta;R_tr)`

with critical load

`M_crit(delta,R_tr) = G_tr(1,delta;R_tr)`.

The support enhancement factor is strictly increasing and exactly invertible:

`zeta_req = (S_req - 1) / [1 + eps (S_req - 2)]`.

For every finite target ratio `R_target > 1`, there exists a stable-side branch point `xi_req in (0,1)` and a corresponding required load `M_req < M_crit`.

Therefore, if the mixed-only branch is below that required load, there is a unique coherent support ratio `zeta_req < zeta_crit < 1/eps` that reaches the target **before** the branch softens out.

So there is no reduced-level support no-go. The remaining question is whether the actual PDE produces a physical `zeta` large enough to meet that exact threshold.

## Further continuation — Stages 32–33

### Stage 32 — Explicit D/N overlap extraction of `zeta`
The coherent support ratio is no longer phenomenological. For the first explicit finite-throat D/N family,

`zeta_n^(phys) = (K_W^(eff) / K_(phi,n)^(eff)) (I_n / I_0)^2`,

with

`I_n = ∫_0^L ds sigma(s) chi_n(s)`,

`chi_n(s) = sqrt(2/L) sin[(n+1/2) pi s / L]`.

For the first uniform local source density `sigma(s)=1`,

`I_n / I_0 = 1 / (2n+1)`.

So the physical coherent support ratio becomes

`zeta_n^(phys) = (K_W^(eff) / K_(phi,n)^(eff)) / (2n+1)^2`.

On the same-operator twin family,

`zeta_n^(twin) = 1 / [ (2n+1)^2 (1 + x n(n+1)) ]`.

The exact lowest-twin value is

`zeta_0^(twin) = 1`.

### Stage 33 — Exact comparison of `zeta_phys` against `zeta_req`
Because Stage 31 reduced support feasibility to `zeta_phys >= zeta_req`, the explicit D/N tower gives exact selection rules.

The lowest symmetric twin lane has

`zeta_0^(twin) = 1`,

so its enhancement factor is exactly

`S_0 = 2`,

independent of `eps`.

Therefore the lowest symmetric twin lane succeeds exactly when

`zeta_req <= 1`,

equivalently

`S_req <= 2`.

For higher D/N harmonics,

`zeta_n^(twin) < 1/(2n+1)^2`,

so they are ruled out immediately if

`zeta_req > 1/(2n+1)^2`.

When they are not ruled out immediately, the exact softness threshold is

`x <= x_max(n; zeta_req)`

with

`x_max(n; zeta_req) = [1 / ((2n+1)^2 zeta_req) - 1] / [n(n+1)]`.

So the explicit coherent support tower is strongly biased toward the lowest twin lane.

## Bottom line after Stages 21–33

The program is now much sharper than when the numbering gap first appeared.

1. Stages 21–23 show that split-`U` physics preserves the scalar placement map but deforms the selected-branch normalization geometry through an exact direction-splitting parameter.
2. Stages 27–29 show that once the actual continuum support data are inserted, the physical wall softening depth is fixed by an exact quadratic, and the first coherent local D/N kernel lands on a tracking branch whose split-`U` deformation worsens the normalization target relative to the old flat branch.
3. Stages 30–31 show there is no reduced-level support no-go: coherent support enhancement can in principle compensate the tracking-branch deficit before softening.
4. Stages 32–33 then make that compensation test microscopic: the lowest symmetric twin D/N lane is an exact universal doubling branch, and every higher support harmonic is strongly overlap- and stiffness-suppressed.

So the next really important question is no longer which closure to use. It is whether the actual moving-throat PDE places the physical coherent support sector on the lowest twin branch with a large enough `zeta_phys` to satisfy the exact threshold.
