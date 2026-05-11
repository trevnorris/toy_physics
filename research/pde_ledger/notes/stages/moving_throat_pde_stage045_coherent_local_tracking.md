# Moving-Throat PDE — Stage 045: Coherent Local D/N Support Kernel and the Exact Tracking Reduction

## Purpose

Stage 27 showed that the first extended continuum kernel generically lands on an exact **intermediate** rank-2 closure, with the source-tied and tracking laws appearing only as special surfaces.

That was the correct generic reduced statement, but it still left one obvious question open:

> what happens if the support field `phi` and the mixed field `W` really are two D/N channels of the **same local throat-support density** rather than two unrelated reduced couplings?

That is the first honest “concrete kernel” question after Stage 27.

This stage answers it.

The main result is that, **within this coherent local-kernel hypothesis**, the first coherent local D/N kernel does **not** land on the generic intermediate closure. It lands exactly on the tracking surface of the reduced Stage-27 branch problem.

So the Stage-27 continuum-selected problem collapses back to a one-parameter selected branch, now with a physically identified total loading

`M_tr = M_mix + M_supp`

and a single common direction factor

`R_tr = R_U = R_phi`.

This is the first explicit evaluation of the Stage-27 residual on a concrete moving-throat kernel family.

---

## 1. Coherent local D/N support kernel

Keep the Stage-22/26 split-`U` continuum operator on the finite throat interval `s in [0,L]`, with

- wall basis `(u_0,u_1)`,
- split internal doublet `U`,
- mixed D/N half-wave `W`,
- support/BdG D/N half-wave `phi`.

Now impose the first concrete local-kernel hypothesis:

> the mixed lane `W` and the support lane `phi` couple to the **same local wall/U support density**.

The minimal local interaction density is

`L_int^(coh)`
`= - int_0^L ds [ lambda_W W(s,t) + lambda_phi phi(s,t) ] [ eta(s,t) - gamma U(s,t) ].`

So the two D/N lanes differ only by their amplitudes `lambda_W` and `lambda_phi`; their wall/U source combination is the same.

After the same modal projection used in Stages 20, 22, and 26, this gives

`L_(eta W)   = - g_W (v.Q) Wbar,`

`L_(UW)      = + g_R (v.U) Wbar,`

`L_(eta phi) = - g_B (v.Q) Phi,`

`L_(Uphi)    = + g_S (v.U) Phi,`

with the exact amplitude pattern

`g_W = lambda_W / sqrt(mu_eta mu_W),`

`g_R = gamma lambda_W / sqrt(mu_U mu_W),`

`g_B = lambda_phi / sqrt(mu_eta mu_phi),`

`g_S = gamma lambda_phi / sqrt(mu_U mu_phi).`

This is more specific than the Stage-26 generic extension: the same local support density forces the wall/U ratio to be identical in the mixed and support lanes.

---

## 2. Exact tracking theorem

Stage 26 proved that support and mixed directions coincide iff

`g_B g_R = g_W g_S`.

For the coherent local kernel above,

`g_B g_R`
`= gamma lambda_W lambda_phi / sqrt(mu_eta mu_phi mu_U mu_W),`

`g_W g_S`
`= gamma lambda_W lambda_phi / sqrt(mu_eta mu_W mu_U mu_phi),`

so

`g_B g_R - g_W g_S = 0`

exactly.

Therefore the Stage-26 codimension-one interference-match surface is now automatic.

### Exact conclusion

> **For the first coherent local D/N support kernel, the reduced Stage-27 branch lands exactly on the tracking surface.**

Equivalently,

`R_phi = R_U`.

So the generic Stage-27 intermediate closure is not the first coherent local-kernel outcome. For this coherent local-kernel model, the outcome is the exact tracking reduction.

---

## 3. Exact common interference ratio and direction factor

Stage 26 defines

`rho_0 = g_R g_U / (K_U g_W),`

`sigma_0 = g_U g_S / (K_U g_B)`.

Under the coherent local kernel,

`rho_0 = sigma_0`.

Define the common value by

`chi_0 := rho_0 = sigma_0.`

Then the exact common direction factor is

`R_tr := R_U = R_phi`
`     = [ 1 + chi_0/(1+delta_U) ] / (1 + chi_0).`

Two exact identities are useful:

`1 - R_tr`
`= chi_0 delta_U / [ (1 + chi_0)(1 + delta_U) ],`

`R_tr - 1/(1+delta_U)`
`= delta_U / [ (1 + chi_0)(1 + delta_U) ].`

So for the constructive branch `chi_0 > 0`, `delta_U > 0`,

`1/(1+delta_U) < R_tr < 1.`

That gives the exact physical range of the tracking factor on this branch.

---

## 4. Exact total loading on the coherent branch

Stage 22 and Stage 26 give the mixed and support baselines

`M_mix`
`= 8 Z_W (1 + chi_0)^2`
`  / [ pi^2 (1 - eps_eta) (1 - eps_W^(split)) ],`

`M_supp`
`= 8 Z_phi (1 + chi_0)^2`
`  / [ pi^2 (1 - eps_eta) (1 - eps_phi^(split)) ].`

Because the branch is exactly tracking, the physical selected-mode problem depends only on the total baseline

`M_tr := M_mix + M_supp,`

so

`M_tr`
`= 8 (1 + chi_0)^2 / [ pi^2 (1 - eps_eta) ]`
`  * [ Z_W / (1 - eps_W^(split)) + Z_phi / (1 - eps_phi^(split)) ].`

So the first concrete local-kernel evaluation of the Stage-27 residual is already much simpler than the generic intermediate closure:

- one common direction factor `R_tr`,
- one total baseline `M_tr`.

---

## 5. Exact collapse of the Stage-27 quadratic branch equation

Stage 27 gave the continuum-selected quadratic branch equation

`xi^2 + B_cont xi + C_cont = 0,`

with mismatch penalty proportional to

`lambda_0 M_mix M_supp (R_U - R_phi)^2.`

On the coherent local-kernel branch,

`R_U - R_phi = 0,`

so the mismatch penalty vanishes identically. The physical branch equation reduces exactly to

`xi^2 + B_tr xi + C_tr = 0,`

with

`B_tr = delta - M_tr ( 1 + lambda_0 R_tr^2 ),`

`C_tr = - delta M_tr.`

Equivalently, the selected branch obeys the exact one-parameter law

`M_tr = G_tr(xi,delta;R_tr)`,

where

`G_tr(xi,delta;R)`
`= 9 xi (xi + delta) / [ 9 delta + (9 + 2 R^2) xi ].`

This is exactly the Stage-23 tracking/source-loading branch, now derived from a concrete local kernel rather than postulated as a special surface.

---

## 6. Exact normalization law on the coherent branch

Because the coherent local kernel lands on the tracking surface, the Stage-27 normalization residual also collapses to the Stage-23 tracking form:

`R_target = F_tr(xi,delta;R_tr)`,

with

`F_tr(xi,delta;R)`
`= [ 9 delta + (9 + 2 R^2) xi ]^2 [ 9 delta + (9 + 2 R) xi ]^2`
`  / [ 81 (1 - xi) ( 9 delta^2 + 18 delta xi + (9 + 2 R^2) xi^2 )^2 ].`

So the full Stage-27 continuum-selected residual now becomes

`R_target - F_tr( xi_phys, delta; R_tr ),`

where `xi_phys` is the physical root of the reduced quadratic above.

This is the first exact “concrete-kernel” form of the normalization test.

---

## 7. Best current theorem statement after Stage 28

### What is now exact

- the first coherent local D/N support kernel,
- the automatic tracking identity `g_B g_R = g_W g_S`,
- the common interference ratio `chi_0`,
- the exact tracking factor `R_tr`,
- the exact total baseline `M_tr`,
- the exact collapse of the Stage-27 quadratic branch equation to the one-parameter tracking law,
- and the exact collapse of the normalization residual to the tracking function `F_tr`.

### What this means physically

The generic Stage-27 intermediate closure is still the correct reduced theorem for the first unrestricted continuum extension. But the first **coherent local D/N kernel** is more special than that: within this reduced coherent model it lands exactly on the tracking surface.

So the next theorem step is no longer to resolve the rank-2 closure ambiguity. That is now done for the first concrete local kernel. The next step is to evaluate how this physical tracking branch compares with the old flat branch and whether the split-`U` deformation helps or hurts the normalization test.
