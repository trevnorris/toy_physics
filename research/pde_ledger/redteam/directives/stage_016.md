---
unit_id: 016
batch: I.2
created_at: 2026-05-21T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-21T13:19:53-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 016

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — missing_verification_script

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage016_parent_throat_action_candidate_mathematica_audit.wl` (CREATE, new file)

**Issue:** The unit is non-checkpoint, non-status-only, but has only a SymPy script. No Mathematica audit exists to cross-verify the parent-throat action's exact Euler-Lagrange equation, the quadratic-fluctuation `K_eta` coefficient, the boundary discharges, or the Y20 angular sector's `+6 T_Omega` stiffness. Create a Mathematica script that independently derives the same eleven claims using Mathematica's own primitives (do not transliterate the SymPy code line-for-line). The script must `Exit[1]` on any failure so the verifier can detect regressions.

**Required change:**

Create the file
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage016_parent_throat_action_candidate_mathematica_audit.wl`
that independently verifies each claim below. Each check must Print a clear pass/fail line and `Exit[1]` on failure. At least one sign-mutation or value-mutation per group must also be exercised (e.g., a flipped sign that the script confirms is nonzero).

Use these Mathematica idioms, not SymPy mirroring:

- For M1, declare `R[t,w,u,v]`, `mu[R[t,w,u,v], w]`, `Tw[R[t,w,u,v], w]`, `TO[R[t,w,u,v], w]`, `USig[R[t,w,u,v], w]` as dependent symbols. Compute the Euler-Lagrange equation via direct `D` chains: `EL = D[L, R[t,w,u,v]] - D[D[L, D[R[t,w,u,v], t]], t] - D[D[L, D[R[t,w,u,v], w]], w] - D[D[L, D[R[t,w,u,v], u]], u] - D[D[L, D[R[t,w,u,v], v]], v]`. Use `Simplify` to compare to the hand-written form.
- For M2-M5, build `L` symbolically in `eps`, then `D[L, eps] /. eps -> 0` for `L1` and `D[L, {eps, 2}] /. eps -> 0` divided by 2 for `L2_raw`. Then construct `K_eta` exactly as in M5 below and verify the IBP-restructured quadratic density matches `(1/2)(mu eta_t^2 - Tw eta_w^2 - TO grad2 - K_eta eta^2)`. The IBP substitution for the cross term `-A eta eta_w` is `-A eta eta_w -> -(1/2) D[A, w] eta^2 + (boundary)`; verify by the product rule: `D[A eta^2, w] = D[A, w] eta^2 + 2 A eta eta_w`, so `2 A eta eta_w = D[A eta^2, w] - D[A, w] eta^2`, hence `-A eta eta_w = -(1/2) D[A eta^2, w] + (1/2) D[A, w] eta^2`; bulk leftover is `+(1/2) D[A, w] eta^2`.
- For M6-M8, use `Limit[expr, w -> Infinity]` and `Limit[expr, w -> -Infinity]` to compute boundary discharges. Use `Integrate[..., {w, -Infinity, Infinity}]` for IBP residuals. Use `Together` to confirm the Lorentzian linear-probe denominator is nontrivial (`1 + w^2`).
- For M9-M10, use `SphericalHarmonicY[2, 0, th, ph]` (or `ExpToTrig` of it). Compute the spherical Laplacian as `Simplify[D[Sin[th] D[Y20, th], th] / Sin[th] + D[Y20, {ph, 2}] / Sin[th]^2]`. Use `Integrate[Sin[th] Y20^2, {ph, 0, 2 Pi}, {th, 0, Pi}]` and `Integrate[Sin[th] (D[Y20, th]^2 + D[Y20, ph]^2 / Sin[th]^2), {ph, 0, 2 Pi}, {th, 0, Pi}]`.
- For M11, declare `q[t, w]`, `mumode[w]`, `Twmode[w]`, `TOmode[w]`, `Kmode[w]` as dependent symbols. Build the closed modal density `(1/2)(mumode q_t^2 - Twmode q_w^2 - (Kmode + 6 TOmode) q^2)`, then derive the modal EL via direct `D` chains and compare to `D[mumode q_t, t] - D[Twmode q_w, w] + (Kmode + 6 TOmode) q`.

**Claim manifest** (for missing-script findings):

- **M1** (Exact nonlinear EL): With `L = (1/2) mu(R,w) R_t^2 - (1/2) Tw(R,w) R_w^2 - (1/2) TO(R,w) (R_u^2 + R_v^2) - USig(R,w)` and `R = R[t,w,u,v]`, verify
  `Simplify[EL_via_D - (D[mu R_t, t] - D[Tw R_w, w] - D[TO R_u, u] - D[TO R_v, v] - (1/2) D[mu, R] R_t^2 + (1/2) D[Tw, R] R_w^2 + (1/2) D[TO, R] (R_u^2 + R_v^2) + D[USig, R])] === 0`.
  Then mutate by flipping the sign of the `D[USig, R]` term (multiplied by 2) and verify `Simplify[...] =!= 0`.

- **M2** (Linear density before IBP): Expand `L` symbolically in `eps` with `R = R0(w) + eps*eta(t,w,Omega)`, `mu = mu0`, `Tw = Tw0 + eps*TwR0*eta + (eps^2/2)*TwRR0*eta^2`, `TO = TO0`, `USig = U0 + eps*UR0*eta + (eps^2/2)*URR0*eta^2`. Compute `L1 = D[L, eps] /. eps -> 0` and verify
  `Simplify[L1 - (-(1/2) TwR0 R0p^2 eta - Tw0 R0p eta_w - UR0 eta)] === 0`.

- **M3** (Linear background coefficient after IBP): Verify the bulk after IBP equals
  `E_bg = d_Tw_R0p - (1/2) TwR0 R0p^2 - UR0`,
  i.e. `Simplify[(L1 + Tw0 R0p eta_w) + d_Tw_R0p eta - E_bg eta] === 0`. The boundary leftover is `[-Tw(R0,w) R0' eta]`.

- **M4** (Raw quadratic density): Compute `L2raw = D[L, {eps, 2}] /. eps -> 0` then divide by 2 and verify
  `Simplify[L2raw - ((1/2) mu0 eta_t^2 - (1/2) Tw0 eta_w^2 - (1/2) TO0 grad2 - TwR0 R0p eta eta_w - (1/4) TwRR0 R0p^2 eta^2 - (1/2) URR0 eta^2)] === 0`.
  Verify the cross-term coefficient explicitly: `Simplify[D[D[L2raw, eta], eta_w] + TwR0 R0p] === 0`.

- **M5** (K_eta after cross-term IBP): With cross-term IBP substitution `-A eta eta_w -> (1/2) d_TwR_R0p eta^2` (bulk only), verify
  `Simplify[L2raw - (-TwR0 R0p eta eta_w) + (1/2) d_TwR_R0p eta^2 - ((1/2) mu0 eta_t^2 - (1/2) Tw0 eta_w^2 - (1/2) TO0 grad2 - (1/2) K_eta eta^2)] === 0`,
  where `K_eta = URR0 - d_TwR_R0p + (1/2) TwRR0 R0p^2`. Then mutate the middle sign (`+ d_TwR_R0p` instead of `-`) and confirm the residual is nonzero.

- **M6** (Concrete boundary discharge): With `Bcon = (1 + w^2) Exp[-w^2]`, `Acon = (1 + w^2/2) Exp[-w^2]`, `etaG = Exp[-w^2/2]`, `etaL = 1/(1 + w^2)`, `Blor = Exp[-w^2]`, verify
  `Limit[-Bcon etaG, w -> Infinity] - Limit[-Bcon etaG, w -> -Infinity] === 0`,
  `Limit[-Blor etaL, w -> Infinity] - Limit[-Blor etaL, w -> -Infinity] === 0`,
  `Limit[-Acon etaG^2 / 2, w -> Infinity] - Limit[-Acon etaG^2 / 2, w -> -Infinity] === 0`,
  `Limit[-Acon etaL^2 / 2, w -> Infinity] - Limit[-Acon etaL^2 / 2, w -> -Infinity] === 0`.
  Also verify the concrete IBP integral identities:
  `Integrate[-Bcon D[etaG, w], {w, -Infinity, Infinity}] === Integrate[D[Bcon, w] etaG, {w, -Infinity, Infinity}]`
  and the analogous quadratic identity.

- **M7** (Lorentzian finite-endpoint discharge): With `Bend = w Sqrt[1 + w^2]`, `etaL = 1/(1 + w^2)`, verify
  `Limit[-Bend etaL, w -> Infinity] - Limit[-Bend etaL, w -> -Infinity] === -2`,
  i.e. the script asserts `... + 2 === 0`.

- **M8** (Nonzero boundary probe): Verify
  `Limit[ArcTan[w], w -> Infinity] - Limit[ArcTan[w], w -> -Infinity] === Pi` (nonzero).
  Also verify `Together[-Exp[-w^2]/(1 + w^2)]` has a nontrivial denominator `(1 + w^2)`.

- **M9** (Y20 spherical-Laplacian eigenvalue): With `Y20 = SphericalHarmonicY[2, 0, th, ph]`, verify
  `Simplify[D[Sin[th] D[Y20, th], th] / Sin[th] + D[Y20, {ph, 2}] / Sin[th]^2 + 6 Y20] === 0`,
  and confirm the same expression with `5 Y20` (not 6) is `=!= 0`.

- **M10** (Angular norm and stiffness): Verify
  `Integrate[Sin[th] Y20^2, {ph, 0, 2 Pi}, {th, 0, Pi}] === 1` and
  `Integrate[Sin[th] (D[Y20, th]^2 + D[Y20, ph]^2 / Sin[th]^2), {ph, 0, 2 Pi}, {th, 0, Pi}] === 6`.

- **M11** (Modal Euler-Lagrange): With `q[t, w]`, `mumode[w]`, `Twmode[w]`, `TOmode[w]`, `Kmode[w]` declared as dependent symbols and closed modal density
  `Lmodal = (1/2)(mumode[w] D[q[t,w], t]^2 - Twmode[w] D[q[t,w], w]^2 - (Kmode[w] + 6 TOmode[w]) q[t,w]^2)`,
  compute `EL = D[Lmodal, q[t,w]] - D[D[Lmodal, D[q[t,w], t]], t] - D[D[Lmodal, D[q[t,w], w]], w]` and verify
  `Simplify[EL + (D[mumode[w] D[q[t,w], t], t] - D[Twmode[w] D[q[t,w], w], w] + (Kmode[w] + 6 TOmode[w]) q[t,w])] === 0`.

The script must structure as:

```
(* Header: paths, claim summary *)
On[Assert];
(* M1: exact EL ... *)
(* M2: linear density ... *)
...
(* M11: modal EL ... *)
Print["STATUS: PASS"];
Exit[0];
```

Each `M*` block:
```
ok = Simplify[<residual>];
If[ok =!= 0, Print["[M*] FAIL: ", ok]; Exit[1]];
Print["[M*] PASS"];
(* mutation check where required *)
muto = Simplify[<mutated residual>];
If[muto === 0, Print["[M* mutation] FAIL: residual collapsed"]; Exit[1]];
Print["[M* mutation] PASS"];
```

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 016` and confirms (a) the new file exists at the path above, (b) the script exits 0, (c) the output transcript contains `[M1] PASS` through `[M11] PASS` plus the corresponding mutation-check passes, and (d) the derived `K_eta` symbolic form and Y20 angular factors agree with the SymPy output transcript lines 38, 44, and 54 of
`/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage016_parent_throat_action_candidate_sympy_audit.txt`.

## Applied: F1

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage016_parent_throat_action_candidate_mathematica_audit.wl`
- summary: Created the Stage 016 Mathematica audit covering M1-M11 with direct residual checks and mutation guards.
- deviation: none
