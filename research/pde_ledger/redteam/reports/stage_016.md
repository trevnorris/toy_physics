---
unit_id: 016
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-21T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 016 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage016_parent_throat_action_candidate_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage016_parent_throat_action_candidate_sympy_audit.txt`
- mathematica output: `(missing)`

## What the script claims to verify

The script claims to verify a minimal nonlinear gauge-fixed parent-throat action
`L_Sigma = (1/2)mu(R,w) R_t^2 - (1/2)T_w(R,w) R_w^2 - (1/2)T_Omega(R,w) |∇_Ω R|^2 - U_Sigma(R,w)`
by (i) deriving the exact nonlinear Euler-Lagrange equation in a local orthonormal angular chart and matching it to a hand-written form; (ii) expanding around a static isotropic background `R = R0(w) + eps*eta` to second order in `eps`, integrating the quadratic cross term `-T_{w,R}(R0,w) R0' eta eta_w` by parts, and asserting the resulting `K_eta(w) = U_{Sigma,RR}(R0,w) - d/dw[T_{w,R}(R0,w) R0'] + (1/2) T_{w,RR}(R0,w) (R0')^2`; (iii) checking generic linear/quadratic IBP identities and concrete boundary-discharge probes on Gaussian/Lorentzian profiles, including a nonzero-endpoint sanity check on `atan(w)`; (iv) computing the `Y_2^0` spherical-Laplacian eigenvalue (-6), the Y20 angular norm (=1), and the Y20 angular gradient norm-squared (=6), and projecting the fluctuation Lagrangian to the closed modal density `(1/2)(mu_mode q_t^2 - Tw_mode q_w^2 - (K_mode + 6 T_Omega) q^2)`; (v) re-deriving the modal Euler-Lagrange equation via sympy's `euler_equations` and matching it to a hand form.

## Assertion inventory

| #  | Script | Line | Form                                                              | Anchored to claim? |
|----|--------|------|-------------------------------------------------------------------|--------------------|
| A1 | sympy  | 53   | `assert_zero(exact_el + exact_el_lhs)`                            | yes                |
| A2 | sympy  | 55   | `assert_nonzero(exact_el + mutated_exact_el_lhs)` (sign mutation) | yes                |
| A3 | sympy  | 65   | `assert_zero(generic linear IBP identity)`                        | yes (product rule) |
| A4 | sympy  | 69   | `assert_zero(generic quadratic IBP identity)`                     | yes (product rule) |
| A5 | sympy  | 73   | `assert_nonzero(mutated quadratic IBP sign)`                      | yes                |
| A6 | sympy  | 84   | `assert_nonzero(atan-w boundary probe)`                           | yes                |
| A7 | sympy  | 88   | `assert_nonzero(lorentz denom)`                                   | yes                |
| A8 | sympy  | 105  | `assert_zero(concrete linear IBP boundary discharge)`             | yes                |
| A9 | sympy  | 106  | `assert_zero(lorentzian linear discharge)`                        | yes                |
| A10| sympy  | 107  | `assert_zero(lorentzian-finite-endpoint discharge + 2)`           | yes                |
| A11| sympy  | 108  | `assert_zero(concrete linear IBP w/ boundary)`                    | yes                |
| A12| sympy  | 109  | `assert_zero(concrete quadratic IBP boundary discharge)`          | yes                |
| A13| sympy  | 110  | `assert_zero(lorentzian quadratic discharge)`                     | yes                |
| A14| sympy  | 111  | `assert_zero(concrete quadratic IBP w/ boundary)`                 | yes                |
| A15| sympy  | 138  | `assert_zero(L1 vs hand form)`                                    | yes                |
| A16| sympy  | 146  | `assert_zero(linear background coefficient after IBP)`            | yes                |
| A17| sympy  | 153  | `assert_zero(raw quadratic cross coefficient)`                    | yes                |
| A18| sympy  | 163  | `assert_zero(quadratic density after IBP = canonical_L2)`         | yes                |
| A19| sympy  | 170  | `assert_nonzero(mutated K_eta sign)`                              | yes                |
| A20| sympy  | 178  | `assert_zero(Y20 spherical-Laplacian + 6 Y20)`                    | yes                |
| A21| sympy  | 179  | `assert_nonzero(Y20 spherical-Laplacian + 5 Y20)`                 | yes                |
| A22| sympy  | 205  | `assert_zero(Y20 modal norm - 1)`                                 | yes                |
| A23| sympy  | 206  | `assert_zero(Y20 modal angular stiffness - 6)`                    | yes                |
| A24| sympy  | 207  | `assert_zero(projected vs closed modal density)`                  | yes                |
| A25| sympy  | 214  | `assert_zero(Y20 fused modal EL equation)`                        | yes                |

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- (proposed) `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage016_parent_throat_action_candidate_mathematica_audit.wl`

**What's wrong:**
Per the unit's manifest entry, stage 016 is neither a checkpoint (`is_checkpoint: false`) nor status-only (`is_status_only_candidate: false`), but only a SymPy script exists. There is no Mathematica audit script anywhere under `/var/projects/toy_physics/research/pde_ledger/mathematica/` for stage 016 (confirmed by directory listing — stage 016 is absent between `stage012` and `stage021`). Without a second independent engine, none of the claims in `assertion inventory` above are cross-checked. The author may have hand-mistranscribed `K_eta = U_RR - d/dw[T_wR R0'] + (1/2)T_wRR (R0')^2`, mis-signed the angular-gradient stiffness (+6 vs -6), or made a sign error in the cross-term IBP rewrite, and the SymPy script alone cannot detect such errors because it derives both sides from the same author's expressions.

**Why this matters:**
The unit derives the `K_eta` coefficient that the entire downstream linear-wall fluctuation theory depends on, plus the angular sector's `+6 T_Omega` Y20 projection that the modal Euler-Lagrange equations feed into. A single-engine verification of these landmark coefficients does not satisfy the second-engine policy for non-checkpoint, non-status-only units.

**Required change:**
Create the file
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage016_parent_throat_action_candidate_mathematica_audit.wl`
that independently re-derives the assertions in the SymPy script using Mathematica primitives (`D`, `Simplify`, `Integrate`, `Limit`, `VariationalMethods` if desired, `SphericalHarmonicY`). The script must `Exit[1]` on any failure (e.g., `If[Simplify[expr] =!= 0, Print["FAIL ..."]; Exit[1]]`) so the verifier can detect regressions. Do not transliterate the SymPy script line-for-line; re-derive each claim using Mathematica idioms (e.g., use `SphericalHarmonicY[2,0,th,ph]` directly; use `Integrate[..., {ph,0,2 Pi}, {th,0,Pi}]`; use `VariationalD` from `VariationalMethods` for the Euler-Lagrange derivation, rather than echoing the hand-written `exact_el_lhs` formula).

**Claim manifest** (M1-M11):
- **M1** Exact nonlinear Euler-Lagrange for `L_Sigma = (1/2) mu(R,w) R_t^2 - (1/2) T_w(R,w) R_w^2 - (1/2) T_Omega(R,w) (R_u^2 + R_v^2) - U_Sigma(R,w)` matches
  `∂_t(mu R_t) - ∂_w(T_w R_w) - ∂_u(T_Omega R_u) - ∂_v(T_Omega R_v) - (1/2) mu_R R_t^2 + (1/2) T_{w,R} R_w^2 + (1/2) T_{Omega,R}(R_u^2 + R_v^2) + U_{Sigma,R} = 0`.
- **M2** Linear-order density (coefficient of `eps^1` at `eps=0` after expanding `L`) equals
  `L1 = -(1/2) T_{w,R}(R0,w) (R0')^2 eta - T_w(R0,w) R0' eta_w - U_{Sigma,R}(R0,w) eta`.
- **M3** Linear bulk coefficient after IBP equals
  `E_bg = d/dw[T_w(R0,w) R0'] - (1/2) T_{w,R}(R0,w) (R0')^2 - U_{Sigma,R}(R0,w)`.
- **M4** Raw quadratic-order density equals
  `L2_raw = (1/2) mu(R0,w) eta_t^2 - (1/2) T_w(R0,w) eta_w^2 - (1/2) T_Omega(R0,w) |∇_Ω eta|^2 - T_{w,R}(R0,w) R0' eta eta_w - (1/4) T_{w,RR}(R0,w) (R0')^2 eta^2 - (1/2) U_{Sigma,RR}(R0,w) eta^2`.
- **M5** After integrating the cross term `-T_{w,R}(R0,w) R0' eta eta_w` by parts, the bulk leftover is `(1/2) d/dw[T_{w,R}(R0,w) R0'] eta^2`, and the resulting canonical quadratic density has stiffness
  `K_eta(w) = U_{Sigma,RR}(R0,w) - d/dw[T_{w,R}(R0,w) R0'] + (1/2) T_{w,RR}(R0,w) (R0')^2`.
- **M6** Concrete-profile boundary-discharge probe for both linear and quadratic boundary terms vanishes on `eta = exp(-w^2/2)` with `B = (1+w^2) exp(-w^2)` and `A = (1 + w^2/2) exp(-w^2)`; vanishes also on `eta = 1/(1+w^2)` with `B = exp(-w^2)`.
- **M7** Lorentzian finite-endpoint probe with `B = w sqrt(1+w^2)` and `eta = 1/(1+w^2)` gives boundary discharge `Limit[-B eta, w->Infinity] - Limit[-B eta, w->-Infinity] = -2`.
- **M8** Nonzero-endpoint sanity probe: `Limit[ArcTan[w], w->Infinity] - Limit[ArcTan[w], w->-Infinity] = Pi` (nonzero).
- **M9** `Y20 = SphericalHarmonicY[2,0,th,ph]` (or its `ExpToTrig`-expanded form) satisfies `(1/Sin[th]) D[Sin[th] D[Y20,th], th] + D[Y20,{ph,2}]/Sin[th]^2 + 6 Y20 = 0` after `Simplify`.
- **M10** Angular integrals: `Integrate[Sin[th] Y20^2, {ph,0,2 Pi}, {th,0,Pi}] = 1` and `Integrate[Sin[th] (D[Y20,th]^2 + D[Y20,ph]^2/Sin[th]^2), {ph,0,2 Pi}, {th,0,Pi}] = 6`.
- **M11** Modal Euler-Lagrange for the closed projected density
  `(1/2)(mu_mode(w) q_t^2 - T_w_mode(w) q_w^2 - (K_mode(w) + 6 T_O_mode(w)) q^2)`
  matches `∂_t(mu_mode q_t) - ∂_w(T_w_mode q_w) + (K_mode + 6 T_O_mode) q = 0`.

**Verification:**
After creation, the verifier runs `redteam exec-mathematica 016`. The new file must exit 0 with explicit prints for each of M1-M11 (e.g., `Print["[M1] PASS"]` after each `If[...; Exit[1]]` guard) and at least one mutation/sanity check per claim (e.g., a sign-flipped variant that triggers `Exit[1]` if reached). The verifier then confirms the resulting `K_eta` and `+6 T_Omega` coefficients agree with the SymPy output transcript at lines 38, 44, and 54 of
`/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage016_parent_throat_action_candidate_sympy_audit.txt`.

## Independent-derivation check (Mathematica)

No `.wl` exists for this unit, so this section is moot for transliteration analysis. See F1 for the missing-script finding. The directive instructs Codex to derive M1-M11 using Mathematica idioms (`SphericalHarmonicY`, `Integrate` over the sphere, `VariationalD` or direct `D` chains, `Limit` for boundary values) rather than line-by-line copying SymPy's `euler_equations` / `sp.diff` choreography.

## Engine cross-check

Only the SymPy engine is present. Cross-check is impossible until the Mathematica script lands.

## Verdict justification

Findings: 1 (missing Mathematica). Attacks I tried on the SymPy script that failed:
- **Cross-term IBP self-consistency**: I expanded `Lexp` by hand to order `eps^2` and matched the script's `L2_raw` term-by-term — every coefficient (including the `-T_{w,R} R0' eta eta_w` cross term, the `-(1/4) T_{w,RR} (R0')^2 eta^2` piece, the `-(1/2) T_Omega grad2` piece, and the `(1/2) mu_0 eta_t^2` kinetic) lines up with the saved output at line 29.
- **`K_eta` definition vs. canonical match**: the assertion at line 163 forces `K_eta = U_{RR} - d/dw[T_{w,R} R0'] + (1/2) T_{w,RR} (R0')^2`. I verified that this exactly cancels the residual `-(1/4) T_{w,RR}(R0')^2 - (1/2) U_{RR} + (1/2) d/dw[T_{w,R} R0']` left over after the cross-term IBP. The hand-set K_eta is a constructive bookkeeping check, but `L2_raw` is derived independently via `sp.diff(Lexp, eps, 2).subs(eps, 0) / 2`, so the assertion is non-tautological.
- **Y20 spherical-Laplacian eigenvalue**: I hand-computed `Δ_Ω(3 cos²θ - 1) = -6 (3 cos²θ - 1)` (the `φ`-derivative vanishes since m=0; the θ-derivative chain gives `-6 sinθ (3 cos²θ - 1)`, dividing by `sinθ` yields `-6 Y20`). The script's assertion `spherical_laplacian_Y20 + 6 Y20 == 0` holds.
- **Angular norm and stiffness**: `∫ sinθ Y20^2 dΩ = 1` (Y_lm normalization) and `∫ sinθ |∇_Ω Y20|^2 dΩ = l(l+1) ∫ sinθ Y20^2 dΩ = 6` (by sphere-IBP). Both confirm via direct integration.
- **Boundary-discharge limits**: I evaluated `-w/sqrt(1+w^2)` at `w → ±∞` to verify the `-2` discharge claim (M7), and confirmed `atan(w)` gives `π`. Both hold.
- **Mutation assertions**: each `assert_nonzero` flips a sign on a non-vanishing piece (e.g., `2*∂USigma/∂R` in M1, `K_eta` mid-term sign, Y20 eigenvalue 5 vs 6); none collapse to zero accidentally.

What didn't hold: the unit lacks the required second engine. Verdict is `findings` (not `clean`) and `stop_cold` is null — the missing-script finding is independently fixable and does not propagate sign errors downstream.

## Self-test notes

Walked the proposed `.wl` claims M1-M11 through the three traps. (1) Variable independence: M1's `D` operator requires `R` to be declared as `R[t,w,u,v]` and the chain rule must come from Mathematica's own dependency tracking, not from an `assert` on a pre-baked formula — the directive flags this explicitly so Codex uses `D[L, R[t,w,u,v]] - D[D[L, D[R[t,w,u,v], t]], t] - ...` or `VariationalD`, not a transliteration of `exact_el_lhs`. (2) Parity: M6's boundary discharges over `R = (-∞, ∞)` rely on the integrand having both `Limit[..., w -> ±Infinity]` go to 0 (Gaussian and Lorentzian decay); M7 deliberately uses a non-decaying probe `w sqrt(1+w^2)` so the limits are finite (-1, +1) and the boundary operator gives `-2`. (3) Trivial-case pre-check: for M9, substituting `th = Pi/2, ph = 0` into `Y20 = (1/4) Sqrt[5/Pi] (3 Cos[th]^2 - 1)` gives `-1/4 Sqrt[5/Pi]` (nonzero) and the Laplacian gives `+6/4 Sqrt[5/Pi]` — consistent with eigenvalue `-6`. Path specification: `.wl` lives under `mathematica/`, not `scripts/`; directive uses the absolute `mathematica/` path.
