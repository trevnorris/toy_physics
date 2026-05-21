---
unit_id: 001
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-20T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 001 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage001_geometry_lift_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage001_geometry_lift_mathematica_audit.txt`

Script mtimes: 2026-04-21 17:04. Output mtimes: 2026-05-11 (sympy 12:38, mathematica 12:40). Outputs are newer than scripts; no `stale_output` finding.

## What the script claims to verify

The Stage 001 scripts verify foundational bookkeeping for the moving-throat PDE skeleton: (i) orthonormality, zero-mean property, and spherical-Laplacian eigenvalues of the real-harmonic basis `{Y00, Y20, Y21c, Y21s, Y22c, Y22s}`; (ii) the mouth-average extraction rule `q00 = 2 sqrt(pi) delta_a` from the truncated expansion; (iii) the chain-rule sign of the linearized confinement variation `delta V_conf = -(V'_wall(Sigma0/ell_c)/ell_c) eta`; (iv) that the densitized modal-wall Lagrangian `L = (1/2) mu_eta q_t^2 - (1/2) T_w q_w^2 - (1/2) K_l q^2` reproduces the expected Euler-Lagrange form, and similarly for its `g(w)`-weighted version; (v) the `K_l = K_eta + ell(ell+1) T_Omega` restoring constant at `ell=0` and `ell=2`; (vi) the sourced-EL with `S_lm + f_ext` forcing; and (vii) two component equations for a representative localized-Maxwell action with `xi` gauge-fixing and source current.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1  | sympy | 90-93   | `integral(Y_lm dOmega) == 0` for {Y20, Y21c, Y21s, Y22c, Y22s} | yes |
| A2  | sympy | 95-98   | `integral(Y00^2 dOmega) - 1 == 0` | yes |
| A3  | sympy | 99-102  | `integral(Y20^2 dOmega) - 1 == 0` | yes |
| A4  | sympy | 103-106 | `integral(Y00 Y20 dOmega) == 0` | yes |
| A5  | sympy | 114     | `mouth_avg - q00/(2 sqrt(pi)) == 0` | yes |
| A6  | sympy | 125     | `-lap_s2(Y00) == 0` | yes |
| A7  | sympy | 126-127 | `-lap_s2(Y_lm) - 6 Y_lm == 0` for ell=2 modes | yes |
| A8  | sympy | 141     | `d/deps Vwall((Sigma0-eps eta)/ell_c)|_0 + eta V'_wall · 1/ell_c == 0` | yes |
| A9  | sympy | 167     | `EL[ldens] - target_dens == 0` | yes |
| A10 | sympy | 177     | `EL[g·ldens] - target_weighted == 0` | yes |
| A11 | sympy | 181-182 | `K_l|ell=0 - K_eta == 0`, `K_l|ell=2 - (K_eta + 6 T_Omega) == 0` | yes |
| A12 | sympy | 192     | `EL[ldens - q·source] - (target_dens - source) == 0` | yes |
| A13 | sympy | 228-229 | localized-Maxwell x- and w-component EL | yes (within "representative" scope) |
| M1  | mathematica | 105-112 | sphere integrals of Y_lm (avg, norm, cross) | yes |
| M2  | mathematica | 121     | `mouth average - q00/(2 sqrt(pi)) == 0` | yes |
| M3  | mathematica | 130-135 | `-lap_s2(Y_lm) - 6 Y_lm == 0` | yes |
| M4  | mathematica | 146     | linearized confinement variation | yes |
| M5  | mathematica | 165     | `EL[ldens] - target_dens == 0` | yes |
| M6  | mathematica | 170     | `EL[g·ldens] - target_weighted == 0` | yes |
| M7  | mathematica | 172-173 | `K_l|ell=0` and `K_l|ell=2` | yes |
| M8  | mathematica | 182     | sourced EL | yes |
| M9  | mathematica | 206-207 | localized-Maxwell components | yes |

Every assertion non-trivially exercises its claim. `target_*` symbols are independently constructed from the canonical PDE form (Newton-like or Maxwell-like) while `el_*` is computed by symbolic Euler-Lagrange machinery (`euler_equations` in SymPy, hand-rolled `eulerLagrange2D` in Mathematica). The two routes could diverge (sign, factor, missed derivative term) — they don't, but the check is real.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl:1-211`
- compared against `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py:60-229`

**What's wrong:**

The `.wl` script mirrors the `.py` script section-by-section and line-by-line, using the same variable choreography, the same intermediate quantities, and the same EL formula. It does not provide an independent derivation of the claims; it ports the SymPy algebra into Mathematica syntax.

Three matched excerpts:

1. **Modal-wall Lagrangian and EL (Section III.1).**

   SymPy (`scripts/.../stage001_..._sympy_audit.py:160-167`):
   ```
   ldens = (
       sp.Rational(1, 2) * mu_eta * sp.diff(q(t, w), t) ** 2
       - sp.Rational(1, 2) * t_w * sp.diff(q(t, w), w) ** 2
       - sp.Rational(1, 2) * k_l * q(t, w) ** 2
   )
   el_dens = euler_equations(ldens, q(t, w), [t, w])[0]
   target_dens = -mu_eta * sp.diff(q(t, w), t, 2) + sp.diff(t_w * sp.diff(q(t, w), w), w) - k_l * q(t, w)
   expect_zero("densitized Euler-Lagrange equation", el_dens.lhs - target_dens)
   ```

   Mathematica (`mathematica/.../stage001_..._mathematica_audit.wl:162-165`):
   ```
   ldens = (1/2) muEta D[qField, t]^2 - (1/2) tW D[qField, w]^2 - (1/2) kEll qField^2;
   elDens = FullSimplify[-eulerLagrange2D[ldens, qField, t, w], Assumptions -> $Assumptions];
   targetDens = -muEta D[qField, {t, 2}] + D[tW D[qField, w], w] - kEll qField;
   expectZero["densitized Euler-Lagrange equation", elDens - targetDens];
   ```

   Same Lagrangian, same `target`, same EL route. Mathematica re-implements `euler_equations` by hand in `eulerLagrange2D` (lines 59-67) and applies a unary minus (line 163) to absorb the sign convention. That manual sign flip is a fingerprint of porting rather than independent derivation.

2. **Confinement chain rule (Section II).**

   SymPy (lines 136-141):
   ```
   expr = vwall((sigma0 - eps * eta) / ell_c)
   first_var = sp.diff(expr, eps).subs(eps, 0)
   target = -eta * sp.diff(vwall(sigma0 / ell_c), sigma0)
   expect_zero("linearized confinement variation", first_var - target)
   ```

   Mathematica (lines 142-146):
   ```
   expr = Vwall[(sigma0 - eps etaMode)/ellc];
   firstVar = FullSimplify[D[expr, eps] /. eps -> 0, Assumptions -> $Assumptions];
   targetVar = -etaMode D[Vwall[sigma0/ellc], sigma0];
   expectZero["linearized confinement variation", firstVar - targetVar];
   ```

   Identical construction. `D[expr, eps] /. eps -> 0` is the Mathematica spelling of `sp.diff(expr, eps).subs(eps, 0)`.

3. **Maxwell EL (Section IV).**

   SymPy (lines 209-226):
   ```
   lmax = (
       sp.Rational(1, 2) * Zloc(w) * Fwx**2
       - sp.Rational(1, 2) * divA**2 / gauge_xi
       + mu0 * (Jx(x, w) * Ax(x, w) + Jw(x, w) * Aw(x, w))
   )
   el_Ax = (
       sp.diff(sp.diff(lmax, sp.diff(Ax(x, w), x)), x)
       + sp.diff(sp.diff(lmax, sp.diff(Ax(x, w), w)), w)
       - sp.diff(lmax, Ax(x, w))
   )
   ...
   target_Ax = sp.diff(Zloc(w) * Fwx, w) - sp.diff(divA, x) / gauge_xi - mu0 * Jx(x, w)
   ```

   Mathematica (lines 195-204):
   ```
   lmax = (1/2) zloc fwx^2 - divA^2/(2 gaugeXi) + mu0 (jxField axField + jwField awField);
   elAx = FullSimplify[
     D[D[lmax, D[axField, x]], x] + D[D[lmax, D[axField, w]], w] - D[lmax, axField],
     ...
   ];
   ...
   targetAx = D[zloc fwx, w] - D[divA, x]/gaugeXi - mu0 jxField;
   ```

   The Mathematica EL is computed by the same explicit `D[D[L, ∂L/∂(field_x)], x] + D[D[L, ∂L/∂(field_w)], w] - D[L, field]` formula as SymPy. No use of Mathematica's `VariationalD` / `VariationalMethods` package, no use of `EulerEquations` from `Variational`, no independent derivation via Hamilton's equations or by direct verification against the canonical Maxwell form `d_mu(Z F^{mu nu}) - (1/xi) d^nu(d.A) = mu0 J^nu`.

   The Mathematica script also does not exploit the symbolic identities Mathematica makes natural — e.g., it could verify the harmonic Laplacian eigenvalues by evaluating `SphericalHarmonicY[l, m, theta, phi]` directly and checking the eigenequation, instead of recomputing `1/sin θ d/dθ (sin θ d/dθ ...)` from scratch (lines 125-128), which exactly mirrors SymPy lines 119-123.

The second-engine policy requires independent derivation. This Mathematica script verifies the same algebra by the same route, so it functions only as a re-run of SymPy in a different language. A bug in the shared derivation pattern (e.g., a sign or factor that appears in both `target` and `EL[L]` because both were constructed from the same mental template) would be missed by both engines.

**Why this matters:**

This is a checkpoint stage (per the manifest) — the bar for second-engine independence is highest. If the two engines share derivation logic, they cannot serve as cross-checks. Downstream Stages 002+ build on these foundational results (modal-wall EL form, source forcing, Maxwell component equations); a shared error here would propagate silently because no truly independent check exists.

**Required change:**

Re-derive at least one substantive claim in each section using a Mathematica-native route that does not echo the SymPy algebra. Specifically:

1. **Section I (harmonic bookkeeping):** Replace the bespoke `y00, y20, ...` polynomial expressions (`.wl:77-82`) with `Sqrt[3] SphericalHarmonicY[2,0,theta,phi]`-style references where possible (or after converting from the complex spherical harmonics to the real basis). Verify the Laplacian eigenvalues via the identity `SphericalHarmonicY[l, m, theta, phi]` is an eigenfunction of the angular Laplacian with eigenvalue `-l(l+1)` — i.e., use Mathematica's built-in knowledge of these functions rather than re-deriving by `(1/sin theta) d/dtheta (sin theta d/dtheta ...)`. Verify orthonormality of the real basis by symbolic integration of pairs that should vanish vs. pairs that should normalize, using `Orthogonalize` or `KroneckerDelta`-based cross-check, not the same `sphereIntegral` pattern as SymPy.

2. **Section III (modal-wall EL):** Compute the Euler-Lagrange equation using Mathematica's `VariationalMethods` package: `Needs["VariationalMethods\`"]` then `EulerEquations[ldens, qField[t,w], {t, w}]`. Compare against `targetDens` constructed the same way as before. This replaces the hand-rolled `eulerLagrange2D` (lines 59-67) with an independent symbolic differentiator and removes the manual sign flip on line 163.

3. **Section IV (Maxwell EL):** Use `VariationalD[lmax, axField, {x, w}]` (with appropriate field-dependence declaration) from `VariationalMethods` instead of the explicit `D[D[L, D[axField, x]], x] + ...` pattern (lines 196-203). Same for `elAw`.

For sections where a Mathematica-native independent route does not exist (e.g., the confinement chain rule in II is so direct that any engine reduces to chain-rule on a single composite function), document inline that this section is intentionally a parallel check and not an independent re-derivation. That single concession should not extend to sections III and IV, which have genuine Mathematica-idiomatic alternatives.

**Verification:**

After the rewrite, the verifier confirms (a) the script still exits 0 with all `expectZero` checks passing, and (b) the patched sections use `EulerEquations` / `VariationalD` / `SphericalHarmonicY` (or equivalent built-in routes) at the spots flagged above. The new output file must show identical numerical residuals (all `= 0`) but the script body must contain calls to those independent operators rather than the hand-rolled mirrors of SymPy.

## Independent-derivation check (Mathematica)

The Mathematica script is a near-line-for-line port of the SymPy script. The same symbol names (modulo camelCase), the same Lagrangian structure, the same `target` constructions, the same explicit EL formula, the same `sphereIntegral` wrapper instead of using Mathematica's built-in spherical-harmonic identities. The unary minus on `eulerLagrange2D` (line 163, 168, 180) is a tell — it reconciles a SymPy/Mathematica sign convention mismatch by post-hoc negation rather than by independent derivation. See F1 for excerpt-level comparison.

## Engine cross-check

Both engines produce identical zero residuals on every shared assertion. Confinement variation: SymPy reports `-η · (∂/∂ξ₁ Vwall(ξ₁))|_{ξ₁=Σ₀/ell_c} / ell_c`; Mathematica reports `-((etaMode·Vwall'[sigma0/ellc])/ellc)` — same expression. All `expect_zero`/`expectZero` assertions report `= 0`. Numerical/symbolic agreement is total — but this agreement is weak evidence given F1: the two engines are running the same calculation in two languages.

## Verdict justification

The math is correct as far as the assertions exercise it. Spherical-harmonic norms, Laplacian eigenvalues, mouth-average extraction, the confinement chain-rule sign, the densitized and weighted EL forms, the source forcing, and the representative Maxwell component equations all hold. No tautological checks, no hardcoded numeric answers, no symbol-assumption errors, no missing branches relative to the stated scope, no stale outputs, no engine disagreement.

The single finding is structural: the Mathematica script is a transliteration of the SymPy script. For a non-checkpoint stage this might be acceptable; for a checkpoint, it undermines the second-engine policy. Verdict: `findings` (one finding, medium severity, no stop-cold).
